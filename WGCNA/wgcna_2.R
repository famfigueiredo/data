dds<- DESeqDataSetFromTximport(txi, coldata, design = ~ batch + Sex + BW)

keep <- rowSums(counts(dds)>=1) >= 30 #perform some prefiltering

dds <- dds[keep,]

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE) #transform while accounting for design 

# directory containing HTSeq count files
directory <-
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/full_dataset/readcounts'

load(file = '~/Documents/PhD/Papers/Paper III/data/RData/sampleTable_heart_group.RData')  # heart sampleTable

# creating dds object
sampleTable_wgcna_subset <-
  sampleTable_heart_group %>% filter(
    treatment %in% c('conu', 'ivld', 'eomes', 'gata3') &
      samplingPoint %in% c('10wpi', '4wpc', '6wpc')
  ) %>% dplyr::select('sample', 'filename', 'n', 'treatment', 'samplingPoint', 'lane')

sampleTable_wgcna_subset <- sampleTable_wgcna_subset %>%
  dplyr::mutate(treatment_sp = paste(treatment, sep = '_', samplingPoint))

dds_heart_wgcna <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_wgcna_subset,
  directory = directory,
  design = ~ treatment + samplingPoint
)

# removing low count genes (<10)
keep <-
  rowSums(counts(dds_heart_wgcna)) >= 10  

dds_heart_wgcna <-
  dds_heart_wgcna[keep, ]


# collapsing replicates (if applicable)
collapsed_heart <- collapseReplicates(dds_heart_wgcna,
                                      groupby = dds_heart_wgcna$n,
                                      run = dds_heart_wgcna$lane)

dds_heart <-
  DESeq(collapsed_heart, parallel = T)

colnames(dds) <- sampleTable_wgcna_subset$sample  # adding treatment names as columns

dds_conu <- dds[, dds$treatment == 'conu']  # subsetting data for 4 wpc

vsd <- varianceStabilizingTransformation(dds_conu, blind = FALSE)
wpn_vsd <- getVarianceStabilizedData(dds_conu)
rv_wpn <- rowVars(wpn_vsd)


q95_wpn <- quantile(rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset

expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

expr_normalized[1:10,1:10]

input_mat = t(expr_normalized)

input_mat[1:5,1:10]           # Look at first 5 rows and 10 columns


picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor     # Return cor function to original namespace

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

write_delim(module_df,
            file = "/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/datagene_modules_conu.txt",
            delim = "\t")


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


# pick out a few modules of interest here
modules_of_interest = c("green", "turquoise", "blue")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


# pull genes from modules of interest
submod = module_df %>% 
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# get normalized expression for the genes in the modules. Calculated with expr_normalized <- wpn_vsd[rv_wpn > q95_wpn,].
# this selects the 5% of the genes with the highest variance. filtering out genes with low variance.
expr_normalized[1:5, 1:10]

subexpr = expr_normalized[submod$gene_id,]  # get the genes associated with the modules of interest

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )  # formatting tibble to include gene_ids, treatment, eigenvalue, and corresponding module

# plotting normalized expression per treatment, per module
submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

# selecting genes associated to the modules of interest
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)


# getting treatment-specific expression data for genes associated with modules of interest
expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]


# only recalculate TOM (topological overlap matrix) for modules of interest. TOM quantifies the similarity between genes based on their gene expression patterns
# The Topological Overlap Matrix (TOM) is a key part of WGCNA that calculates
# the similarity between gene pairs based on both their direct co-expression and
# their relationships with other genes in the network. It is used to identify
# gene modules and better understand the underlying biology of gene
# co-expression patterns.
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)


edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

edge_list %>% na.omit()


# Filter edge list by TOM value threshold (e.g., TOM > 0.1)
threshold <- 0.1

# Keep only rows where correlation is above the threshold
filtered_edge_list <- edge_list %>% filter(correlation > threshold)

# Check the filtered edge list
filtered_edge_list %>% na.omit()


filtered_edge_list <- filtered_edge_list %>%
  filter(correlation > 0.3) %>%
  na.omit()


# Create a graph from the filtered edge list
g <- graph_from_data_frame(filtered_edge_list, directed = FALSE)

# Plot the network with the filtered edges
plot(g, vertex.size = 5, vertex.label.cex = 0.7, vertex.color = "lightblue", edge.width = 0.5)


# Extract the unique genes involved in the network
genes_in_network <- unique(c(filtered_edge_list$gene1, filtered_edge_list$gene2))

# Print the genes to check
print(genes_in_network)




# 5. Associating Modules and Phenotypes

nGenes = ncol(input_mat)
nSamples = nrow(input_mat)


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)


