#### Loading packages ####
suppressPackageStartupMessages({
  library('tidyverse')
  library('DESeq2')
  library('BiocParallel')
  library('org.Hs.eg.db')
  library('WGCNA')
  library('gridExtra')
  library('CorLevelPlot')
  library('gprofiler2')
  library('clusterProfiler')
  library('datapasta')
  register(MulticoreParam(10))
})

# clone repo in your Terminal (no '')
'git clone https://github.com/famfigueiredo/data.git'

# point to readcount directory
directory <-
  '~/readcounts/full_dataset'

# Heart ----
# load sampleTable
load(file = '~/RData/sampleTables/sampleTable_heart_group.RData')  # heart sampleTable

# creating dds object by selecting treatments of interest, sampling point, and relevant columns
sampleTable_4wpc <-
  sampleTable_heart_group %>% filter(
    treatment %in% c('conu', 'ivld', 'eomes', 'gata3') &
      samplingPoint %in% c('4wpc')
  ) %>% dplyr::select('sample', 'filename', 'n', 'treatment', 'samplingPoint', 'lane')


sampleTable_4wpc$treatment <- droplevels(sampleTable_4wpc$treatment)  # dropping factor levels

levels(sampleTable_4wpc$treatment) <- c('conu', 'ivld', 'eomes', 'gata3')  # re-setting treatment levels

sampleTable_4wpc$sample <- paste0("s", sampleTable_4wpc$sample)  # adding 's' to sample name


# Creating dds object from sampleTable
dds_4wpc <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_4wpc,
  directory = directory,
  design = ~ 1
)

counts(dds_4wpc)
#######
# # remove all genes with counts < 15 in more than 75% of samples (24*.75 = 18)
# ## suggested by WGCNA on RNAseq FAQ
# dds75 <- dds_4wpc[rowSums(counts(dds_4wpc) >= 15) >= 18,]
# nrow(dds75)
# 
# # perform variance stabilization
# dds_norm <- vst(dds75)
# assay(dds_norm) %>% 
#   head()

#######
# Extracting count matrix to filter out low count genes
counts_matrix <- counts(dds_4wpc)
df <- counts_matrix %>%
  as.data.frame() %>%
  dplyr::filter(rowSums(.) >= 50)

# Re-formatting into dds
dds <- DESeqDataSetFromMatrix(
  countData = df,
  colData = sampleTable_4wpc,
  design = ~ 1
)
nrow(dds)
# VST
dds_norm <- vst(dds)
# Extract normalized counts from dds object

normalized_counts <- assay(dds_norm) %>% 
  t()

normalized_counts[1:10, 1:10]

# Picking soft thresholding power for constructing scale-free network
sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed",
                         verbose = 3
)

# Calculate soft-thresholding powers
powers <- c(1:20) # Range of powers to test
sft <- pickSoftThreshold(normalized_counts, powerVector = powers, verbose = 5)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

a1 <- ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()


a2 <- ggplot(sft_df, aes(Power, mean.k., label = Power)) +
  geom_text(nudge_y = .1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

soft_threshold_heart4wpc <- arrangeGrob(a1, a2, nrow = 2)
plot(soft_threshold_heart4wpc)

normalized_counts[1:10, 1:10]

ggsave('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/WGCNA/samplingPoint_test/soft_threshold_heart4wpc.png', 
       soft_threshold_heart4wpc, width = 12, height = 10, dpi = 'retina')

# Network construction. Finding gene co-expression modules in WGCNA
bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 15000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          minModuleSize = 30, # minimum module size. Default is 30
                          power = 10, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use number module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation so we should set a seed
                          verbose = 3
)

save(bwnet, file = 'bwnet_4wpc.RData')  # saving bwnet for later
# load('bwnet_4wpc.RData')

table(bwnet$colors)  # checking module assignments
length(unique(bwnet$colors))  # counting unique modules

# merged <- mergeCloseModules(
#   normalized_counts, 
#   bwnet$colors, 
#   cutHeight = 0.25,  # Lower this value (e.g., 0.15 or 0.2) for more merging
#   verbose = 3
# )
# 
# # Replace original modules with merged results
# modules$colors <- merged$colors
# modules$MEs <- merged$newMEs

# Extracting eigengenes
module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

all.equal(sampleTable_4wpc$sample, rownames(module_eigengenes))  # checking if sample names match eigengenes

# Create the design matrix from the `treatment` variable
des_mat_treat <- model.matrix(~ sampleTable_4wpc$treatment)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat_treat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df)

sig_modules <- stats_df %>% filter(adj.P.Val < 0.05)
print(sig_modules)

binary_treatment <- model.matrix(~ treatment - 1, data = sampleTable_4wpc)  # binarizing phenotype data
rownames(binary_treatment) <- sampleTable_4wpc$sample

# define numbers of genes and samples
nSamples <- nrow(normalized_counts)
nGenes <- ncol(normalized_counts)

# Calculate module-trait relationships
module.trait.corr <- cor(module_eigengenes, binary_treatment, use = "p")  # MEs = module eigengenes
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

heatmap.data <- merge(module_eigengenes, binary_treatment, by = 'row.names')

heatmap.data <- 
  heatmap.data %>% column_to_rownames(var = 'Row.names')

names(heatmap.data)
names(heatmap.data)[39:42] <- c('CONU', 'IV-LD', 'EOMES', 'GATA-3')  # Renaming variables for plotting

corrplot_heart4wpc <- CorLevelPlot(
  heatmap.data,
  x = names(heatmap.data)[39:42],
  y = names(heatmap.data)[1:38],
  col = c('blue', 'lightblue', 'white', 'darkorange', 'red'),
  main = 'Module-trait relationship\nHeart, 4wpc',
  cexMain = 2,
  titleX = 'Treatments',
  titleY = 'Modules',
  rotTitleY = 90,
  cexTitleX = 2,
  cexTitleY = 2,
  fontLabX = 1,
  fontLabY = 1,
  fontTitleY = 1,
  fontTitleX = 1,
  fontMain = 1
)

lattice::trellis.par.set(
  fontsize = list(text = 12),  # Set font size
  fontfamily = "Times"        # Set font family (e.g., Times)
)  # changing font didn't work. size did, though.

# Saving a trellis object: Render using graphics device and explicitly call print() to display the plot within the device
png('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/WGCNA/samplingPoint_test/corrplot_heart4wpc.png', 
    width = 12, height = 10, units = "in", res = 300)
print(corrplot_heart4wpc) # Render the trellis plot
dev.off() # Close the device

stats_df[1:10, 1:5]

# By plotting the module-trait correlations in a heatmap, we can identify which modules are significant for each treatment.
# Now we can look into the interesting modules, to see how they perform across treatments.
# Afterwards, we can extract the genes belonging to said modules, for enrichment.

# Sample assignment to module eigengenes
ME_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(sampleTable_4wpc %>%
                      dplyr::select(sample, treatment),
                    by = c('sample')
  )

ggplot(ME_df,
       aes(x = treatment,
           y = ME10,
           color = treatment)) +  # Closing parenthesis for ggplot()
  # A boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic() +
  labs(x = 'Treatment Group',
       y = 'Module Eigengene',
       title = 'ME10 Analysis plot') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_x_discrete(labels = c('CONU', 'IV-LD', 'EOMES', 'GATA-3')) +
  # Center title and remove legend
  theme(
    legend.position = "none",
    # Removes legend
    axis.text = element_text(size = 12),
    # Axis text size
    axis.title = element_text(size = 14),
    # Axis label size
    plot.title = element_text(# Center title
      size = 16, hjust = 0.5, face = "bold")  # hjust = 0.5 centers the title
    )

ggsave('~/Documents/PhD/Papers/Paper III/data/WGCNA/samplingPoint_test/ME21_analysis_plot.png',
       width = 10, height = 7, dpi = 'retina')


ME_df %>% 
  filter(treatment == 'gata3') %>% 
  dplyr::select(ME21) %>% 
  summarise(median_ME21 = median(ME21))

?chooseTopHubInEachModule

normalized_counts[] <- sapply(normalized_counts, as.numeric)

# Identify hub genes in a module
## Filter genes based on module membership (kME)
### Select the module of interest (e.g., "21")
MEs <- bwnet$MEs
module <- "ME21"
# moduleGenes <- (moduleColors == module)

# Calculate module membership (kME)
MM <- as.data.frame(cor(normalized_counts, MEs[, paste0("ME", module)], use = "p"))
colnames(MM) <- "kME"

# Order genes by kME values (highest to lowest)
hubGenes <- MM[order(-MM$kME), , drop = FALSE]
head(hubGenes, 10)  # Top 10 hub genes


gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear
  dplyr::mutate(module = paste0("ME", module))

# Extracting genes from ME21
ME21 <- gene_module_key %>%
  dplyr::filter(module == 'ME21')

# Saving ME21 genes
readr::write_tsv(ME21,
                 file = file.path('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/WGCNA/samplingPoint_test/ME21.tsv')
)

# Converting  genes from ssalar to hsapiens orthologs
ME21_orth <- gorth(
  query = ME21$gene,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

ME21_orth[1:20, 1:7]

result_ME21 <- gost(ME21_orth$ortholog_name, organism = 'hsapiens')

p <- gostplot(result_ME21, capped = TRUE, interactive = FALSE)

# Filter for GO:BP terms
biological_processes <- as_tibble(result_ME21$result[result_ME21$result$source == 'GO:BP', ])
reactome_processes <- as_tibble(result_ME21$result[result_ME21$result$source == 'REAC', ])


as_tibble(result_ME21$result)
# View the filtered results
head(biological_processes)
reactome_processes

# Selecting terms of interest for Manhattan plot
highlight_terms = c('REAC:R-HSA-977606', 'REAC:R-HSA-174577', 'REAC:R-HSA-166658', 'REAC:R-HSA-173736', 'GO:0006959', 'GO:0002455', 'GO:0006958', 'GO:0006956', 'GO:0006957', 'GO:0016064', 'GO:0019724', 'GO:0002449', 'GO:0016064', 'GO:0050776', 'GO:0002253', 'GO:0002252')

pp <- publish_gostplot(p, highlight_terms = highlight_terms,
                       width = NA, height = NA, filename = NULL )


pp + ggtitle('Manhattan plot') +
  theme(
  legend.position = "none",
  # Removes legend
  axis.text = element_text(size = 12),
  # Axis text size
  axis.title = element_text(size = 14),
  # Axis label size
  plot.title = element_text(# Center title
    size = 16, hjust = 0.5, face = "bold")  # hjust = 0.5 centers the title
)

# Rename query for better annotation (optional)
pathway_results <- gost(query = ME21_orth$ortholog_name, organism = "hsapiens", sources = c('GO:BP', 'REAC'))
pathway_results$result <- pathway_results$result %>% 
  mutate(query = str_replace(query, 'query_1', 'ME21 functional analysis'))
p <- gostplot(pathway_results, capped = TRUE, interactive = FALSE)
pp <- publish_gostplot(p, highlight_terms = highlight_terms,
                       width = NA, height = NA, filename = NULL)

# Saving Manhattan plot
ggsave('~/Documents/PhD/Papers/Paper III/data/WGCNA/samplingPoint_test/manhattan_plot_ME21.png',
       width = 15, height = 12, dpi = 'retina')


# ORA ----
ME21 <- read.table('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/WGCNA/samplingPoint_test/ME21_heart4wpc.tsv', header = T)

ME21_orth <- gorth(
  query = ME21$gene,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

go_enrich <- enrichGO(gene = ME21_orth$ortholog_name,
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

dotplot(go_enrich,
        x = 'Count',
        label_format = 40)


library(enrichplot)
upsetplot(go_enrich)

install.packages('wordcloud')
library(wordcloud)

wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Set1"), max.words = 25)

install.packages('wesanderson')
library(wesanderson)

wordcloud(
  words = wcdf$term, 
  freq = wcdf$V1, 
  col=terrain.colors(length(wcdf$term) , alpha=0.9) , rot.per=0.3 
)

go_df <- as_tibble(go_enrich) %>%
  arrange((p.adjust)) %>%
  top_n(20, wt = p.adjust) %>%
  mutate(Count = sapply(strsplit(as.character(geneID), '/'), length))

dotplot(
  go_enrich,
  showCategory = 10,
  size = NULL,
  split = NULL,
  font.size = 12,
  title = "",
  label_format = 30)


# GSEA ---- 
source('~/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/functions.R')
gsea_formatting(heart_res_gata3_vs_conu_4wpc, 'heart', 'gata3', '4wpc')
heart_gsea_simplified_results_gata3_4wpc <-
  simplify(heart_gsea_results_gata3_4wpc)

heart_entrez_gene_list_gata3_4wpc <- entrez_gene_list

top10_high_nes <-
  as_tibble(heart_gsea_simplified_results_gata3_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(20, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_results_gata3_4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(20, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_gata3_4wpc <-
  bind_rows(top10_high_nes, bottom10_low_nes)

as_tibble(heart_gsea_simplified_results_gata3_4wpc) %>% pull(ID) %>% intersect(., highlight_terms)  # no intersected terms
as_tibble(heart_gsea_results_gata3_4wpc) %>% pull(ID) %>% intersect(., highlight_terms)  # no intersected terms

low_high_nes_gata3_4wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(
    aes(size = setSize),
    shape = 1,
    stroke = 0.2,
    color = 'red'
  ) +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 300)) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_gata3_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_gata3_4wpc$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'GATA3, 4WPC, heart tissue') +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = 'Times New Roman'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm'),
    strip.text = element_text(size = 24),
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(
      color = 'black',
      linewidth = .05,
      linetype = 2
    )
  ) + guides(color = guide_legend(override.aes = list(size = 5)),
             # increase point size in gene count legend
             size = guide_legend(override.aes = list(
               shape = 1,
               fill = NA,
               stroke = .5,
               color = 'red'
             ))) +
  facet_grid(. ~ Regulation)

# library(ReactomePA)
# y_gata3_4wpc <-
#   gsePathway(
#     heart_entrez_gene_list_gata3_4wpc,
#     # the gsea_formatting function removes the duplicates from this object
#     pvalueCutoff = .2,
#     pAdjustMethod = 'BH',
#     eps = 1e-300,
#     nPermSimple = 100000,
#     verbose = F
#   )
# 
# as_tibble(y_gata3_4wpc) %>% arrange(-NES) %>% print(n = 100)
# 
# gata3_4wpc_pathways <- as_tibble(y_gata3_4wpc) %>% 
#   arrange(NES) %>% 
#   dplyr::select(., Description, NES) %>%
#   mutate(NES = sprintf('%.3f', NES))  # format NES to 3 decimal places


gsea_formatting(heart_res_gata3_vs_eomes_4wpc, 'heart', 'GATA3', '4wpc')
heart_gsea_simplified_results_GATA3_4wpc <-
  simplify(heart_gsea_results_GATA3_4wpc)

# No significantly regulated terms between GATA3 and EOMES

############################################################################################################


# lightyellow_df <- module_eigengenes %>%
#   tibble::rownames_to_column("sample") %>%
#   # Here we are performing an inner join with a subset of metadata
#   dplyr::inner_join(sampleTable_4wpc %>%
#                       dplyr::select(sample, treatment),
#                     by = c('sample')
#   )
# 
# darkgreen <- gene_module_key %>%
#   dplyr::filter(module == "MEdarkgreen")
# 
# readr::write_tsv(darkgreen,
#                  file = file.path('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/WGCNA/samplingPoint_test/darkgreen.tsv')
# )
# 
# 
# darkgreen_orth <- gorth(
#   query = darkgreen$gene,
#   source_organism = 'ssalar',
#   target_organism = 'hsapiens',
#   mthreshold = 1,
#   filter_na = T
# )
# 
# readr::write_tsv(darkgreen_orth,
#                  file = file.path('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/WGCNA/samplingPoint_test/darkgreen_orth.tsv')
# )
# 
# ggplot(
#   lightyellow_df,
#   aes(
#     x = treatment,
#     y = MElightyellow,
#     color = treatment
#   )
# ) +
#   # a boxplot with outlier points hidden (they will be in the sina plot)
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   # A sina plot to show all of the individual data points
#   ggforce::geom_sina(maxwidth = 0.3) +
#   theme_classic()
# 
# 
# # MEdarkgreen is higher in gata3, while MElightyellow is higher in EOMES
# 
# normalized_counts[] <- sapply(normalized_counts, as.numeric)
# module.gene.mapping <- as.data.frame(bwnet$colors)
# darkred.module.genes <- module.gene.mapping %>%
#   filter(`bwnet$colors` == 'darkred') %>%
#   rownames()
# 
# # Identify hub genes in a module
# ## Filter genes based on module membership (kME)\
# 
# ### Select the module of interest (e.g., "green")
# MEs <- bwnet$MEs
# module <- "darkred"
# # moduleGenes <- (moduleColors == module)
# 
# # Calculate module membership (kME)
# MM <- as.data.frame(cor(normalized_counts, MEs[, paste0("ME", module)], use = "p"))
# colnames(MM) <- "kME"
# 
# MM[1:20, ]
# 
# # Order genes by kME values (highest to lowest)
# hubGenes <- MM[order(-MM$kME), , drop = FALSE]
# head(hubGenes, 10)  # Top 10 hub genes
# 
# library(igraph)
# 
# # Extract adjacency matrix for the selected module
# adjacency <- adjacency(normalized_counts[, darkred.module.genes], power = 10)
# 
# # Convert adjacency matrix to graph object
# graph <- graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE)
# 
# # Filter for top hub genes (e.g., top 20)
# topHubGenes <- rownames(head(hubGenes, 50))
# 
# 
# subgraph <- induced_subgraph(graph, vids = topHubGenes)
# 
# 
# # Plot the network
# plot(subgraph,
#      vertex.size = 10,
#      vertex.label.cex = 0.7,
#      edge.width = E(subgraph)$weight * 2,
#      main = paste("Top 50 Hub Genes in", module, "Module"))
# 
# 
# 
# 
# ############################
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("CEMiTool")
# library('CEMiTool')
# 
# # Extracting count matrix to filter out low count genes
# counts_matrix <- counts(dds_4wpc)
# df <- counts_matrix %>% 
#   as.data.frame() %>% 
#   dplyr::filter(rowSums(.) >= 50)
# 
# # Re-formatting into dds
# dds <- DESeqDataSetFromMatrix(
#   countData = df,
#   colData = sampleTable_4wpc,
#   design = ~ 1
# )
# 
# # VST
# dds_norm <- vst(dds)
# 
# # Extract normalized counts from dds object
# normalized_counts <- assay(dds_norm)
# 
# head(normalized_counts[, 1:10])
# 
# cem <- cemitool(as.data.frame(normalized_counts))
# 
# cem
# 
# nmodules(cem)
# head(module_genes(cem))
# 
# hubs <- get_hubs(cem, 10)
# 
# generate_report(cem)
# 
# save_plots(cem, 'all')
# sample_annot <- sampleTable_4wpc %>% dplyr::select(., sample, treatment)
# sample_annot <- sample_annot %>% dplyr::rename(., sampleName = sample, Class = treatment)
# 
# cem <- cemitool(as.data.frame(normalized_counts), sample_annot)
# 
# 
# cem <- cemitool(
#   as.data.frame(normalized_counts),
#   sample_annot,
#   sample_name_column = "sampleName"
# )
# 
# cem <- mod_gsea(cem)
# cem <- plot_gsea(cem)
# 
# show_plot(cem, 'gsea')
# 
# 
# # plot gene expression within each module
# cem <- plot_profile(cem)
# plots <- show_plot(cem, "profile")
# plots[9]







