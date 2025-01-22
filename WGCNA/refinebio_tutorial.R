# Tutorial
# (https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html)

directory <-
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/full_dataset/readcounts'
# Heart ----
load(file = '~/Documents/PhD/Papers/Paper III/data/RData/sampleTable_heart_group.RData')  # heart sampleTable

# creating dds object
sampleTable_wgcna <-
  sampleTable_heart_group %>% filter(
    treatment %in% c('conu', 'ivld', 'eomes', 'gata3') &
      samplingPoint %in% c('10wpi', '4wpc', '6wpc')
  ) %>% dplyr::select('sample', 'filename', 'n', 'treatment', 'samplingPoint', 'lane')

head(sampleTable_wgcna)

sampleTable_wgcna$treatment <- droplevels(sampleTable_wgcna$treatment)
levels(sampleTable_wgcna$treatment) <- c('conu', 'ivld', 'eomes', 'gata3')
sampleTable_wgcna$samplingPoint <- droplevels(sampleTable_wgcna$samplingPoint)
levels(sampleTable_wgcna$samplingPoint) <- c('10wpi', '4wpc', '6wpc')
levels(sampleTable_wgcna$samplingPoint)


dds_heart_wgcna <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_wgcna,
  directory = directory,
  design = ~ 1
)

counts_matrix <- counts(dds_heart_wgcna)
df <- counts_matrix %>% 
  as.data.frame() %>% 
  dplyr::filter(rowSums(.) >= 50)

dds <- DESeqDataSetFromMatrix(
  countData = df,
  colData = sampleTable_wgcna,
  design = ~ 1
)

dds_norm <- vst(dds)


normalized_counts <- assay(dds_norm) %>% 
  t()

normalized_counts[1:10, 1:10]

sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)



sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
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


bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 10000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = 9, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          verbose = 3
                          # so we should set a seed
)


module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

all.equal(sampleTable_wgcna$sample, rownames(module_eigengenes))

# Create the design matrix from the `samplingPoint` variable
des_mat_sp <- model.matrix(~ sampleTable_wgcna$samplingPoint)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat_sp)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df)


module_7_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(sampleTable_wgcna %>%
                      dplyr::select(sample, samplingPoint),
                    by = c('sample')
  )


ggplot(
  module_7_df,
  aes(
    x = samplingPoint,
    y = ME7,
    color = samplingPoint
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


# The expression of module 7 decreases from pre- to post- challenge
# What genes are part of module 7?

gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

gene_module_key %>%
  dplyr::filter(module == "ME7")


make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = sampleTable_wgcna,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sample")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(sample, samplingPoint, n) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sample") %>%
    # Arrange by patient and time point
    dplyr::arrange(samplingPoint, n) %>%
    # Store sample
    tibble::column_to_rownames("sample")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    samplingPoint = col_annot_df$samplingPoint,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(samplingPoint = c("10wpi" = "#f1a340", "4wpc" = "#998ec3", '6wpc' = 'red'))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
}



mod_7_heatmap <- make_module_heatmap(module_name = "ME7")

# Print out the plot
mod_7_heatmap







# Create a design matrix from the `treatment` variable
des_mat_treat <- model.matrix(~ sampleTable_wgcna$treatment)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat_treat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df)





module_27_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(sampleTable_wgcna %>%
                      dplyr::select(sample, treatment),
                    by = c('sample')
  )


ggplot(
  module_27_df,
  aes(
    x = treatment,
    y = ME27,
    color = treatment
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()



module_21_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(sampleTable_wgcna %>%
                      dplyr::select(sample, treatment),
                    by = c('sample')
  )


ggplot(
  module_21_df,
  aes(
    x = treatment,
    y = ME21,
    color = treatment
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()

