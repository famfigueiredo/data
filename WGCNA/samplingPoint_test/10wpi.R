directory <-
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/full_dataset/readcounts'

# Heart ----
load(file = '~/Documents/PhD/Papers/Paper III/data/RData/sampleTable_heart_group.RData')  # heart sampleTable

# creating dds object
sampleTable_10wpi <-
  sampleTable_heart_group %>% filter(
    treatment %in% c('conu', 'ivld', 'eomes', 'gata3') &
      samplingPoint %in% c('10wpi')
  ) %>% dplyr::select('sample', 'filename', 'n', 'treatment', 'samplingPoint', 'lane')

head(sampleTable_10wpi)

sampleTable_10wpi$treatment <- droplevels(sampleTable_10wpi$treatment)
levels(sampleTable_10wpi$treatment) <- c('conu', 'ivld', 'eomes', 'gata3')


dds_10wpi <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_10wpi,
  directory = directory,
  design = ~ 1
)

counts_matrix <- counts(dds_10wpi)
df <- counts_matrix %>% 
  as.data.frame() %>% 
  dplyr::filter(rowSums(.) >= 50)

dds <- DESeqDataSetFromMatrix(
  countData = df,
  colData = sampleTable_10wpi,
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
                          power = 12, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          verbose = 3
                          # so we should set a seed
)


module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)


all.equal(sampleTable_10wpi$sample, rownames(module_eigengenes))

# Create the design matrix from the `treatment` variable
des_mat_treat <- model.matrix(~ sampleTable_10wpi$treatment)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat_treat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df)


module_6_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(sampleTable_10wpi %>%
                      dplyr::select(sample, treatment),
                    by = c('sample')
  )


ggplot(
  module_6_df,
  aes(
    x = treatment,
    y = ME6,
    color = treatment
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()




















