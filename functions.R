# Helper function to create a Markdown table
markdown_table <- function(data) {
  # Ensure data is a data frame
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }
  
  # Limit data to 20 rows
  data <- head(data, 100)
  
  # Get the header
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Format rows
  rows <- apply(data, 1, function(row) {
    # Truncate long strings for readability
    formatted_row <- sapply(row, function(cell) {
      if (is.character(cell) && nchar(cell) > 50) {
        paste0(substr(cell, 1, 47), "...")
      } else {
        cell
      }
    })
    paste("|", paste(formatted_row, collapse = " | "), "|")
  })
  
  # Combine header, separator, and rows
  c(header, separator, rows)
}


# Improved helper function to select top and bottom 10
get_top_bottom_pathways <- function(data, n = 10, database) {
  # Get the top N pathways with positive NES
  top_n_pathways <- data %>%
    as_tibble() %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    top_n(n, wt = NES) %>%
    dplyr::select(Description, NES)
  
  # Get the bottom N pathways with negative NES
  bottom_n_pathways <- data %>%
    as_tibble() %>%
    filter(NES < 0) %>%
    arrange(-NES) %>%
    top_n(n, wt = -NES) %>%
    dplyr::select(Description, NES)
  
  # Combine the top and bottom N pathways
  combined_pathways <- bind_rows(top_n_pathways, bottom_n_pathways)
  
  # Create the object name based on the database parameter
  object_name <- paste0(database, "_pathways")
  
  # Assign the combined pathways tibble to a global environment object
  assign(object_name, combined_pathways, envir = .GlobalEnv)
  
  return(combined_pathways)
}

# running gsea starting from a DESeq results table. uses hsapiens orthologs
gsea_formatting <-
  function(results_table, tissue, treatment, sampling_point) {
    # Install and load required packages
    required_packages <-
      c('dplyr', 'gprofiler2', 'clusterProfiler', 'org.Hs.eg.db')
    installed_packages <- rownames(installed.packages())
    
    for (pkg in required_packages) {
      if (!(pkg %in% installed_packages)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
    
    # Convert rownames to column 'ensembl'
    results_df <-
      tibble::rownames_to_column(as.data.frame(results_table), var = 'ensembl')
    
    # Convert salmon genes to human orthologs
    orthologs <- gorth(
      query = rownames(results_table),
      source_organism = 'ssalar',
      target_organism = 'hsapiens',
      mthreshold = 1,
      filter_na = TRUE
    )
    
    # Select relevant variables and join with ortholog data
    merged_df <- results_df %>%
      left_join(orthologs, by = c('ensembl' = 'input')) %>%
      dplyr::select(ensembl,
                    ortholog_name,
                    ortholog_ensg,
                    log2FoldChange,
                    padj,
                    description) %>%
      na.omit()
    
    # Order genes by fold change
    ordered_df <- merged_df[order(-merged_df$log2FoldChange),]
    
    # Prepare matrix for GSEA
    gene_list <- ordered_df$log2FoldChange
    names(gene_list) <- ordered_df$ortholog_name
    
    # Prepare matrix for gsePathway
    ordered_entrez <-
      bitr(ordered_df$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)
    entrez_genes <-
      ordered_df %>% left_join(
        ordered_entrez,
        by = c('ortholog_name' = 'SYMBOL'),
        relationship = 'many-to-many'
      ) %>% dplyr::select(ENTREZID, log2FoldChange)
    distinct_genes <-
      entrez_genes %>% distinct(ENTREZID, .keep_all = T)
    entrez_gene_list <<- distinct_genes$log2FoldChange
    names(entrez_gene_list) <<- distinct_genes$ENTREZID
    
    # Run GSEA
    gsea_results <- gseGO(
      gene_list,
      keyType = 'SYMBOL',
      OrgDb = org.Hs.eg.db,
      ont = 'BP',
      pvalueCutoff = 0.05,
      pAdjustMethod = 'BH',
      verbose = T,
      eps = 1e-300,
      nPermSimple = 10000
    )
    
    # Assign the results to a variable including treatment and sampling_point in the name
    results_name <-
      paste0(tissue, '_', 'gsea_results_', treatment, '_', sampling_point)
    assign(results_name, gsea_results, envir = .GlobalEnv)
    
    return(gsea_results)
  }

# Significant genes grabs a DESeq2 result table, subsets the genes with padj < 0.05, and selects the ID, log2FoldChange, and padj columns. Also arranges log2FC in a descending manner.
significant_genes <- function(results_files) {
  b <- as.data.frame(subset(results_files, padj < 0.1)) %>%
    rownames_to_column(var = 'ID') %>%
    as_tibble()
  
  sig_genes <- b %>%
    dplyr::select(ID, log2FC = log2FoldChange, adjusted_p.val = padj, pvalue) %>%
    dplyr::arrange(desc(log2FC))
  
  return(sig_genes)
}

# Significant genes metrics creates a table with information about number of up/down significantly regulated genes
sig_genes_metrics <- function(significant_genes) {
  as.data.frame(significant_genes) %>% filter(adjusted_p.val < 0.1) %>%
    dplyr::mutate(Regulation = ifelse(
      log2FC < 0,
      'downregulated',
      ifelse(log2FC > 0, 'upregulated', NA)
    )) %>%
    filter(!is.na(Regulation)) %>%
    dplyr::count(Regulation)
}

# improved data wrangling function selects significantly regulated genes, converts ssalar IDs to hsapiens orthologs, and joins both
improved_data_wrangling <-
  function(results_table, treatment, sampling_point) {
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      install.packages("dplyr")
    }
    
    if (!requireNamespace("gprofiler2", quietly = TRUE)) {
      install.packages("gprofiler2")
    }
    
    # Load the required packages
    library(dplyr)
    library(gprofiler2)
    
    # select significant genes from DESeq2 results table
    a <-
      significant_genes(results_table)
    
    # convert ssalar gene IDs to human orthologs
    orth_hs <- gorth(
      query = a$ID,
      source_organism = 'ssalar',
      target_organism = 'hsapiens',
      mthreshold = 1,
      filter_na = T
    )
    
    # join significant genes table with human ortholog names
    results <- a %>% left_join(orth_hs, by = c('ID' = 'input')) %>%
      dplyr::select(.,
                    ID,
                    ortholog_name,
                    log2FC,
                    adjusted_p.val,
                    pvalue,
                    ortholog_ensg,
                    description)
    
    # getting gene lists for ORA
    # ora_up <<-
    #   results %>% drop_na() %>% dplyr::filter(log2FC > 0) %>%  pull(ortholog_ensg)
    # 
    # ora_down <<-
    #   results %>% drop_na() %>% dplyr::filter(log2FC < 0) %>%  pull(ortholog_ensg)
    
    # create results name
    results_name <-
      paste('results', treatment, sampling_point, sep = '_')
    
    # Add treatment column
    results <- results %>%
      mutate(treatment = treatment)
    
    # assign the results data frame to the dynamic name
    assign(results_name, results, envir = .GlobalEnv)
    
  }