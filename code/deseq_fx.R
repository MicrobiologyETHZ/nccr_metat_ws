# data_dir: salmon file directory
require(tidyverse)
require(data.table)
require(DESeq2)

##############################

align_samples_and_counts <- function(count_data, sample_data, conditions, gene_col, sample_col) {
  count_data <- count_data %>% as.data.frame()
  rownames(count_data) <- count_data[[gene_col]]
  sample_data <- sample_data %>% filter(sample_id %in% colnames(count_data))
  conditions <- unlist(strsplit(conditions, ","))
  sample_data <- sample_data %>% unite("group", all_of(conditions), remove = FALSE)
  count_data <- count_data[, sample_data[[sample_col]]]
  return(list("count_data" = count_data, "sample_data" = sample_data))
}

deseq_norm_mat <- function(raw_data) {

  # Adapted from "How to normalize metatranscriptomic
  # count data for differential expression analysis
  # raw_data holds $count_data -> raw count data for 1 species
  # and $sample_data -> sample data matching the samples in count_data 
  count_data <- raw_data$count_data
  col_data <- raw_data$sample_data
  count_data_row <- nrow(count_data)
  count_data <- round(count_data) # Makes sure it is all integers
  if (sum(rowSums(count_data ==0) == 0) != 0) {

  
    dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~group)
    colData(dds)$group = factor(colData(dds)$group, levels=unique(colData(dds)$group))
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    dds <- DESeq(dds, quiet = FALSE)
    return(dds)}
    #normalize the data
    #norm_count <- count_data/rep(dds@colData@listData$sizeFactor, each = (count_data_row))
    #return(norm_count)}

  else {
    return(data.frame())
  }
}

deseq_norm_by_taxon <- function(count_data, sample_data, conditions, gene_col, sample_id_col, taxon_col){
  taxa <- count_data[[taxon_col]] %>% unique()
  print(paste0("Normalizing within each taxon. Number of taxa: ", length(taxa)))
  df_list <- list()
  for (i in seq_along(taxa)){
    print(taxa[[i]])
    taxa_data <- count_data %>% filter(!!sym(taxon_col) == taxa[[i]]) %>% select(-!!sym(taxon_col))
    raw_data <- align_samples_and_counts(taxa_data, sample_data, conditions, gene_col, sample_id_col)
    df_list[[i]] <- deseq_norm_mat(raw_data) %>% rownames_to_column(gene_col)
  }
  print("Finished normalization")
  return (rbindlist(df_list, fill=TRUE))
}

deseq_de_on_full <- function(fdf, sample_data, conditions, gene_col, sample_id_col) {

    #fdf: concatenated dataaframes that were normalized by taxa
    rownames(fdf) <- fdf[[gene_col]]
    fdf <- fdf %>% select(-!!sym(gene_col)) %>% round %>% rownames_to_column(gene_col)
    raw_data <- align_samples_and_counts(fdf, sample_data, conditions, gene_col, sample_id_col)
    dds <- DESeqDataSetFromMatrix(countData = raw_data$count_data, colData = raw_data$sample_data, design = ~group)
    colData(dds)$group = factor(colData(dds)$group, levels=unique(colData(dds)$group))
    #stop DESeq2 from performing additional normalization
    normFactors <- matrix(1, ncol =ncol(raw_data$count_data), nrow = nrow(raw_data$count_data))
    normalizationFactors(dds) <- normFactors
    print("Running DE analysis")
    dds <- DESeq(dds)
    return(dds)
    }


deseq_on_metat_taxon <- function(count_data, sample_data, conditions, gene_col, sample_id_col, taxon_col){
    
    # Step 1: Normalize by taxon
    fdf <- deseq_norm_by_taxon(count_data, sample_data, conditions, gene_col, sample_id_col, taxon_col)
    # Run DESeq
    res <- deseq_de_on_full(fdf, sample_data, conditions, gene_col, sample_id_col)
    results <- list("norm_counts"=fdf, "dds"=res)
    return(results)
}

###################################


run_deseq_on_featcnts <- function(sample_data_file, count_file,
                                  conditions, outDir, prefix, rundeseq = TRUE, filter_value=10) {

  # Validating the data
  print('Validationg data')
  sample_data <- read.csv(sample_data_file, header = TRUE)
  conditions <- unlist(strsplit(conditions, ","))
  sample_data <- sample_data %>% unite("group", conditions, remove = FALSE)

  # Importing the data
  cnts <- read_csv(count_file)
  cnts <- cnts %>% as.data.frame()
  rownames(cnts) <- cnts$ID
  sample_data <- sample_data %>% filter(sample_id %in% colnames(cnts))
  print(dim(cnts))
  # Differential analysis
  print('Diff Analysis')
  cnts <- cnts[, sample_data$sample_id]
  rownames(sample_data) <- sample_data$sample_id
  if (rundeseq == TRUE) {
    dds <- DESeqDataSetFromMatrix(
      countData = cnts,
      colData = sample_data,
      design = ~group
    )
    print('Filtering')
    keep <- rowSums(counts(dds)) >= filter_value
    dds <- dds[keep, ]
    dds <- estimateSizeFactors(dds)
    norm_cnts <- counts(dds, normalized = TRUE)
    dds <- DESeq(dds)
    vsd <- vst(dds)
    write.csv(assay(vsd), file.path(outDir, paste0(prefix, "featCnts-vsd.csv")))
    write.csv(norm_cnts, file.path(outDir, paste0(prefix, "norm_cnts.csv")))
    return(list("dds" = dds, "sample_data" = sample_data))
  }
  return(list())
}


run_deseq_on_salmon <- function(sample_data_file, data_dir, tx2gene_file,
                                conditions, outDir, prefix, rundeseq = TRUE) {

  # Validating the data

  sample_data <- read.csv(sample_data_file, header = TRUE, row.names = 1)
  samples <- rownames(sample_data)
  sample_files <- file.path(data_dir, paste0(samples, "_quant"), "quant.sf")
  names(sample_files) <- samples
  tx2gene <- read.csv(tx2gene_file)
  conditions <- unlist(strsplit(conditions, ","))

  sample_data <- sample_data %>% unite("group", conditions, remove = FALSE)

  # Importing the data
  txi <- tximport(sample_files, type = "salmon", tx2gene = tx2gene)
  print(dim(txi$counts))
  write.csv(txi$abundance, file.path(outDir, paste0(prefix, "salmon-gene-tpms.csv")))

  # Differential analysis

  if (rundeseq == TRUE) {
    dds <- DESeqDataSetFromTximport(txi, sample_data, ~group)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    dds <- DESeq(dds)
    vsd <- vst(dds)

    write.csv(assay(vsd), file.path(outDir, paste0(prefix, "salmon-vsd.csv")))
    return(list("dds" = dds, "sample_data" = sample_data))
  }
  return(list())
}


get_contrasts <- function(sample_data, pattern1, pattern2) {
  first_contrasts <- unique(sample_data$group[grepl(pattern1, sample_data$group)])
  second_contrasts <- unique(sample_data$group[grepl(pattern2, sample_data$group)])
  first_contrasts <- setdiff(first_contrasts, second_contrasts)
  expand.grid(first_contrasts, second_contrasts)
}



get_results <- function(dds, c1, c2, outDir, lfct, alpha, prefix) {
  res <- results(dds, contrast = c("group", c1, c2), alpha = alpha, lfcThreshold = lfct)
  res$contrast <- paste0(c1, "_vs_", c2)
  res <- as.data.frame(res) %>% rownames_to_column("ID")
  write.csv(res, file.path(outDir, paste0(prefix, c1, "_vs_", c2, "_l", lfct, "a", alpha, "_results.csv")),
    row.names = FALSE
  )
  return(res)
}


get_all_results <- function(sample_data, pattern1, pattern2, dds, outDir,lfct, alpha, prefix) {
  contrasts <- get_contrasts(sample_data, pattern1, pattern2)
  f1 <- contrasts$Var1
  f2 <- contrasts$Var2
  for (var in 1:length(f1)) {
    print(paste0(as.character(f1[var]), " vs. ", as.character(f2[var])))
    res <- get_results(dds, as.character(f1[var]), as.character(f2[var]), outDir,lfct, alpha, prefix)
  }
}
