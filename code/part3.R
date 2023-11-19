
# Load the libraries

library(tidyverse)
library(stringr)
library(stats)
library(stats4)
library(caret)
library(ggplot2)

# Read data

# METAT: metatranscriptomic data
# METAG: metagenomic data

# For models M1 and M2 we are only using METAT data. For model M4, we will need both METAT and METAG data

metat <- read.table("datasets/part3/true-exp/true-exp.mtx_abunds.tsv.gz", header = TRUE,
                 sep = "\t", row.names=1)

metag <- read.table("datasets/part3/true-exp/true-exp.mgx_abunds.tsv.gz", header = TRUE,
                    sep="\t", row.names=1)

# Look at the data
# metat table contains Phenotype as well as readcounts for 100 samples
# each `gene name` is composed of 'species' name (ex. BUG0001), and 'orthologous group' name (ex. GROUP000001)
metat[0:5, 0:5]

# Ground truth data contains a list of `genes` that are differentially expressed between phenotypes

ground_truth <- read.table("datasets/part3/true-exp/true-exp.mtx_spiked.tsv.gz", 
                           header = FALSE,sep = "\t", 
                           col.names = c('gene_name', 'GT'))
head(ground_truth)
dim(ground_truth)


# Preprocess the data

# Create metatdata table (i.e. phenotype  for each sample)
metat_pheno <- metat[1:2, ] %>%
  t() %>%
  data.frame() 

# Remove metadata from METAT column
metat <- metat[-c(1, 2), ]

# Save METAG sequencing depth, we will need it later
metag_seqdepth <- metag[2,] %>% t()
colnames(metag_seqdepth) <- c('metag_seqdepth')

# Remove metadata from METAG table
metag <- metag[-c(1, 2), ]


# Assign each 'transcript' to a 'species'
metat$bug <- str_split_fixed(rownames(df), "_", 2)[, 1]

# Create table with 'taxonomic' abundances based on METAT data 
# We will need this for model M2.

bug_metat <- metat %>%
  group_by(bug) %>%
  summarize(across(all_of(sample_cols), sum)) %>% 
  pivot_longer(cols = -bug, names_to = "sample_id", values_to = "bug_abundance") 

head(bug_metat)

# Get sample names
sample_cols <- grep("SAMPLE", names(metat), value = TRUE)


# Go through example of building all the models for one gene 

example_gene = "BUG0013_GROUP001489"
example_gene_bug <- str_split_fixed(example_gene, "_", 2)[, 1]

# Get METAT data for the gene
gene_df <- metat[example_gene, sample_cols] %>% t() %>% as.data.frame()

# Get METAG data for the gene
gene_metag <- metag[example_gene, sample_cols] %>% t() %>% as.data.frame() 
colnames(gene_metag) <- c("metag")

#gene_df[[example_gene]] %>% hist()

# Get species abundance
gene_bug_abundance <- bug_metat %>% filter(bug == example_gene_bug) %>% 
                      column_to_rownames('sample_id')

# Put all of this together with phenotype data
gene_df <- cbind(metat_pheno, gene_df, gene_metag, gene_bug_abundance, metag_seqdepth)

head(gene_df)

# Filtering

# Semi-strict: Remove sample if both gene abundance and bug abundance are 0 (for models M1 and M2)
gene_df_semi <- gene_df %>% filter((!!as.symbol(example_gene)) > 0 | bug_abundance > 0) # 63 samples

# For M4, we filter based on METAG abundance of that gene
gene_df_semistrict_metag<- gene_df %>% filter((!!as.symbol(example_gene)) > 0 & metag > 0)

# Strict: Remove sample if either gene abundance and species abundance are 0 (for models M1 and M2)
gene_df_strict <- gene_df %>% filter((!!as.symbol(example_gene)) > 0 & bug_abundance > 0) # 47 samples

# For M4, we filter based on METAG abundance of that gene
gene_df_strict_metag<- gene_df %>% filter((!!as.symbol(example_gene)) > 0 & metag > 0)

# Normalization / Models
# We will do this example using the 'strict' filtered dataframe

##  Model 1: Total sum scaling

M1_col <- paste0(example_gene, "_M1")

# `SeqDepth` is METAT total sequencing depth 
gene_df_strict['tss'] <- gene_df_strict[example_gene]/gene_df_strict['SeqDepth']

# Add pseudo count and log transform
pseudo_count <- min(gene_df_strict[gene_df_strict[, 'tss'] > 0, 'tss']) / 2
gene_df_strict[M1_col] <- log10(gene_df_strict[, 'tss'] + pseudo_count)


## Model 2: Within taxon normalisation

M2_col <- paste0(example_gene, "_M2")
# Normalising by species abundance
gene_df_strict[M2_col] <- gene_df_strict[example_gene]/gene_df_strict['bug_abundance']
# Add pseudo count and log transform
pseudo_count <- min(gene_df_strict[gene_df_strict[, M2_col] > 0, M2_col], na.rm=TRUE) / 2
gene_df_strict[M2_col] <- log10(gene_df_strict[, M2_col] + pseudo_count)


## Model 4: RNA/DNA ratio

M4_col <- paste0(example_gene, "_M4")
# Note using a different dataframe here, because filtering was done differently
# First perform total sum scaling for both METAT and METAG data
gene_df_strict_metag['tss'] <- gene_df_strict_metag[example_gene]/gene_df_strict_metag['SeqDepth']
gene_df_strict_metag['metag_tss'] = gene_df_strict_metag['metag']/gene_df_strict_metag['metag_seqdepth']

# Calculate ratio of METAT/METAG abundance
gene_df_strict_metag[M4_col] <- gene_df_strict_metag['tss']/gene_df_strict_metag['metag_tss']
# Add pseudocount and log transform
pseudo_count <- min(gene_df_strict_metag[gene_df_strict_metag[, M4_col] > 0, M4_col], na.rm=TRUE) / 2
gene_df_strict_metag[M4_col] <- log10(gene_df_strict_metag[, M4_col] + pseudo_count)



# Run linear model

fixed_effects <- c("Phenotype")

# This function gets the coefficients, std error and p-values from the fit linear model
summary_function <- function(fit) {
  lm_summary <- summary(fit)$coefficients
  para <- as.data.frame(lm_summary)[-1, -3]
  para$name <- rownames(lm_summary)[-1]
  return(para)
}

# For model M1 we are modelling log trnasformed TSS expression values (M1_col) against the phenotype

# For model M2 are modelling log transformed taxon normalised expression values (M2_col) against the phenotype

model_cols <- c(M1_col, M2_col)
models <- list()
for (model_col in model_cols){
    formula <- paste0(model_col, " ~ ", paste(fixed_effects, collapse = " + "))
    print(formula)
    md <- glm(formula = formula, data = gene_df_strict, family = gaussian()) %>%
      suppressWarnings(suppressMessages(update(., .~.)))
    models[[as.symbol(model_col)]] <- summary_function(md)
}

# For model M4 we are modelling RNA/DNA ratio (M4_col) against the phenotype

formula <- paste0(M4_col, "~", "Phenotype")
print(formula)
md <- glm(formula = formula, data = gene_df_strict_metag, family = gaussian()) %>%
  suppressWarnings(suppressMessages(update(., .~.)))
models[[as.symbol(M4_col)]] <- summary_function(md)
summary <- do.call(rbind, models)



# To run this on all of the genes, we're putting this together into a function

# Run the model for each gene

model_gene <- function(gene_df, expr_col, filter_strategy='strict', fixed_effects){
  
  # Filter 
  if (filter_strategy == 'strict'){
    gene_df <- gene_df %>% filter((!!as.symbol(expr_col)) > 0 & bug_abundance > 0)
    gene_df_metag <- gene_df %>% filter((!!as.symbol(expr_col)) > 0 & metag > 0)
  }
  else if (filter_strategy == 'semistrict'){
    gene_df <- gene_df %>% filter((!!as.symbol(expr_col)) > 0 | bug_abundance > 0)
    gene_df_metag <- gene_df %>% filter((!!as.symbol(expr_col)) > 0 | metag > 0) 
  }
  else {print("No filtering applied")}
  
  # Normalize
  # M1: divide by sequencing depth
  M1_col <- paste0(expr_col, "_M1")
  gene_df[M1_col] <- gene_df[expr_col]/gene_df['SeqDepth']
  pseudo_count <- min(gene_df[gene_df[, M1_col] > 0, M1_col]) / 2
  gene_df[M1_col] <- log10(gene_df[, M1_col] + pseudo_count)
  
  #M2: divide by bug abundance
  M2_col <- paste0(expr_col, "_M2")
  gene_df[M2_col] <- gene_df[expr_col]/gene_df['bug_abundance']
  pseudo_count <- min(gene_df[gene_df[, M2_col] > 0, M2_col], na.rm=TRUE) / 2
  gene_df[M2_col] <- log10(gene_df[, M2_col] + pseudo_count)
  
  
  M4_col <- paste0(expr_col, "_M4")

  gene_df_metag['tss'] <- gene_df_metag[expr_col]/gene_df_metag['SeqDepth']
  pseudo_count <- min(gene_df_metag[gene_df_metag[, 'tss'] > 0, 'tss'], na.rm=TRUE) / 2
  gene_df_metag['tss'] <- gene_df_metag['tss'] + pseudo_count
  gene_df_metag['metag_tss'] = gene_df_metag['metag']/gene_df_metag['metag_seqdepth']
  pseudo_count <- min(gene_df_metag[gene_df_metag[, 'metag_tss'] > 0, 'metag_tss'], na.rm=TRUE) / 2
  gene_df_metag['metag_tss'] <-  gene_df_metag['metag_tss'] + pseudo_count
  gene_df_metag[M4_col] <- gene_df_metag['tss']/gene_df_metag['metag_tss']
  gene_df_metag[M4_col] <- log10(gene_df_metag[, M4_col])
  
  
  # Run the model
  expr_cols <- c(M1_col, M2_col)
  if (nrow(gene_df) > 0 && sum(gene_df[[expr_col]]) > 0) {
    models <- list()
    for (model_col in expr_cols){
      formula <- paste0(model_col, " ~ ", paste(fixed_effects, collapse = " + "))
      print(formula)
      md <- glm(formula = formula, data = gene_df, family = gaussian()) %>%
      suppressWarnings(suppressMessages(update(., .~.)))
      models[[as.symbol(model_col)]] <- summary_function(md)
    }
    
    formula <- paste0(M4_col, " ~ ", "Phenotype")
    print(formula)
    if (nrow(gene_df_metag) > 0&& sum(gene_df_metag[[expr_col]]) > 0){
    md <- glm(formula = formula, data = gene_df_metag, family = gaussian()) %>%
      suppressWarnings(suppressMessages(update(., .~.)))
    models[[as.symbol(M4_col)]] <- summary_function(md)}
    return (models)
  }
  
}


# To make the code run faster, we're only going to test a subset of genes
# Feel free to run this on the full dataframe
test_df <- metat %>% select(all_of(sample_cols)) %>% sample_n(1000, replace=FALSE)

outputs <- list()
for (i in 1:nrow(test_df)) {
  #Subset to gene dataframe and transform
  gene_df <- test_df[i, ] %>% t() %>% as.data.frame()
  expr_col <- colnames(gene_df)
  gene_metag <- metag[expr_col, sample_cols] %>% t() %>% as.data.frame() 
  colnames(gene_metag) <- c("metag")
  
  current_bug <- str_split_fixed(expr_col, "_", 2)[, 1]
  # Get information about the bug abundance for that gene across samples
  gene_bug_abundance <- bug_metat %>% filter(bug == current_bug) %>% column_to_rownames('sample_id')
  gene_df <- cbind(metat_pheno, gene_df, gene_metag, gene_bug_abundance, metag_seqdepth)
  
  outputs <- append(outputs, model_gene(gene_df, expr_col, 'semistrict', fixed_effects))
}

outputs <- do.call(rbind, outputs)
outputs$model <- str_split_fixed(rownames(outputs), "_", 3)[, 3]
outputs$gene_name <- str_split_fixed(rownames(outputs), "_M", 2)[, 1]
res <- outputs%>% select(gene_name, "Pr(>|t|)", model) %>% pivot_wider(names_from = model, values_from = "Pr(>|t|)")

# Multiple testing correction 
res$M1_qvals <- as.numeric(p.adjust(res[["M1"]], method = "BH"))
res$M2_qvals <- as.numeric(p.adjust(res[["M2"]], method = "BH"))
res$M4_qvals <- as.numeric(p.adjust(res[["M4"]], method = "BH"))


# Evaluate results

res <- res %>% 
        full_join(ground_truth, by='gene_name')

# To calculate FPR and TPR we assign all 'real' differentially expressed genes 1 and the rest 0
# For model predictions, we assign 1 to all genes with q-val < 0.01, 0 otherwise
res <- res %>% mutate(GT = abs(coalesce(as.numeric(GT), 0))) %>% 
      mutate(M1_predicted = as.numeric(M1_qvals < 0.05), M2_predicted = as.numeric(M2_qvals < 0.05),
             M4_predicted = as.numeric(M4_qvals < 0.05))


# 

calculate_rates <- function(predicted, reference){
  u <- union(predicted, reference)
  xtab <- table(factor(predicted, u), factor(reference, u))
  cm <- caret::confusionMatrix(xtab, positive = "1")
  
}


reference <- res$GT

metrics <- list()
cm_M1 <- calculate_rates(res$M1_predicted, reference)
metrics$M1 <- c(cm_M1$byClass[["Sensitivity"]], 1-cm_M1$byClass[['Specificity']])
names(metrics$M1) <- c('TPR', 'FRP')

cm_M2 <- calculate_rates(res$M2_predicted, reference)
metrics$M2 <- c(cm_M2$byClass[["Sensitivity"]], 1-cm_M2$byClass[['Specificity']])
names(metrics$M2) <- c('TPR', 'FRP')

cm_M4 <- calculate_rates(res$M4_predicted, reference)
metrics$M4 <- c(cm_M4$byClass[["Sensitivity"]], 1-cm_M4$byClass[['Specificity']])
names(metrics$M4) <- c('TPR', 'FRP')


final_df <- data.frame(metrics) %>% rownames_to_column('metric') %>%
  pivot_longer(cols=-c(metric), names_to='Model', values_to = 'Proportion')

# TPR (Sensetivity)

# FPR (1 - Specificity)


# Basic barplot

p<-ggplot(data=final_df, aes(x=Model, y=Proportion)) +
  geom_bar(stat="identity")+ facet_wrap('metric')
p


# The end
