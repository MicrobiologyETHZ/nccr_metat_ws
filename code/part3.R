# Load the libraries
library(tidyverse)
library(stringr)
library(stats)
library(stats4)
library(caret)
library(ggplot2)

# Read data

# For models M1 through M3 we are only using transcriptomic data

metat <- read.table("data/part3/true-exp.mtx_abunds.tsv", header = TRUE,
                 sep = "\t", row.names=1)

metag <- read.table("data/part3/true-exp.mgx_abunds.tsv", header = TRUE,
                    sep="\t", row.names=1)

# Ground truth data

ground_truth <- read.table("data/part3/true-exp.mtx_spiked.tsv", header = FALSE,
                           sep = "\t",  col.names = c('gene_name', 'GT'))


# todo gz the file

# For models M4 through M6, we are using transcriptomic and genomic data

# todo load dataset

# metat_metag <- read.table(...)


# Preprocess the data (might do this before the workshop)
# Create metatdata table (i.e. phenotype for each sample)
metat_pheno <- metat[1:2, ] %>%
  t() %>%
  data.frame() 
metat <- metat[-c(1, 2), ]

metag_seqdepth <- metag[2,] %>% t()
colnames(metag_seqdepth) <- c('metag_seqdepth')
metag <- metag[-c(1, 2), ]
# Exclude metadata from the table


# Assign each 'transcript' to a 'taxonomy'
metat$bug <- str_split_fixed(rownames(df), "_", 2)[, 1]

# Create table with 'taxonomic' abundances 
bug_metat <- metat %>%
  group_by(bug) %>%
  summarize(across(all_of(sample_cols), sum)) %>% 
  pivot_longer(cols = -bug, names_to = "sample_id", values_to = "bug_abundance") 

sample_cols <- grep("SAMPLE", names(df), value = TRUE)


# Go through example of one gene 

example_gene = "BUG0013_GROUP001489"

example_gene = "BUG0077_GROUP001555"
gene_df <- metat[example_gene, sample_cols] %>% t() %>% as.data.frame()
gene_metag <- metag[example_gene, sample_cols] %>% t() %>% as.data.frame() 
colnames(gene_metag) <- c("metag")
gene_df[[example_gene]] %>% hist()

# Add bug abundance and phenotype data
current_bug <- str_split_fixed(example_gene, "_", 2)[, 1]
gene_bug_abundance <- bug_metat %>% filter(bug == current_bug) %>% column_to_rownames('sample_id')
gene_df <- cbind(metat_pheno, gene_df, gene_metag, gene_bug_abundance, metag_seqdepth)

# Filtering

## Semi-strict: Remove sample if both gene abundance and bug abundance are 0
gene_df_semi <- gene_df %>% filter((!!as.symbol(example_gene)) > 0 | bug_abundance > 0) # 63 smaples

## Strict: Remove sample if either gene abundance and bug abundance are 0
gene_df_strict <- gene_df %>% filter((!!as.symbol(example_gene)) > 0 & bug_abundance > 0) # 47 samples
gene_df_strict_metag<- gene_df %>% filter((!!as.symbol(example_gene)) > 0 & metag > 0)

# Normalization /Models

##  Model 1: Total sum scaling
M1_col <- paste0(example_gene, "_M1")
gene_df_strict['tss'] <- gene_df_strict[example_gene]/gene_df_strict['SeqDepth']
# Add pseudo count and log transform
pseudo_count <- min(gene_df_strict[gene_df_strict[, 'tss'] > 0, 'tss']) / 2
gene_df_strict[M1_col] <- log10(gene_df_strict[, 'tss'] + pseudo_count)


## Model 2: Within taxon normalisation

M2_col <- paste0(example_gene, "_M2")
gene_df_strict[M2_col] <- gene_df_strict[example_gene]/gene_df_strict['bug_abundance']
pseudo_count <- min(gene_df_strict[gene_df_strict[, M2_col] > 0, M2_col], na.rm=TRUE) / 2
gene_df_strict[M2_col] <- log10(gene_df_strict[, M2_col] + pseudo_count)


## Model 4: RNA/DNA ration
M4_col <- paste0(example_gene, "_M4")
gene_df_strict_metag['tss'] <- gene_df_strict_metag[example_gene]/gene_df_strict_metag['SeqDepth']
gene_df_strict_metag['metag_tss'] = gene_df_strict_metag['metag']/gene_df_strict_metag['metag_seqdepth']
gene_df_strict_metag[M4_col] <- gene_df_strict_metag['tss']/gene_df_strict_metag['metag_tss']
pseudo_count <- min(gene_df_strict_metag[gene_df_strict_metag[, M4_col] > 0, M4_col], na.rm=TRUE) / 2
gene_df_strict_metag[M4_col] <- log10(gene_df_strict_metag[, M4_col] + pseudo_count)




gene_df$BUG0034_GROUP001402_M1 %>% hist()
gene_df$BUG0034_GROUP001402_M2 %>% hist()

# Run linear model

fixed_effects <- c("Phenotype")

summary_function <- function(fit) {
  lm_summary <- summary(fit)$coefficients
  para <- as.data.frame(lm_summary)[-1, -3]
  para$name <- rownames(lm_summary)[-1]
  return(para)
}


model_cols <- c(M1_col, M2_col)
models <- list()
for (model_col in model_cols){
    formula <- paste0(model_col, " ~ ", paste(fixed_effects, collapse = " + "))
    print(formula)
    md <- glm(formula = formula, data = gene_df_strict, family = gaussian()) %>%
      suppressWarnings(suppressMessages(update(., .~.)))
    models[[as.symbol(model_col)]] <- summary_function(md)
}

formula <- paste0(M4_col, "~", "Phenotype")
print(formula)
md <- glm(formula = formula, data = gene_df_strict_metag, family = gaussian()) %>%
  suppressWarnings(suppressMessages(update(., .~.)))
models[[as.symbol(M4_col)]] <- summary_function(md)
summary <- do.call(rbind, models)



# Putting this together into a function

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


res$M1_qvals <- as.numeric(p.adjust(res[["M1"]], method = "BH"))
res$M2_qvals <- as.numeric(p.adjust(res[["M2"]], method = "BH"))
res$M4_qvals <- as.numeric(p.adjust(res[["M4"]], method = "BH"))


# Evaluate results

res <- res %>% 
        full_join(ground_truth, by='gene_name')

res <- res %>% mutate(GT = abs(coalesce(as.numeric(GT), 0))) %>% 
      mutate(M1_predicted = as.numeric(M1_qvals < 0.01), M2_predicted = as.numeric(M2_qvals < 0.01),
             M4_predicted = as.numeric(M4_qvals < 0.01))


# 

predicted_M1 <- res$M1_predicted
predicted_M2 <- res$M2_predicted
predicted_M4 <- res$M4_predicted
reference <- res$GT

metrics <- list()


u <- union(predicted_M1, reference)
xtab <- table(factor(predicted_M1, u), factor(reference, u))
cm <- caret::confusionMatrix(xtab, positive = "1")


metrics$M1 <- c(cm$byClass[["Sensitivity"]], 1-cm$byClass[['Specificity']])
names(metrics$M1) <- c('TPR', 'FRP')


u <- union(predicted_M2, reference)
xtab <- table(factor(predicted_M2, u), factor(reference, u))
cm <- caret::confusionMatrix(xtab, positive = "1")
metrics$M2 <- c(cm$byClass[["Sensitivity"]], 1-cm$byClass[['Specificity']])
names(metrics$M2) <- c('TPR', 'FRP')


u <- union(predicted_M4, reference)
xtab <- table(factor(predicted_M4, u), factor(reference, u))
cm <- caret::confusionMatrix(xtab, positive = "1")
metrics$M4 <- c(cm$byClass[["Sensitivity"]], 1-cm$byClass[['Specificity']])
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
#............................
      t() %>%
      as.data.frame()%>%
      rownames_to_column(var = "sample_id") %>% 
      pivot_longer(cols = -c(sample_id, Phenotype, SeqDepth), names_to = "gene_name", values_to = "raw_cnt")


write.csv(data, "data/part3/true-exp.mtx_abunds.pivot.tsv")


process_sample <-function(df, names){
    print('starting')
    names(df) <- names
    df <- df %>%
      t() %>%
      as.data.frame()%>%
      rownames_to_column(var = "sample_id") %>% 
    pivot_longer(cols = -c(sample_id, Phenotype, SeqDepth), names_to = "gene_name", values_to = "raw_cnt")
    df$bug <- str_split_fixed(df$gene_name, "_", 2)[,1]
  # Calculate 'tss' column
    df$tss <-  df$raw_cnt / df$SeqDepth * 1e6

  # Calculate 'bug_cnt' and merge it into the original dataframe
    bug_cnt_df <- df %>%
      group_by(sample_id, bug) %>%
      summarise(bug_cnt = sum(raw_cnt)) 
    df <- df %>%
      left_join(bug_cnt_df, by = c("sample_id", "bug"))
    # Calculate 'within_bug' and 'bug_perc'
    df <- df %>%
      mutate(within_bug = raw_cnt / bug_cnt * 1e6,
             bug_perc = bug_cnt / SeqDepth, sample_id='x')

    return (df)
}


res <- sapply(data, process_sample, names=rownames(data)) %>% as.data.frame() 

# Reshape the data
df <- df %>%
  pivot_longer(cols = -c(sample_id, Phenotype, SeqDepth), names_to = "gene_name", values_to = "raw_cnt")
library(stringr)
df$bug <- str_split_fixed(df$gene_name, "_", 2)[,1]
# Calculate 'tss' column
df$tss <-  df$raw_cnt / df$SeqDepth * 1e6
duckdb_register(con, "df", df)
# Calculate 'bug_cnt' and merge it into the original dataframe
bug_cnt_df <- tbl(con, 'df') %>%
  group_by(sample_id, bug) %>%
  summarise(bug_cnt = sum(raw_cnt)) %>%
   collect()

dbWriteTable(con, "bug_cnt_df", bug_cnt_df)


# Calculate 'bug_cnt' and merge it into the original dataframe
bug_cnt_df <- df %>%
  group_by(sample_id, bug) %>%
  summarise(bug_cnt = sum(raw_cnt)) 


df <- df %>%
  left_join(bug_cnt_df, by = c("sample_id", "bug"))

# Calculate 'within_bug' and 'bug_perc'
df <- df %>%
  mutate(within_bug = raw_cnt / bug_cnt * 1e6,
         bug_perc = bug_cnt / SeqDepth)

# Create 'semi_df' and 'strict_df'
semi_df <- df %>%
  filter(raw_cnt > 0 | bug > 0)

strict_df <- df %>%
  filter(raw_cnt > 0 & bug > 0)





#Filter

gene_df <- gene_df %>% filter((!!as.symbol(expr_col)) > 0 & bug_abundance > 0)
pseudo_count <- min(gene_df[gene_df[, expr_col] > 0, expr_col]) / 2
gene_df[paste0(expr_col, "_log")] <- log10(gene_df[, expr_col] + pseudo_count)
if (!is.null(gene_df) && any(fixed_effects %in% colnames(gene_df))) {
  formula <- paste0(expr_col, "_log ~ ", paste(fixed_effects, collapse = " + "))
  print(formula)
  md <- glm(formula = formula, data = gene_df, family = gaussian()) %>%
    suppressWarnings(suppressMessages(update(., .~.)))
  outputs[[expr_col]] <-  summary_function(md)
}

}
outputs <- do.call(rbind, outputs)

#outputs$qval <- as.numeric(p.adjust(outputs[["Pr(>|t|)"]], method = "BH"))
