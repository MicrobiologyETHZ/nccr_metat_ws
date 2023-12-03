library(DESeq2)
library(tidyverse)
source("code/deseq_fx.R")
library(data.table)


date <- Sys.Date()

# Set hypothesis test params
lfct <- 0.0
alpha <- 0.01


# Testing synthetic data
metat <- read.table("datasets/part3/true-exp/true-exp.mtx_abunds.tsv.gz", header = TRUE,
                    sep = "\t", row.names=1)

ground_truth <- read.table("datasets/part3/true-exp/true-exp.mtx_spiked.tsv.gz", 
                           header = FALSE,sep = "\t", 
                           col.names = c('gene_name', 'GT'))


metat_pheno <- metat[1:2, ] %>%
  t() %>%
  data.frame() 

metat_pheno <- metat_pheno %>% rownames_to_column('sample_id')
metat_pheno$Phenotype <- as.factor(metat_pheno$Phenotype)

# Have to add 1, DESeq can't handle data with so many 0
metat <- metat[-c(1, 2), ] +1


# Example of global normalisation 

s5 <- c("SAMPLE0001", "SAMPLE0002",
        "SAMPLE0003", "SAMPLE0004", "SAMPLE0005")

cnts <- metat[0:20, s5]
sd <- metat_pheno %>% filter(sample_id %in% s5)
rownames(sd) <- sd$sample_id
ex_dds <- DESeqDataSetFromMatrix(
  countData = cnts,
  colData = sd,
  design = ~Phenotype
)
ex_dds <- estimateSizeFactors(ex_dds)
norm_cnts <- counts(ex_dds, normalized = TRUE)
cnts <- cnts %>% rownames_to_column('gene_name') %>%  pivot_longer(cols=-gene_name, names_to = 'sample_id', values_to='Before DESeq2')
norm_cnts <- norm_cnts %>%as.data.frame() %>%  rownames_to_column('gene_name') %>%  pivot_longer(cols=-gene_name, names_to = 'sample_id', values_to='After DESeq2')
test_graph <- cnts %>% full_join(norm_cnts, by=c('gene_name', 'sample_id')) %>% pivot_longer(cols=-c(gene_name, sample_id), names_to = 'method', values_to = 'Counts')


test_graph$method = factor(test_graph$method, levels=c('Before DESeq2','After DESeq2'))
ggplot(test_graph, aes(fill=gene_name, y=Counts, x=sample_id)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~method)+
  theme_bw()

test_graph2 <- test_graph %>% filter(gene_name == 'BUG0001_GROUP000001')
ggplot(test_graph2, aes(fill=gene_name, y=Counts, x=sample_id)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~method)+
  theme_bw()


# Now running within taxon normalisation
metat$bug <- str_split_fixed(rownames(metat), "_", 2)[, 1]
conditions <- "Phenotype,"

taxa <- c("BUG0001", 'BUG0002')
taxon_col <- 'bug'
gene_col <- 'gene_col'
sample_id_col <- 'sample_id'
print(paste0("Normalizing within each taxon. Number of taxa: ", length(taxa)))
df_list <- list()
for (i in seq_along(taxa)){
  print(taxa[[i]])
  taxa_data <- metat %>% filter(!!sym(taxon_col) == taxa[[i]]) %>% 
    select(-!!sym(taxon_col))%>%
    rownames_to_column('gene_col')
  raw_data <- align_samples_and_counts(taxa_data, metat_pheno, conditions, 
                                       gene_col, sample_id_col)
  df_list[[i]] <- deseq_norm_mat(raw_data) #%>% rownames_to_column(gene_col)
}


test <- metat %>% filter(bug %in% taxa) %>%
  select(c("SAMPLE0001", "SAMPLE0002",
           "SAMPLE0003", "SAMPLE0004", "SAMPLE0005", 'bug')) %>%
  rownames_to_column('gene') 

test_graph <- test %>% 
  pivot_longer(cols = -c(bug, gene), names_to = "sample_id", values_to = "Raw")
sf0 <- dds$sizeFactor[0:5]%>% as.data.frame()
colnames(sf0) <- c('global')
sf0 <- sf0 %>% rownames_to_column('sample_id') 
sf1 <- df_list[[1]]$sizeFactor[0:5] %>% as.data.frame()
colnames(sf1) <- c('BUG0001')
sf2 <- df_list[[2]]$sizeFactor[0:5] %>% as.data.frame()
colnames(sf2) <- c('BUG0002')
sf <- cbind(sf1, sf2) %>% rownames_to_column('sample_id') %>%
  pivot_longer(cols = -c(sample_id), names_to='bug',values_to='sf' )
test_graph <- sf %>% full_join(test_graph, by=c('bug', 'sample_id'))
test_graph$Taxon_normalisation = test_graph$Raw/test_graph$sf
test_graph2 <- sf0 %>% full_join(test_graph, by=c('sample_id'))
test_graph2$Global_normalisation = test_graph2$Raw/test_graph2$global


test_graph2 <- test_graph2 %>% pivot_longer(cols=-c('sample_id', 'global', 'bug', 'sf', 'gene'), names_to='method', values_to='Counts')

test_graph2$method <- factor(test_graph2$method, levels=c('Raw', 'Global_normalisation', 'Taxon_normalisation'))

ggplot(test_graph2, aes(fill=bug, y=Counts, x=sample_id)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap('method')+
  theme_bw()



# Analysis using global normalisation

cnts <- metat[, metat_pheno$sample_id]
rownames(metat_pheno) <- metat_pheno$sample_id
dds <- DESeqDataSetFromMatrix(
  countData = cnts,
  colData = metat_pheno,
  design = ~Phenotype
)
colData(dds)$Phenotype = factor(colData(dds)$Phenotype, levels=unique(colData(dds)$Phenotype))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds, quiet = FALSE)
res_glo <- results(dds, contrast = c("Phenotype", "0", "1"), alpha = alpha, lfcThreshold = lfct)







# First running DESeq2 on the full count data without normalisation. 
# Count data had +1 added to it. Otherwise DESeq2 fails (due to low coverage)
deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix, TRUE, 10)

dds <- deseq_output$dds
sd <- deseq_output$sample_data

get_all_results(sd, pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)







ggplot(test_graph, aes(fill=bug, y=norm_ab, x=sample_id)) + 
  geom_bar(position="stack", stat="identity")





deseq_output <- deseq_on_metat_taxon(count_data, sample_data, conditions, "ID", "sample_id", "genome")
dds <- deseq_output$dds
get_all_results(sd, pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)


# Now running on each individual bug one by one. 
prefix <- paste0(date, "_true-exp-deseq-taxon-one-by-one-")
bugs <- count_data$genome %>% unique()
df_list <- list()
for (i in seq_along(bugs)){
    data <- count_data %>% filter(genome == bugs[[i]])
    print(bugs[[i]])
    raw_data <- align_samples_and_counts(data, sample_data, conditions, "ID", "sample_id")
    dds <- DESeqDataSetFromMatrix(countData = raw_data$count_data, colData = raw_data$sample_data, design = ~group)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    dds <- DESeq(dds, quiet = FALSE)
    df_list[[i]] <- results(dds, contrast = c("group", pattern1, pattern2), alpha = alpha, lfcThreshold = lfct) %>% 
                        as.data.frame() %>% mutate(bug=bugs[[i]]) %>%
                      rownames_to_column("ID")
}

final_result <- rbindlist(df_list, fill=TRUE)
write.csv(final_result, file.path(out_dir, paste0(prefix, pattern1, "_vs_", pattern2, "_l", lfct, "a", alpha, "_results.csv")))

# Now applying this to OLIGO and LCM data
out_dir <- file.path(root, configs$output_dir)

count_data_file <- file.path(root, configs$deseq_count_file)
sample_data_file <- file.path(root, configs$sample_data_file)

count_data <- read.csv(count_data_file)
sample_data <- read.csv(sample_data_file)
conditions <- "Mouse,Treatment"
prefix <- paste0(date, "_oligo-lcm-within-taxon-")
deseq_output <- deseq_on_metat_taxon(count_data, sample_data, conditions, "ID", "sample_id", "genome")
dds <- deseq_output$dds


pattern1 <- "LCM_D"
pattern2 <- "LCM_PBS_D1"
get_all_results(colData(dds), pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)


pattern1 <- "Oligo_PBS"
pattern2 <- "LCM_PBS_D1"
get_all_results(colData(dds), pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)
write.csv(deseq_output$norm_counts, file.path(out_dir, paste0(prefix, 'norm_cnts.csv')))




# Need to analyze OLIGO and LCM experiments on their own
oligo <- sample_data %>% filter(Mouse == 'Oligo')
oligo_counts <- count_data %>% select(c(c('ID', 'genome'), oligo$sample_id))
oligo_counts <- oligo_counts %>% filter(genome != "ASF457")
conditions <- 'Treatment'
prefix <- paste0(date, "_oligo-alone-within-taxon-")
oligo_output <- deseq_on_metat_taxon(oligo_counts, oligo, conditions, "ID", "sample_id", "genome")
dds <- oligo_output$dds
pattern1 <- "LPS"
pattern2 <- "PBS"
get_all_results(colData(dds), pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)

write.csv(oligo_output$norm_counts, file.path(out_dir, paste0(prefix, 'norm_cnts.csv')))

# Now LCM
lcm <- sample_data %>% filter(Mouse == 'LCM')
lcm_counts <- count_data %>% select(c(c('ID', 'genome'), lcm$sample_id))
lcm_counts <- lcm_counts %>% filter(genome != "ASF457")
conditions <- 'Treatment'
prefix <- paste0(date, "_lcm-alone-within-taxon-")
lcm_output <- deseq_on_metat_taxon(lcm_counts, lcm, conditions, "ID", "sample_id", "genome")
dds <- lcm_output$dds
pattern1 <- "^D"
pattern2 <- "PBS_D1"
get_all_results(colData(dds), pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)

write.csv(lcm_output$norm_counts, file.path(out_dir, paste0(prefix, 'norm_cnts.csv')))



# have to set the following env variables: export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
# Settings set to active terminal
configs <- config::get(file="03_2023_sal_transcriptomics/config.yaml")
root <- configs$root
sample_data_file <- file.path(root, configs$sample_data_file)

# OLIGO EXPERIMENT
conditions <- "Treatment"
out_dir <- file.path(root, configs$output_dir)

lfct <- 0.0
alpha <- 0.01
date <- Sys.Date()

# Analysing strains together
prefix <- paste0(date, "_oligo-metat-fc-deseq-")
count_file <- file.path(root, configs$oligo_fc_raw)
deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix, TRUE)

dds <- deseq_output$dds
sd <- deseq_output$sample_data

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)

# Analysing different strains seperately
# I48

prefix <- paste0(date, "_oligo-i48-fc-deseq-")
count_file <- file.path(root, configs$oligo_i48_fc_raw)
deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sample_data

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)

# YL27

prefix <- paste0(date, "_oligo-yl27-fc-deseq-")
count_file <- file.path(root, configs$oligo_yl27_fc_raw)

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sample_data

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)


#-----------------------------------------------------------------------


# YL32

prefix <- "08-06-23-YL32-high-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL32_high_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)



# YL32 - all

prefix <- "08-06-23-YL32-all-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL32_all_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)




# YL27

prefix <- "08-06-23-YL27-all-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL27_all_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)



# YL58

prefix <- "08-06-23-YL58-med-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL58_med_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)



# YL58

prefix <- "08-06-23-YL58-all-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL58_all_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)

