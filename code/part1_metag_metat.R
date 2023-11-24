
library(data.table)
library(tidyverse)
library(R.utils)
library(patchwork)
library(GGally)


profile <- fread("datasets/part1/MGS_K03040_K03043_tara.tsv.gz",sep="\t", 
                 header=T,data.table = F,tmpdir=".")
sample_info <- fread("datasets/part1/sample_info.csv",sep=",",
                header=T, data.table = F, tmpdir=".")


# Assign -1 fraction median gene length

if (length(which(profile$length < 0)) > 0){
  med_length = median(profile$length[which(profile$length > 0)])
  profile$length[which(profile$length < 0)]<- med_length
}



profile_lengthnorm <- profile %>% 
                      mutate_if(is.numeric, function(x){x / profile$length})



profile_lengthnorm<-profile[,1:4]
for (i in 5:ncol(profile)){
  cat("Normalizing by gene length: sample",colnames(profile)[i],"\n")
  tmp<-profile[,i]/profile$length %>%
    as.data.frame()
  colnames(tmp)<-colnames(profile)[i]
  profile_lengthnorm<-profile_lengthnorm %>%
    bind_cols(tmp)
}

mgs <- c("K06942", "K01889", "K01887", "K01875", "K01883", 
         "K01869", "K01873", "K01409", "K03106", "K03110")



profile_lengthnorm_mgnorm <- profile_lengthnorm[, 1:4]

for (i in 5:ncol(profile_lengthnorm)){
  cat("Normalizing by 10 MGs: sample",colnames(profile_lengthnorm)[i],"\n")
  mg_median<-profile_lengthnorm %>%
    select(KO,abundance=all_of(colnames(profile_lengthnorm)[i])) %>%
    filter(KO %in% mgs) %>%
    group_by(KO) %>% summarise(abundance=sum(abundance)) %>%
    ungroup() %>% summarise(mg_median=median(abundance)) %>%
    pull()
  tmp<-profile_lengthnorm[,i]/mg_median
  tmp<-tmp %>% as.data.frame()
  colnames(tmp)<-colnames(profile_lengthnorm)[i]
  profile_lengthnorm_mgnorm<-profile_lengthnorm_mgnorm %>%
    bind_cols(tmp)
}



rp_ab <- profile %>%
  select(-reference, -length, -Description) %>%
  filter(KO %in% c("K03040", "K03043")) %>%
  pivot_longer(-KO, names_to = "sample", values_to = "inserts") %>% 
  filter(grepl('METAG', sample)) %>% 
  group_by(KO, sample) %>% summarize(inserts = sum(inserts)) %>%
  pivot_wider(names_from = "KO", values_from = "inserts")
  
  
  
rp_ab_lengthnorm <- profile_lengthnorm %>%
  select(-reference, -length, -Description) %>%
  filter(KO %in% c("K03040", "K03043")) %>%
  pivot_longer(-KO, names_to = "sample", values_to = "inserts_lengthnorm") %>%
  filter(grepl('METAG', sample)) %>% 
  group_by(KO, sample) %>% summarise(inserts_lengthnorm = sum(inserts_lengthnorm)) %>%
  pivot_wider(names_from = "KO", values_from = "inserts_lengthnorm")

g1 <- ggplot(data = rp_ab, aes(x = K03040, y = K03043)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = (1342/329)) +
  geom_abline(linetype = 2) +
  xlim(range(rp_ab$K03040, rp_ab$K03043)) +
  ylim(range(rp_ab$K03040, rp_ab$K03043)) +
  xlab("K03040: rpoA\n(DNA-directed RNA polymerase subunit alpha)") +
  ylab("K03043: rpoB\n(DNA-directed RNA polymerase subunit beta)") +
  labs(title = "Insert counts", subtitle = "Slope ~ 4 which corresponds to the ratio of gene lengths\n(K03040: 1,342 aa; K03043: 329 aa in E. coli K-12)") +
  coord_fixed() +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 7))


g2 <- ggplot(data = rp_ab_lengthnorm, aes(x = K03040, y = K03043)) +
  geom_point(alpha = 0.5) +
  geom_abline(linetype = 2) +
  xlim(range(rp_ab_lengthnorm$K03040, rp_ab_lengthnorm$K03043)) +
  ylim(range(rp_ab_lengthnorm$K03040, rp_ab_lengthnorm$K03043)) +
  xlab("K03040: rpoA\n(DNA-directed RNA polymerase subunit alpha)") +
  ylab("K03043: rpoB\n(DNA-directed RNA polymerase subunit beta)") +
  labs(title = "Gene-length normalized insert counts", subtitle = "Slope ~ 1 once insert counts are corrected for differences\nin gene lengths") +
  coord_fixed() +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 7))

g1|g2


mgs_ab_lengthnorm <- profile_lengthnorm %>%
  select(-reference, -Description, -length) %>%
  filter(KO %in% mgs) %>%
  pivot_longer(-KO, names_to = "sample", values_to = "inserts_lengthnorm") %>%
  group_by(KO, sample) %>% summarise(inserts_lengthnorm = sum(inserts_lengthnorm)) %>%
  ungroup() %>% group_by(sample) %>% summarise(median_mgs = median(inserts_lengthnorm)) %>%
  inner_join(sample_info, by = c("sample" = "sample_metag"))


g3 <- ggplot(data = mgs_ab_lengthnorm, aes(x = sample_metag_nreads, y = median_mgs)) +
  geom_point(alpha = 0.5) +
  #geom_smooth(method = "lm") +
  #scale_x_log10() +
  #scale_y_log10() +
  xlab("Sequencing depth (number of reads)") +
  ylab("Median abundance of the 10 universal\nand single-copy marker genes") +
  theme_bw() +
  theme(legend.title = element_blank())


g3


mgs_ab_lengthnorm <- profile_lengthnorm %>%
  select(-reference, -Description, -length) %>%
  filter(KO %in% mgs) %>%
  pivot_longer(-KO, names_to = "sample", values_to = "inserts_lengthnorm") %>%
  group_by(KO, sample) %>% summarise(inserts_lengthnorm = sum(inserts_lengthnorm)) %>%
  inner_join(sample_info, by = c("sample" = "sample_metag")) %>%
  select(KO, sample, inserts_lengthnorm) %>%
  pivot_wider(names_from = "KO", values_from = "inserts_lengthnorm")

g4 <- ggpairs(data = mgs_ab_lengthnorm %>% column_to_rownames("sample")) +
  scale_x_log10() +
  scale_y_log10()

g4

# Combining Metagenomic and metatranscriptomic data

gc_profile <- fread("datasets/part1/K03704_tara_lengthnorm_percell.tsv.gz", sep = "\t", header = T, data.table = F, tmpdir = ".")
#sample_info <- fread("sample_info_tara.tsv", sep = "\t", header = T, data.table = F, tmpdir = ".")


ko_profile <- gc_profile %>%
  group_by(KO) %>% summarise(across(starts_with("TARA"), sum)) %>%
  as.data.frame()


# Compute the gene abundance, transcript abundance and expression for the pairs of metaG-metaT samples
# The expression is just the ratio of transcript_abundance to gene_abundance
tmp_sample_info <- sample_info %>%
  select(sample_metag, sample_metat) %>%
  mutate(sample_pair = paste(sample_metag, sample_metat, sep = "-"))
tmp_metag <- ko_profile %>%
  select(KO, all_of(tmp_sample_info$sample_metag)) %>%
  pivot_longer(-KO, names_to = "sample_metag", values_to = "gene_abundance")
tmp_metat <- ko_profile %>%
  select(KO, all_of(tmp_sample_info$sample_metat)) %>%
  pivot_longer(-KO, names_to = "sample_metat", values_to = "transcript_abundance")
final_profile <- tmp_sample_info %>%
  left_join(tmp_metag, by = "sample_metag") %>%
  left_join(tmp_metat, by = c("KO", "sample_metat")) %>%
  mutate(expression = transcript_abundance/gene_abundance)

toplot <- final_profile %>%
  filter(KO == "K03704") %>%
  left_join(sample_info, by = c("sample_metag","sample_metat"))



g_metat <- ggplot(data = toplot, aes(y = transcript_abundance, x = Temperature, color =  polar)) +
  geom_point() +
  geom_smooth(method = "gam", se = T, formula = y ~ s(x, bs = "cs", k=5)) +
  scale_color_manual(values = c("darkgreen", "darkblue"))+
  ylab("Transcript abundance") +
  theme_bw() +
  theme(legend.position = "none")


g_metag <- ggplot(data = toplot, aes(y = gene_abundance, x = Temperature, color = polar)) +
  geom_point() +
  geom_smooth(method = "gam", se = T, formula = y ~ s(x, bs = "cs", k=5)) +
  #scale_y_log10() +
  #coord_flip() +
  scale_color_manual(values = c("darkgreen", "darkblue")) +
  ylab("Gene abundance") +
  theme_bw() +
  theme(legend.position = "none")


g_exp <- ggplot(data = toplot, aes(y = expression, x = Temperature, color = polar)) +
  geom_point() +
  geom_smooth(method = "gam", se = T, formula = y ~ s(x, bs = "cs", k=5)) +
  #scale_y_log10() +
  #coord_flip() +
  scale_color_manual(values = c("darkgreen", "darkblue")) +
  ylab("Gene expression") +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank())
g <- g_metag | g_exp | g_metat
g


