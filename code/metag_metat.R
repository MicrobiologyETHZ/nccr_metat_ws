
library(data.table)
library(tidyverse)
library(R.utils)

library(patchwork)
library(GGally)


profile <- fread("data/metat_metag_test_profile.csv.gz",sep=",", 
                 header=T,data.table = F,tmpdir=".")
sample_info <- fread("data/metat_metag_test_meta.csv.gz",sep=",",
                     header=T, data.table = F, tmpdir=".")

profile[1:4]%>% tail(1)

if (length(which(profile$length < 0)) > 0){
  med_length = median(profile$length[which(profile$length > 0)])
  profile$length[which(profile$length < 0)]<- med_length
}

profile[1:4]%>% tail(1)

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

mg_median <- profile_lengthnorm %>% select(-length) %>% 
            filter(KO %in% mgs) %>% 
            group_by(KO) %>% summarise(across(contains('Sample'), sum)) %>%
            summarise(across(contains('Sample'), median))


profile_lengthnorm_mgnorm<- profile_lengthnorm[,1:4]
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
  pivot_longer(-KO, names_to = "sample", values_to = "inserts") %>% head()
  filter(sample,  grepl('METAG')) %>% 
  group_by(KO, sample) %>% summarise(inserts = sum(inserts)) %>%
  pivot_wider(names_from = "KO", values_from = "inserts")
  
  
  
rp_ab_lengthnorm <- profile_lengthnorm %>%
  select(-reference, -length, -Description) %>%
  filter(KO %in% c("K03040", "K03043")) %>%
  pivot_longer(-KO, names_to = "sample", values_to = "inserts_lengthnorm") %>% head() %>%
  filter(sample,  grepl('METAG')) %>% head()
  #group_by(KO, sample) %>% summarise(inserts_lengthnorm = sum(inserts_lengthnorm)) %>%
  #pivot_wider(names_from = "KO", values_from = "inserts_lengthnorm")
