#install.packages('arrow')
install.packages('duckdb')

library(arrow)
library(duckdb)
library(tidyverse)

library(dplyr, warn.conflicts = FALSE)



con <- dbConnect(duckdb())



df <- read.table("data/gene_profile_tara_lengthnorm_percell.tsv.gz", stringsAsFactors = FALSE, 
                 sep = "\t")

duckdb_register(con, "df", df)

df_tbl <- to_duckdb(df, con, "df")
test <- df %>%
  head() 
  #filter(KO == "K03704") %>%
  



df$tss <-  df$raw_cnt / df$SeqDepth * 1e6
duckdb_register(con, "df", df)
# Calculate 'bug_cnt' and merge it into the original dataframe
bug_cnt_df <- tbl(con, 'df') %>%
  group_by(sample_id, bug) %>%
  summarise(bug_cnt = sum(raw_cnt)) %>%
  collect()

csv_ds <- open_dataset("original_csv", 
                       format = "csv",
                       partitioning = c("year", "month"))

## this reads each csv file in the csv_ds dataset and converts it to a .parquet file
write_dataset(csv_ds, 
              "converted_parquet", 
              format = "parquet",
              partitioning = c("year", "month"))
}


data_dir <- "data/part3/"
null_ds <- open_dataset("data/part3/tru-exp",
                         schema = schema(geneid=string(),
                                         BUG=string(),
                                         sample_id=string(),
                                         mgx_abunds=int64(),
                                         mtx_abunds=int64(),
                                         Phenotype=int64()
                                         ))

con <- DBI::dbConnect(duckdb::duckdb())
null_tbl <- to_duckdb(null_ds, con, "null")

null_tbl %>%
  group_by(sample_id, BUG) %>%
  summarize(bug_mtx = sum(mtx_abunds, na.rm=TRUE), bug_mgx = sum(mgx_abunds, 
                                                                 na.rm=TRUE)) %>%
  print()


