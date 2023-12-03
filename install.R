
install.packages(c("tidyverse", "stats", "stats4", "caret", "ggplot2"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.16")
BiocManager::install("DESeq2")