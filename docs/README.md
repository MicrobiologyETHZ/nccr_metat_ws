
Binder link: [![Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MicrobiologyETHZ/nccr_metat_ws/HEAD)

# Part 1: Analysis using matching metatranscriptomic/metagenomic data 
![metat_metag](metat_metag.png)

# Part 2: Analysis using metatranscriptomic data alone (defined communities) 
![metat_comm](/docs/images/metat_community.png)

# Part 3: Differential expression analysis with linear models

In the part of the workshop, we will reproduce some of the analysis from the paper [Statistical approaches for differential expression analysis in metatranscriptomics](https://doi.org/10.1093/bioinformatics/btab327). In this paper, the authors use different synthetic metatranscriptomic dataset to test a number of filtering, normalisation and modeling strategies to identify the best analysis method. 


## Models 

We will reproduce 3 models from the M1, M2, and M4 as they are described in the paper. 

## Filtering strategies



## Download the data and code

The easiest way to follow along is to clone this repository and navigate to it in RStudio.

```
git clone https://github.com/MicrobiologyETHZ/nccr_metat_ws
```

All of the synthetic data generated in the paper can be found [here](). We will be working with only 1 of the datasets, you can download the data from [here](../data/part3/true-exp/). 


This dataset contains 3 files:

- `true-exp.mtx_abunds.tsv.gz` contains simulated transcriptomic counts for each gene  as well as **phenotype** information (1 or 0) across 100 samples.
- `true-exp.mgx_abunds.tsv.gz` contains metagenomic counts for each gene across 100 samples
- `true-exp.mtx_spiked.tsv.gz` contains 'ground truth' information about with genes are diffrentially expressed between the phenotypic groups

All of the analysis code can be found in the file [part3.R](../code/part3.R).



## Homework

Using what you've learned today, reproduce the remaining models (M3, M5, and M6).
