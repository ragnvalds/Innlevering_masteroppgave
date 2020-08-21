###################Correlation figure qpcr rna_seq#########


## Fold change analysis, comparisons between RNA-seq and qpcr

source("./R/lib_fun.R")

#We’ll use the ggpubr R package for an easy ggplot2-based data visualization

#Install the latest version from GitHub as follow (recommended):
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")

#install.packages("ggpubr")
library("ggpubr")

##load data qpcr and seq

correlation <- readRDS("./derivedData/cor_seq_qpcr.RDS")



ggscatter(correlation, x = "seq", y = "qpcr", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "RNAseq", ylab = "qPCR")







