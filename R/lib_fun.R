## RNA-seq pipeline #######
# Libraries and functions #

 #BiocManager::install("clusterProfiler")

#BiocManager::install("rnaseqcomp")

#if (!requireNamespace('BiocManager', quietly = TRUE))
  #install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')

## Bioconductor packages
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(edgeR)
library(limma)
library(Glimma)
library(variancePartition)
library(qvalue)
library(preprocessCore)

library(topGO)

library(DHARMa)
library(rnaseqcomp)
library(biomaRt)
library(Rnmr1D)

# Cran packages

library(RColorBrewer)

library(tidyverse)
library(readxl)
library(knitr)
library(dplyr)
library(nlme)
library(emmeans)
library(feather)
library(glmmTMB)
library(lme4)
library(mgcv)

library(svMisc)
library(foreach)

library(doSNOW)
library(parallel)
library(doParallel)
#install.packages("huxtable")
library(huxtable)
#install.packages("magick")
#install.packages("webshot")
#webshot::install_phantomjs()
library(magick)
library(webshot)


### For plotting 
library(cowplot)


# github installations

# devtools::install_github("DarwinAwardWinner/rctutils")
# devtools::install_github("dhammarstrom/publR")
# devtools::install_github("dhammarstrom/qpcrpal")
library(qpcR)
library(qpcrpal)
library(publR)
#library(rctutils)


######### Functions #######################
# Extracts last letter of a string
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## Extract which row is quantile nnn

which.quantile <- function (x, probs = 0.5, na.rm = FALSE){
  if (! na.rm & any (is.na (x)))
    return (rep (NA_integer_, length (probs)))
  
  o <- order (x)
  n <- sum (! is.na (x))
  o <- o [seq_len (n)]
  
  nppm <- n * probs - 0.5
  j <- floor(nppm)
  h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
  j <- j + h
  
  j [j == 0] <- 1
  o[j]
}
