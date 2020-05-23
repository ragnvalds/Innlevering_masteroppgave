

## Fold change analysis, comparisons between RNA-seq and qpcr

source("./R/lib_fun.R")


              ######load qpcr data#####

##load data
qpcr_results <- readRDS("./derivedData/qpcr_results.RDS")


# Ensemble IDs from file
linc_ensembl <- read_excel("./data/name_LINCs.xlsx") %>%
  dplyr::select(lincs, gene) %>%
  filter(!is.na(gene)) %>%
  mutate(Gene = lincs, 
         ensembl = gene, 
         gene = Gene) %>%
  dplyr::select(gene, ensembl) %>%
  data.frame()




#####RNA seq results#####

#####Response at rest pre-post#####


### with timepoint only at rest, model:tissue_offset_lib_size_normalized

time_rest <- readRDS("./derivedData/DE/mixedmodel2_timemodel.RDS")

rnaseq_w2pre_t<- time_rest %>%
  filter(gene %in% linc_ensembl$ensembl, 
         coef == "timew12", 
         model == "tissue_offset_lib_size_normalized") %>%
  dplyr::select(gene, estimate) %>%
  mutate(type = "rnaseq") %>%
  print()

qpcr_results %>%
  filter(coef == "timepointw12") %>%
  dplyr::select(gene, estimate = Value) %>%
  mutate(type = "qpcr") %>%
  rbind(rnaseq_w2pre_t) %>%
  pivot_wider(names_from = type, values_from = estimate, values_fn = list(estimate = mean)) %>%
  ggplot(aes(qpcr, rnaseq)) + geom_point() + 
  xlab("qPCR") +
  ylab("RNAseq") +
  labs(title ="Figure 5: Timepointw12 fold change, timepoint only", subtitle = "Model: Tissue offset lib size normalized") +
  
  labs(caption = "(Based on data from qpcr and gene sequencing)")

  


### with timepoint and sets at rest, model:tissue_offset_lib_size_normalized

full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")

rnaseq_w2pre_full_t<- full_rest %>%
  filter(gene %in% linc_ensembl$ensembl, 
         coef == "timew12", 
         model == "tissue_offset_lib_size_normalized") %>%
  dplyr::select(gene, estimate) %>%
  mutate(type = "rnaseq") %>%
  print()

qpcr_results %>%
  filter(coef == "timepointw12") %>%
  dplyr::select(gene, estimate = Value) %>%
  mutate(type = "qpcr") %>%
  rbind(rnaseq_w2pre_full_t) %>%
  pivot_wider(names_from = type, values_from = estimate, values_fn = list(estimate = mean)) %>%
  ggplot(aes(qpcr, rnaseq)) + geom_point() + 
  xlab("qPCR") +
  ylab("RNAseq") +
  labs(title = "Figure 6: Timepointw12 fold change, timepoint and sets", subtitle = "Model: Tissue offset lib size normalized") +
  
  labs(caption = "(Based on data from qpcr and gene sequencing)")



#####Acute respons. only lib size normalized#####


#### only time model, timepoint acute

time_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel_timeonly.RDS")

rnaseq_acute<- time_acute %>%
  filter(gene %in% linc_ensembl$ensembl, 
         coef == "timew2post", 
         model == "lib_size_normalized") %>%
  dplyr::select(gene, estimate) %>%
  mutate(type = "rnaseq") %>%
  print()

qpcr_results %>%
  filter(coef == "timepointw2post") %>%
  dplyr::select(gene, estimate = Value) %>%
  mutate(type = "qpcr") %>%
  rbind(rnaseq_acute) %>%
  pivot_wider(names_from = type, values_from = estimate, values_fn = list(estimate = mean)) %>%
  ggplot(aes(qpcr, rnaseq)) + geom_point() + 
  xlab("qPCR") +
  ylab("RNAseq") +
  labs(title = "Figure 7: Timepoint acute fold change, timepoint only", subtitle = "Lib size normalized") +
  
  labs(caption = "(Based on data from qpcr and gene sequencing)")



#####timepoint and sets accute

full_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel.RDS")

rnaseq_acute_full<- time_acute %>%
  filter(gene %in% linc_ensembl$ensembl, 
         coef == "timew2post", 
         model == "lib_size_normalized") %>%
  dplyr::select(gene, estimate) %>%
  mutate(type = "rnaseq") %>%
  print()

qpcr_results %>%
  filter(coef == "timepointw2post") %>%
  dplyr::select(gene, estimate = Value) %>%
  mutate(type = "qpcr") %>%
  rbind(rnaseq_acute_full) %>%
  pivot_wider(names_from = type, values_from = estimate, values_fn = list(estimate = mean)) %>%
  ggplot(aes(qpcr, rnaseq)) + geom_point() + 
  xlab("qPCR") +
  ylab("RNAseq") +
  labs(title = "Figure 8: Timepointw2post acute fold change, timepoint and sets", subtitle = "Lib size normalized") +
  
  labs(caption = "(Based on data from qpcr and gene sequencing)")


