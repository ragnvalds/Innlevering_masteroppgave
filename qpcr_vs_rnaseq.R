

## Fold change analysis, comparisons between RNA-seq and qpcr

source("./R/lib_fun.R")


######load qpcr data#####

##load data
qpcr_results <- readRDS("./derivedData/qpcr_all_info.RDS")


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
  labs(title = "Correlation rested state, single sets") 




#####timepoint and sets accute

full_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel.RDS")

rnaseq_acute_full<- full_acute %>%
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
  labs(title = "Correlation acute response, single sets")




#####################multiple sets################

#####Response at rest pre-post#####

### with timepoint and sets at rest, model:tissue_offset_lib_size_normalized

full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")

rnaseq_w2pre_full_t_m<- full_rest %>%
  filter(gene %in% linc_ensembl$ensembl, 
         coef == "timew12:setsmultiple", 
         model == "tissue_offset_lib_size_normalized") %>%
  dplyr::select(gene, estimate) %>%
  mutate(type = "rnaseq") %>%
  print()

qpcr_results %>%
  filter(coef == "timepointw12",
         sets == "multiple") %>%
  dplyr::select(gene, estimate = Value) %>%
  mutate(type = "qpcr") %>%
  rbind(rnaseq_w2pre_full_t_m) %>%
  pivot_wider(names_from = type, values_from = estimate, values_fn = list(estimate = mean)) %>%
  ggplot(aes(qpcr, rnaseq)) + geom_point() + 
  xlab("qPCR") +
  ylab("RNAseq") +
  labs(title = "Correlation rested state, multiple sets") 


  

