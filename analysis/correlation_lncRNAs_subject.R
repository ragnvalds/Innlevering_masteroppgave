######Correlation lnCRNAs and subjects######


source("./R/lib_fun.R")

  

#####RNAseq data#####


seq_counts <- readRDS("./derivedData/lncRNA_counts_all_info.RDS")


######load qpcr data#####

##load data
qpcr <- readRDS("./derivedData/qpcr_results.RDS")

qpcr_data <- read_excel("./data/qpcr_data.xlsx") %>% 
  mutate(hgnc = target)



# Ensemble IDs from file
linc_ensembl <- read_excel("./data/name_LINCs.xlsx") %>%
  dplyr::select(lincs, gene) %>%
  filter(!is.na(gene)) %>%
  mutate(Gene = lincs, 
         ensembl = gene, 
         gene = Gene) %>%
  dplyr::select(gene, ensembl) %>%
  data.frame()

 

###inner join linc_ensemble and qpcr_results

linc_ensembl <- linc_ensembl %>% 
  mutate(gene = factor(gene, levels = c("GAS5_F1R1",
                                        "LincAKO17368_F2R2",
                                        "LincMD1_F3R3", 
                                        "LNC1405_F3R3",
                                        "LNC310con1_F2R2",
                                        "LNC310con2_F5R5",
                                        "Parrot_F1R1",
                                        "PVT1_F1R1",                                          
                                        "RMRP_F1R1"),
                            labels = c("GAS5 F1R1",
                                       "LincAKO17368 F2R2",
                                       "LincMD1 F3R3", 
                                       "LNC1405 F3R3",
                                       "LNC310con1 F2R2",
                                       "LNC310con2 F5R5",
                                       "Parrot F1R1",
                                       "PVT1 F1R1",
                                       "RMRP F1R1"))) %>% 
  mutate(hgnc = gene)

  
         
        
         


###join qpcr data and name LINCs to create a file with gene and ensemble ID

qpcr_all_info <- qpcr_data %>% 
  inner_join(linc_ensembl)


#####innerjoin qpcr and seq data   #HOW Do THIS WORK?

#qpcr_seq_matrix <- seq_counts %>%
  #cbind(qpcr_all_info)

####Make correlation plot with seq data and qpcr data.

  seq_counts %>% 
  #mutate(type = "rnaseq") %>% 
  group_by(time, sets) %>%
  ggplot(aes(time, count), group = paste(subject, sets), color = as.factor(subject))+ 
  geom_point()+
    geom_line()+
    facet_grid(gene ~ ., scales = "free" )


  #mutate(type = "rnaseq") %>% 
  group_by(time, sets) %>%
  ggplot(aes(time, count)) + #group = paste(subject, sets), color = as.factor(subject))) 
    geom_point() + 
  geom_line()+
  #geom_point(position = position_dodge(width = 0.4)) +
  #geom_line(position = position_dodge(width = 0.4)) +
  labs(caption = "Figure 10: ")+
  xlab("Timepoint") +
  ylab("count")

  



qpcr_all_info %>% 
  mutate(type = "qpcr") %>%
  rbind(seq) %>% 
  pivot_wider(names_from = type, values_from = count, values_fn = list(count = mean)) %>%
    print()
  
  
  



    
    
    
    
    
  
  


   #######maybe use some of this?


### with timepoint and sets at rest, model:tisse_offset_lib_size_normalized

full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")


full_seq_rest <- full_rest %>%
  filter(gene %in% linc_ensembl$ensembl, 
         model == "tissue_offset_lib_size_normalized") %>%
  print()
  #dplyr::select(gene, estimate) %>%
  mutate(type = "rnaseq") %>%
  mutate(time = factor(time, levels = c("timew0", "timew2pre", "timew2post", "timew12"))) %>%
  ggplot(aes(timepoint, emmean, color = sets, group = paste(lincs, sets))) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  facet_grid(lincs ~ ., scales = "free" ) +
  geom_pointrange(data = qpcr_est, aes(timepoint, emmean,
                                       ymin = lower.CL,
                                       ymax = upper.CL), 
                  position = position_dodge(width = 0.5), 
                  shape = 21)
  
  
  print()

qpcr %>%
  filter(coef == "timepointw12") %>%
  dplyr::select(gene, estimate = Value) %>%
  mutate(type = "qpcr") %>%
  rbind(full_seq_rest) %>%
  pivot_wider(names_from = type, values_from = estimate, values_fn = list(estimate = mean)) %>%
  print()
  
  

mutate(timepoint = factor(timepoint, levels = c("w0", "w2pre", "w2post", "w12"))) %>%
  ggplot(aes(timepoint, emmean, color = sets, group = paste(lincs, sets))) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  facet_grid(lincs ~ ., scales = "free" ) +
  geom_pointrange(data = qpcr_est, aes(timepoint, emmean,
                                       ymin = lower.CL,
                                       ymax = upper.CL), 
                  position = position_dodge(width = 0.5), 
                  shape = 21)







  
  
  
  filter(sets == "single") %>%
  group_by(subject, timepoint, sets) %>%
  summarise(rna = mean(rna.weight, na.rm = TRUE)) %>%
  ggplot(aes(timepoint, rna, group = paste(subject, sets), color = as.factor(subject))) + geom_point() + 
  geom_line()+
  geom_point(position = position_dodge(width = 0.4)) +
  geom_line(position = position_dodge(width = 0.4)) +
  labs(caption = "Figure 9: RNA.per.mg in RNA seq data, single sets")+
  xlab("Timepoint") +
  ylab("RNA.per.mg")



# Remove objects from environment
rm(list = ls())