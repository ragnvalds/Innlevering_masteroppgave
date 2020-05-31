
# Source libraries

source("./R/lib_fun.R")



# Load data
# Mixed-effects negative binomial count models were fitted and saved in 
# ./R/dge_list_models.R. Results saved in RDS files (see below).
# no need to run these models again.


time_rest <- readRDS("./derivedData/DE/mixedmodel2_timemodel.RDS")
full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")

time_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel_timeonly.RDS")
full_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel.RDS")


# load list of all genes with ensemble info
all_genes <- readRDS("./derivedData/all_genes.RDS")

#load lncRNA

lncRNA <- readRDS("./derivedData/lncRNA.RDS")



##### Filter and adjust p-values for differentially expression #####



####rest w2pre #####


  

restw2 <- time_rest %>%
  filter(model %in% c("tissue_offset_lib_size_normalized"),
         coef %in% c("timew2pre")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()



restw2 %>%
  mutate(p.value = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
   
  #filter(pt == "sig" & fc == "sig") %>% 
  #mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  ggplot(aes(estimate, -log10(p.val), color = p.value)) + geom_point() + 
  facet_grid( ~ coef)+
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.text = element_text(size = 10), 
  legend.title = element_text(size = 12))+
  ###slå sammen pt og fc. fremvise alle med log 2 over 0,5 eller under -0,5. tabell 1+ lincer med lavest p verdi.
  xlab("log2FC pre-post") %>%    
    print()


  
  
  
                          ######Count lnc genes DE##############
  
DE_rest_w2 <- restw2 %>%
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig" & 
         inc_dec == "increase" &
         #fc == "sig" & 
         fc == "sig") %>%
  NROW() %>% 
  print()


   

# DE timew2pre 169
# DE increase 164
# DE decrease 5



                ##############time w12#################


restw12 <- time_rest %>%
  filter(model %in% c("tissue_offset_lib_size_normalized"),
         coef %in% c("timew12")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


restw12 %>%
  mutate(p.value = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  
  #filter(pt == "sig" & fc == "sig") %>% 
  #mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  ggplot(aes(estimate, -log10(p.val), color = p.value)) + geom_point() + 
  facet_grid( ~ coef)+
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))+###slå sammen pt og fc. fremvise alle med log 2 over 0,5 eller under -0,5. tabell 1+ lincer med lavest p verdi.
  xlab("log2FC pre-post ")%>% 
  print()


                         ###############count DE lnc##############


DE_rest_w12 <- restw12 %>%
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig" & 
           inc_dec == "decrease" &
           #fc == "sig" & 
           fc == "sig") %>%
  NROW() %>% 
  print()


# DE lnc timew12 64
# 64 increase





                  #####acute model timepoint w2post #####

  acute_w2post <- time_acute %>%
  filter(model %in% c("lib_size_normalized"),
         coef %in% c("timew2post")) %>%
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


acute_w2post%>%
  mutate(p.value = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  ggplot(aes(estimate, -log10(p.val), color = p.value)) + geom_point() +
  facet_grid( ~ coef)+
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))+
  xlab("log2FC acute") 
                     
                                  ######Count lnc genes DE#############


DE_acute_w2post <- acute_w2post %>%
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig" & 
           inc_dec == "decrease" &
           #fc == "sig" & 
           fc == "sig") %>%
  NROW() %>% 
  print()



#DE w2post 102
#increase 40
# decrease 62








 
 

 
 
 
                        ############ multiple sets###############



####rest model timepoint and sets multiple sets#####


full_rest_w2pre <- full_rest %>%
  filter(model %in% c("tissue_offset_lib_size_normalized"),
         coef %in% c("timew2pre:setsmultiple")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_rest_w2pre %>%
  mutate(p.value = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  ggplot(aes(estimate, -log10(p.val), color = p.value)) + geom_point() + 
  facet_grid( ~ coef)+
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))+
  xlab("log2FC pre-post") 


 
                    ######Count lnc genes DE##############

  DE_sets_w2pre <- full_rest_w2pre %>%
    mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
           fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
    mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
    filter(pt == "sig" ) %>%  
             #inc_dec == "increase" &
             #fc == "sig" & 
             #fc == "sig") %>%
    #NROW() %>% 
    print()


# 0 DE lncs at rest w2pre multiple sets





#######################timepoint w12################

full_rest_w12 <- full_rest %>%
  filter(model %in% c("tissue_offset_lib_size_normalized"),
         coef %in% c("timew12:setsmultiple")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_rest_w12 %>%
  mutate(p.value = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  ggplot(aes(estimate, -log10(p.val), color = p.value)) + geom_point() + 
  facet_grid( ~ coef)+
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))+
  xlab("log2FC pre-post") 



######Count lnc genes DE##############

DE_sets_w12 <- full_rest_w12 %>%
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig" ) %>%  
  #inc_dec == "increase" &
  #fc == "sig" & 
  #fc == "sig") %>%
  #NROW() %>% 
  print()

# 0 DE w12 multiple sets

  
  


#####acute model timepoint and sets multiple sets#####

full_acute_w2post <- full_acute %>%
  filter(model %in% c("lib_size_normalized"),
         coef %in% c("timew2post:setsmultiple")) %>%
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()




full_acute_w2post %>%
  mutate(p.value = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  ggplot(aes(estimate, -log10(p.val), color = p.value)) + geom_point() + 
  facet_grid( ~ coef)+
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))+
  xlab("log2FC acute") 


                   ######Count lnc genes DE##############

DE_sets_w2post <- full_acute_w2post %>%
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(fc == "sig") %>% 
           #inc_dec == "decrease" &
           #fc == "sig" & 
           #fc == "sig") %>%
  #NROW() %>% 
  print()

# 0 DE lncs acute w2post multiple sets











                        















                          ########### from qpcr###################




####rest model timepoint and sets single sets#####


full_rest_lnc_pcr <- full_rest %>%
  filter(model %in% c("lib_size_normalized", "tissue_offset_lib_size_normalized"),
         coef %in% c("timew2pre", "timew12")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only qpcr lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  filter(gene %in% goi) %>% 
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_rest_lnc_pcr %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  ggplot(aes(estimate, -log10(p.val), color = paste(pt,inc_dec))) + geom_point() + 
  facet_grid(model ~ coef)+
  xlab("log2FC") +
  labs(title = "Rested state qpcr lncs, single sets")


               ######Count lnc genes DE##############

DE_qpcr_rest_single <- full_rest_lnc_pcr %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
filter(pt == "sig",
       model == "tissue_offset_lib_size_normalized",
       coef == "timew12") %>%
  #inc_dec == "increase") %>%
  #NROW() %>% 
  print()

# 5 DE increased rest single



#####acute model timepoint and sets single sets#####

full_acute_lnc_pcr <- full_acute %>%
  filter(model %in% c("lib_size_normalized"),
         coef %in% c("timew2post")) %>%
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only qpcr lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  filter(gene %in% goi) %>% 
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_acute_lnc_pcr %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  ggplot(aes(estimate, -log10(p.val), color = paste(pt,inc_dec))) + geom_point() +
  facet_grid(model ~ coef)+
  xlab("log2FC") +
  labs(title = "Acute respons qpcr lncs, single sets")

                ######Count lnc genes DE##############

DE_qpcr_acute_single <- full_acute_lnc_pcr %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
filter(pt == "sig",
       coef == "timew2post") %>%
  dplyr::select(gene) %>% 
  #inc_dec == "increase") %>%
  #NROW() %>% 
  print()

# 3 DE decreased acute w2post single sets.



############ multiple sets###############



####rest model timepoint and sets multiple sets#####


full_rest_lnc_m_pcr <- full_rest %>%
  filter(model %in% c("lib_size_normalized", "tissue_offset_lib_size_normalized"),
         coef %in% c("timew2pre:setsmultiple", "timew12:setsmultiple")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only qpcr lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  filter(gene %in% goi) %>% 
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_rest_lnc_m_pcr %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  ggplot(aes(estimate, -log10(p.val), color = paste( pt, inc_dec))) + geom_point() + 
  facet_grid(model ~ coef)+
  xlab("log2FC") +
  labs(title = "Rested state qpcr lncs, multiple sets")


                ######Count lnc genes DE##############

DE_qpcr_rest_multiple <- full_rest_lnc_m_pcr %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig",
    model == "tissue_offset_lib_size_normalized",
       coef == "timew12:setsmultiple") %>% 
  #inc_dec == "increase") %>%
  #NROW() %>% 
  print()

# 0 DE lnc at rest multiple sets.



#####acute model timepoint and sets multiple sets#####

full_acute_lnc_m_pcr <- full_acute %>%
  filter(model %in% c("lib_size_normalized"),
         coef %in% c("timew2post:setsmultiple")) %>%
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  filter(gene %in% goi) %>% 
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_acute_lnc_m_pcr %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  ggplot(aes(estimate, -log10(p.val), color = paste(pt,inc_dec))) + geom_point() +
  facet_grid(model ~ coef)+
  xlab("log2FC") +
  labs(title = "Acute response qpcr lncs, multiple sets")




                  ######Count lnc genes DE##############

DE_qpcr_acute_multiple <- full_acute_lnc_m_pcr %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
filter(pt == "ns",
       coef == "timew2post:setsmultiple") %>% 
  #inc_dec == "increase") %>%
  #NROW() %>% 
  print()


# 0 De lncs at acute w2post multiple sets.









###########################testing volcano plot###################

library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels


input <-  #convert the rownames to a column
volc = ggplot(input, aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance
  scale_color_manual(values=c("black", "red")) + 
  ggtitle("Your title here") #e.g. 'Volcanoplot DESeq2'
volc+geom_text_repel(data=head(input, 20), aes(label=gene))






######rest model time only #####

#time_rest_lnc <- time_rest %>%
#filter(model %in% c("lib_size_normalized", "tissue_offset_lib_size_normalized"),
# coef %in% c("timew2pre", "timew12")) %>%

#dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%

# filter out only lncRNA

#filter(gene %in% lncRNA$ensembl_gene_id) %>%
# P-value adjustments
#group_by(model, coef) %>%
# mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
#print()


#time_rest_lnc %>%
# mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
#ggplot(aes(estimate, -log10(p.val), color = pt)) + geom_point() + 
# facet_grid(model ~ coef)+
# xlab("log2FC") +
# labs(caption = "Figure 1: Differential expressed lncRNAs at rest, timepoint only model")





###### Acute model time only #####

#time_acute_lnc <- time_acute %>%
#filter(model %in% c("lib_size_normalized"),
# coef %in% c("timew2post")) %>%
#dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%

# filter out only lncRNA

#filter(gene %in% lncRNA$ensembl_gene_id) %>%
# P-value adjustments
#group_by(model, coef) %>%
#mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
#print()


#time_acute_lnc %>%
#mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
#mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
#ggplot(aes(estimate, -log10(p.val), color = paste(pt,inc_dec))) + geom_point() +
#facet_grid(model ~ coef)+
#xlab("log2FC") +
#labs(caption = "Figure 3: Differential expressed lncRNAs acute respons, timepoint only model")



#######mutate(timepoint = factor(Timepoint, levels = c("PreSupp",
                                                #"PreExc",
                                               # "ThreeW", 
                                                #"PostExc")))


#mutate(timepoint = factor(timepoint, levels = c("pre",
                                                #"meso1",
                                               ## "meso2", 
                                               # "meso3"),
                          #labels = c("Pre-\ntraining",
                                     ##"Meso-\ncycle 1",
                                     #"Meso-\ncycle 2", 
                                    # "Meso-\ncycle 3")),
       #group = factor(group, levels = c("MIX", "DECR", "INCR"),
                      #labels = c("Mix", "Decrease", "Increase")), 
       #VO2max.kg = VO2.max/weight.T1) %>%
  
 # group_by(timepoint, group) %>%
 # summarise(m = mean(VO2max.kg, na.rm = TRUE), 
           # s = sd(VO2max.kg, na.rm = TRUE)) %>%
  
  
  #ggplot(aes(timepoint, m, color = group, group = group)) + 
  #geom_errorbar(aes(ymin = m-s, ymax = m+s),
                #width = 0.05,
               # position = position_dodge(width = -0.2)) +
  #geom_point(position = position_dodge(width = -0.2)) + 
  #geom_line(position = position_dodge(width = -0.2)) + 
  #xlab("Time-point") + ylab(expression(VO["2max"]~~L~min^-1)) +
  #scale_color_manual(values = c("orange", "yellow", "blue"), 
                     #name = "Training Strategy") +
  #theme_classic() +
  #theme(legend.position = c(1, 0.50),
        #legend.text = element_text(size = 6), 
       # legend.title = element_text(size = 8),
        ##legend.background = element_rect(fill = "lightgray"),
        #legend.key = element_rect(fill = "dark gray")) +
  
  #scale_y_continuous(limit = c(30, 80)) +
  #labs(caption = "Figure 1: Development of VO2max over time in the cycling study") 


#forandre bakgrunn punkter ilegend
#Forandre punkter i legend
#scale_x_discrete(labels=c("pre" = "Pre-training", 
                         # "meso1" = "Meso-cycle 1",
                         # "meso2" = "Meso-cycle 2", 
                         # "meso3" = "Meso-cycle 3"))














