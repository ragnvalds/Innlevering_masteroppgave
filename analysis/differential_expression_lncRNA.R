
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




### Find lncRNA from ensamble #########

# Makes a vector of all transcripts (after filtering)
all_genes <- unique(time_rest$gene) 



ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# List attributes from biomart
listAttributes(ensembl)

# Make a data frame of all genes with symbol, biotype etc...
all_genes <- getBM(attributes=c('ensembl_gene_id', 'gene_biotype', 'entrezgene_id', 'hgnc_symbol'), 
              filters = 'ensembl_gene_id',
      values = all_genes, 
      mart = ensembl)







## Genes that are also analyzed by qPCR
goi <- c("ENSG00000234741", "ENSG00000268518",  "ENSG00000225613",  "ENSG00000185847",
         "ENSG00000249515",  "ENSG00000223784",  "ENSG00000249859",  "ENSG00000269900")

all_genes %>%
  filter(ensembl_gene_id %in% goi) %>%
  print()


#####lncRNA identified from seq data#####


#### make a data frame of lncRNA 

lncRNA <- all_genes %>%
  filter(gene_biotype == "lncRNA") %>%
  print()




###pull number of lncRNAs idetified


lnc_count <- lncRNA %>% 
  nrow() %>%
  print()




##### Filter and adjust p-values for differentially expression #####


####rest model timepoint and sets single sets#####


  
full_rest_lnc <- full_rest %>%
  filter(model %in% c("lib_size_normalized", "tissue_offset_lib_size_normalized"),
         coef %in% c("timew2pre", "timew12")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_rest_lnc %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
 ggplot(aes(estimate, -log10(p.val), color = paste(pt,inc_dec))) + geom_point() + 
  facet_grid(model ~ coef)+
  xlab("log2FC")+
  labs(title = "Rested state,single sets")


  
  



#####acute model timepoint and sets single sets#####

  full_acute_lnc <- full_acute %>%
  filter(model %in% c("lib_size_normalized"),
         coef %in% c("timew2post")) %>%
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_acute_lnc %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  ggplot(aes(estimate, -log10(p.val), color = paste(pt,inc_dec))) + geom_point() +
  facet_grid(model ~ coef)+
  xlab("log2FC") +
  labs(title = "Acute respons,single sets")



############ multiple sets###############



####rest model timepoint and sets multiple sets#####


full_rest_lnc_m <- full_rest %>%
  filter(model %in% c("lib_size_normalized", "tissue_offset_lib_size_normalized"),
         coef %in% c("timew2pre:setsmultiple", "timew12:setsmultiple")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()


full_rest_lnc_m %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  ggplot(aes(estimate, -log10(p.val), color = paste(pt,inc_dec))) + geom_point() + 
  facet_grid(model ~ coef)+
  xlab("log2FC") +
  labs(title = "Rested state,multiple sets")






#####acute model timepoint and sets multiple sets#####

full_acute_lnc_m <- full_acute %>%
  filter(model %in% c("lib_size_normalized"),
         coef %in% c("timew2post:setsmultiple")) %>%
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()




full_acute_lnc_m %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  ggplot(aes(estimate, -log10(p.val), color = paste(pt,inc_dec))) + geom_point(aes(label=gene)) +
  facet_grid(model ~ coef)+
  xlab("log2FC") +
  labs(title = "Acute respons, multiple sets")







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














