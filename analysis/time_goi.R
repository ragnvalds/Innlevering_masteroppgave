##########find p-value and direction for the time_goi###########

#source libraries
source("./R/lib_fun.R")

#load data
time_goi <- readRDS("./derivedData/time_goi.RDS")

time_rest <- readRDS("./derivedData/DE/mixedmodel2_timemodel.RDS")

lncRNA <- readRDS("./derivedData/lncRNA.RDS") 

#filter goi in time_rest

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

superlnc <- restw2 %>%
  ungroup() %>% 
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig" & 
           #inc_dec == "increase" &
           #fc == "sig" & 
           fc == "sig") %>%
  filter(gene %in% time_goi) %>% 
  dplyr::select(gene, estimate, p.adj, inc_dec) %>% 
  print()


###arrange dataset on p.adjust


superlnc <- arrange(superlnc,p.adj )


  
  

superlnc1 <- superlnc %>% 
mutate(as.character(superlnc$p.adj)) %>%
  print()
  
  superlnc1 <- superlnc %>% 
  dplyr:: select(gene, estimate, inc_dec, "as.character(superlnc$p.adj)") %>% 
  print()

library(kableExtra)

superlnc1 %>% 
  kable(format = "html",col.names = c("Ensemble gene ID", 
                                      "log2FC",
                                      "Gene expression",
                                      "P.adjust"),
        align = c("l", "l", "l", "l"), results = "asis") %>% 
        #caption = "DE LncRNAs all timepoints ") %>% 
  kable_styling(bootstrap_options = "striped", font_size = 21) %>% 
  print()



