#####analysis high and low volume#####

source("./R/lib_fun.R")


# Load data
# Mixed-effects negative binomial count models were fitted and saved in 
# ./R/dge_list_models.R. Results saved in RDS files (see below).
# no need to run these models again.



full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")
full_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel.RDS")




##### Find lncRNA from ensamble #########

# Makes a vector of all transcripts (after filtering)
all_genes <- unique(full_rest$gene) 


ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# List attributes from biomart
listAttributes(ensembl)

# Make a data frame of all genes with symbol, biotype etc...
all_genes <- getBM(attributes=c('ensembl_gene_id', 'gene_biotype', 'entrezgene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id',
                   values = all_genes, 
                   mart = ensembl)

#filter lncRNAs
lncRNA <- all_genes %>%
  filter(gene_biotype == "lncRNA") %>%
  print()




##### Filter and adjust p-values for differentially expression #####

######rest model full #####

lnc_rest <- full_rest %>%
  filter(model %in% c("tissue_offset_lib_size_normalized")) %>% 
         #coef %in% c("(intercept)", "timew2pre", "timew12", "timew0:setsmultiple", "timew2pre:setsmultiple", "timew12:setsmultiple")) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(gene, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>% 
  print()


lnc_rest %>%
  ggplot(aes(estimate, -log10(p.val), color = pt)) + geom_point() + 
  facet_grid(model ~ coef)+
  xlab("Estimate") +
  labs(caption = "")



#######acute model full#######


lnc_acute <- full_acute %>%
  filter(model %in% c("lib_size_normalized")) %>% 
  #coef %in% c("(intercept)", "timew2post", "setsmultiple", "timew2post:setsmultiple" )) %>%
  
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(gene, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>% 
  print()


lnc_acute %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  ggplot(aes(estimate, -log10(p.val), color = pt)) + geom_point() + 
  facet_grid(model ~ coef)+
  xlab("Estimate") +
  labs(caption = "")










######moderate volume####

###	**av disse LNCs endret 114 (2 weeks)**
  
qpcr_t_mod <-qpcr_results %>%
  filter(coef %in% "timew2pre:setsmultiple") %>% 
  print()

two_rest_mod <- qpcr_t_mod %>% 
  filter(pt == "sig") %>%
  NROW() %>% 
  print()


#increase/decrease

inc_dec_tw_mod <- two_w_mod_lnc %>% 
  filter(pt == "sig") %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  print()

inc_two_mod <- inc_dec_tw_mod %>% 
  filter(inc_dec=="increase") %>% 
  NROW() %>% 
  print()

#94 increase


dec_two_mod <- inc_dec_tw_mod %>% 
  filter(inc_dec=="decrease") %>% 
  NROW() %>% 
  print()
#20 decrease



###76 (2weeks acute)

two_w_acute_mod <- qpcr_results %>%
  filter(model %in% c("lib_size_normalized"), 
         coef %in% c("timew2post:setsmultiple" )) %>%
  print()

two_acut_mod <- two_w_acute_mod%>% 
  filter(pt == "sig") %>%
  NROW() %>% 
  print()


#increase/decrease

inc_dec_tw_a_mod <- two_w_acute_mod %>% 
  filter(pt == "sig") %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  print()


inc_two_a_mod <- inc_dec_tw_a_mod %>% 
  filter(inc_dec=="increase") %>% 
  NROW() %>% 
  print()

# 31 increase


dec_two_a_mod <- inc_dec_tw_a_mod %>% 
  filter(inc_dec=="decrease") %>% 
  NROW() %>% 
  print()
# 45 decrease




###82 DE 12 weeeks med moderat treningsvolum 

  twelve_w_rest_mod <- lnc_rest %>%
  filter(model %in% c("tissue_offset_lib_size_normalized"),
         coef %in% c("timew12:setsmultiple")) %>% 
           print()

  
  twelve_rest_mod <- twelve_w_rest_mod%>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()

  
  
  #increase/decrease
  
  
  inc_dec_twelve_mod <- twelve_w_rest_mod %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
    print()
  
  
  inc_twelve_mod <- inc_dec_twelve_mod %>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  # 51 increase
  
  
  dec_twelve_mod <- inc_dec_twelve_mod %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  
  # 31 decrease 
  
  
  
  
  
  ####low training volume#####

##2 weeks##
  two_w_low_lnc <-lnc_rest %>%
    filter(model %in% c("tissue_offset_lib_size_normalized"),
           coef %in% "timew2pre") %>% 
    print()
  
  two_rest_low <- two_w_low_lnc%>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  #686 lnc changed
  
  
  #increase/decrease
  
  inc_dec_tw_low <- two_w_low_lnc %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
    print()
  
  
  inc_two_low <- inc_dec_tw_low %>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  # 634 increase
  
  
  dec_two_low <- inc_dec_tw_low %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  
  # 52 decrease 
  
  
  

###2weeks acute####
  two_w_acute_low <- lnc_acute %>%
    filter(model %in% c("lib_size_normalized"), 
           coef %in% c("timew2post" )) %>%
    print()
  
  
  two_acute_low <- two_w_acute_low%>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  #544 lnc changed
  
  
  
  #increase/decrease
  inc_dec_tw_a_low <- two_w_acute_low %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
    print()
  
  
  inc_tw_a_low <- inc_dec_tw_a_low %>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  # 196 increase
  
  
  dec_tw_a_low <- inc_dec_tw_a_low %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  
  # 348 decrease
  
  
  
  
###12 weeks####
  twelve_w_rest_low <- lnc_rest %>%
    filter(model %in% c("tissue_offset_lib_size_normalized"),
           coef %in% c("timew12")) %>% 
    print()
  
  
  twelve_rest_low <- twelve_w_rest_low%>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  #432 lnc changed
  
  
  
  #increase/decrease
  
  inc_dec_twelve_low <- twelve_w_rest_low %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
    print()
  
  
  inc_twelve_low <- inc_dec_twelve_low %>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  # 428 increase
  
  
  dec_twelve_low <- inc_dec_twelve_low %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  
  # 4 decrease
  
  

  

