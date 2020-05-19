# Source libraries

source("./R/lib_fun.R")

# Load data
# Mixed-effects negative binomial count models were fitted and saved in 
# ./R/dge_list_models.R. Results saved in RDS files (see below).
# no need to run these models again.

full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")
full_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel.RDS")
### Find lncRNA from ensamble #########

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










#####lncRNA identified from seq data#####


#### make a data frame of lncRNA 

lncRNA <- all_genes %>%
  filter(gene_biotype == "lncRNA") %>%
  print()


##### Filter and adjust p-values for differentially expression #####

######rest model time only #####

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


######DE lncs count at rest#####

#####Total DE lncs at rest######

DE_count <- full_rest_lnc %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  filter(model %in% "lib_size_normalized",
         coef %in% "timew2pre",
         pt %in% "sig") %>%
 print()

DE_count_lnc_rest <- DE_count %>% 
  NROW() %>% 
  print()
#NROW gives number of rows in a column.
#323 genes DE

#####count lnc DE increase######
pos_neg_lnc_rest <- DE_count %>%
    mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
           print()

inc_rest <- pos_neg_lnc_rest %>% 
  filter(inc_dec=="increase") %>% 
  NROW() %>% 
  print()

#223 DE increase
  

##### count lnc DE decrease######
dec_rest <- pos_neg_lnc_rest %>% 
  filter(inc_dec=="decrease") %>% 
  NROW() %>% 
  print()
 
#100 DE decrease
 

##### count DE qPCR lncs#####
## Genes that are also analyzed by qPCR
goi <- c("ENSG00000234741", "ENSG00000268518",  "ENSG00000225613",  "ENSG00000185847",
         "ENSG00000249515",  "ENSG00000223784",  "ENSG00000249859",  "ENSG00000269900")
DE_count_lnc_q <- DE_count %>%
  filter(gene %in% goi) %>%
  print()

DE_qpcr_rest <- DE_count_lnc_q %>%  pull(gene)

#"ENSG00000225613" and "ENSG00000234741 is DE



######DE lncs count acute#####

#####acute model timepoint and sets#####

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

#####Total DE lncs acute ######

DE_count_acute <- full_acute_lnc %>%
  mutate(pt = if_else(p.adj < 0.05, "sig", "ns")) %>%
  filter(model %in% "lib_size_normalized",
         coef %in% "timew2post",
         pt %in% "sig") %>%
  print()

DE_count_lnc_acute <- DE_count_acute %>% 
  NROW() %>% 
  print()
#NROW gives number of rows in a column.
#385 genes DE


#####count lnc DE acute increase######

pos_neg_lnc_acute <- DE_count_acute %>%
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>% 
  print()

inc_acute <- pos_neg_lnc_acute %>% 
  filter(inc_dec=="increase") %>% 
  NROW() %>% 
  print()

#146 DE increase

##### count lnc DE acute decrease######


dec_acute <- pos_neg_lnc_acute %>% 
  filter(inc_dec=="decrease") %>% 
  NROW() %>% 
  print()

#239 DE decrease




##### count DE acute qPCR lncs#####
## Genes that are also analyzed by qPCR
goi <- c("ENSG00000234741", "ENSG00000268518",  "ENSG00000225613",  "ENSG00000185847",
         "ENSG00000249515",  "ENSG00000223784",  "ENSG00000249859",  "ENSG00000269900")
DE_count_lnc_acute_qpcr <- DE_count_acute %>%
  filter(gene %in% goi) %>%
  print()

DE_qpcr_acute <- DE_count_lnc_acute_qpcr %>%  pull(gene)


