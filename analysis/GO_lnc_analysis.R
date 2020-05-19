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


###### Acute model time only #####

full_rest_lnc <- full_rest %>%
  filter(model %in% c("lib_size_normalized"),
         coef %in% c("timew12")) %>%
  dplyr::select(gene, model, coef, estimate, se, z.val, p.val) %>%
  
  # filter out only lncRNA
  
  filter(gene %in% lncRNA$ensembl_gene_id) %>%
  # P-value adjustments
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  print()

###### Go analysis #######
library(clusterProfiler); library(org.Hs.eg.db)

temp <- full_rest_lnc %>%
  filter(p.adj < 0.05) %>%
  filter(coef == "timew12") %>%
  inner_join(lncRNA %>% mutate(gene = ensembl_gene_id)) %>%
  pull(gene)

go <- enrichGO(gene = temp,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "MF",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               universe,
               qvalueCutoff = 0.2,
               minGSSize = 10,
               maxGSSize = 500,
               readable = FALSE,
               pool = TRUE)



go_numbers <- getBM(attributes=c('ensembl_gene_id', 'gene_biotype', 'entrezgene_id'), 
                    filters = 'ensembl_gene_id',
                    values = temp, 
                    mart = ensembl)

go_numbers
