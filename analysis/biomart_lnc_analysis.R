#####lncRNA information from biomart#####


library(tidyverse)
library(readxl)
library(knitr)
library(dplyr)
library(biomaRt)
library(kableExtra)
### Find lncRNA from ensamble #########


time_rest <- readRDS("./derivedData/DE/mixedmodel2_timemodel.RDS") #na="NA" tells read_excel how missing values are coded

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

lncRNA <- data.frame (all_genes, na = "NA") %>%
  filter(gene_biotype == "lncRNA") %>%
  print()

###make a table of all lncRNAs identified

lnc_table <- lncRNA %>% 
dplyr::select(ensembl_gene_id, gene_biotype, entrezgene_id, hgnc_symbol) %>%
  kable(format = "html", col.names = c("Ensemble gene ID", 
                                       "Gene biotype",
                                       "Entrezgene ID", 
                                       "Hgnc symbol"), 
        caption = "Table 3: All lncRNAs  identified") %>% 
  kable_styling(bootstrap_options = "striped", font_size = 12) %>%
  print()
  
  




###pull number of lncRNAs identified


lnc_count <- lncRNA %>% 
  nrow() %>%
  print()


###lncRNAs entrezgene_id


e_id <- lncRNA %>% 
  pull(entrezgene_id)

#NROW gives you number of rows. na.omit removes NA.

ent_id_count <- NROW(na.omit(e_id)) %>%
  print()


remove_na<- data.frame (lncRNA) %>% 
  print()


###make table of lncs with entrezgene ID

#na.omit removes NA
entrezgene_id_lnc <- na.omit(remove_na) 


#use select to select wanted variables
entrezgene_id_lnc %>% 
  dplyr::select(ensembl_gene_id, gene_biotype, entrezgene_id, hgnc_symbol) %>%
kable(format = "html", col.names = c("Ensemble gene ID", 
                                     "Gene biotype",
                                     "Entrezgene ID", 
                                     "Hgnc symbol"), 
      caption = "Table 3: All lncRNAs  identified with entrexgene ID") %>% 
  kable_styling(bootstrap_options = "striped", font_size = 12) %>%  
  save_kable(file = "figures/all_entrezgene_id.png")



#does not work because NA is not recognized.
##lncRNAs hgnc symbol
count_hgnc <- lncRNA %>% pull(hgnc_symbol) %>% 
  print()

NROW(na.omit(count_hgnc)) %>%
  print()





###lncRNAs used in qPCR analysis found in seq data

goi <- c("ENSG00000234741", "ENSG00000268518",  "ENSG00000225613",  "ENSG00000185847",
         "ENSG00000249515",  "ENSG00000223784",  "ENSG00000249859",  "ENSG00000269900")

qPCR_Lnc <- all_genes %>%
  filter(ensembl_gene_id %in% goi) %>%
  print()

qPCR_Lncs_identf <- qPCR_Lnc %>% 
  nrow() %>%
  print()