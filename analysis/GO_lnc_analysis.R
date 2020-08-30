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

qpcr_lncs <- all_genes %>%
  filter(ensembl_gene_id %in% goi) %>%
  print()




#####lncRNA identified from seq data#####


#### make a data frame of lncRNA 

lncRNA <- all_genes %>%
  filter(gene_biotype == "lncRNA") %>%
  print()




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





##################### Go analysis ###################################

library(clusterProfiler); library(org.Hs.eg.db)


#####################try with all genes############
#use all_genes


temp_all <- all_genes %>%
  #filter(p.adj < 0.05) %>%
  #filter(coef == "timew12") %>%
  #inner_join(lncRNA %>% mutate(gene = ensembl_gene_id)) %>%
  pull(entrezgene_id)

       
go_all <- enrichGO(gene = temp_all,
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



show(go_all)



###################make plot########################

library(enrichplot)

x <- go_all


#enrichplot::cnetplot(x)    To much info, overplotting


enrichplot::goplot(
  x,
  showCategory = 10,
  color = "p.adjust",
  layout = "sugiyama",
  geom = "text"
)

#####################with all lncRNAs##########

#use lncRNA



temp_lnc <- lncRNA %>%
  #filter(p.adj < 0.05) %>%
  #filter(coef == "timew12") %>%
  #inner_join(lncRNA %>% mutate(gene = ensembl_gene_id)) %>%
  pull(entrezgene_id)


go_lnc <- enrichGO(gene = temp_lnc,
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



show(go_lnc)



###################make plot########################

library(enrichplot)

y <- go_lnc


enrichplot::cnetplot(y)


enrichplot::goplot(
  y,
  showCategory = 10,
  color = "p.adjust",
  layout = "sugiyama",
  geom = "text"
)

#################try with full_rest lncRNAs ###################

temp <- full_rest_lnc %>%
  filter(p.adj < 0.05) %>%
  #filter(coef == "timew12") %>%
  inner_join(lncRNA %>% mutate(gene = ensembl_gene_id)) %>%
  pull(gene)


###########check if genes can be mapped

bitr_kegg(temp_lnc, fromType = "kegg", toType = "Path", organism = "hsa")
bitr_kegg(temp, fromType = "kegg", toType = "Module", organism = "hsa")


###########they can't be mapped. 100 %!!!!!!!!!!!!!!


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


###################make plot########################

library(enrichplot)

z <- go


enrichplot::cnetplot(z)


enrichplot::goplot(
  z,
  showCategory = 10,
  color = "p.adjust",
  layout = "sugiyama",
  geom = "text"
)


