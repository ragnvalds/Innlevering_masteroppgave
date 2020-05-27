
source("./R/lib_fun.R")



time_rest <- readRDS("./derivedData/DE/mixedmodel2_timemodel.RDS")
full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")

time_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel_timeonly.RDS")
full_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel.RDS")


# load list of all genes with ensemble info
all_genes <- readRDS("./derivedData/all_genes.RDS")

#load lncRNA

lncRNA <- readRDS("./derivedData/lncRNA.RDS")



#####################make tables for results. 40 highest p.value##############


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


DE_rest_w2 <- restw2 %>%
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig" & 
           #inc_dec == "increase" &
           #fc == "sig" & 
           fc == "sig") %>%
  print()

# DE timew2pre 169
# DE increase 164
# DE decrease 5




###arrange dataset on p.adjust

w2_arranged <- arrange(DE_rest_w2,p.adj )

######sort out 40 best p adj values. Try to use top_n#############

w2_p_table <- w2_arranged %>% 
  top_n(40, -p.adj) %>% 
  print()


##innerjoin with lncRNA to get hgnc symbol####
w2_p_table <- w2_p_table %>% inner_join(lncRNA)


####write excel spreadsheet####
library(writexl)
write_xlsx(w2_p_table, "./figures/w2_p_table.xlsx")

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




DE_rest_w12 <- restw12 %>%
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig" & 
           #inc_dec == "decrease" &
           #fc == "sig" & 
           fc == "sig") %>%
  print()


# DE lnc timew12 64
# 64 increase




###arrange dataset on p.adjust

w12_arranged <- arrange(DE_rest_w12,p.adj )

######sort out 40 best p adj values. Try to use top_n#############

w12_p_table <- w12_arranged %>% 
  top_n(40, -p.adj) %>% 
  print()


##innerjoin with lncRNA to get hgnc symbol####
w12_p_table <- w12_p_table %>% inner_join(lncRNA)


####write excel spreadsheet####
library(writexl)
write_xlsx(w12_p_table, "./figures/w12_p_table.xlsx")


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




DE_acute_w2post <- acute_w2post %>%
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns"),
         fc = if_else(estimate < -0.5 | estimate > 0.5 , "sig", "ns" )) %>% #significant if over 0,5 and -0,5
  mutate(inc_dec = if_else(estimate >0.00,"increase", "decrease")) %>%
  filter(pt == "sig" & 
           #inc_dec == "decrease" &
           #fc == "sig" & 
           fc == "sig") %>%
  print()



#DE w2post 102
#increase 40
# decrease 62


###arrange dataset on p.adjust

w2post_arranged <- arrange(DE_acute_w2post,p.adj )

######sort out 40 best p adj values. Try to use top_n#############

w2post_p_table <- w2post_arranged %>% 
  top_n(40, -p.adj) %>% 
  print()


##innerjoin with lncRNA to get hgnc symbol####
w2post_p_table <- w2post_p_table %>% inner_join(lncRNA)

####write excel spreadsheet####
library(writexl)
write_xlsx(w2post_p_table, "./figures/w2post_p_table.xlsx")


# Libraries
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)

# Load datasets
data1 <- DE_rest_w2 
data2 <- DE_acute_w2post 
data3 <- DE_rest_w12 

##trying to rbind datasets. It works. Now we have a dataset with all the sig lnc across different timepoints.

# Add datasets vertically, remember to ungroup.
all_data <- rbind(data1, data2, data3) %>% 
  ungroup %>% 
  dplyr::select(gene, coef)



# library
library(VennDiagram)

#make x

x <- list(
  all_data %>% filter(coef =="timew2pre") %>% dplyr:: select(gene) %>% unlist(), 
  all_data %>% filter(coef =="timew2post")  %>% dplyr:: select(gene) %>% unlist(), 
  all_data %>% filter(coef =="timew12") %>% dplyr:: select(gene) %>% unlist()) 

#cMake the plot
venn.diagram(x,
             category.names = c("Timew2pre(169)" , "Timew2post (102)" , "Timew12 (64)"),
             output = TRUE,
             imagetype="png") %>% 
  print()
          
height = 480 , 
width = 480 , 
resolution = 300,
compression = "lzw",
lwd = 1,
col=c("#440154ff", '#21908dff', '#fde725ff'),
fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
cex = 0.5,
fontfamily = "sans",
cat.cex = 0.3,
cat.default.pos = "outer",
cat.pos = c(-27, 27, 135),
cat.dist = c(0.055, 0.055, 0.085),
cat.fontfamily = "sans",
cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
rotation = 1
)




###################use venndiagram from bioconductor#####################


             
             
             
             
             
             
           