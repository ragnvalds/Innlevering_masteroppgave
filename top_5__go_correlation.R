# Source libraries

source("./R/lib_fun.R")

#Use library(DGCA)?



#import top 5 lncRNAs to compare profile with all the genes in the seq data.

#top5_w2pre <- readRDS("./derivedData/top5_w2pre.RDS")
#top5_w2post <- readRDS("./derivedData/top5_w2post.RDS")
#top5_w12 <- readRDS("./derivedData/top5_w12.RDS")


#import all seq data with info
all_seq <- readRDS("./derivedData/seq_counts_all_info.RDS")

##########calculate changescore seq_data

seq_change <- all_seq %>%
  #filter(sets == "single") %>%
  spread(time, gene) %>%
  mutate(change = (w12 / w2pre)) %>%
  print()
  #summarise(m =mean(change, na.rm =TRUE), s =sd(change, na.rm =TRUE)) 

###############changescore qpcr data###########
##pcr_change <- readRDS("./derivedData/qpcr_results_estimated_means.RDS")


#############calculate changescore selected seq_data

#top 5 genes based on significance level. Need to have top five based on change score. 
w2pre_sig <- c("ENSG00000214548", 	
           "ENSG00000172965", 
           "ENSG00000204054", 	
           "ENSG00000230630", 
           "ENSG00000130600")

w2post_sig <- c("ENSG00000259820", 
            "ENSG00000221817", 
            "ENSG00000270605", 
            "ENSG00000242902", 
            "ENSG00000212719")

w12_sig <- c("ENSG00000286214", 
         "ENSG00000250208", 
         "ENSG00000260807", 
         "ENSG00000286191", 
         "ENSG00000272168")


#####top 5 based on change score.
seq_change %>% 
  mutate(pt = if_else(p.adj < 0.01, "sig", "ns") %>% 
           print()


pre <- arrange(seq_change,change )

######sort out 5 best p change scores. Try to#############

  
pre5 <- pre %>% 
  top_n(40, -change) %>% 
  print()





<- seq_change %>%
  filter(gene %in% w2pre) %>%
  print()

