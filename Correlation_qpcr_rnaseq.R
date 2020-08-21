################Correlation qpcr vs RNAseq

source("./R/lib_fun.R")

library("ggpubr")

##load data
qpcr_results <- readRDS("./derivedData/qpcr_results.RDS")

qpcr_results_norm <- readRDS("./derivedData/qpcr_results_norm.RDS")




#first normalize qpcr to rna pr mg
rna_dat <- read_excel("./data/RNAamount.xlsx") %>%
   mutate(RNA.per.mg = (conc*elution.volume) / prot.mrna1) %>% 
   filter(include == "incl") %>%
   filter(timepoint != "w2post",
          sets == "single") %>%
   dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
   mutate(tissue = 1000 / RNA.per.mg, 
          tl = log(tissue),
          rl = log(RNA.per.mg),
          rna.weight = RNA.per.mg/tissue)  %>%
   print()
  


rtest<- rna_dat %>% 
   summarise(mean(rna.weight))



#normalizing factor 2.874087

#log factor 1.041577

#normalizing on RNA.pr,mg.factor 357.6457
#log 5.866178


# Ensemble IDs from file
linc_ensembl <- read_excel("./data/name_LINCs.xlsx") %>%
   dplyr::select(lincs, gene) %>%
   filter(!is.na(gene)) %>%
   mutate(Gene = lincs, 
          ensembl = gene, 
          gene = Gene) %>%
   dplyr::select(gene, ensembl) %>%
   data.frame()


time_rest <- readRDS("./derivedData/DE/mixedmodel2_timemodel.RDS")

seq<- time_rest %>%
   filter(gene %in% linc_ensembl$ensembl, 
          coef == "timew12", 
          model == "tissue_offset_lib_size_normalized") %>%
   dplyr::select(gene, estimate) %>%
   mutate(type = "rnaseq") %>%
   #mutate(norm = estimate/1.041577) %>% 
   print()


qpcr <- qpcr_results %>%
   filter(coef == "timepointw12",
          sets == "setssingle") %>%
   dplyr::select(gene, estimate = Value) %>%
   mutate(type = "qpcr") %>%
   #mutate(norm = estimate/	131.2127) %>% 
   print()

cor_qs <- merge(seq, qpcr, by="gene", all = T) %>% 
   print()

saveRDS(cor_qs, "./derivedData/cor_qs.RDS")
 qs <- readRDS("./derivedData/cor_qs.RDS")

 gout <- c("ENSG00000249515", "ENSG00000249515", "ENSG00000249859")
 
 qs <- cor_qs %>% 
  filter(!gene %in% gout)
 
 qs <- qs %>% 
    mutate(seq = estimate.x,
           qpcr = estimate.y)
 
   
 

ggscatter(qs, x = "qpcr" , y = "seq", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "RNAseq", ylab = "qPCR")
#Is the covariation linear? Yes, form the plot above, the relationship is linear. In the situation where the scatter plots show curved patterns, we are dealing with nonlinear association between the two variables.

#Are the data from each of the 2 variables (x, y) follow a normal distribution?
   #Use Shapiro-Wilk normality test –> R function: shapiro.test()
#and look at the normality plot —> R function: ggpubr::ggqqplot()
#hapiro-Wilk test can be performed as follow:
   #Null hypothesis: the data are normally distributed
#Alternative hypothesis: the data are not normally distributed

#From the output, the two p-values are greater than the significance level 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.

#Visual inspection of the data normality using Q-Q plots (quantile-quantile plots). Q-Q plot draws the correlation between a given sample and the normal distribution.

# Shapiro-Wilk normality test for seq
shapiro.test(qs$seq) # => p = 0.7093
# Shapiro-Wilk normality test for qpcr
shapiro.test(qs$qpcr) # => p = 0.6963


library("ggpubr")
# seq
ggqqplot(qs$seq, ylab = "seq")
# qpcr
ggqqplot(qs$qpcr, ylab = "qpcr")

res <- cor.test(qs$seq, qs$qpcr, 
                method = "pearson")
res

#data:  qs$seq and qs$qpcr

#t = -3.2196, df = 4,
#p-value = 0.03229
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
   #-0.9832071 -0.1219751
#sample estimates:
   #cor 
#-0.8494492 




#In the result above :
   
  # t is the t-test statistic value (t = -3.2196),
#df is the degrees of freedom (df= 4),
#p-value is the significance level of the t-test (p-value = 0.03229).
#conf.int is the confidence interval of the correlation coefficient at 95% (conf.int = [-0.9832071, -0.1219751]);
#sample estimates is the correlation coefficient (Cor.coeff = -0.8494492).


#Interpretation of the result
#The p-value of the test is 0.03229, which is less than the significance level alpha = 0.05. We can conclude that 
#seq and qpcr are significantly correlated with a correlation coefficient of -0.85 and p-value of 0.03229 .

#Access to the values returned by cor.test() function
#The function cor.test() returns a list containing the following components:
   
   #p.value: the p-value of the test
#estimate: the correlation coefficient
# Extract the p.value
