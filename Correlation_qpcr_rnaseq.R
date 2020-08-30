################Correlation qpcr vs RNAseq

source("./R/lib_fun.R")









qpcr_results_norm <- readRDS("./derivedData/qpcr_results_norm.RDS")








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
full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")

seq<- time_rest %>%
   filter(gene %in% linc_ensembl$ensembl, 
          #coef == "timew12", 
          model == "tissue_offset_lib_size_normalized") %>%
   dplyr::select(gene, estimate, coef) %>%
   #mutate(type = "rnaseq") %>%
 print()

#sort out "nf" in coef

seqt <- seq %>% 
   filter(!(coef == "nf")) %>% 
   filter(!(coef == "(Intercept)")) %>%
   print()


seqtime <- seqt %>% 
   sort(coef) %>% 
   print()
   


######sort out genes thats not in seq data
gout <- c("ENSG00000249515", "ENSG00000249515", "ENSG00000249859")

q <- qpcr_results_norm %>% 
   filter(!gene %in% gout)


r <- q %>%
   dplyr::select(gene, estimate = Value, coef) %>%
   mutate(labels=c("Intercept" = "(Intercept)",
                   "Intercept" = "(Intercept)",
                   "Intercept" = "(Intercept)",
                   "Intercept" = "(Intercept)",
                   "Intercept" = "(Intercept)",
                   "Intercept" = "(Intercept)",
            "timepointw2pre" = "timew2pre",
            "timepointw2pre" = "timew2pre",
            "timepointw2pre" = "timew2pre",
            "timepointw2pre" = "timew2pre",
            "timepointw2pre" = "timew2pre",
            "timepointw2pre" = "timew2pre",
            "	timepointw12" = "timew12",
            "	timepointw12" = "timew12",
            "	timepointw12" = "timew12",
            "	timepointw12" = "timew12",
            "	timepointw12" = "timew12",
            "	timepointw12" = "timew12",
            "timepointw0:setsmultiple" = "timew0:setsmultiple",
            "timepointw0:setsmultiple" = "timew0:setsmultiple",
            "timepointw0:setsmultiple" = "timew0:setsmultiple",
            "timepointw0:setsmultiple" = "timew0:setsmultiple",
            "timepointw0:setsmultiple" = "timew0:setsmultiple",
            "timepointw0:setsmultiple" = "timew0:setsmultiple",
            "timepointw2pre:setsmultiple" = "timew2pre:setsmultiple",
            "timepointw2pre:setsmultiple" = "timew2pre:setsmultiple",
            "timepointw2pre:setsmultiple" = "timew2pre:setsmultiple",
            "timepointw2pre:setsmultiple" = "timew2pre:setsmultiple",
            "timepointw2pre:setsmultiple" = "timew2pre:setsmultiple",
            "timepointw2pre:setsmultiple" = "timew2pre:setsmultiple",
            "timepointw12:setsmultiple" = "timew12:setsmultiple",
            "timepointw12:setsmultiple" = "timew12:setsmultiple",
            "timepointw12:setsmultiple" = "timew12:setsmultiple",
            "timepointw12:setsmultiple" = "timew12:setsmultiple",
            "timepointw12:setsmultiple" = "timew12:setsmultiple",
            "timepointw12:setsmultiple" = "timew12:setsmultiple")) %>% 
   #mutate(type = "qpcr") %>%
   print()

qpcr1 <- r %>% 
   mutate(coef = labels, labels = coef) %>% 
   filter(!(coef =="timew0:setsmultiple")) %>% 
   print()
   
qpcr2 <- qpcr1 %>% 
   filter(!(coef =="timew2pre:setsmultiple")) %>%
   dplyr::select(!(labels)) %>% 
   print()

qpcr3 <-qpcr2 %>% 
   filter(!(coef =="timew12:setsmultiple")) %>% 
   print()

qpcr <-qpcr3 %>% 
   filter(!(coef =="(Intercept)")) %>% 
   print()


cor_qs <- merge(seqtime, qpcr,by = "gene", all = T) %>%
  
   
   print()
 

###saveRDS(cor_qs, "./derivedData/cor_qs.RDS")
 ####qs <- readRDS("./derivedData/cor_qs.RDS")

qs <- cor_qs %>% 
   mutate(seq = estimate.x,
          qpcr = estimate.y) %>%
  print()
 
 
   
 

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
shapiro.test(qs$qpcr) # => p = 0.1333


library("ggpubr")
# seq
ggqqplot(qs$seq, ylab = "seq")
# qpcr
ggqqplot(qs$qpcr, ylab = "qpcr")

res <- cor.test(qs$seq, qs$qpcr, 
                method = "pearson")
res

#data:  qs$seq and qs$qpcr

#t = 2.2439, df = 4,
#p-value = 0.08824
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
   #-0.1650172  0.9702552
#sample estimates:
   #cor 
#0.7465191  




#In the result above :
   
  # t is the t-test statistic value (t = 2.2439),
#df is the degrees of freedom (df= 4),
#p-value is the significance level of the t-test (p-value = 0.08824).
#conf.int is the confidence interval of the correlation coefficient at 95% (conf.int = [-0.1650172  0.9702552]);
#sample estimates is the correlation coefficient (Cor.coeff = 0.7465191).


#Interpretation of the result
#The p-value of the test is 0.03229, which is less than the significance level alpha = 0.05. We can conclude that 
#seq and qpcr are significantly correlated with a correlation coefficient of -0.85 and p-value of 0.03229 .

#Access to the values returned by cor.test() function
#The function cor.test() returns a list containing the following components:
   
   #p.value: the p-value of the test
#estimate: the correlation coefficient
# Extract the p.value
