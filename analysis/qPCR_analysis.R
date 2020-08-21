######qPCR analysis######

source("./R/lib_fun.R")

### Model qpcr data 


qpcr.dat <- readRDS("./derivedData/qpcr.replicates_LINCs.Rds")


qdat <- qpcr.dat %>%
  mutate(Ra = eff^-cq, 
         technical = paste0(subject, timepoint, leg, cdna), 
         biological = paste0(subject, timepoint, leg, sets)) %>%
  print()



qdat_norm <- qpcr.dat %>%
  mutate(Ra = eff^-cq, 
         technical = paste0(subject, timepoint, leg, cdna), 
         biological = paste0(subject, timepoint, leg, sets)) %>%
  mutate(Ra = Ra/1.041577) %>% 
print()





# make model

m1 <- lme(log(Ra) ~ 0 + target + target:timepoint + target:timepoint:sets,
          random = list(subject = ~ 1, technical = ~ 1), 
          data = qdat_norm, 
          control=list(msMaxIter=120,
                       opt = "nloptwrap", msVerbose=TRUE),
          method = "REML", na.action = na.exclude)

##### per total-RNA models ######

# Fixed and random effects formulas for the first step of model building
# Compared to Method paper (Hammarstr��t al 2018), the fixed effects are reuced
# to only contain gene-specific time + time:sets



## m1 is the model assuming homoscedastic errors.
#m1 <- lme(fixed, random = random, data = qdat.full,
#          control=list(msMaxIter=120,
#                       opt = "nloptwrap", msVerbose=TRUE),
#          method = "REML", na.action = na.exclude) # This allows to check progress
#

#### Models allowing for heteroscedasticity in the residuals. ####
# The variance is believed to be different from gene to gene due to different expression, primer design etc.
# One can also expect that residual variance to be related to cq-values as higher cq-values 
# implies higher meassurement error due to stochastic variation in the beginning of amplification.

# Variance functions:
varfun1 <- varIdent(form = ~1|target) # Heterogeneity per gene
varfun2 <- varPower(form = ~cq|target) # Power of the variance covariate cq, per gene stratum
varfun3 <- varFixed(~cq) # Variance increases as a function of cq
varfun4 <- varExp(form =~cq|target) # Variance changes as a exponentil of cq-values
varfun5 <- varConstPower(form = ~cq|target)

m2 <- update(m1, weights = varfun1)
# m3 <- update(m1, weights = varfun2)
m4 <- update(m1, weights = varfun3)
m5 <- update(m1, weights = varfun4)# m6 <- update(m1, weights = varfun5)

### Test models
anova(m1, m2, m5)


intervals(m5)



# Ensemble IDs from file
linc_ensembl <- read_excel("./data/name_LINCs.xlsx") %>%
  dplyr::select(lincs, gene) %>%
  filter(!is.na(gene)) %>%
  mutate(Gene = lincs, 
         ensembl = gene, 
         gene = Gene) %>%
  dplyr::select(gene, ensembl) %>%
  data.frame()



qpcr_results_norm <- coef(summary(m5)) %>%
data.frame() %>%
  mutate(gene = rownames(.)) %>%
  separate(gene, into = c("gene", "coef", "sets"), sep = ":") %>% 
  mutate(coef = if_else(is.na(coef), "Intercept", coef), 
         gene = as.character(gsub("target", "", gene)), 
         gene = gsub(" ", "_", gene)) %>%
  inner_join(linc_ensembl) %>%
  mutate(gene = ensembl) %>%
  mutate(sets = if_else(is.na(sets), "Multiple", sets)) %>% 
  print()

##save data as RDS

saveRDS(qpcr_results_norm, file = "./derivedData/qpcr_results_norm.RDS")

                      #######lnc from qPCR change on different timepoint########
  
  

                          ######moderate training volume######
  
##load data
qpcr_results <- readRDS("./derivedData/qpcr_results.RDS")

###	**av disse LNCs endret 3 (2 weeks)**
  
  
  qpcr_t_mod <-qpcr_results %>%
    mutate(pt = if_else(p.value < 0.05, "sig", "ns")) %>%
    filter(coef %in% "timepointw2pre",
           sets %in% "Multiple") %>% 
    print()
  
  
  qpcr_two_mod <- qpcr_t_mod  %>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  

  
  #increase/decrease
  
  q_inc_dec_two_mod <- qpcr_t_mod %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(Value >0.00,"increase", "decrease")) %>% 
    print()
  
  q_inc_two_mod <- q_inc_dec_two_mod%>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  #2 increase
  
  
  q_dec_two_mod <- q_inc_dec_two_mod %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  #1 decrease
  
  
  
  
  
  
  ### 4 changed (2weeks acute)
  
  
  qpcr_ta_mod <-qpcr_results %>%
    mutate(pt = if_else(p.value < 0.05, "sig", "ns")) %>%
    filter(coef %in% "timepointw2post",
           sets %in% "Multiple") %>% 
    print()
  
  
  
  qpcr_two_a_mod <- qpcr_ta_mod  %>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  #4 changed
  
  
  
  
  #increase/decrease
  
  q_inc_dec_two_a_mod <- qpcr_ta_mod %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(Value >0.00,"increase", "decrease")) %>% 
    print()
  
  q_inc_two_a_mod <- q_inc_dec_two_a_mod%>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  #2 increase
  
  
  q_dec_two_a_mod <- q_inc_dec_two_a_mod %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  #2 decrease
  
  
  
  
  
  
  
  
  ### 3 changed at 12 weeeks med moderat treningsvolum 
  
  
  qpcr_twel_mod <-qpcr_results %>%
    mutate(pt = if_else(p.value < 0.05, "sig", "ns")) %>%
    filter(coef %in% "timepointw12",
           sets %in% "Multiple") %>% 
    print()
  
  
  qpcr_twelve_mod <- qpcr_twel_mod  %>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  #3 changed
  
  
  
  
  #increase/decrease
  
  q_inc_dec_twelve_mod <- qpcr_twel_mod %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(Value >0.00,"increase", "decrease")) %>% 
    print()
  
  q_inc_two_a_mod <- q_inc_dec_twelve_mod%>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  #3 increase
  
  
  q_dec_twelve_mod <-q_inc_dec_twelve_mod %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  
  #0 decrease
  
  
  
                         
  
  
  
                               ####low training volume#####
  
 
  
  ###	**av disse LNCs endret 1  at(2 weeks)**
  
  
  qpcr_t_low <-qpcr_results %>%
    mutate(pt = if_else(p.value < 0.05, "sig", "ns")) %>%
    filter(coef %in% "timepointw2pre",
           sets %in% "setssingle") %>% 
    print()
  
  
  qpcr_two_low <- qpcr_t_low  %>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  
  
  
  #increase/decrease
  
  q_inc_dec_two_low <- qpcr_t_low %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(Value >0.00,"increase", "decrease")) %>% 
    print()
  
  q_inc_two_low <- q_inc_dec_two_low%>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  #1 increase
  
  
  q_dec_two_low <- q_inc_dec_two_low %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  #0 decrease
  
  
  
  
  
  
  ### 1 changed (2weeks acute)
  
  
  qpcr_ta_low <-qpcr_results %>%
    mutate(pt = if_else(p.value < 0.05, "sig", "ns")) %>%
    filter(coef %in% "timepointw2post",
           sets %in% "setssingle") %>% 
    print()
  
  
  
  qpcr_two_a_low <- qpcr_ta_low  %>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  #1 changed
  
  
  
  
  #increase/decrease
  
  q_inc_dec_two_a_low <- qpcr_ta_low %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(Value >0.00,"increase", "decrease")) %>% 
    print()
  
  q_inc_two_a_low <- q_inc_dec_two_a_low%>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  #0 increase
  
  
  q_dec_two_a_low <- q_inc_dec_two_a_low %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  #1 decrease
  
  
  
  
  
  
  
  
  ### 0 changed at 12 weeeks med moderat treningsvolum 
  
  
  qpcr_twel_low <-qpcr_results %>%
    mutate(pt = if_else(p.value < 0.05, "sig", "ns")) %>%
    filter(coef %in% "timepointw12",
           sets %in% "setssingle") %>% 
    print()
  
  
  qpcr_twelve_low <- qpcr_twel_low  %>% 
    filter(pt == "sig") %>%
    NROW() %>% 
    print()
  #0 changed
  
  
  
  
  #increase/decrease
  
  q_inc_dec_twelve_low <- qpcr_twel_low %>% 
    filter(pt == "sig") %>%
    mutate(inc_dec = if_else(Value >0.00,"increase", "decrease")) %>% 
    print()
  
  q_inc_two_a_low <- q_inc_dec_twelve_low%>% 
    filter(inc_dec=="increase") %>% 
    NROW() %>% 
    print()
  
  #0 increase
  
  
  q_dec_twelve_low <-q_inc_dec_twelve_low %>% 
    filter(inc_dec=="decrease") %>% 
    NROW() %>% 
    print()
  #0 decrease
  
  







