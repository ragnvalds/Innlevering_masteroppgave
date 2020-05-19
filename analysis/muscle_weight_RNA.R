#######per muscle weight RNA#######

source("./R/lib_fun.R")
#use estimates from seq and qpcr to compare


##load data
qpcr_results <- readRDS("./derivedData/qpcr_results.RDS")
  
  
time_rest <- readRDS("./derivedData/DE/mixedmodel2_timemodel.RDS")
full_rest <- readRDS("./derivedData/DE/mixedmodel2_fullmodel.RDS")

time_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel_timeonly.RDS")
full_acute <- readRDS("./derivedData/DE/mixedmodel2_acutemodel.RDS")


# Ensemble IDs from file
linc_ensembl <- read_excel("./data/name_LINCs.xlsx") %>%
  dplyr::select(lincs, gene) %>%
  filter(!is.na(gene)) %>%
  mutate(Gene = lincs, 
         ensembl = gene, 
         gene = Gene) %>%
  dplyr::select(gene, ensembl) %>%
  data.frame()
  
  
##MUSCLE weight cDNA synthesis

# Muscle weights in cDNA synthesis are loaded ##
rna_dat <- read_excel("./data/RNAamount.xlsx") %>%
  mutate(RNA.per.mg = (conc*elution.volume) / prot.mrna1) %>% 
  filter(include == "incl") %>%
  #filter(timepoint != "w2post") %>%
  dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
  mutate(tissue = 1000 / RNA.per.mg, 
         tl = log(tissue))  %>%
  print()
  
  
  # Visualize muscle weight to RNA relationship
  rna_dat %>%
    ggplot(aes(log(tissue), log(RNA.per.mg))) + geom_point() + geom_smooth(method = "lm")
  
  # There might be sample loss, use the below code to create an outlier variable
  
  # Create linear model 
  mod <- with(rna_dat, lm(log(RNA.per.mg) ~ log(tissue)))
  
  
  # Calculates standardized residuals
  rna_dat$resid<- resid(mod)/sd(resid(mod))
  # store predicted values for diagnostic plotting
  rna_dat$pred<- predict(mod)
  
  rna_dat %>%
    ggplot(aes(pred, resid)) + geom_point() +
    scale_y_continuous(breaks = c(-20, -10, -5, -4, 0, 4, 5), 
                       limits = c(-20, 5))
  
  ## Adding outlier variable
  rna_dat %>%
    mutate(outlier = if_else(resid < -4 | resid > 4 , "out", "in")) %>%# 4 units deviance from zero
    ggplot(aes(pred, resid, color = outlier)) + geom_point() +
    scale_y_continuous(breaks = c(-20, -10, -5, -4, 0, 4, 5), 
                       limits = c(-20, 5))
  
  rna <- rna_dat %>%
    filter(sets %in% c("multiple", "single")) %>%
    mutate(outlier = if_else(resid < -4 | resid > 4 , "out", "in"), # 4 units deviance from zer
           timepoint = factor(time, levels = c("w0","w2pre", "w2post", "w12")),
           sample = paste0(subject, sets, time), 
          rna.weight = RNA.per.mg/tissue) %>%
    print()
  
  ####### Modelling ###########
  
  
  rna %>%
    filter(sets == "single") %>%
    group_by(subject, timepoint, sets) %>%
    summarise(rna = mean(rna.weight, na.rm = TRUE)) %>%
    ggplot(aes(timepoint, rna, group = paste(subject, sets), color = as.factor(subject))) + geom_point() + 
    geom_line()+
    geom_point(position = position_dodge(width = 0.4)) +
    geom_line(position = position_dodge(width = 0.4)) +
    labs(caption = "Figure 9: RNA.per.mg in RNA seq data, single sets")+
    xlab("Timepoint") +
    ylab("RNA.per.mg")
  
  
  
  rna %>%
    filter(sets == "multiple") %>%
    group_by(subject, timepoint, sets) %>%
    summarise(rna = mean(rna.weight, na.rm = TRUE)) %>%
    ggplot(aes(timepoint, rna, group = paste(subject, sets), color = as.factor(subject))) + geom_point() + 
    geom_line()+
    geom_point(position = position_dodge(width = 0.4)) +
    geom_line(position = position_dodge(width = 0.4)) +
    labs(caption = "Figure 10: RNA.per.mg in RNA seq data, multiple sets")+
    xlab("Timepoint") +
    ylab("RNA.per.mg")
    
  
    
    #ggplot(aes(timepoint, m, color = group, group = group)) + 
    #geom_errorbar(aes(ymin = m-s, ymax = m+s),
    #width = 0.05,
    # position = position_dodge(width = -0.2)) +
    #geom_point(position = position_dodge(width = -0.2)) + 
    #geom_line(position = position_dodge(width = -0.2)) + 
    #xlab("Time-point") + ylab(expression(VO["2max"]~~L~min^-1)) +
    #scale_color_manual(values = c("orange", "yellow", "blue"), 
  #name = "Training Strategy") +
  #theme_classic() +
  #theme(legend.position = c(1, 0.50),
  #legend.text = element_text(size = 6), 
  # legend.title = element_text(size = 8),
  ##legend.background = element_rect(fill = "lightgray"),
  #legend.key = element_rect(fill = "dark gray")) +
  
  #scale_y_continuous(limit = c(30, 80)) +
  #labs(caption = "Figure 1: Development of VO2max over time in the cycling study") 
  
  
  #forandre bakgrunn punkter ilegend
  #Forandre punkter i legend
  #scale_x_discrete(labels=c("pre" = "Pre-training", 
  # "meso1" = "Meso-cycle 1",
  # "meso2" = "Meso-cycle 2", 
