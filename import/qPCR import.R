


## qPCR import
 library(qpcrpal); library(dplyr); library(qpcR); library(readxl); library(ggplot2); library(stringr); library(tidyr)

### Prepare a batch of imported data 
 batch <- prepare_batch("./data/exports/", equipment = "quant", skip = 21) %>%
   model_qpcr()
 
 
 ## Perform model comparisons
model.tests <- test_models(batch)


# Plot best models
data.frame(table(model.tests$best.model, model.tests$target)) %>%
  ggplot(aes(Var2, Freq, fill = Var1)) + geom_bar(stat = "identity") + coord_flip()
  

# Best model per target are stored for modeling with best model
# Plot best model per target
data.frame(table(model.tests$best.model, model.tests$target)) %>%
  ggplot(aes(Var2, Freq, fill = Var1)) + geom_bar(stat="identity") + coord_flip()

# Best model per target are stored for modeling with best model
best.models <- data.frame(target = names(apply(table(model.tests$best.model, model.tests$target), 2, which.max)),
                          model = as.character(row.names(table(model.tests$best.model, model.tests$target))[apply(table(model.tests$best.model, model.tests$target), 2, which.max)]))

# Remove targets that are not to be used (bad primers)
best.models <- best.models %>%
  filter(!(target %in% c("NEAT1 F2R2", "MTs16 F3R3", 
                         "MTs12 F2R2", "MALAT1 F2R2", 
                         "Lnc31 F3R3", "LincRAM F2R2", 
                         "LincPINT F2R2", "LincP21 F3R3", 
                         "H19 F1R1", "DRRRNA F3R3"))) %>%
  print()

## load data with best model
qpcrbatch <- prepare_batch("./data/exports/", equipment = "quant", skip = 21) 

results <- list()

# Loop through all targets in best.models data frame
for(i in 1:nrow(best.models)){
  
  results[[i]] <- qpcrbatch %>%
    filter(target == best.models[i,1]) %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>% # use the best model in each model_qpcr
    analyze_models() # analyze models for cpD2
  
}

#GO to this point if you need to rerun data


# combine all results and str split id variables
qpcrdat <- rbind_all(results) 
id.var <- str_split_fixed(qpcrdat$ID, "_", 4) 
colnames(id.var) <- c("subject", "timepoint", "leg", "target")  
qpcrdat <- cbind(id.var, qpcrdat[,-1])



## estimate efficiencies ##
efficiencies <- list()

# use the same loop to analyze efficiencies
for(i in  1:nrow(best.models)){
  
  efficiencies[[i]] <- qpcrbatch %>%
    filter(target == best.models[i,1]) %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>%
    analyze_efficiency(method = "cpD2", model = "linexp", offset = -3)
  
}



# combine results and use str split to extract id variables

efficiencies <- rbind_all(efficiencies) 
id.var <- str_split_fixed(efficiencies$ID, "_", 4) 
colnames(id.var) <- c("subject", "timepoint", "leg", "target")  
efficiencies <- cbind(id.var, efficiencies[,-1])

efficiencies %>%
  filter(eff > 1.5 & eff < 2.5)%>% # remove outliers from efficiency estimation
  separate(target, into = c("target", "cDNA"), sep = "_") %>%
  group_by(target)%>%
  summarise(efficiency = mean(eff, na.rm = TRUE),
            max.eff = max(eff, na.rm = TRUE),
            min.eff = min(eff, na.rm = TRUE),
            sd.eff = sd(eff, na.rm = TRUE))%>%
  ggplot(aes(target, efficiency)) + geom_point() + 
  coord_flip() 


effs <- efficiencies %>%
  filter(eff > 1.5 & eff < 2.5)%>% # remove outliers from efficiency estimation
  separate(target, into = c("target", "cdna"), sep = "_") %>%
  group_by(target)%>%
  summarise(eff = mean(eff, na.rm = TRUE))



## Extract information on replicates

replicates <- data.frame(str_split(qpcrdat$target, "_", simplify = TRUE))

colnames(replicates) <- c("target", "cdna")


qpcrdat <- cbind(replicates, qpcrdat[, -4])


qpcrdat

## Combine all qPCR parameters in the study to a data.frame containing all replicates
qpcrdat.replicates <- qpcrdat %>%

  inner_join(effs, by = "target") %>%
  dplyr::select(subject, leg, timepoint, target, cdna, cpD2, eff.y) %>%
  mutate(cq = cpD2,
         eff = eff.y) %>%
  ungroup() %>%
  data.frame()


#load sheet with FP and sets info

info <- read.csv("./data/oneThreeSetLeg.csv", sep = ";") %>%
  pivot_longer(names_to = "sets", 
               values_to = "leg", multiple:single) %>%
  print()

qpcrdat.replicates <- qpcrdat.replicates %>%
  inner_join(info) %>%
  print()
### removed this code and then it worked:
#separate(target, into = c("target", "cdna"), sep = "_") %>%



### Save data ####

saveRDS(qpcrdat.replicates, "./derivedData/qpcr.replicates_LINCs.Rds")






# Remove objects from environment
rm(list = ls())
