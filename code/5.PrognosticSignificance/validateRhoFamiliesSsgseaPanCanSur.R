
library(dplyr)
library(reshape2)
library(tibble)
library(survival)

######################
tcgaSsgseaRes <- readRDS(file = '/data/rhoDiffValida/rhoFam.GPL570Sur.rds')

##########
tumorSsgseaScores <- subset(tcgaSsgseaRes, !is.na(OS_status) & !is.na(OS_month)) %>% 
  dcast(SampleIDs~Groups, value.var = "Scores") 

panCanRhoEnrichScoreCliData <- tcgaSsgseaRes %>% select(SampleIDs, cancer_type, GEO_number, OS_status, OS_month) %>% 
  distinct(SampleIDs, .keep_all = TRUE) %>% inner_join(tumorSsgseaScores, by = join_by(SampleIDs)) %>% 
  dplyr::rename(patientId = SampleIDs)

# breast   colon  kidney    lung stomach 
# 934    1053      22    1087     618 
##########

##############Rho scores---OS
osCanTypes <- c('breast', 'colon', 'kidney', 'lung', 'stomach')
rhoFimilies <- c('DH_DOCKER_dom', 'RhoGAP_dom', 'Small_GTPase_Rho')


panCanEnrichScoreSurPvalue <- lapply(osCanTypes, function(cancer){
  
  
  canNrExpCliData <- subset(panCanRhoEnrichScoreCliData, cancer_type == cancer)
  
  surPvalue <- lapply(rhoFimilies, function(gene){
    
    
    naSum <- sum(is.na(canNrExpCliData[, gene]) | canNrExpCliData[, gene] == 0, na.rm = TRUE)
    
    if(naSum/nrow(canNrExpCliData) > 0.2){
      
      pValue <- setNames(c(NA, NA, NA), c('HR', 'coxPvalue', 'logrankPvalue'))
      
    }else{
      
      
      nrExpCliData <- canNrExpCliData[, c('OS_month', 'OS_status', gene)]
      colnames(nrExpCliData) <- c('time', 'event', 'geneExp')
      
      nrExpCliData <- subset(nrExpCliData, !is.na(geneExp)) #  & geneExp > 0
      nrExpCliData <- nrExpCliData %>% mutate(expGroup = geneExp > median(geneExp))
      
      
      univ.model <- coxph(Surv(time, event)~geneExp, data = nrExpCliData)
      univ.p <- summary(univ.model)$coefficients[1, 5]
      hr <- summary(univ.model)$coefficients[1, 2]
      
      
      lr.test <- survdiff(Surv(time, event)~expGroup, data = nrExpCliData)
      lr.p <- 1 - pchisq(lr.test$chisq, length(lr.test$n) - 1)
      
      pValue <- setNames(c(hr, univ.p, lr.p), c('HR', 'coxPvalue', 'logrankPvalue'))
      
    }
    
    
    return(pValue)
  })
  
  surPvalue <- as.data.frame(do.call(rbind, surPvalue)) %>% 
    mutate(disease = cancer, geneSymbol = rhoFimilies, coxQvalue = p.adjust(coxPvalue, method = 'BH'), 
           logrankQvalue = p.adjust(logrankPvalue, method = 'BH'))
  
  
  return(surPvalue)
})

panCanEnrichScoreSurPvalue <- bind_rows(panCanEnrichScoreSurPvalue)


######################

panCanEnrichScoreSurPvalue <- sapply(osCanTypes, function(cancer){
  
  
  canNrExpCliData <- subset(panCanRhoEnrichScoreCliData, cancer_type == cancer)
  
  canNrExpCliData <- split.data.frame(canNrExpCliData, f = ~ GEO_number)
  
  eachDataSurPvalue <- sapply(names(canNrExpCliData), function(GEO_number){
   eachGeoData <- canNrExpCliData[[GEO_number]]
   
   surPvalue <- sapply(rhoFimilies, function(gene){
     
     
     naSum <- sum(is.na(eachGeoData[, gene]) | eachGeoData[, gene] == 0, na.rm = TRUE)
     
     if(naSum/nrow(eachGeoData) > 0.2){
       
       pValue <- setNames(c(NA, NA, NA), c('HR', 'coxPvalue', 'logrankPvalue'))
       
     }else{
       
       
       nrExpCliData <- eachGeoData[, c('OS_month', 'OS_status', gene)]
       colnames(nrExpCliData) <- c('time', 'event', 'geneExp')
       
       nrExpCliData <- subset(nrExpCliData, !is.na(geneExp)) #  & geneExp > 0
       nrExpCliData <- nrExpCliData %>% mutate(expGroup = geneExp > median(geneExp))
       
       
       univ.model <- coxph(Surv(time, event)~geneExp, data = nrExpCliData)
       univ.p <- summary(univ.model)$coefficients[1, 5]
       hr <- summary(univ.model)$coefficients[1, 2]
       
       
       lr.test <- survdiff(Surv(time, event)~expGroup, data = nrExpCliData)
       lr.p <- 1 - pchisq(lr.test$chisq, length(lr.test$n) - 1)
       
       pValue <- setNames(c(hr, univ.p, lr.p), c('HR', 'coxPvalue', 'logrankPvalue'))
       
     }
     
     return(pValue)
   }, simplify = F) 
    
   surPvalue <- bind_rows(surPvalue, .id = 'Group')
   return(surPvalue)
   
  }, simplify = F)

  eachDataSurPvalue <- bind_rows(eachDataSurPvalue, .id = 'GEO_number') %>% 
    mutate(coxQvalue = p.adjust(coxPvalue, method = 'BH'), logrankQvalue = p.adjust(logrankPvalue, method = 'BH'))
  
  return(eachDataSurPvalue)
}, simplify = F)

panCanEnrichScoreSurPvalue <- bind_rows(panCanEnrichScoreSurPvalue, .id = 'disease')


# saveRDS(panCanEnrichScoreSurPvalue, file = '/data/panCanRhoEnrichScoreSurPvalue.rds')

