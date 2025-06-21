
library(dplyr)
library(reshape2)
library(tibble)
library(survival)

######################
tcgaSsgseaRes <- readRDS(file = '/data/tcgaRhoFamiliesSsgsea2.0Res.rds')
load(file = '/data/tcgaExpSams.RData')
tcgaPanCanCliData <- readRDS('/data/tcgaPanCanCliData.rds')

##########
tumorSsgseaScores <- subset(tcgaSsgseaRes, SampleIDs %in% tumorExpSams) %>% 
  dcast(PATIENT_BARCODE~Groups, value.var = "Scores") %>% rename(patientId = PATIENT_BARCODE)

tcgaPanCanCliData <- tcgaPanCanCliData %>% rownames_to_column(var = 'patientId')
panCanRhoEnrichScoreCliData <- merge.data.frame(tumorSsgseaScores, tcgaPanCanCliData, by = 'patientId')

##########

##############Rho scores---OS
osCanTypes <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
                'LIHC','LUAD','LUSC','MESO','OV','PAAD','PRAD','READ','SARC','SKCM','STAD','THCA','UCEC','UCS','UVM')

rhoFimilies <- c('DH_DOCKER_dom', 'RhoGAP_dom', 'Small_GTPase_Rho')


panCanEnrichScoreSurPvalue <- lapply(osCanTypes, function(cancer){
  
  
  canNrExpCliData <- subset(panCanRhoEnrichScoreCliData, cancer_type == cancer & !is.na(os) & !is.na(os_time))
  
  surPvalue <- lapply(rhoFimilies, function(gene){
    
    
    naSum <- sum(is.na(canNrExpCliData[, gene]) | canNrExpCliData[, gene] == 0, na.rm = TRUE)
    
    if(naSum/nrow(canNrExpCliData) > 0.2){
      
      pValue <- setNames(c(NA, NA, NA), c('HR', 'coxPvalue', 'logrankPvalue'))
      
    }else{
      
      
      nrExpCliData <- canNrExpCliData[, c('os_time', 'os', gene)]
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

panCanEnrichScoreSurPvalue <- do.call(rbind.data.frame, panCanEnrichScoreSurPvalue)


saveRDS(panCanEnrichScoreSurPvalue, file = '/data/panCanRhoEnrichScoreSurPvalue.rds')


##############Rho EnrichGroup---OS

tumorRhoEnrichGroup <- subset(tcgaSsgseaRes, SampleIDs %in% tumorExpSams) %>% 
  mutate(EnrichGroup = ifelse(FDR >= 0.25, 'Neutral', ifelse(Scores > 0, 'Up', 'Down'))) %>% subset(EnrichGroup != 'Neutral') %>% 
  dcast(PATIENT_BARCODE~Groups, value.var = "EnrichGroup") %>% rename(patientId = PATIENT_BARCODE)


panCanRhoEnrichGroupCliData <- merge.data.frame(tumorRhoEnrichGroup, tcgaPanCanCliData, by = 'patientId')


panCanEnrichGroupSurPvalue <- lapply(osCanTypes, function(cancer){
  
  
  canNrExpCliData <- subset(panCanRhoEnrichGroupCliData, cancer_type == cancer & !is.na(os) & !is.na(os_time))
  
  surPvalue <- lapply(rhoFimilies, function(gene){
    
    
    # naSum <- sum(is.na(canNrExpCliData[, gene]) | canNrExpCliData[, gene] == 0, na.rm = TRUE)
    # 
    # if(naSum/nrow(canNrExpCliData) > 0.2){
    #   
    #   pValue <- setNames(c(NA, NA, NA), c('HR', 'coxPvalue', 'logrankPvalue'))
    #   
    # }else{
    #   
      
      nrExpCliData <- canNrExpCliData[, c('os_time', 'os', gene)]
      colnames(nrExpCliData) <- c('time', 'event', 'geneExp')
      
      nrExpCliData <- subset(nrExpCliData, !is.na(geneExp)) #  & geneExp > 0

     
      univ.model <- coxph(Surv(time, event)~geneExp, data = nrExpCliData)
      univ.p <- summary(univ.model)$coefficients[1, 5]
      hr <- summary(univ.model)$coefficients[1, 2]
      
      
      lr.test <- survdiff(Surv(time, event)~geneExp, data = nrExpCliData)
      lr.p <- 1 - pchisq(lr.test$chisq, length(lr.test$n) - 1)
      
      pValue <- setNames(c(hr, univ.p, lr.p), c('HR', 'coxPvalue', 'logrankPvalue'))
      
    # }
    
    
    return(pValue)
  })
  
  surPvalue <- as.data.frame(do.call(rbind, surPvalue)) %>% 
    mutate(disease = cancer, geneSymbol = rhoFimilies, coxQvalue = p.adjust(coxPvalue, method = 'BH'), 
           logrankQvalue = p.adjust(logrankPvalue, method = 'BH'))
  
  
  return(surPvalue)
})

panCanEnrichGroupSurPvalue <- do.call(rbind.data.frame, panCanEnrichGroupSurPvalue)


saveRDS(panCanEnrichGroupSurPvalue, file = '/data/panCanRhoEnrichGroupSurPvalue.rds')



