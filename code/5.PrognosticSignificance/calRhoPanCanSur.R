
library(tibble)
library(survival)

#########################
tcgaPanCanCliData <- readRDS('/data/tcgaPanCanCliData.rds')
tcgaExpData <- readRDS(file = '/data/tcgaFirehoseExpData.rds')
load(file = '/data/tcgaExpSams.RData')
# tumorExpSams, normalExpSams, pairedSamsPats
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

##############
rhoFamilies <- rhoFamilies %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID)) %>% 
  subset(NCBI.Gene.ID %in% rownames(tcgaExpData$ACC$TPM) & !duplicated(NCBI.Gene.ID))


disTpmExp <- lapply(tcgaExpData, function(expData){
  
  expData <- expData$TPM[rhoFamilies$NCBI.Gene.ID, ]
  rownames(expData) <- rhoFamilies$Approved.symbol
  
  return(expData)
})

disTpmExp <- cbind.data.frame(disTpmExp)
colnames(disTpmExp) <- stringr::str_split(colnames(disTpmExp), "\\.", simplify = TRUE)[, 2]

disTpmExp <- disTpmExp[, intersect(colnames(disTpmExp), tumorExpSams)]
colnames(disTpmExp) <- substr(colnames(disTpmExp), 1, 12)

##############
comSams <- intersect(colnames(disTpmExp), rownames(tcgaPanCanCliData))
disTpmExp <- disTpmExp[, comSams]
tcgaPanCanCliData <- tcgaPanCanCliData[comSams, ]

tcgaPanCanCliData <- tcgaPanCanCliData %>% rownames_to_column(var = 'patientId')
disTpmExp <- as.data.frame(t(disTpmExp)) %>% rownames_to_column(var = 'patientId')

panCanRhoExpCliData <- merge.data.frame(disTpmExp, tcgaPanCanCliData, by = 'patientId')

##############OS
osCanTypes <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
                'LIHC','LUAD','LUSC','MESO','OV','PAAD','PRAD','READ','SARC','SKCM','STAD','THCA','UCEC','UCS','UVM')


rhoGenes <- rhoFamilies$Approved.symbol


# univariate Cox and Log-rank test
panCanSurUniPvalue <- lapply(osCanTypes, function(cancer){
  
  
  canNrExpCliData <- subset(panCanRhoExpCliData, cancer_type == cancer & !is.na(os) & !is.na(os_time))
  
  surPvalue <- lapply(rhoGenes, function(gene){
    
    # gene <- 'HTR2B'
    
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
    mutate(disease = cancer, geneSymbol = rhoGenes, coxQvalue = p.adjust(coxPvalue, method = 'BH'), 
           logrankQvalue = p.adjust(logrankPvalue, method = 'BH'))
  

  return(surPvalue)
})

panCanSurUniPvalue <- do.call(rbind.data.frame, panCanSurUniPvalue)


# saveRDS(panCanSurUniPvalue, file = '/data/panCanRhoSurPvalue.rds')



# multivariate Cox
panCanSurMultiPvalue <- lapply(osCanTypes, function(cancer){
  
  # cancer = 'KIRC'
  canNrExpCliData <- subset(panCanRhoExpCliData, cancer_type == cancer & !is.na(os) & !is.na(os_time))
  
  surPvalue <- lapply(rhoGenes, function(gene){
    
    # gene <- 'RHOH'
    
    naSum <- sum(is.na(canNrExpCliData[, gene]) | canNrExpCliData[, gene] == 0, na.rm = TRUE)
    
    if(naSum/nrow(canNrExpCliData) > 0.2){
      
      pValue <- setNames(c(NA, NA), c('mulHR', 'mulcoxPvalue'))
      
    }else{
      
      
      nrExpCliData <- canNrExpCliData[, c('os_time', 'os', gene, 'age', 'gender', 'race', 'ajcc_stage', 'clinical_stage', 'histological_grade')]
      colnames(nrExpCliData) <- c('time', 'event', 'geneExp', 'age', 'gender', 'race', 'ajcc_stage', 'clinical_stage', 'histological_grade')
      
      nrExpCliData <- subset(nrExpCliData, !is.na(geneExp)) #  & geneExp > 0
      
      na_counts <- nrExpCliData %>% select(age, gender, race, ajcc_stage, clinical_stage, histological_grade) %>% summarise_all(~sum(is.na(.)))
      na_counts <- na_counts/nrow(nrExpCliData)
      vars <- names(na_counts)[na_counts < 0.3]
      n_distinct_counts <- nrExpCliData[, vars] %>% summarise(across(everything(), n_distinct))
      vars <- names(n_distinct_counts)[n_distinct_counts > 1]
      
      formula <- reformulate(termlabels = c(vars, 'geneExp'), response = "Surv(time, event)")
      
      nrExpCliData <- na.omit(nrExpCliData[, c('time', 'event', 'geneExp', vars)])
      multiv.model <- coxph(formula, data = nrExpCliData)
      
      multiv.p <- summary(multiv.model)$coefficients['geneExp', 5]
      hr <- summary(multiv.model)$coefficients['geneExp', 2]
      
      pValue <- setNames(c(hr, multiv.p), c('mulHR', 'mulcoxPvalue'))
    }
    
    
    return(pValue)
  })
  
  surPvalue <- as.data.frame(do.call(rbind, surPvalue)) %>%
    mutate(disease = cancer, geneSymbol = rhoGenes, mulcoxQvalue = p.adjust(mulcoxPvalue, method = 'BH'))
  
  print(cancer)
  return(surPvalue)
})

panCanSurMultiPvalue <- do.call(rbind.data.frame, panCanSurMultiPvalue)


panCanSurUniPvalue <- readRDS(file = '/data/panCanRhoSurPvalue.rds')

panCanSurPvalue <- merge(panCanSurUniPvalue, panCanSurMultiPvalue, by = c('disease', 'geneSymbol'))

# table(subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05)$mulcoxPvalue < 0.05)
# FALSE  TRUE 
# 173   418 
# > 418/(173+418)
# [1] 0.7072758

saveRDS(panCanSurPvalue, file = '/data/panCanRhoSurUniMulLogPvalue.rds')








