
library(dplyr)
library(reshape2)
library(tibble)
library(survival)

######################
tcgaSsgseaRes <- readRDS(file = '/data/rhoDiffValida/rhoFam.GPL570Sur.rds')

tumorSsgseaScores <- subset(tcgaSsgseaRes, !is.na(OS_status) & !is.na(OS_month)) %>% 
  dcast(SampleIDs~Groups, value.var = "Scores")

tcgaPanCanCliData <- tcgaSsgseaRes %>% select(SampleIDs, OS_status, OS_month, cancer_type) %>% distinct(SampleIDs, .keep_all = TRUE)

panCanRhoEnrichScoreCliData <- merge.data.frame(tumorSsgseaScores, tcgaPanCanCliData, by = 'SampleIDs')
##########

##############Rho scores---OS
osCanTypes <- c('breast', 'colon', 'lung', 'stomach')

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

panCanEnrichScoreSurPvalue <- do.call(rbind.data.frame, panCanEnrichScoreSurPvalue)

#           HR    coxPvalue logrankPvalue disease       geneSymbol    coxQvalue logrankQvalue
# 1  0.9436687 0.0921209303   0.132875630  breast    DH_DOCKER_dom 0.2763627908    0.19931345
# 2  0.9968394 0.9196956141   0.407382410  breast       RhoGAP_dom 0.9598972975    0.40738241
# 3  1.0022449 0.9598972975   0.082143511  breast Small_GTPase_Rho 0.9598972975    0.19931345
# 4  1.0353038 0.0821790456   0.054196405   colon    DH_DOCKER_dom 0.1232685685    0.16258922
# 5  1.0095133 0.7350862599   0.960526952   colon       RhoGAP_dom 0.7350862599    0.96052695
# 6  1.0487826 0.0499271285   0.132837162   colon Small_GTPase_Rho 0.1232685685    0.19925574
# 7  0.9376190 0.0001496975   0.001835037    lung    DH_DOCKER_dom 0.0004249113    0.00550511
# 8  0.9360331 0.0002832742   0.007541435    lung       RhoGAP_dom 0.0004249113    0.01131215
# 9  0.9954040 0.8327835520   0.980955523    lung Small_GTPase_Rho 0.8327835520    0.98095552
# 10 1.0431146 0.0463058101   0.185249259 stomach    DH_DOCKER_dom 0.0497609256    0.18524926
# 11 1.0492467 0.0389062065   0.006582627 stomach       RhoGAP_dom 0.0497609256    0.01974788
# 12 1.0568765 0.0497609256   0.041056698 stomach Small_GTPase_Rho 0.0497609256    0.06158505
# saveRDS(panCanEnrichScoreSurPvalue, file = '/data/panCanRhoEnrichScoreSurPvalue.rds')


##############Rho EnrichGroup---OS
tumorRhoEnrichGroup <- subset(tcgaSsgseaRes, !is.na(OS_status) & !is.na(OS_month)) %>% 
  mutate(EnrichGroup = ifelse(FDR >= 0.25, 'Neutral', ifelse(Scores > 0, 'Up', 'Down'))) %>% subset(EnrichGroup != 'Neutral') %>% 
  dcast(SampleIDs~Groups, value.var = "EnrichGroup")


panCanRhoEnrichGroupCliData <- merge.data.frame(tumorRhoEnrichGroup, tcgaPanCanCliData, by = 'SampleIDs')


panCanEnrichGroupSurPvalue <- lapply(osCanTypes, function(cancer){
  
  
  canNrExpCliData <- subset(panCanRhoEnrichGroupCliData, cancer_type == cancer)
  
  surPvalue <- lapply(rhoFimilies, function(gene){
    
    
    # naSum <- sum(is.na(canNrExpCliData[, gene]) | canNrExpCliData[, gene] == 0, na.rm = TRUE)
    # 
    # if(naSum/nrow(canNrExpCliData) > 0.2){
    #   
    #   pValue <- setNames(c(NA, NA, NA), c('HR', 'coxPvalue', 'logrankPvalue'))
    #   
    # }else{
    #   
    
    nrExpCliData <- canNrExpCliData[, c('OS_month', 'OS_status', gene)]
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

#           HR    coxPvalue logrankPvalue disease       geneSymbol    coxQvalue logrankQvalue
# 1  0.8378659 0.4079979938  0.4069099872  breast    DH_DOCKER_dom 0.6931851379  0.6908352100
# 2  0.9837328 0.9423431776  0.9416413364  breast       RhoGAP_dom 0.9423431776  0.9416413364
# 3  0.8137811 0.4621234253  0.4605568066  breast Small_GTPase_Rho 0.6931851379  0.6908352100
# 4  1.2759203 0.0861096957  0.0863824994   colon    DH_DOCKER_dom 0.1355067639  0.1338229233
# 5  1.1197669 0.5167739595  0.5178230165   colon       RhoGAP_dom 0.5167739595  0.5178230165
# 6  1.2877669 0.0903378426  0.0892152822   colon Small_GTPase_Rho 0.1355067639  0.1338229233
# 7  0.6482621 0.0002701157  0.0002419654    lung    DH_DOCKER_dom 0.0008103471  0.0007258963
# 8  0.6720862 0.0012904055  0.0011887073    lung       RhoGAP_dom 0.0019356082  0.0017830609
# 9  0.9212660 0.5120721301  0.5101332606    lung Small_GTPase_Rho 0.5120721301  0.5101332606
# 10 1.4766652 0.0126985147  0.0121126423 stomach    DH_DOCKER_dom 0.0380955442  0.0363379270
# 11 1.4219488 0.0300121125  0.0290833010 stomach       RhoGAP_dom 0.0450181688  0.0436249515
# 12 1.2990805 0.1242103028  0.1232045589 stomach Small_GTPase_Rho 0.1242103028  0.1232045589
# saveRDS(panCanEnrichGroupSurPvalue, file = '/data/panCanRhoEnrichGroupSurPvalue.rds')



