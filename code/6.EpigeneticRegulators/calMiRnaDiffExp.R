
library('limma')
library('dplyr')
library('tibble')

###############

load(file = '/data/panCanMiRnaExpData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')



pairdSams <- substr(colnames(panCanPairdTurMiRnaExp), 1, 12)
pairdSams <- tcgaPanCanSamples %>% subset(PATIENT_BARCODE %in% pairdSams)


samSta <- pairdSams %>% group_by(DISEASE) %>% count() %>% subset(n >=10)
pairdSams <- subset(pairdSams, DISEASE %in% samSta$DISEASE)


tumorSams <- pairdSams %>% mutate(samBarcode = SAMPLE_BARCODE, samType = 'Tumor')
normalSams <- pairdSams %>% mutate(samBarcode = paste0(PATIENT_BARCODE, '-11'), samType = 'Normal')


pairdSams <- rbind.data.frame(tumorSams, normalSams)
pairdSams <- split.data.frame(pairdSams, f = pairdSams$DISEASE)


# differently expressed miRNA
diseases <- names(pairdSams)

mirDiffExp <- lapply(diseases, function(disease){
  
  
  diseaseSams <- pairdSams[[disease]]
  
  tM <- panCanPairdTurMiRnaExp[, subset(diseaseSams, samType == 'Tumor')$samBarcode]
  nM <- panCanPairdNormMiRnaExp[, subset(diseaseSams, samType == 'Normal')$samBarcode]
  
  
  diffM <- lapply(rownames(tM), function(geneName){
    
    tnv <- data.frame(tMv = as.numeric(tM[geneName, ]), nMv = as.numeric(nM[geneName, ]))
    tnv <- subset(tnv, !is.na(tMv) & !is.na(nMv))
    
    
    if(nrow(tnv)/ncol(tM) > 0.2 & nrow(tnv) >=10){
      pValue <- wilcox.test(Pair(tMv, nMv) ~ 1, data = tnv)$p.value
      logFC <- log2(mean(tnv$tMv)/mean(tnv$nMv))
      
      return(data.frame(Approved.symbol = geneName, logFC, pValue, Disease = disease, numPat = nrow(tnv)))      
      
    }else{

      return(data.frame(Approved.symbol = geneName, logFC = NA, pValue = NA, Disease = disease, numPat = NA))     
    }
  })
  
  diffM <- do.call(rbind, diffM)
  diffM <- subset(diffM, !is.na(numPat))
  diffM$adjPVal <- p.adjust(diffM$pValue, method = 'fdr')
  
  return(diffM)
})

mirDiffExp <- do.call(rbind.data.frame, mirDiffExp)



saveRDS(mirDiffExp, file = '/data/panCanMiRnaDiffExp.rds')


