

library('dplyr')
library('tibble')

##############
load(file = '/data/panCanMethyData.RData')
# panCanTurMethy, panCanPairdTurMethy, panCanPairdNormMethy

anno450k <- readRDS(file = '/data/methAnno450k.rds')
rhoFamilie <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilie <- rhoFamilie %>% distinct(Approved.symbol, .keep_all = TRUE)


nrAnno450k <- subset(anno450k, RefGene_Name %in% rhoFamilie$Approved.symbol)
nrpromoterAnno <- subset(nrAnno450k, RefGene_Group %in% c("5'UTR", "TSS1500", "TSS200"))


panCanPairdTurMethy <- rownames_to_column(panCanPairdTurMethy, var = 'Name')
panCanPairdNormMethy <- rownames_to_column(panCanPairdNormMethy, var = 'Name')


panCanPairdTurMethy <- panCanPairdTurMethy %>% inner_join(nrpromoterAnno[, c('Name', 'RefGene_Name')], by = 'Name')
panCanPairdNormMethy <- panCanPairdNormMethy %>% inner_join(nrpromoterAnno[, c('Name', 'RefGene_Name')], by = 'Name')


##############
panCanPairdTurMethy <- panCanPairdTurMethy %>% group_by(RefGene_Name) %>% select(-Name, -RefGene_Name) %>%
  summarize(across(everything(), mean, na.rm = TRUE)) %>% column_to_rownames(var = 'RefGene_Name')   


panCanPairdNormMethy <- panCanPairdNormMethy %>% group_by(RefGene_Name) %>% select(-Name, -RefGene_Name) %>%
  summarize(across(everything(), mean, na.rm = TRUE)) %>% column_to_rownames(var = 'RefGene_Name')   


###############
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')


pairdSams <- substr(colnames(panCanPairdTurMethy), 1, 12)
pairdSams <- tcgaPanCanSamples %>% subset(PATIENT_BARCODE %in% pairdSams)


samSta <- pairdSams %>% group_by(DISEASE) %>% count() %>% subset(n >=10)
pairdSams <- subset(pairdSams, DISEASE %in% samSta$DISEASE)


tumorSams <- pairdSams %>% mutate(samBarcode = SAMPLE_BARCODE, samType = 'Tumor')
normalSams <- pairdSams %>% mutate(samBarcode = paste0(PATIENT_BARCODE, '-11'), samType = 'Normal')


pairdSams <- rbind.data.frame(tumorSams, normalSams)
pairdSams <- split.data.frame(pairdSams, f = pairdSams$DISEASE)



# "BLCA" "BRCA" "COAD" "ESCA" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA" "UCEC"
# differently methylation
diseases <- names(pairdSams)

diffMethy <- lapply(diseases, function(disease){
  
  
  diseaseSams <- pairdSams[[disease]]
  
  tM <- panCanPairdTurMethy[, subset(diseaseSams, samType == 'Tumor')$samBarcode]
  nM <- panCanPairdNormMethy[, subset(diseaseSams, samType == 'Normal')$samBarcode]
  
  
  diffM <- lapply(rownames(tM), function(geneName){
    
    tnv <- data.frame(tMv = as.numeric(tM[geneName, ]), nMv = as.numeric(nM[geneName, ]))
    tnv <- subset(tnv, !is.na(tMv) & !is.na(nMv))
    
    pValue <- wilcox.test(Pair(tMv, nMv) ~ 1, data = tnv)$p.value
    logFC <- log2(mean(tnv$tMv)/mean(tnv$nMv))
    
    return(data.frame(Approved.symbol = geneName, logFC, pValue, Disease = disease))
  })
  
  diffM <- do.call(rbind, diffM)
  diffM$adjPVal <- p.adjust(diffM$pValue, method = 'fdr')
  
  return(diffM)
})

diffMethy <- do.call(rbind.data.frame, diffMethy)


saveRDS(diffMethy, file = '/data/panCanDiffMethy.rds')

