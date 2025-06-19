
library('dplyr')

tcgaExpData <- readRDS(file = '/data/tcgaFirehoseExpData.rds')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')


tcgaPanCanExpSamples <- unlist(lapply(tcgaExpData, function(Exp){
  
  expSams <- colnames(Exp$TPM)
  
  return(expSams)
  
}), use.names = FALSE)


tumorExpSams <- intersect(tcgaPanCanSamples$SAMPLE_BARCODE, tcgaPanCanExpSamples)
normalExpSams <- intersect(paste0(tcgaPanCanSamples$PATIENT_BARCODE, '-11'), tcgaPanCanExpSamples)

pairedSamsPats <- intersect(substr(tumorExpSams, 1, 12), substr(normalExpSams, 1, 12))

# n >10
pairedSamsDis <- subset(tcgaPanCanSamples, PATIENT_BARCODE %in% pairedSamsPats) %>% count(DISEASE) %>% subset(n>=10)

pairedSamsPats <- subset(tcgaPanCanSamples, PATIENT_BARCODE %in% pairedSamsPats & DISEASE %in% pairedSamsDis$DISEASE)$PATIENT_BARCODE


save(tumorExpSams, normalExpSams, pairedSamsPats, file = '/data/tcgaExpSams.RData')


