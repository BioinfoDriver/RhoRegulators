library('dplyr')
library('tibble')
library('Hmisc')
library('data.table')

###########
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

load(file = '/data/panCanMiRnaExpData.RData')
# panCanTurMiRnaExp, panCanPairdTurMiRnaExp, panCanPairdNormMiRnaExp


# geneInfo, panCanTurGeneExp, panCanPairdTurGeneExp, panCanPairdNormGeneExp
tcgaExpData <- readRDS(file = '/data/tcgaFirehoseExpData.rds')


load(file = '/data/tcgaExpSams.RData')
# tumorExpSams, normalExpSams, pairedSamsPats


panCanTurGeneExp <- lapply(tcgaExpData, function(expData){
  
  expData <- expData$TPM
  expData <- expData[, intersect(colnames(expData), tumorExpSams)]
  
  return(expData)
})
panCanTurGeneExp <- bind_cols(panCanTurGeneExp)


###########
rhoFamilie <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilie <- rhoFamilie %>% distinct(Approved.symbol, .keep_all = TRUE) %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID))


comSmas <- intersect(colnames(panCanTurMiRnaExp), colnames(panCanTurGeneExp))


turMiRnaExp <- panCanTurMiRnaExp[, comSmas]
turGeneExp <- panCanTurGeneExp[rhoFamilie$NCBI.Gene.ID, comSmas]
rownames(turGeneExp) <- rhoFamilie$Approved.symbol


tcgaPanCanSamples <- subset(tcgaPanCanSamples, SAMPLE_BARCODE %in% comSmas)
diseaseSams <- split.data.frame(tcgaPanCanSamples, f = ~ DISEASE)



###########

setwd('/data/TargetScan')

targetScan <- read.table(file = 'Predicted_Targets_Context_Scores.default_predictions.txt', header = T, sep = '\t', stringsAsFactors = F)
targetScan <- subset(targetScan, Gene.Tax.ID == 9606)
targetScan <- targetScan[!duplicated(targetScan[, c('Gene.Symbol', 'miRNA')]), ]


###########

miRnaTarByDisease <- lapply(diseaseSams, function(diseaseSam){
  
  # diseaseSam <- diseaseSams[[1]]
  # Matrix of Correlations and P-values
  mirExp <- t(turMiRnaExp[, diseaseSam$SAMPLE_BARCODE])
  geneExp <- t(turGeneExp[, diseaseSam$SAMPLE_BARCODE])
  
  # filtering-new add
  mirExp <- mirExp[, apply(mirExp, 2, function(x) (sum(is.na(x)) + sum(x == 0, na.rm = T))/nrow(mirExp) < 0.1)]
  geneExp <- geneExp[, apply(geneExp, 2, function(x) (sum(is.na(x)) + sum(x == 0, na.rm = T))/nrow(geneExp) < 0.1)]
  
  
  rcorrRes <- rcorr(mirExp, geneExp, type = "spearman")
  
  corRho <- rcorrRes$r[colnames(mirExp), colnames(geneExp)]
  corPvalue <- rcorrRes$P[colnames(mirExp), colnames(geneExp)]
  
  corFdr <- matrix(p.adjust(corPvalue, method = 'fdr'), nrow = nrow(corPvalue), byrow = FALSE)
  dimnames(corFdr) <- dimnames(corPvalue)
  
  
  sifCor <- corRho < -0.25 & corFdr < 1.0e-05
  sifCor <- sifCor %>% melt() %>% rename(miRNA = Var1, TargetGene = Var2, sig = value) %>% subset(sig == TRUE)
  
  
  # miRnaTar <- sifCor %>% inner_join(hasMiRnaTar, by = c('miRNA', 'TargetGene'))
  miRnaTar <- sifCor %>% inner_join(targetScan, by = join_by(miRNA == miRNA, TargetGene == Gene.Symbol))
  
  
  if(nrow(miRnaTar) > 0){
    
    miRnaTar$DISEASE <- diseaseSam$DISEASE[1]
    return(miRnaTar)
    
  }else{
    return(NULL)
    
  }
})


miRnaTarByDisease <- do.call(rbind.data.frame, miRnaTarByDisease) %>% remove_rownames()

# saveRDS(miRnaTarByDisease, file = '/data/miRnaTarByDisease.rds')

saveRDS(miRnaTarByDisease, file = '/data/miRnaTarByDiseaseTargetScan.rds')

