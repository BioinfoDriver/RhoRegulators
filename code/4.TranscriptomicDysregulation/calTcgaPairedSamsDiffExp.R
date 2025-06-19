
library('limma')
library('dplyr')
library('tibble')
library("DESeq2")

#########################
tcgaExpData <- readRDS(file = '/data/tcgaFirehoseExpData.rds')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

load(file = '/data/tcgaExpSams.RData')
# tumorExpSams, normalExpSams, pairedSamsPats


# paired samples
disSamPairs <- subset.data.frame(tcgaPanCanSamples, PATIENT_BARCODE %in% pairedSamsPats) 

tumorSams <- disSamPairs %>% mutate(samBarcode = SAMPLE_BARCODE, samType = 'Tumor')
normalSams <- disSamPairs %>% mutate(samBarcode = paste0(PATIENT_BARCODE, '-11'), samType = 'Normal')
pairdSams <- rbind.data.frame(tumorSams, normalSams) %>% split.data.frame(f = ~ DISEASE)


# differently expressed genes
# limma
pairedSamsLimmaDiffExp <- lapply(names(pairdSams), function(disease){
  
  
  diseaseSams <- pairdSams[[disease]]
  diseaseExpData <- tcgaExpData[[disease]][['TPM']]
  
  
  diseaseSams$PATIENT_BARCODE <- factor(diseaseSams$PATIENT_BARCODE)
  diseaseSams$samType <- factor(diseaseSams$samType, levels = c('Normal', 'Tumor'))
  
  design <- model.matrix(~PATIENT_BARCODE+samType, data = diseaseSams)
  
  
  eset <- cbind.data.frame(diseaseExpData[, subset(diseaseSams, samType == 'Tumor')$samBarcode], 
                           diseaseExpData[, subset(diseaseSams, samType == 'Normal')$samBarcode])
  
  
  # filter
  filterIndex <- apply(eset, 1, function(x) sum(x == 0, na.rm = TRUE) + sum(is.na(x)))/ncol(eset) > 0.2
  print(table(filterIndex))
  
  eset <- eset[!filterIndex, ]
  
  
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)
  difRes <- topTable(fit, coef="samTypeTumor", adjust.method = "BH", n = Inf)
  
  difRes <- difRes %>% mutate(DISEASE = disease) %>% rownames_to_column(var = "geneID")
  
  return(difRes)
})

pairedSamsLimmaDiffExp <- do.call(rbind.data.frame, pairedSamsLimmaDiffExp)


# Deseq2
RunDESeq2 <- function(countMatrix, pData){
  
  library(DESeq2);
  dds <- DESeqDataSetFromMatrix(countData=countMatrix, colData=pData, design = ~ PATIENT_BARCODE + samType);
  
  dds <- DESeq(dds, parallel=T);
  dds <- replaceOutliersWithTrimmedMean(dds);
  
  res <- results(dds, alpha = 0.1, cooksCutoff=FALSE);
  res <- res[order(res$padj), ];
  
  return(res);
}

pairedSamsDeseq2DiffExp <- lapply(names(pairdSams), function(disease){
  
  
  diseaseSams <- pairdSams[[disease]]
  diseaseExpData <- tcgaExpData[[disease]][['COUNT']]
  
  
  diseaseSams$PATIENT_BARCODE <- factor(diseaseSams$PATIENT_BARCODE)
  diseaseSams$samType <- factor(diseaseSams$samType, levels = c('Normal', 'Tumor'))
  

  expMatrix <- cbind.data.frame(diseaseExpData[, subset(diseaseSams, samType == 'Tumor')$samBarcode], 
                           diseaseExpData[, subset(diseaseSams, samType == 'Normal')$samBarcode])
  
  
  # filter
  filterIndex <- apply(expMatrix, 1, function(x) sum(x == 0, na.rm = TRUE) + sum(is.na(x)))/ncol(expMatrix) > 0.2
  print(table(filterIndex))
  
  expMatrix <- expMatrix[!filterIndex, ]
  
  
  difRes <- RunDESeq2(countMatrix = expMatrix, pData = diseaseSams)
  
  difRes <- as.data.frame(difRes)
  difRes <- difRes %>% mutate(DISEASE = disease) %>% rownames_to_column(var = "geneID")
  
  return(difRes)
})

pairedSamsDeseq2DiffExp <- do.call(rbind.data.frame, pairedSamsDeseq2DiffExp)


########
geneInfo <- read.csv(file= '/data/Homo_sapiens.gene_info', sep = '\t', header = TRUE, stringsAsFactors = FALSE)


pairedSamsDeseq2DiffExp <- pairedSamsDeseq2DiffExp %>% mutate(geneID = as.integer(geneID)) %>% 
  left_join(geneInfo[, c('GeneID', 'Symbol')], by = join_by(geneID == GeneID))


pairedSamsLimmaDiffExp <- pairedSamsLimmaDiffExp %>% mutate(geneID = as.integer(geneID)) %>% 
  left_join(geneInfo[, c('GeneID', 'Symbol')], by = join_by(geneID == GeneID))

save(pairedSamsDeseq2DiffExp, pairedSamsLimmaDiffExp, file = '/data/tcgaPairedSamsDiffExp.RData')





