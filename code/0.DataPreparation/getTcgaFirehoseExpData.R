
setwd('/data/Firehose/RNASEQ/RSEMExp')
files <- list.files(path = ".", pattern = 'Level_3__RSEM_genes__data.data', recursive = TRUE)

tcgaExpData <- lapply(files, function(fileNames){
  
  expData<- read.csv(file = fileNames, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  
  expCount <- expData[-1, seq(1, ncol(expData), 3)]
  expTpm <- expData[-1, seq(2, ncol(expData), 3)]
  
  colnames(expCount) <- substr(gsub('\\.', '-', colnames(expCount)), 1, 15)
  colnames(expTpm) <- substr(gsub('\\.', '-', colnames(expTpm)), 1, 15)
  
  # Raw count and scaled estimate (TPM) (PMID:34430923)
  # We extracted the gene expression data from "illuminahiseq_rnaseqv2-RSEM_genes" files.
  # From these data, we used "raw_count" values as counts and 
  # we calculated transcripts per million (TPM) from "scaled_estimate" values multiplied by 1,000,000. 
  
  library(dplyr)
  expCount <- expCount %>% mutate(across(where(is.character), as.numeric))
  expCount <- round(expCount)

  expTpm <- expTpm %>% mutate(across(where(is.character), as.numeric))
  expTpm <- expTpm * 10^6
  expTpm <- log2(expTpm + 1.0)
  
  rowNames <- do.call(rbind, strsplit(rownames(expCount), split = '\\|')) 

  rownames(expCount) <- rowNames[, 2]
  rownames(expTpm) <- rowNames[, 2]
  
  
  expZscore <- as.data.frame(t(scale(t(expTpm), center = TRUE, scale = TRUE)))
  
  return(list(COUNT = expCount, TPM = expTpm, Zscore = expZscore))
})

names(tcgaExpData) <- gsub('org_', '', do.call(rbind, strsplit(files, split = '\\.'))[, 3])

saveRDS(tcgaExpData, file = '/data/tcgaFirehoseExpData.rds')


