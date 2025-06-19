
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

setwd('/data/Bailey_Cell_2018_CancerDriverGenes/MutSig2CV/')

list.files <- list.files()
list.files <- list.files[!grepl("^(README|UCEChyb)\\.txt$", list.files)]


# 'TMEM133' <- 'ARHGAP42'
# 'SRGAP2' <- 'SRGAP2', missing
# 'ARHGAP40' <- 'ARHGAP40', missing
# 'HMHA1' <- 'ARHGAP45'
# 'ARHGAP23' <- 'ARHGAP23', missing
# 'ARHGEF33' <- 'ARHGEF33', missing



mutSigGeneQvalue <- lapply(list.files, function(file.name){
  
  mutSigRes <- read.csv(file = file.name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  mutSigRes$gene[mutSigRes$gene == 'TMEM133'] <- 'ARHGAP42'
  mutSigRes$gene[mutSigRes$gene == 'HMHA1'] <- 'ARHGAP45'
  
  mutSigRes <- subset(mutSigRes, gene %in% rhoFamilies$Approved.symbol)
  
  mutSigRes <- mutSigRes[, c('gene', 'qvalue')]
  colnames(mutSigRes) <- c('gene', paste0(gsub('.txt', '', file.name), c('qvalue')))
  
  
  return(mutSigRes)
})

mutSigGeneQvalue <- Reduce(x = mutSigGeneQvalue, function(x, y) merge(x, y, by = 'gene'))


mutSigGenePvalue <- lapply(list.files, function(file.name){
  
  mutSigRes <- read.csv(file = file.name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  mutSigRes$gene[mutSigRes$gene == 'TMEM133'] <- 'ARHGAP42'
  mutSigRes$gene[mutSigRes$gene == 'HMHA1'] <- 'ARHGAP45'
  
  mutSigRes <- subset(mutSigRes, gene %in% rhoFamilies$Approved.symbol)
  
  mutSigRes <- mutSigRes[, c('gene', 'pvalue')]
  colnames(mutSigRes) <- c('gene', paste0(gsub('.txt', '', file.name), c('pvalue')))
  
  
  return(mutSigRes)
})

mutSigGenePvalue <- Reduce(x = mutSigGenePvalue, function(x, y) merge(x, y, by = 'gene'))

save(mutSigGenePvalue, mutSigGeneQvalue, file = '/data/rhoFamiliesMutSig2CV.RData')

