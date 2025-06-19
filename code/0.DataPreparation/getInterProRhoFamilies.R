
setwd('/data/InterPro/Domain')
listFiles <- list.files()[-1]

rhoFamilies <- sapply(listFiles, function(fileName){
  
  fimilyGenes <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  return(fimilyGenes)
}, simplify = FALSE)


rhoFamilies <- do.call(rbind, rhoFamilies)
rhoFamilies <- tibble::remove_rownames(rhoFamilies)


saveRDS(rhoFamilies, file = '/data/InterProRhoFamilies.rds')
