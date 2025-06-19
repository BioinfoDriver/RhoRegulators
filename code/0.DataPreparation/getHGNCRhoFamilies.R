

setwd('/data/HGNC')
listFiles <- list.files()[-1]


rhoFamilies <- sapply(listFiles, function(fileName){
  
  fimilyGenes <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  fimilyGenes$classesOfRho <- gsub(pattern = '.txt', replacement = '', x = fileName)
  
  return(fimilyGenes)
}, simplify = FALSE)


rhoFamilies <- do.call(rbind, rhoFamilies)


# rownames(rhoFamilies) <- rhoFamilies$Approved.symbol
# rhoFamilies <- subset(rhoFamilies, Locus.type != 'pseudogene')

rhoFamilies <- tibble::remove_rownames(rhoFamilies)

saveRDS(rhoFamilies, file = '/data/HGNCRhoFamilies.rds')

