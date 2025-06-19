library('dplyr')
library('tibble')

###########################
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
geneInfo <- read.csv(file= '/data/gene_with_protein_product.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

###########################
setwd('/data/FireBro')
fileNames <- dir(recursive = TRUE, pattern = 'all_thresholded.by_genes.txt$')


panCanGeneCNV <- lapply(fileNames, function(fileName){
  
  geneCNV <- read.csv(file = fileName, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

  geneCNV <- geneCNV %>% select(-Locus.ID, -Cytoband) %>% column_to_rownames(var = 'Gene.Symbol') %>% 
    rename_all(~substr(., 1, 15)) %>% rename_all(~gsub("\\.", "-", .)) %>% mutate_all(~replace(., . %in% c(1, -1), 0))
  
  return(geneCNV)
})


names(panCanGeneCNV) <- stringr::str_split(fileNames, "_", simplify = TRUE)[, 1]
panCanGeneCNV <- panCanGeneCNV[!names(panCanGeneCNV) %in% c('COADREAD', 'GBMLGG', 'KIPAN', 'STES')]


mergeF <- function(x, y){
  
  z <- merge(x, y, by = 'row.names') %>% column_to_rownames(var = 'Row.names')
  return(z)
}

panCanGeneCNV <- Reduce(function(x, y) mergeF(x, y), panCanGeneCNV)

##############
panCanGeneCNV <- panCanGeneCNV[, intersect(tcgaPanCanSamples$SAMPLE_BARCODE, colnames(panCanGeneCNV))]

# setdiff(rhoFamilies$Approved.symbol, rownames(panCanGeneCNV))
# "ARHGAP45"

rownames(panCanGeneCNV)[rownames(panCanGeneCNV) == 'HMHA1'] <- 'ARHGAP45'


##############

rhoCopyNumData <- panCanGeneCNV[unique(rhoFamilies$Approved.symbol), ]
proteinGeneCopyNumData <- panCanGeneCNV[intersect(rownames(panCanGeneCNV), geneInfo$symbol), ]


save(rhoCopyNumData, proteinGeneCopyNumData, panCanGeneCNV, file = '/data/panCanCopyNumData.RData')

