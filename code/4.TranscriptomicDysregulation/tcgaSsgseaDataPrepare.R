library('dplyr')
library('cmapR')
library('tibble')

#############GCT
tcgaExpData <- readRDS(file = '/data/tcgaFirehoseExpData.rds')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')


geneInfo <- read.csv(file= '/data/Homo_sapiens.gene_info', sep = '\t', 
                     header = TRUE, stringsAsFactors = FALSE)
geneInfo <- subset(geneInfo, type_of_gene  == 'protein-coding') %>% mutate(GeneID = as.character(GeneID))


setwd('/data/ssGSEAdata')
# ssGSEAPre <- lapply(split.data.frame(tcgaPanCanSamples, f = ~ DISEASE), function(cancerType){
#   
#   
#   expZscore <- tcgaExpData[[cancerType$DISEASE[1]]]$Zscore
#   expZscore <- expZscore[, !is.na(match(substr(colnames(expZscore), 1, 12), cancerType$PATIENT_BARCODE))]
#   
#   
#   indices <- apply(expZscore, 1, function(x) all(is.na(x)))
#   expZscore <- expZscore[!indices, ]
#   
#   expZscore <- tibble::rownames_to_column(expZscore, var = "GeneID")
#   
#   expZscore <- geneInfo[, c('GeneID', 'Symbol', 'description')] %>% 
#     inner_join(expZscore, join_by(GeneID)) %>% mutate(GeneID = NULL)
#  
#   
#   # GCT
#   gctHeader <- c("#1.2", paste(dim(expZscore), collapse="\t"),
#                  paste(c("NAME", "Description", colnames(expZscore)[-c(1:2)]), collapse="\t"))
#   
#   gctContent <- c(gctHeader, 
#                   apply(expZscore, 1, function(x) paste(x, collapse="\t")))
#   
# 
#   write(gctContent, paste0(cancerType$DISEASE[1], "exp.gct"), sep="\n")  
#   
#   return(NULL)
# })
# 



ssGSEAPre <- lapply(split.data.frame(tcgaPanCanSamples, f = ~ DISEASE), function(cancerType){
  
  
  expZscore <- tcgaExpData[[cancerType$DISEASE[1]]]$Zscore
  expZscore <- expZscore[, !is.na(match(substr(colnames(expZscore), 1, 12), cancerType$PATIENT_BARCODE))]
  
  
  indices <- apply(expZscore, 1, function(x) all(is.na(x)))
  expZscore <- expZscore[!indices, ]
  
  expZscore <- tibble::rownames_to_column(expZscore, var = "GeneID")
  
  expZscore <- geneInfo[, c('GeneID', 'Symbol')] %>% 
    inner_join(expZscore, join_by(GeneID)) %>% mutate(GeneID = NULL) %>% column_to_rownames(var = 'Symbol')
  
  
  # GCT
  expZscore <- new("GCT", mat=as.matrix(expZscore))
  write_gct(ds = expZscore, ofile = cancerType$DISEASE[1], appenddim = TRUE, ver = 3)
  
  return(NULL)
})




#############GMT
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

rhoFamilies <- rhoFamilies %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID), 
                                      rhoGroup = ifelse(Class %in% c('DH_domain', 'DOCKER_dom'), 'DH_DOCKER_dom', Class))

genesets <- split(rhoFamilies$Approved.symbol, rhoFamilies$rhoGroup)
gmtContent <- sapply(seq(length(genesets)), function(i) paste(c(names(genesets)[i], 'NA', genesets[[i]]), collapse="\t"))

write(gmtContent, 'rhoFamilies.gmt', sep="\n")


