

library('tibble')
library('ConsensusClusterPlus')
library('dplyr')

#########################
tcgaExpData <- readRDS(file = '/data/tcgaFirehoseExpData.rds')
load(file = '/data/tcgaExpSams.RData')
# tumorExpSams, normalExpSams, pairedSamsPats
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

##############
rhoFamilies <- rhoFamilies %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID)) %>% 
  subset(NCBI.Gene.ID %in% rownames(tcgaExpData$ACC$Zscore) & !duplicated(NCBI.Gene.ID))

# regulon
rhoFamilies <- rhoFamilies %>% subset(Class != 'Small_GTPase_Rho')



# Zscore
disZscoreExp <- lapply(tcgaExpData, function(expData){
  
  expData <- expData$Zscore[rhoFamilies$NCBI.Gene.ID, ]
  rownames(expData) <- rhoFamilies$Approved.symbol
  
  return(expData)
})

disZscoreExp <- cbind.data.frame(disZscoreExp)
colnames(disZscoreExp) <- stringr::str_split(colnames(disZscoreExp), "\\.", simplify = TRUE)[, 2]

disZscoreExp <- disZscoreExp[, intersect(colnames(disZscoreExp), tumorExpSams)]
colnames(disZscoreExp) <- substr(colnames(disZscoreExp), 1, 12)



# TPM
disTpmExp <- lapply(tcgaExpData, function(expData){
  
  expData <- expData$TPM[rhoFamilies$NCBI.Gene.ID, ]
  rownames(expData) <- rhoFamilies$Approved.symbol
  
  return(expData)
})

disTpmExp <- cbind.data.frame(disTpmExp)
colnames(disTpmExp) <- stringr::str_split(colnames(disTpmExp), "\\.", simplify = TRUE)[, 2]

disTpmExp <- disTpmExp[, intersect(colnames(disTpmExp), tumorExpSams)]
colnames(disTpmExp) <- substr(colnames(disTpmExp), 1, 12)

# Filtering
indices <- apply(disTpmExp, 1, function(x) {sum(x == 0)/ length(x) < 0.2})
disTpmExp <- disTpmExp[indices, ]



# Top 50 most variable genes
# reg
topVarGenes <- names(sort(apply(disTpmExp, 1, sd), T)[1:50])

# DH_domain DOCKER_dom    Rho_GDI RhoGAP_dom 
# 26          5          1         18 

disTopVarGenesZscoreExp <- as.matrix(disZscoreExp[topVarGenes, ])


##############
# Consensus Cluster

setwd('/result/Section5/')

ccRes <- ConsensusClusterPlus(
  d=disTopVarGenesZscoreExp, maxK = 9, reps=100, pItem=0.8, pFeature=1, clusterAlg="hc", 
  title="panCanCluster", innerLinkage="ward.D2", finalLinkage="ward.D2", distance="pearson", ml=NULL,
  tmyPal=NULL, seed=1024, plot='png', writeTable=TRUE,
  weightsItem=NULL, weightsFeature=NULL, verbose=TRUE, corUse="na.or.complete")

saveRDS(ccRes, file = 'panCanCluster_reg_top50.rds')



