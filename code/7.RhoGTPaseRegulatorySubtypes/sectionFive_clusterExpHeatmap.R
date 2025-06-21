
library('dplyr')
library('tibble')
library('RColorBrewer')
library('pheatmap')
####################
tcgaExpData <- readRDS(file = '/data/tcgaFirehoseExpData.rds')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
difExpQvalue <- readRDS(file = '/data/ccClusterDiffExpGene.rds')
ccClusterMajorGenes <- readRDS(file = '/data/ccClusterMajorGenes.rds')


topGenes <- rhoFamilies %>% subset.data.frame(Approved.symbol %in% names(difExpQvalue)) %>% remove_rownames() %>% 
  column_to_rownames(var = 'NCBI.Gene.ID') %>% 
  mutate(Class = recode(Class, 'Rho_GDI' = 'GDI', 'RhoGAP_dom' = 'GAP',
                        'Small_GTPase_Rho' = 'Rho', 'DH_domain' = 'GEF', 'DOCKER_dom' = 'GEF'), Locus.type = NULL)


####################
setwd('/result/Section5/panCanCluster')

panCanCluster <- read.csv(file='panCanCluster.k=6.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
panCanCluster <- panCanCluster %>% dplyr::rename(PATIENT_BARCODE = V1, Clusters = V2) %>% 
  left_join(tcgaPanCanSamples[, c('PATIENT_BARCODE', 'SAMPLE_BARCODE', 'DISEASE')], by = 'PATIENT_BARCODE') %>% 
  arrange(Clusters, DISEASE) %>% mutate(Clusters = paste0('Cluster', Clusters), PATIENT_BARCODE = NULL) %>% 
  column_to_rownames(var = 'SAMPLE_BARCODE')


#################### TPM
disZscoreExp <- sapply(tcgaExpData, function(expData){
  
  expData <- expData$Zscore
  
  return(expData)
}, simplify = FALSE)

disZscoreExp <- bind_cols(disZscoreExp)
disZscoreExp <- disZscoreExp[rownames(topGenes), rownames(panCanCluster)]
rownames(disZscoreExp) <- topGenes$Approved.symbol

disZscoreExp[disZscoreExp > 2] <- 2
disZscoreExp[disZscoreExp < -2] <- -2


####################
rowAnno <- topGenes %>% remove_rownames() %>% column_to_rownames(var = 'Approved.symbol')

mycol <- brewer.pal(6, "Set3")
annColors <- list(Clusters = c("Cluster1" = mycol[1], "Cluster2" = mycol[2], "Cluster3" = mycol[3], 
                               "Cluster4" = mycol[4], "Cluster5" = mycol[5], "Cluster6" = mycol[6]))


pdf(file = '/result/Section5/panCanClusterExp.pdf', height = 7)
pheatmap(mat = disZscoreExp[ccClusterMajorGenes$Approved.symbol, ], cluster_cols = F, cluster_rows = F, 
         treeheight_row = 10, fontsize_row = 6, fontsize = 8, 
         color = colorRampPalette((c("#3878C1", "white","#AB221F")))(300), border_color = NA,
         annotation_col = panCanCluster, annotation_colors = annColors, annotation_row = rowAnno,
         show_colnames = F, show_rownames = T, na_col = "black")

dev.off()




####################
setwd('/result/Section5/panCanCluster')

panCanCluster <- read.csv(file='panCanCluster.k=6.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
panCanCluster <- panCanCluster %>% rename(PATIENT_BARCODE = V1, Clusters = V2) %>% 
  left_join(tcgaPanCanSamples[, c('PATIENT_BARCODE', 'DISEASE', 'SUBTYPE')], by = 'PATIENT_BARCODE') %>% arrange(Clusters, DISEASE) %>% 
  column_to_rownames(var = 'PATIENT_BARCODE') %>%  mutate(Clusters = paste0('Cluster', Clusters))


panCanCluster <- subset(panCanCluster, DISEASE %in% c('BRCA', 'CESC', 'COAD', 'ESCA', 'GBM', 'HNSC', 'LGG', 'READ', 'SARC', 'STAD', 'TGCT', 'UCEC'))
panCanCluster <- subset(panCanCluster, SUBTYPE != 'NA')

#######
fisherF <- function(data) {
  
  print(table(data[, c('Clusters', 'SUBTYPE')]))
  res <- chisq.test(table(data[, c('Clusters', 'SUBTYPE')]))
  
  return(data.frame(pvalue = res$p.value))
  
}

enrichPvalue <- panCanCluster %>% group_by(DISEASE) %>% do(fisherF(.))
enrichPvalue$fdr <- p.adjust(enrichPvalue$pvalue, method = 'fdr')


# DISEASE   pvalue      fdr
# <chr>      <dbl>    <dbl>
# 1 BRCA    5.68e-44 6.82e-43
# 2 CESC    1.02e-18 2.46e-18
# 3 COAD    1.02e- 6 1.53e- 6
# 4 ESCA    1.94e- 3 2.58e- 3
# 5 GBM     9.10e- 1 9.10e- 1
# 6 HNSC    2.07e-12 3.55e-12
# 7 LGG     1.52e-41 9.14e-41
# 8 READ    1.89e- 2 2.26e- 2
# 9 SARC    4.30e-15 8.59e-15
# 10 STAD    6.77e-22 2.71e-21
# 11 TGCT    9.01e-20 2.70e-19
# 12 UCEC    1.35e- 1 1.47e- 1

#######
panCanCluster <- panCanCluster %>% mutate(Clusters = gsub('Cluster', '', Clusters))


pdf(file = '/result/Section5/ccClusterAssoWithMolSubtype.pdf', width = 3, height = 3)
sapply(split.data.frame(panCanCluster, f = ~DISEASE), function(data){
  
  
  data <- data %>% arrange(desc(SUBTYPE))
  data$Clusters <- as.numeric(data$Clusters)
  
  pheatmap(mat = data[, 'Clusters', F], cluster_cols = F, cluster_rows = F, treeheight_row = 10, fontsize_row = 6, fontsize = 8, 
           color = brewer.pal(6, "Set3"), border_color = NA, cellwidth = 10,
           annotation_row = data[, 'SUBTYPE', F], show_colnames = F, show_rownames = F, na_col = "black", main = data$DISEASE[1])
})

dev.off()





