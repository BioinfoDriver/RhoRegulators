library('dplyr')
library('ggplot2')
library('ComplexHeatmap')
library('circlize')

######################
setwd('/result/Section5/panCanCluster')
panCanCluster <- read.csv(file='panCanCluster.k=6.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
panCanCluster <- panCanCluster %>% dplyr::rename(PATIENT_BARCODE = V1, Clusters = V2)  %>% 
  mutate(Clusters = paste0('Cluster', Clusters))


tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
panCanCluster <- panCanCluster %>% left_join(tcgaPanCanSamples, by = 'PATIENT_BARCODE')


clusterNum <- table(panCanCluster$DISEASE, panCanCluster$Clusters)
clusterPec <- as.data.frame(t(apply(clusterNum, 1, function(x) x/sum(x))))


######################
colors <- colorRamp2(breaks = c(0, 1), colors = c("white","red"))


DiagFunc <- function(num, pec){
  function(j, i, x, y, width, height, fill){
    
    if(pec[i, j] > 0){
      
      grid.text(num[i, j], x, y, gp = gpar(fontsize = 7))
    }
  }
}


p1 <- Heatmap(clusterPec, col = colors, 
              show_heatmap_legend = F, cluster_rows = F, cluster_columns = F, 
              cell_fun = DiagFunc(num = clusterNum, pec = clusterPec))

lgd <- list(Legend(title = "Percent", col_fun = colors, 
                   at = c(0,0.25, 0.5, 0.75, 1),  direction = "horizontal"))


pdf(file = '/result/Section5/clusterNumHeatmap.pdf', width = 6)

draw(p1, annotation_legend_list = lgd,
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)

dev.off()



