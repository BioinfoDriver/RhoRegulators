
library('dplyr')
library('tibble')
library('reshape2')
library('ComplexHeatmap')
library('circlize')

########################
load(file = '/data/genesGisticPeakQvalue.RData')
nrGisticPeakQvalue <- rhoGisticPeakQvalue

load(file = '/data/tcgaPairedSamsDiffExp.RData')

nRDiffExp <- pairedSamsDeseq2DiffExp
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% distinct(Approved.symbol, .keep_all = TRUE) %>% 
  mutate(Class = recode(Class, DH_domain = 1, DOCKER_dom = 1, Rho_GDI = 2, RhoGAP_dom = 3, Small_GTPase_Rho = 4))


nRDiffExp <- subset(nRDiffExp, Symbol %in% rhoFamilies$Approved.symbol)

panCanSurPvalue <- readRDS(file = '/data/panCanRhoSurPvalue.rds')



#############
nrGisticPeakQvalue <- nrGisticPeakQvalue %>% as.data.frame() %>% select(Approved.symbol, Descriptor, Wide.Peak.Limits, q.values, Disease, Direction)
nrGisticPeakQvalue <- nrGisticPeakQvalue[!duplicated(nrGisticPeakQvalue), ]


nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(log2FoldChange > 1 & padj < 0.05, 'Up', ifelse(log2FoldChange < -1 & padj < 0.05, 'Down', 'Neutral')))
nRDiffExp <- nRDiffExp %>% inner_join(nrGisticPeakQvalue,  by = join_by(DISEASE == Disease, Symbol == Approved.symbol))



panCanSurPvalue <- subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05)
panCanSurPvalue <- panCanSurPvalue %>% mutate(OS_status = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival')) %>% 
  dplyr::rename(DISEASE = disease, Approved.symbol = geneSymbol)


SCNAE_assoRes <- nRDiffExp %>% left_join(panCanSurPvalue, by = join_by(DISEASE, Symbol == Approved.symbol)) %>% 
  subset((Direction == 'Deletion' & Estatus == 'Down') | (Direction == 'Amplification' & Estatus == 'Up'))


#####

# write.table(SCNAE_assoRes, file = '/result/Section4/SCNAE_assoRes.txt',
#             col.names = T, row.names = F, sep = '\t', quote = F)

#####




########################
# SCNAE_assoRes <- SCNAE_assoRes %>% mutate(geneDisease = paste(Approved.symbol, DISEASE, sep = '_'))


geneStatus <- dcast(SCNAE_assoRes, Symbol ~ DISEASE, value.var = "log2FoldChange") %>% column_to_rownames(var = 'Symbol')
scnaStatus <- dcast(SCNAE_assoRes, Symbol ~ DISEASE, value.var = "Direction") %>% column_to_rownames(var = 'Symbol')
osStatus <- dcast(SCNAE_assoRes, Symbol ~ DISEASE, value.var = "HR") %>% column_to_rownames(var = 'Symbol')


rhoFamilies <- subset(rhoFamilies, Approved.symbol %in% rownames(geneStatus)) %>% arrange(Class)
geneStatus <- geneStatus[rhoFamilies$Approved.symbol, ]
scnaStatus <- scnaStatus[rhoFamilies$Approved.symbol, ]
osStatus <- osStatus[rhoFamilies$Approved.symbol, ]


gcolFun  <- colorRamp2(breaks = c(-6, 0, 4.0), colors = c("#3878C1", "white", "#AB221F"))
# mcolFun  <- colorRamp2(breaks = c(-1, 0, 3.5), colors = c( "#9FBA95", "white","#B696B6"))

# gcolFun  <- colorRamp2(breaks = c(-5.5, 0, 4.5), colors = c("#FFC107", "white", "#E91E63"))
# mcolFun  <- colorRamp2(breaks = c(-1, 0, 3.5), colors = c( "#2196F3", "white","#4CAF50"))




DiagFunc <- function(gStatus, mStatus, surStatus){
  function(j, i, x, y, width, height, fill){
    
    if(!is.na(mStatus[i, j])){
      if(mStatus[i, j] == 'Deletion'){
        grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                     unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                     gp = gpar(fill = "#9FBA95", col = "black"))        
        
      }else{
        grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                     unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                     gp = gpar(fill = "#B696B6", col = "black"))       
        
      }
    }
    
    if(!is.na(gStatus[i, j]))
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = gpar(fill = gcolFun(gStatus[i, j]), col = "black"))    
    
    
    
    if(is.na(surStatus[i, j])){
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(col = "white", fill = NA))
      
      
    }else if(surStatus[i, j] > 1){
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(col = '#FB8C62', fill = NA))
      
    }else{
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(col = "#8C9FCA", fill = NA))  
      
    }
  }
}



polt1 <- Heatmap(geneStatus, column_title = "Methylation status across cancer types",
                 rect_gp = gpar(type = "none"),
                 show_heatmap_legend = F,
                 cluster_rows = F,
                 cluster_columns = F, 
                 cell_fun = DiagFunc(gStatus = geneStatus, mStatus = scnaStatus, surStatus = osStatus))


lgd <- list(Legend(title = "Exp log2Fc", 
                   col_fun = gcolFun, 
                   at = c(-6.0, -3.0, 0, 2.0, 4.0), 
                   direction = "horizontal"))



pdf(file = '/result/Section4/scnaExpressionOS.pdf')

draw(polt1, annotation_legend_list = lgd,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)

dev.off()

