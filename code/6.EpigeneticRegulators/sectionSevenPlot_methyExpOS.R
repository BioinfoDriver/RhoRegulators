library('dplyr')
library('tibble')
library('reshape2')
library('ComplexHeatmap')
library('circlize')

#############
diffMethy <- readRDS(file = '/data/panCanDiffMethy.rds')

load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp

rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% distinct(Approved.symbol, .keep_all = TRUE) %>% 
  mutate(Class = recode(Class, DH_domain = 1, DOCKER_dom = 1, Rho_GDI = 2, RhoGAP_dom = 3, Small_GTPase_Rho = 4))

nRDiffExp <- subset(nRDiffExp, Symbol %in% rhoFamilies$Approved.symbol)

panCanSurPvalue <- readRDS(file = '/data/panCanRhoSurPvalue.rds')


#############
# diffMethy <- diffMethy %>% mutate(Mstatus = ifelse(logFC > 0 & adjPVal < 0.05, 'Hyper', ifelse(logFC < 0 & adjPVal < 0.05, 'Hypo', 'Neutral'))) %>% 
#   subset(Mstatus != 'Neutral') %>% dplyr::rename(MlogFC = logFC, MpValue = pValue, DISEASE = Disease, MadjPVal = adjPVal)


diffMethy <- diffMethy %>% mutate(Mstatus = ifelse(logFC > 0.3 & adjPVal < 0.05, 'Hyper', ifelse(logFC < -0.3 & adjPVal < 0.05, 'Hypo', 'Neutral'))) %>% 
  subset(Mstatus != 'Neutral') %>% dplyr::rename(MlogFC = logFC, MpValue = pValue, DISEASE = Disease, MadjPVal = adjPVal)


nRDiffExp <- nRDiffExp %>% subset(abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  mutate(status = ifelse(log2FoldChange > 0,  'Up', 'Down')) %>% dplyr::rename(Approved.symbol = Symbol)


ME_assoRes <- nRDiffExp %>% inner_join(diffMethy, by = c('DISEASE', 'Approved.symbol')) %>% 
  subset((Mstatus == 'Hyper' & status == 'Down') | (Mstatus == 'Hypo' & status == 'Up'))


panCanSurPvalue <- subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05)
panCanSurPvalue <- panCanSurPvalue %>% mutate(Direction = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival')) %>% 
  dplyr::rename(DISEASE = disease, Approved.symbol = geneSymbol)


ME_assoRes <- ME_assoRes %>% left_join(panCanSurPvalue, by = c('DISEASE', 'Approved.symbol'))


#####
# write.table(ME_assoRes, file = '/result/Section4/ME_assoRes.txt',
#             col.names = T, row.names = F, sep = '\t', quote = F)
#####


#####
# write.table(ME_assoRes, file = '/result/Section4/ME_Exp_explain.txt',
#             col.names = T, row.names = F, sep = '\t', quote = F)
#####



########################

geneStatus <- dcast(ME_assoRes, Approved.symbol ~ DISEASE, value.var = "log2FoldChange") %>% column_to_rownames(var = 'Approved.symbol')
methyStatus <- dcast(ME_assoRes, Approved.symbol ~ DISEASE, value.var = "MlogFC") %>% column_to_rownames(var = 'Approved.symbol')
osStatus <- dcast(ME_assoRes, Approved.symbol ~ DISEASE, value.var = "HR") %>% column_to_rownames(var = 'Approved.symbol')


rhoFamilies <- subset(rhoFamilies, Approved.symbol %in% rownames(geneStatus)) %>% arrange(Class)
geneStatus <- geneStatus[rhoFamilies$Approved.symbol, ]
methyStatus <- methyStatus[rhoFamilies$Approved.symbol, ]
osStatus <- osStatus[rhoFamilies$Approved.symbol, ]


gcolFun  <- colorRamp2(breaks = c(-6.0, 0, 3.5), colors = c("#3878C1", "white", "#AB221F"))
mcolFun  <- colorRamp2(breaks = c(-2, 0, 2.5), colors = c( "#9FBA95", "white","#B696B6"))

# gcolFun  <- colorRamp2(breaks = c(-5.5, 0, 4.5), colors = c("#FFC107", "white", "#E91E63"))
# mcolFun  <- colorRamp2(breaks = c(-1, 0, 3.5), colors = c( "#2196F3", "white","#4CAF50"))




DiagFunc <- function(gStatus, mStatus, surStatus){
  function(j, i, x, y, width, height, fill){
    
    if(!is.na(mStatus[i, j]))
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                   unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                   gp = gpar(fill = mcolFun(mStatus[i, j]), col = "black"))
    
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
                 cell_fun = DiagFunc(gStatus = geneStatus, mStatus = methyStatus, surStatus = osStatus))


lgd <- list(Legend(title = "Exp log2Fc", 
                   col_fun = gcolFun, 
                   at = c(-6.0, -3.0, 0, 2.0, 3.5), 
                   direction = "horizontal"), 
            Legend(title = "Methy log2Fc", 
                   col_fun = mcolFun, 
                   at = c(-2, -1, 0, 1.5, 2.5), 
                   direction = "horizontal"))


pdf(file = '/result/Section4/methylationExpressionOS.pdf')


draw(polt1, annotation_legend_list = lgd,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)

dev.off()