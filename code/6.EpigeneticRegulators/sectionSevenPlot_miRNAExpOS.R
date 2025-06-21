
library('dplyr')
library('tibble')
library('reshape2')
library('ComplexHeatmap')
library('circlize')
library('tibble')
#############
diffMirRna <- readRDS(file = '/data/panCanMiRnaDiffExp.rds')

load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp

rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% distinct(Approved.symbol, .keep_all = TRUE) %>% 
  mutate(Class = recode(Class, DH_domain = 1, DOCKER_dom = 1, Rho_GDI = 2, RhoGAP_dom = 3, Small_GTPase_Rho = 4))

nRDiffExp <- subset(nRDiffExp, Symbol %in% rhoFamilies$Approved.symbol)

panCanSurPvalue <- readRDS(file = '/data/panCanRhoSurPvalue.rds')

miRnaTarByDisease <- readRDS(file = '/data/miRnaTarByDiseaseTargetScan.rds')


#############
# diffMirRna <- diffMirRna %>% mutate(MirStatus = ifelse(logFC > 0 & adjPVal < 0.05, 'Up', ifelse(logFC < 0 & adjPVal < 0.05, 'Down', 'Neutral'))) %>% 
#   subset(MirStatus != 'Neutral') 


diffMirRna <- diffMirRna %>% mutate(MirStatus = ifelse(logFC > 0.585 & adjPVal < 0.05, 'Up', ifelse(logFC < -0.585 & adjPVal < 0.05, 'Down', 'Neutral'))) %>% 
  subset(MirStatus != 'Neutral') 
colnames(diffMirRna) <- c('miRNA', 'MirlogFC', 'MirP.Value', 'DISEASE', 'pairedPats', 'Miradj.P.Val', 'MirStatus')


diffMirRna <- diffMirRna %>% inner_join(miRnaTarByDisease, by = c('DISEASE', 'miRNA'))


nRDiffExp <- nRDiffExp %>% subset(abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  mutate(status = ifelse(log2FoldChange > 0,  'Up', 'Down'))



MirE_assoRes <- nRDiffExp %>% inner_join(diffMirRna,  by = join_by(DISEASE == DISEASE, Symbol == TargetGene)) %>% 
  subset((MirStatus == 'Up' & status == 'Down') | (MirStatus == 'Down' & status == 'Up'))


panCanSurPvalue <- subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05) #  
panCanSurPvalue <- panCanSurPvalue %>% mutate(Direction = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival')) %>% 
  dplyr::rename(DISEASE = disease, Symbol = geneSymbol)


MirE_assoRes <- MirE_assoRes %>% left_join(panCanSurPvalue, by = c('DISEASE', 'Symbol'))


#####

# write.table(MirE_assoRes, file = '/result/Section4/MirE_assoRes.txt',
#             col.names = T, row.names = F, sep = '\t', quote = F)

#####


# write.table(MirE_assoRes, file = '/result/Section4/Mir_Exp_explain.txt',
#             col.names = T, row.names = F, sep = '\t', quote = F)



#############
MirE_assoRes <- MirE_assoRes %>% mutate(geneDisease = paste(Symbol, DISEASE, sep = '_'))

geneStatus <- dcast(MirE_assoRes, miRNA ~ geneDisease, value.var = "log2FoldChange") %>% column_to_rownames(var = 'miRNA')
mirStatus <- dcast(MirE_assoRes, miRNA ~ geneDisease, value.var = "MirlogFC") %>% column_to_rownames(var = 'miRNA')
osStatus <- dcast(MirE_assoRes, miRNA ~ geneDisease, value.var = "HR") %>% column_to_rownames(var = 'miRNA')


rhoFamilies <- stringr::str_split(colnames(geneStatus), pattern = '_', simplify = T) %>% as.data.frame() %>% dplyr::rename(Approved.symbol = V1, Disease = V2) %>% 
  left_join(rhoFamilies, by = join_by(Approved.symbol)) %>% arrange(Class) %>% mutate(Approved.symbol = paste(Approved.symbol, Disease, sep = '_'))

geneStatus <- geneStatus[, rhoFamilies$Approved.symbol]
mirStatus <- mirStatus[, rhoFamilies$Approved.symbol]
osStatus <- osStatus[, rhoFamilies$Approved.symbol]



gcolFun  <- colorRamp2(breaks = c(-4.5, 0, 1.5), colors = c("#3878C1", "white", "#AB221F"))
mcolFun  <- colorRamp2(breaks = c(-1.0, 0, 3.0), colors = c( "#9FBA95", "white","#B696B6"))

# gcolFun  <- colorRamp2(breaks = c(-5.5, 0, 2.5), colors = c("#FFC107", "white", "#E91E63"))
# mcolFun  <- colorRamp2(breaks = c(-2.5, 0, 3.5), colors = c( "#2196F3", "white","#4CAF50"))



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



polt1 <- Heatmap(geneStatus, column_title = "miRNA status across cancer types",
                 rect_gp = gpar(type = "none"),
                 show_heatmap_legend = F,
                 cluster_rows = F,
                 cluster_columns = F, 
                 cell_fun = DiagFunc(gStatus = geneStatus, mStatus = mirStatus, surStatus = osStatus))


lgd <- list(Legend(title = "Exp log2Fc", 
                   col_fun = gcolFun, 
                   at = c(-4.5, -2.0, 0, 1.5), 
                   direction = "horizontal"), 
            Legend(title = "miRNA log2Fc", 
                   col_fun = mcolFun, 
                   at = c(-1.0, 0, 1.5, 3.0), 
                   direction = "horizontal"))


pdf(file = '/result/Section4/miRNAExpressionOS.pdf')


draw(polt1, annotation_legend_list = lgd,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)

dev.off()



