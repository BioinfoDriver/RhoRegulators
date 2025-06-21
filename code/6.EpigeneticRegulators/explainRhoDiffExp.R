library('dplyr')
library('tibble')
library('ggplot2')

########################
load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp

rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% select(Approved.symbol) %>% distinct() %>% mutate(Direction = 'Neutral')


nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(log2FoldChange > 1 & padj < 0.05, 'Up', ifelse(log2FoldChange < -1 & padj < 0.05, 'Down', 'Neutral')))
nRDiffExp <- nRDiffExp %>% subset(Symbol %in% rhoFamilies$Approved.symbol)


#############
downGene <- nRDiffExp %>% subset.data.frame(Estatus == 'Down') %>% select(DISEASE, Symbol) %>% split.data.frame(., f=~DISEASE)
upGene <- nRDiffExp %>% subset.data.frame(Estatus == 'Up') %>% select(DISEASE, Symbol) %>% split.data.frame(., f=~DISEASE)

downGene <- sapply(downGene, function(x) return(x$Symbol))
upGene <- sapply(upGene, function(x) return(x$Symbol))






########################SCNA
load(file = '/data/genesGisticPeakQvalue.RData')
nrGisticPeakQvalue <- rhoGisticPeakQvalue
nrGisticPeakQvalue <- nrGisticPeakQvalue %>% as.data.frame() %>% select(Approved.symbol, Disease, Direction) %>% distinct()


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% select(Approved.symbol) %>% distinct() %>% mutate(Direction = 'Neutral')


nrGisticPeakQvalue <- split.data.frame(nrGisticPeakQvalue, f = ~ Disease)
nrGisticPeakQvalue <- lapply(nrGisticPeakQvalue, function(data){
  
  tmp <- rhoFamilies %>% subset(!(Approved.symbol %in% data$Approved.symbol)) %>% mutate(Disease = data$Disease[1])
  
  data <- bind_rows(data, tmp)
  return(data)
})
nrGisticPeakQvalue <- bind_rows(nrGisticPeakQvalue)


#############
load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp


nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(log2FoldChange > 1 & padj < 0.05, 'Up', ifelse(log2FoldChange < -1 & padj < 0.05, 'Down', 'Neutral')))
nRDiffExp <- nRDiffExp %>% subset(Symbol %in% rhoFamilies$Approved.symbol)

nRDiffExp <- nRDiffExp %>% left_join(nrGisticPeakQvalue,  by = join_by(DISEASE == Disease, Symbol == Approved.symbol)) %>% 
  mutate(Direction = factor(Direction, levels = c('Deletion', 'Neutral', 'Amplification')))


#############
downDelExpla <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% arrange(Direction) %>% distinct(Symbol, .keep_all = T) %>% ungroup() %>% 
    subset(Estatus == 'Down' & Direction == 'Deletion') %>% select(DISEASE, Symbol) %>% split.data.frame(., f=~DISEASE)
  
upAmpExpla <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% arrange(desc(Direction)) %>% distinct(Symbol, .keep_all = T) %>% ungroup() %>% 
  subset(Estatus == 'Up' & Direction == 'Amplification') %>% select(DISEASE, Symbol) %>% split.data.frame(., f=~DISEASE)


downDelExpla <- sapply(downDelExpla, function(x) return(x$Symbol))
upAmpExpla <- sapply(upAmpExpla, function(x) return(x$Symbol))





########################Methy
diffMethy <- readRDS(file = '/data/panCanDiffMethy.rds')
diffMethy <- diffMethy %>% mutate(Mstatus = ifelse(logFC > 0 & adjPVal < 0.05, 'Hyper', ifelse(logFC < 0 & adjPVal < 0.05, 'Hypo', 'Neutral')))


load(file = '/ata/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp
nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(log2FoldChange > 1 & padj < 0.05, 'Up', ifelse(log2FoldChange < -1 & padj < 0.05, 'Down', 'Neutral')))


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
nRDiffExp <- subset(nRDiffExp, Symbol %in% rhoFamilies$Approved.symbol)


nRDiffExp <- nRDiffExp %>% left_join(diffMethy[, c('Disease', 'Approved.symbol', 'Mstatus')], 
                                     by = join_by(DISEASE == Disease, Symbol == Approved.symbol)) %>% 
  mutate(Mstatus = ifelse(is.na(Mstatus), 'Other', Mstatus))


#############
downHyperExpla <- nRDiffExp %>% subset(Estatus == 'Down' & Mstatus == 'Hyper') %>% select(DISEASE, Symbol) %>% split.data.frame(., f=~DISEASE)

upHypoExpla <- nRDiffExp %>% subset(Estatus == 'Up' & Mstatus == 'Hypo') %>% select(DISEASE, Symbol) %>% split.data.frame(., f=~DISEASE)

downHyperExpla <- sapply(downHyperExpla, function(x) return(x$Symbol))
upHypoExpla <- sapply(upHypoExpla, function(x) return(x$Symbol))





########################MIRna
diffMirRna <- readRDS(file = '/data/panCanMiRnaDiffExp.rds')
load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
miRnaTarByDisease <- readRDS(file = '/data/miRnaTarByDiseaseTargetScan.rds')

#############
diffMirRna <- diffMirRna %>% mutate(MirStatus = ifelse(logFC > 0 & adjPVal < 0.05, 'Up', ifelse(logFC < -0 & adjPVal < 0.05, 'Down', 'Neutral'))) %>% 
  select(Approved.symbol, Disease, MirStatus) %>% dplyr::rename(miRNA = Approved.symbol, DISEASE = Disease)

diffMirRna <- diffMirRna %>% inner_join(miRnaTarByDisease, by = c('DISEASE', 'miRNA')) %>% select(miRNA, DISEASE, MirStatus, TargetGene)

nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(log2FoldChange > 1 & padj < 0.05, 'Up', ifelse(log2FoldChange < -1 & padj < 0.05, 'Down', 'Neutral')))
nRDiffExp <- subset(nRDiffExp, Symbol %in% rhoFamilies$Approved.symbol)

nRDiffExp <- nRDiffExp %>% left_join(diffMirRna,  by = join_by(DISEASE == DISEASE, Symbol == TargetGene)) %>% 
  mutate(MirStatus = ifelse(is.na(MirStatus), 'Other', MirStatus))



#############
downMirnaUpExpla <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% distinct(Symbol, MirStatus, .keep_all = T) %>% 
  arrange(desc(MirStatus)) %>% distinct(Symbol, .keep_all = T) %>% ungroup() %>% 
  subset(Estatus == 'Down' & MirStatus == 'Up') %>% select(DISEASE, Symbol) %>% split.data.frame(., f=~DISEASE)


upMirnaDownExpla <- nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
  arrange(desc(MirStatus)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% ungroup() %>% 
  subset(Estatus == 'Up' & MirStatus == 'Down') %>% select(DISEASE, Symbol) %>% split.data.frame(., f=~DISEASE)


downMirnaUpExpla <- sapply(downMirnaUpExpla, function(x) return(x$Symbol))
upMirnaDownExpla <- sapply(upMirnaDownExpla, function(x) return(x$Symbol))




########################
downGeneExpla <- sapply(names(downGene), function(Disease){
  
  dataList <- list(Del = downDelExpla[[Disease]], Hyper = downHyperExpla[[Disease]], miRNAUp = downMirnaUpExpla[[Disease]], 
       Other = setdiff(downGene[[Disease]], c(downDelExpla[[Disease]], downHyperExpla[[Disease]], downMirnaUpExpla[[Disease]])))
  
  dataList <- Filter(Negate(is.null), dataList)
  
  return(dataList)
}, simplify = F)


upGeneExpla <- sapply(names(upGene), function(Disease){
  
  dataList <- list(Amp = upAmpExpla[[Disease]], Hypo = upHypoExpla[[Disease]], miRNADown = upMirnaDownExpla[[Disease]], 
                   Other = setdiff(upGene[[Disease]], c(upAmpExpla[[Disease]], upHypoExpla[[Disease]], upMirnaDownExpla[[Disease]])))
  
  dataList <- Filter(Negate(is.null), dataList)
  
  return(dataList)
}, simplify = F)


#######################Plot
library(ggVennDiagram)

downvennPlot <- lapply(names(downGeneExpla), function(disease){
  
  vennPlot <- try(ggVennDiagram(downGeneExpla[[disease]], label_size = 2, set_size = 3) + ggtitle(disease) + theme(legend.position = "none"))
  
  return(vennPlot)
})



upvennPlot <- lapply(names(upGeneExpla), function(disease){
  
  vennPlot <- try(ggVennDiagram(upGeneExpla[[disease]], label_size = 2, set_size = 3) + ggtitle(disease) + theme(legend.position = "none"))
  
  return(vennPlot)
})



pdf('/result/Section4/expPlain.pdf')
cowplot::plot_grid(ncol = 5, plotlist = downvennPlot[c(1:5, 7:14)])
cowplot::plot_grid(ncol = 5, plotlist = upvennPlot[c(1, 5:14)])
dev.off()

# down
# $KICH
# $KICH$Other
# [1] "DOCK11"   "ARHGAP29" "ARHGEF33" "NET1"     "RHOB"     "CHN2"     "STARD13"  "MCF2"     "PLEKHG4"  "SYDE1"    "NGEF"     "RND2"     "TIAM1"    "FARP1"   
# [15] "STARD8"   "PLEKHG7"  "DOCK1"    "ARHGAP24" "ARHGAP10" "PIK3R1"   "FGD5"     "PLEKHG5"  "SYDE2"    "ARHGEF25" "PLEKHG1"  "SRGAP3"   "PLEKHG4B" "ARHGEF4" 
# [29] "ARHGEF19" "ARHGAP23" "RHOJ"     "RHOC"     "DLC1"     "PLEKHG6"  "PREX2" 

# up
# $BRCA
# $BRCA$Other
# [1] "DEPDC1B"   "ECT2"      "RACGAP1"   "ARHGAP11A" "PIK3R2"    "ARHGEF39"  "ARHGAP39"  "ARHGAP11B" "RND1"      "FGD3"      "RHOF"      "VAV3"     
# 
# 
# $COAD
# $COAD$Other
# [1] "ECT2"      "PLEKHG4"   "ARHGEF19"  "DEPDC1B"   "ARHGEF28"  "ARHGEF38"  "ARHGAP11B"
# 
# 
# $ESCA
# $ESCA$Other
# [1] "ECT2"      "ARHGAP11A" "DEPDC1B"   "SH3BP1"    "VAV2"      "ARHGAP11B" "RAC2"      "RHOV"      "PLEKHG4"   "RACGAP1"   "PLEKHG5"   "FGD6"      "ARHGAP8"  
# [14] "PLEKHG2"   "ARHGAP39"  "RHOH"      "PLEKHG6"   "CHN1"      "ARHGAP27"

