
library('tibble')
library('dplyr')
library('ComplexHeatmap')

####################
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
tcgaSsgseaRes <- readRDS(file = '/data/tcgaRhoFamiliesSsgsea2.0Res.rds')
load(file = '/data/tcgaExpSams.RData')


#######
ssGseaGroup <- tcgaSsgseaRes %>% subset(SampleIDs %in% tumorExpSams) %>% group_by(DISEASE, Groups) %>% 
  mutate(Class = ifelse(Scores > median(Scores), 'Up', 'Down')) %>% 
  ungroup() %>% merge.data.frame(tcgaPanCanSamples[, c('PATIENT_BARCODE', 'SUBTYPE')], by = 'PATIENT_BARCODE') %>% 
  subset(SUBTYPE != 'Not_Applicable' & SUBTYPE != 'NA')


#######
fisherF <- function(data) {
  
  # print(table(data[, c('Class', 'SUBTYPE')]))
  res <- chisq.test(table(data[, c('Class', 'SUBTYPE')]))
  
  return(data.frame(pvalue = res$p.value))
  
}

enrichPvalue <- ssGseaGroup %>% group_by(DISEASE, Groups) %>% do(fisherF(.))
enrichPvalue$fdr <- p.adjust(enrichPvalue$pvalue, method = 'fdr')

setwd('/result/Section5/subtypeAsso/')

write.table(enrichPvalue, sep = '\t', quote = F, row.names = F, col.names = T, file = 'subtypeAsso.txt')

#######
for(data in split.data.frame(ssGseaGroup, f = ~ DISEASE)){
  pdf(width = 4, file = paste0(data$DISEASE[1], '.pdf'))

  annoData <- data %>% select(PATIENT_BARCODE, SUBTYPE) %>% distinct() %>% column_to_rownames(var = 'PATIENT_BARCODE')
  classData <- data %>% reshape2::dcast(PATIENT_BARCODE ~ Groups, value.var = "Class") %>% column_to_rownames(var = 'PATIENT_BARCODE')

  annoData <- merge(annoData, classData, by= 'row.names') %>% column_to_rownames(var = 'Row.names') %>% arrange(SUBTYPE)
  
  ha <- HeatmapAnnotation(df = annoData)
  
  draw(ha)
  
  dev.off()
}


