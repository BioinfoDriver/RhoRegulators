library('dplyr')
library('tibble')
library('ggplot2')
########################

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
fisherT <- function(data) {
  
  
  data <- data %>% distinct(Symbol, Estatus, MirStatus, .keep_all = T)
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(MirStatus = ifelse(MirStatus == 'Up', 'Up', 'Other'))
  downData <- downData %>% arrange(desc(MirStatus)) %>% distinct(Symbol, .keep_all = T)
  
  downtable <- table(downData$Estatus, downData$MirStatus)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  
  
  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(MirStatus = ifelse(MirStatus == 'Down', 'Down', 'Other'))
  downData <- downData %>% arrange(MirStatus) %>% distinct(Symbol, .keep_all = T)
  
  
  uptable <- table(upData$Estatus, upData$MirStatus)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}


fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.))

#    DISEASE      pValue        OR status
# 1     BLCA 0.012665254 0.3452513   Down
# 2     BLCA          NA        NA     Up
# 3     BRCA 0.643661805 0.7790729   Down
# 4     BRCA 0.072184132 0.1436829     Up
# 5     COAD          NA        NA   Down
# 6     COAD          NA        NA     Up
# 7     ESCA 0.108857624 0.2811917   Down
# 8     ESCA          NA        NA     Up
# 9     HNSC 0.183987437 0.5050974   Down
# 10    HNSC 0.564635041 0.7431481     Up
# 11    KICH          NA        NA   Down
# 12    KICH          NA        NA     Up
# 13    KIRC 0.061100987 0.3964930   Down
# 14    KIRC 0.016887834 8.0966823     Up
# 15    KIRP 0.364818691 2.3305853   Down
# 16    KIRP 0.128637528       Inf     Up
# 17    LIHC 1.000000000       Inf   Down
# 18    LIHC 0.575238695 0.7516823     Up
# 19    LUAD 0.029101245 0.3679685   Down
# 20    LUAD 0.120985680 0.2128317     Up
# 21    LUSC 0.080533128 0.4847227   Down
# 22    LUSC 0.342726478       Inf     Up
# 23    PRAD 0.787034897 0.7641133   Down
# 24    PRAD 0.002628697 0.0000000     Up
# 25    STAD 0.020153363 0.2443807   Down
# 26    STAD 0.389823828 0.4058112     Up
# 27    THCA 0.012625821 0.2189696   Down
# 28    THCA 1.000000000 0.8950683     Up



# Down
# nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>%
#   arrange(desc(MirStatus)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% group_by(Estatus, MirStatus) %>% count()

# > fisher.test(matrix(c(101, 280, 303, 1403), 2, byrow = T))$p.value
# [1] 0.0001810133
# > fisher.test(matrix(c(101, 280, 303, 1403), 2, byrow = T))$estimate
# odds ratio 
# 1.669745 


# Up
# nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
#   arrange(MirStatus) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% group_by(Estatus, MirStatus) %>% count()

# > fisher.test(matrix(c(13, 224, 95, 1611), 2, byrow = T))$p.value
# [1] 1
# > fisher.test(matrix(c(13, 224, 95, 1611), 2, byrow = T))$estimate
# odds ratio 
# 0.9841706 


########################
######Down
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% distinct(Symbol, MirStatus, .keep_all = T) %>% 
  arrange(desc(MirStatus)) %>% distinct(Symbol, .keep_all = T) %>% count(MirStatus)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, 
                                                               MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
  arrange(desc(MirStatus)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(MirStatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, 
                                                     MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


######Up
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% distinct(Symbol, MirStatus, .keep_all = T) %>% 
  arrange(MirStatus) %>% distinct(Symbol, .keep_all = T) %>% count(MirStatus)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, 
                                                               MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot3 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
  arrange(MirStatus) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(MirStatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, 
                                                     MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot4 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


pdf('/result/Section4/miRNAEnrich.pdf')

plot1
plot2
plot3
plot4

dev.off()




###################################regulon
diffMirRna <- readRDS(file = '/data/panCanMiRnaDiffExp.rds')
load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% subset(Class != 'Small_GTPase_Rho') %>% distinct(Approved.symbol, .keep_all = T)
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
fisherT <- function(data) {
  
  
  data <- data %>% distinct(Symbol, Estatus, MirStatus, .keep_all = T)
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(MirStatus = ifelse(MirStatus == 'Up', 'Up', 'Other'))
  downData <- downData %>% arrange(desc(MirStatus)) %>% distinct(Symbol, .keep_all = T)
  
  downtable <- table(downData$Estatus, downData$MirStatus)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  
  
  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(MirStatus = ifelse(MirStatus == 'Down', 'Down', 'Other'))
  downData <- downData %>% arrange(MirStatus) %>% distinct(Symbol, .keep_all = T)
  
  
  uptable <- table(upData$Estatus, upData$MirStatus)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}


fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.))

#    DISEASE      pValue        OR status
# 1     BLCA 0.009912363 0.2900518   Down
# 2     BLCA          NA        NA     Up
# 3     BRCA 0.457516545 0.6983883   Down
# 4     BRCA 0.068973954 0.1379868     Up
# 5     COAD          NA        NA   Down
# 6     COAD          NA        NA     Up
# 7     ESCA 0.116812037 0.2895338   Down
# 8     ESCA          NA        NA     Up
# 9     HNSC 0.391362232 0.5811373   Down
# 10    HNSC 0.542193092 0.6861528     Up
# 11    KICH          NA        NA   Down
# 12    KICH          NA        NA     Up
# 13    KIRC 0.048293671 0.3500124   Down
# 14    KIRC 0.003958333       Inf     Up
# 15    KIRP 0.359219027 2.4147499   Down
# 16    KIRP 0.218608578       Inf     Up
# 17    LIHC 1.000000000       Inf   Down
# 18    LIHC 1.000000000 0.8582574     Up
# 19    LUAD 0.041580473 0.3668422   Down
# 20    LUAD 0.120317320 0.2088729     Up
# 21    LUSC 0.037403897 0.4129736   Down
# 22    LUSC 0.597606795       Inf     Up
# 23    PRAD 1.000000000 1.1109316   Down
# 24    PRAD 0.002505593 0.0000000     Up
# 25    STAD 0.068442241 0.3232586   Down
# 26    STAD 0.350769802 0.3396460     Up
# 27    THCA 0.010101508 0.1994665   Down
# 28    THCA 1.000000000 0.9016738     Up



# Down
# nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>%
#   arrange(desc(MirStatus)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% group_by(Estatus, MirStatus) %>% count()

# > fisher.test(matrix(c(89, 244, 272, 1226), 2, byrow = T))$p.value
# [1] 0.0005869607
# > fisher.test(matrix(c(89, 244, 272, 1226), 2, byrow = T))$estimate
# odds ratio 
# 1.643561 


# Up
# nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
#   arrange(MirStatus) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% group_by(Estatus, MirStatus) %>% count()

# > fisher.test(matrix(c(12, 201, 88, 1410), 2, byrow = T))$p.value
# [1] 1
# > fisher.test(matrix(c(12, 201, 88, 1410), 2, byrow = T))$estimate
# odds ratio 
# 0.9565979 



########################
######Down
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% distinct(Symbol, MirStatus, .keep_all = T) %>% 
  arrange(desc(MirStatus)) %>% distinct(Symbol, .keep_all = T) %>% count(MirStatus)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, 
                                                               MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
  arrange(desc(MirStatus)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(MirStatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, 
                                                     MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


######Up
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% distinct(Symbol, MirStatus, .keep_all = T) %>% 
  arrange(MirStatus) %>% distinct(Symbol, .keep_all = T) %>% count(MirStatus)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, 
                                                               MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot3 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
  arrange(MirStatus) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(MirStatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, 
                                                     MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot4 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



pdf('/result/Section4/regulonmiRNAEnrich.pdf')

plot1
plot2
plot3
plot4

dev.off()






#################################Rho
diffMirRna <- readRDS(file = '/data/panCanMiRnaDiffExp.rds')
load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% subset(Class == 'Small_GTPase_Rho')

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
fisherT <- function(data) {
  
  
  data <- data %>% distinct(Symbol, Estatus, MirStatus, .keep_all = T)
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(MirStatus = ifelse(MirStatus == 'Up', 'Up', 'Other'))
  downData <- downData %>% arrange(desc(MirStatus)) %>% distinct(Symbol, .keep_all = T)
  
  downtable <- table(downData$Estatus, downData$MirStatus)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  
  
  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(MirStatus = ifelse(MirStatus == 'Down', 'Down', 'Other'))
  downData <- downData %>% arrange(MirStatus) %>% distinct(Symbol, .keep_all = T)
  
  
  uptable <- table(upData$Estatus, upData$MirStatus)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}


fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.))

#    DISEASE    pValue        OR status
# 1     BLCA 1.0000000 1.0000000   Down
# 2     BLCA        NA        NA     Up
# 3     BRCA 1.0000000 1.6224103   Down
# 4     BRCA        NA        NA     Up
# 5     COAD        NA        NA   Down
# 6     COAD        NA        NA     Up
# 7     ESCA        NA        NA   Down
# 8     ESCA        NA        NA     Up
# 9     HNSC 0.2982456 0.1581138   Down
# 10    HNSC 1.0000000       Inf     Up
# 11    KICH        NA        NA   Down
# 12    KICH        NA        NA     Up
# 13    KIRC 1.0000000       Inf   Down
# 14    KIRC 1.0000000 0.5223288     Up
# 15    KIRP        NA        NA   Down
# 16    KIRP 1.0000000       Inf     Up
# 17    LIHC        NA        NA   Down
# 18    LIHC        NA        NA     Up
# 19    LUAD 0.4052288 0.2567763   Down
# 20    LUAD        NA        NA     Up
# 21    LUSC 1.0000000 1.3097484   Down
# 22    LUSC        NA        NA     Up
# 23    PRAD 0.2621259 0.1843261   Down
# 24    PRAD        NA        NA     Up
# 25    STAD 0.0877193 0.0000000   Down
# 26    STAD 1.0000000       Inf     Up
# 27    THCA        NA        NA   Down
# 28    THCA 1.0000000       Inf     Up


# Down
# nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>%
#   arrange(desc(MirStatus)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% group_by(Estatus, MirStatus) %>% count()

# > fisher.test(matrix(c(12, 36, 31, 177), 2, byrow = T))$p.value
# [1] 0.1313104
# > fisher.test(matrix(c(12, 36, 31, 177), 2, byrow = T))$estimate
# odds ratio 
# 1.897845 


# Up
nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
  arrange(MirStatus) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% group_by(Estatus, MirStatus) %>% count()

# > fisher.test(matrix(c(1, 23, 7, 201), 2, byrow = T))$p.value
# [1] 0.5884587
# > fisher.test(matrix(c(1, 23, 7, 201), 2, byrow = T))$estimate
# odds ratio 
# 1.247151 



########################
######Down
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% distinct(Symbol, MirStatus, .keep_all = T) %>% 
  arrange(desc(MirStatus)) %>% distinct(Symbol, .keep_all = T) %>% count(MirStatus)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, 
                                                               MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
  arrange(desc(MirStatus)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(MirStatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, 
                                                     MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


######Up
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% distinct(Symbol, MirStatus, .keep_all = T) %>% 
  arrange(MirStatus) %>% distinct(Symbol, .keep_all = T) %>% count(MirStatus)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, 
                                                               MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot3 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% distinct(DISEASE, Symbol, MirStatus, .keep_all = T) %>% 
  arrange(MirStatus) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(MirStatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, 
                                                     MirStatus = factor(MirStatus, levels = c('Up', 'Neutral', 'Other', 'Down')))

plot4 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = MirStatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



pdf('/result/Section4/rhomiRNAEnrich.pdf')

plot1
plot2
plot3
plot4

dev.off()
