library('dplyr')
library('tibble')
library('ggplot2')
########################
diffMethy <- readRDS(file = '/data/panCanDiffMethy.rds')
diffMethy <- diffMethy %>% mutate(Mstatus = ifelse(logFC > 0 & adjPVal < 0.05, 'Hyper', ifelse(logFC < 0 & adjPVal < 0.05, 'Hypo', 'Neutral')))


load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp
nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(log2FoldChange > 1 & padj < 0.05, 'Up', ifelse(log2FoldChange < -1 & padj < 0.05, 'Down', 'Neutral')))


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
nRDiffExp <- subset(nRDiffExp, Symbol %in% rhoFamilies$Approved.symbol)


nRDiffExp <- nRDiffExp %>% left_join(diffMethy[, c('Disease', 'Approved.symbol', 'Mstatus')], 
                                      by = join_by(DISEASE == Disease, Symbol == Approved.symbol)) %>% 
  mutate(Mstatus = ifelse(is.na(Mstatus), 'Other', Mstatus))


fisherT <- function(data) {
  
  
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other'))
  
  
  downtable <- table(downData$Estatus, downData$Mstatus)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  
  
  
  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(Mstatus = ifelse(Mstatus == 'Hypo', 'Hypo', 'Other'))
  
  uptable <- table(upData$Estatus, upData$Mstatus)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}

fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.)) 
fisherPvalue$FDR <- p.adjust(fisherPvalue$pValue, method = 'fdr')

#    DISEASE     pValue        OR status       FDR
# 1     BLCA 1.00000000 0.0000000   Down 1.0000000
# 2     BLCA 0.39856226 0.5180363     Up 0.7401871
# 3     BRCA 0.02546163 2.6063474   Down 0.1880533
# 4     BRCA 0.36132138       Inf     Up 0.7226428
# 5     COAD 0.02893128 3.3383399   Down 0.1880533
# 6     COAD 1.00000000       Inf     Up 1.0000000
# 7     ESCA 0.68067804 1.2257097   Down 0.9314542
# 8     ESCA 1.00000000       Inf     Up 1.0000000
# 9     HNSC 0.13308777 2.9931813   Down 0.3460282
# 10    HNSC 1.00000000 1.3713692     Up 1.0000000
# 11    KICH         NA        NA   Down        NA
# 12    KICH         NA        NA     Up        NA
# 13    KIRC 0.11174835 2.2538023   Down 0.3460282
# 14    KIRC 0.53672803 0.6633206     Up 0.8232736
# 15    KIRP 0.58984250 1.3154151   Down 0.8519947
# 16    KIRP 0.10288268 0.3399966     Up 0.3460282
# 17    LIHC 0.13186448 3.1242819   Down 0.3460282
# 18    LIHC 0.12902612       Inf     Up 0.3460282
# 19    LUAD 0.02631519 3.6165058   Down 0.1880533
# 20    LUAD 0.01214482 0.2057023     Up 0.1880533
# 21    LUSC 1.00000000 0.9470255   Down 1.0000000
# 22    LUSC 0.53829426 0.6835895     Up 0.8232736
# 23    PRAD 0.48354899 1.6722091   Down 0.8232736
# 24    PRAD 0.35157558 0.3514954     Up 0.7226428
# 25    STAD 1.00000000 0.6342847   Down 1.0000000
# 26    STAD 1.00000000 1.6810235     Up 1.0000000
# 27    THCA 0.30472338 4.0329291   Down 0.7202553
# 28    THCA 0.08895079 0.3283908     Up 0.3460282


# > table(nRDiffExp$Estatus, nRDiffExp$Mstatus)
#         Hyper Hypo Neutral Other
# Down       68   33     113   167
# Neutral   179  164     706   657
# Up         15   30      75   117
# Down
# > fisher.test(matrix(c(68, 313, 179, 1527), 2, byrow = T))$p.value
# [1] 0.0001453857
# > fisher.test(matrix(c(68, 313, 179, 1527), 2, byrow = T))$estimate
# odds ratio 
# 1.852775 

# Up
# > fisher.test(matrix(c(30, 207, 164, 1542), 2, byrow = T))$p.value
# [1] 0.1642053
# > fisher.test(matrix(c(30, 207, 164, 1542), 2, byrow = T))$estimate
# odds ratio 
# 1.362433


########################
# dataStat <- nRDiffExp %>% subset(Estatus != 'Up') %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other')) %>%
#   group_by(DISEASE, Estatus) %>% count(Mstatus)

dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% count(Mstatus)


dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Mstatus = factor(Mstatus, levels = c('Hyper', 'Neutral', 'Other', 'Hypo')))


plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Mstatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


# totalStat <- nRDiffExp %>% subset(Estatus != 'Up') %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other')) %>% 
#   group_by(Estatus) %>% count(Mstatus)

totalStat <- nRDiffExp %>% group_by(Estatus) %>% count(Mstatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Mstatus = factor(Mstatus, levels = c('Hyper', 'Neutral', 'Other', 'Hypo')))


plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Mstatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


pdf('/result/Section4/methyEnrich.pdf')

plot1
plot2

dev.off()




#####################################regulon
diffMethy <- readRDS(file = '/data/panCanDiffMethy.rds')
diffMethy <- diffMethy %>% mutate(Mstatus = ifelse(logFC > 0 & adjPVal < 0.05, 'Hyper', ifelse(logFC < 0 & adjPVal < 0.05, 'Hypo', 'Neutral')))


load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp
nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(log2FoldChange > 1 & padj < 0.05, 'Up', ifelse(log2FoldChange < -1 & padj < 0.05, 'Down', 'Neutral')))


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% subset(Class != 'Small_GTPase_Rho') %>% distinct(Approved.symbol,.keep_all = TRUE)
nRDiffExp <- subset(nRDiffExp, Symbol %in% rhoFamilies$Approved.symbol)


nRDiffExp <- nRDiffExp %>% left_join(diffMethy[, c('Disease', 'Approved.symbol', 'Mstatus')], 
                                      by = join_by(DISEASE == Disease, Symbol == Approved.symbol)) %>% 
  mutate(Mstatus = ifelse(is.na(Mstatus), 'Other', Mstatus))


fisherT <- function(data) {
  
  
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other'))
  
  
  downtable <- table(downData$Estatus, downData$Mstatus)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  
  
  
  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(Mstatus = ifelse(Mstatus == 'Hypo', 'Hypo', 'Other'))
  
  uptable <- table(upData$Estatus, upData$Mstatus)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}

fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.)) 
fisherPvalue$FDR <- p.adjust(fisherPvalue$pValue, method = 'fdr')

#    DISEASE     pValue        OR status       FDR
# 1     BLCA 1.00000000 0.0000000   Down 1.0000000
# 2     BLCA 0.63915461 0.6620461     Up 0.9775306
# 3     BRCA 0.02892305 2.7458327   Down 0.3135344
# 4     BRCA 0.35941359       Inf     Up 0.7188272
# 5     COAD 0.03662392 3.3618539   Down 0.3135344
# 6     COAD 1.00000000       Inf     Up 1.0000000
# 7     ESCA 0.68957595 1.1796721   Down 0.9960541
# 8     ESCA 1.00000000       Inf     Up 1.0000000
# 9     HNSC 0.10818209 3.4223120   Down 0.3515918
# 10    HNSC 0.60726989       Inf     Up 0.9775306
# 11    KICH         NA        NA   Down        NA
# 12    KICH         NA        NA     Up        NA
# 13    KIRC 0.09746785 2.3494206   Down 0.3515918
# 14    KIRC 0.32334028 0.5317764     Up 0.7005706
# 15    KIRP 0.77492770 1.1748824   Down 1.0000000
# 16    KIRP 0.05128723 0.2409372     Up 0.3135344
# 17    LIHC 0.07194235 4.4611838   Down 0.3135344
# 18    LIHC 0.12486737       Inf     Up 0.3607280
# 19    LUAD 0.02457903 3.8046496   Down 0.3135344
# 20    LUAD 0.07235409 0.2835642     Up 0.3135344
# 21    LUSC 0.80640606 1.2286652   Down 1.0000000
# 22    LUSC 0.51003851 0.6846911     Up 0.8840668
# 23    PRAD 0.40939147 1.8483774   Down 0.7602984
# 24    PRAD 1.00000000       Inf     Up 1.0000000
# 25    STAD 1.00000000 0.6387657   Down 1.0000000
# 26    STAD 1.00000000 1.5521065     Up 1.0000000
# 27    THCA 0.24357017 6.8816611   Down 0.5757113
# 28    THCA 0.14672768 0.3467688     Up 0.3814920


# > table(nRDiffExp$Estatus, nRDiffExp$Mstatus)
#         Hyper Hypo Neutral Other
# Down       61   29      80   163
# Neutral   155  140     574   629
# Up          9   23      65   116
# Down
# > fisher.test(matrix(c(61, 272, 155, 1343), 2, byrow = T))$p.value
# [1] 0.0001050999
# > fisher.test(matrix(c(61, 272, 155, 1343), 2, byrow = T))$estimate
# odds ratio 
# 1.942372

# Up
# > fisher.test(matrix(c(23, 190, 140, 1358), 2, byrow = T))$p.value
# [1] 0.5322511
# > fisher.test(matrix(c(23, 190, 140, 1358), 2, byrow = T))$estimate
# odds ratio 
# 1.174061


########################
# dataStat <- nRDiffExp %>% subset(Estatus != 'Up') %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other')) %>%
#   group_by(DISEASE, Estatus) %>% count(Mstatus)

dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% count(Mstatus)


dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Mstatus = factor(Mstatus, levels = c('Hyper', 'Neutral', 'Other', 'Hypo')))


plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Mstatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


# totalStat <- nRDiffExp %>% subset(Estatus != 'Up') %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other')) %>% 
#   group_by(Estatus) %>% count(Mstatus)

totalStat <- nRDiffExp %>% group_by(Estatus) %>% count(Mstatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Mstatus = factor(Mstatus, levels = c('Hyper', 'Neutral', 'Other', 'Hypo')))


plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Mstatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


pdf('/result/Section4/regulonMethyEnrich.pdf')

plot1
plot2

dev.off()




#####################################Rho
diffMethy <- readRDS(file = '/data/panCanDiffMethy.rds')
diffMethy <- diffMethy %>% mutate(Mstatus = ifelse(logFC > 0 & adjPVal < 0.05, 'Hyper', ifelse(logFC < 0 & adjPVal < 0.05, 'Hypo', 'Neutral')))


load(file = '/data/tcgaPairedSamsDiffExp.RData')
nRDiffExp <- pairedSamsDeseq2DiffExp
nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(log2FoldChange > 1 & padj < 0.05, 'Up', ifelse(log2FoldChange < -1 & padj < 0.05, 'Down', 'Neutral')))


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% subset(Class == 'Small_GTPase_Rho') %>% distinct(Approved.symbol,.keep_all = TRUE)
nRDiffExp <- subset(nRDiffExp, Symbol %in% rhoFamilies$Approved.symbol)


nRDiffExp <- nRDiffExp %>% left_join(diffMethy[, c('Disease', 'Approved.symbol', 'Mstatus')], 
                                      by = join_by(DISEASE == Disease, Symbol == Approved.symbol)) %>% 
  mutate(Mstatus = ifelse(is.na(Mstatus), 'Other', Mstatus))


fisherT <- function(data) {
  
  
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other'))
  
  
  downtable <- table(downData$Estatus, downData$Mstatus)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  
  
  
  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(Mstatus = ifelse(Mstatus == 'Hypo', 'Hypo', 'Other'))
  
  uptable <- table(upData$Estatus, upData$Mstatus)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}

fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.)) 
fisherPvalue$FDR <- p.adjust(fisherPvalue$pValue, method = 'fdr')

#    DISEASE     pValue        OR status       FDR
# 1     BLCA         NA        NA   Down        NA
# 2     BLCA 0.39560440 0.2335493     Up 1.0000000
# 3     BRCA 1.00000000 1.7400212   Down 1.0000000
# 4     BRCA 1.00000000       Inf     Up 1.0000000
# 5     COAD 0.55979102 2.5220116   Down 1.0000000
# 6     COAD         NA        NA     Up        NA
# 7     ESCA         NA        NA   Down        NA
# 8     ESCA         NA        NA     Up        NA
# 9     HNSC 1.00000000 0.0000000   Down 1.0000000
# 10    HNSC 0.11111111 0.0000000     Up 0.7916667
# 11    KICH         NA        NA   Down        NA
# 12    KICH         NA        NA     Up        NA
# 13    KIRC 1.00000000 1.6099403   Down 1.0000000
# 14    KIRC 1.00000000       Inf     Up 1.0000000
# 15    KIRP 0.46470588 2.7673628   Down 1.0000000
# 16    KIRP 1.00000000       Inf     Up 1.0000000
# 17    LIHC 1.00000000 0.0000000   Down 1.0000000
# 18    LIHC         NA        NA     Up        NA
# 19    LUAD 1.00000000 0.0000000   Down 1.0000000
# 20    LUAD 0.06535948 0.0000000     Up 0.7916667
# 21    LUSC 0.25000000 0.0000000   Down 0.9500000
# 22    LUSC 1.00000000 0.7661432     Up 1.0000000
# 23    PRAD 1.00000000 0.9208139   Down 1.0000000
# 24    PRAD 0.12500000 0.0000000     Up 0.7916667
# 25    STAD         NA        NA   Down        NA
# 26    STAD 1.00000000       Inf     Up 1.0000000
# 27    THCA         NA        NA   Down        NA
# 28    THCA 0.25000000 0.0000000     Up 0.9500000


# > table(nRDiffExp$Estatus, nRDiffExp$Mstatus)
#         Hyper Hypo Neutral Other
# Down        7    4      33     4
# Neutral    24   24     132    28
# Up          6    7      10     1
# Down
# > fisher.test(matrix(c(7, 41, 24, 184), 2, byrow = T))$p.value
# [1] 0.6231211
# > fisher.test(matrix(c(7, 41, 24, 184), 2, byrow = T))$estimate
# odds ratio 
# 1.307478

# Up
# > fisher.test(matrix(c(7, 17, 24, 184), 2, byrow = T))$p.value
# [1] 0.0255189
# > fisher.test(matrix(c(7, 17, 24, 184), 2, byrow = T))$estimate
# odds ratio 
# 3.13569 

########################
# dataStat <- nRDiffExp %>% subset(Estatus != 'Up') %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other')) %>%
#   group_by(DISEASE, Estatus) %>% count(Mstatus)

dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% count(Mstatus)


dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Mstatus = factor(Mstatus, levels = c('Hyper', 'Neutral', 'Other', 'Hypo')))


plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Mstatus)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


# totalStat <- nRDiffExp %>% subset(Estatus != 'Up') %>% mutate(Mstatus = ifelse(Mstatus == 'Hyper', 'Hyper', 'Other')) %>% 
#   group_by(Estatus) %>% count(Mstatus)

totalStat <- nRDiffExp %>% group_by(Estatus) %>% count(Mstatus)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Mstatus = factor(Mstatus, levels = c('Hyper', 'Neutral', 'Other', 'Hypo')))


plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Mstatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


pdf('/result/Section4/rhoMethyEnrich.pdf')

plot1
plot2

dev.off()