

library('dplyr')
library('tibble')
library('ggplot2')
########################
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
fisherT <- function(data) {
  
  
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(Direction = ifelse(Direction == 'Deletion', 'Deletion', 'Other'))
  downData <- downData %>% arrange(Direction) %>% distinct(Symbol, .keep_all = T)
  
  downtable <- table(downData$Estatus, downData$Direction)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  
  
  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(Direction = ifelse(Direction == 'Amplification', 'Amplification', 'Other'))
  upData <- upData %>% arrange(desc(Direction)) %>% distinct(Symbol, .keep_all = T)
  
  uptable <- table(upData$Estatus, upData$Direction)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}


fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.))

#    DISEASE     pValue         OR status
# 1     BLCA 0.76743572 1.27980148   Down
# 2     BLCA 0.11571965 0.19498374     Up
# 3     BRCA 1.00000000 0.81531863   Down
# 4     BRCA 1.00000000        Inf     Up
# 5     COAD 1.00000000 0.87073307   Down
# 6     COAD         NA         NA     Up
# 7     ESCA 0.68067804 1.22570974   Down
# 8     ESCA 0.21716039        Inf     Up
# 9     HNSC 0.67421551 1.25948044   Down
# 10    HNSC 0.04372522 0.09272265     Up
# 11    KICH         NA         NA   Down
# 12    KICH 0.07946305 0.10086499     Up
# 13    KIRC 0.76400009 0.68596616   Down
# 14    KIRC 1.00000000 1.29327661     Up
# 15    KIRP 1.00000000 0.94394709   Down
# 16    KIRP 0.41103716 0.37528431     Up
# 17    LIHC 0.31168249 1.76973078   Down
# 18    LIHC 1.00000000        Inf     Up
# 19    LUAD 0.31864417 0.50890731   Down
# 20    LUAD 1.00000000        Inf     Up
# 21    LUSC 0.79651897 0.77223420   Down
# 22    LUSC 0.13403826 0.14445489     Up
# 23    PRAD 1.00000000 0.79865613   Down
# 24    PRAD 1.00000000        Inf     Up
# 25    STAD 1.00000000 0.58743995   Down
# 26    STAD 0.30504587 0.27053395     Up
# 27    THCA 0.24533367 2.02748862   Down
# 28    THCA 1.00000000        Inf     Up
# 



# Down
# nRDiffExp %>% group_by(Estatus) %>% arrange(Direction) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% 
#   group_by(Estatus, Direction) %>% count()
# 

# > fisher.test(matrix(c(47, 334, 218, 1488), 2, byrow = T))$p.value
# [1] 0.8650802
# > fisher.test(matrix(c(47, 334, 218, 1488), 2, byrow = T))$estimate
# odds ratio 
# 0.9604942


# Up
# nRDiffExp %>% group_by(Estatus) %>% arrange(desc(Direction)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>%
#   group_by(Estatus, Direction) %>% count()

# > fisher.test(matrix(c(14, 223, 74, 1632), 2, byrow = T))$p.value
# [1] 0.3147577
# > fisher.test(matrix(c(14, 223, 74, 1632), 2, byrow = T))$estimate
# odds ratio 
# 1.384298


########################
####down
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% arrange(Direction) %>% distinct(Symbol, .keep_all = T) %>% count(Direction)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% arrange(Direction) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(Direction)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


####up
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% arrange(desc(Direction)) %>% distinct(Symbol, .keep_all = T) %>% count(Direction)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot3 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% arrange(desc(Direction)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(Direction)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot4 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



pdf('/result/Section4/scnaEnrich.pdf')

plot1
plot2
plot3
plot4

dev.off()



################################################regulon
########################
load(file = '/data/genesGisticPeakQvalue.RData')
nrGisticPeakQvalue <- rhoGisticPeakQvalue
nrGisticPeakQvalue <- nrGisticPeakQvalue %>% as.data.frame() %>% subset(Class != 'Small_GTPase_Rho') %>% 
  select(Approved.symbol, Disease, Direction) %>% distinct()


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% subset(Class != 'Small_GTPase_Rho') %>% 
  select(Approved.symbol) %>% distinct() %>% mutate(Direction = 'Neutral')


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

nRDiffExp <- nRDiffExp %>% left_join(nrGisticPeakQvalue, by = join_by(DISEASE == Disease, Symbol == Approved.symbol)) %>% 
  mutate(Direction = factor(Direction, levels = c('Deletion', 'Neutral', 'Amplification')))


#############
fisherT <- function(data) {
  
  
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(Direction = ifelse(Direction == 'Deletion', 'Deletion', 'Other'))
  downData <- downData %>% arrange(Direction) %>% distinct(Symbol, .keep_all = T)
  
  downtable <- table(downData$Estatus, downData$Direction)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  

  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(Direction = ifelse(Direction == 'Amplification', 'Amplification', 'Other'))
  upData <- upData %>% arrange(desc(Direction)) %>% distinct(Symbol, .keep_all = T)
  
  uptable <- table(upData$Estatus, upData$Direction)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}


fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.))


#    DISEASE     pValue         OR status
# 1     BLCA 0.35021170 1.69677219   Down
# 2     BLCA 0.11146001 0.18823379     Up
# 3     BRCA 0.68173350 1.22425173   Down
# 4     BRCA 1.00000000        Inf     Up
# 5     COAD 1.00000000 0.76328625   Down
# 6     COAD         NA         NA     Up
# 7     ESCA 0.63606079 1.54302344   Down
# 8     ESCA 0.35640697        Inf     Up
# 9     HNSC 0.31033319 2.11929969   Down
# 10    HNSC 0.04819973 0.09760869     Up
# 11    KICH         NA         NA   Down
# 12    KICH 0.03859649 0.00000000     Up
# 13    KIRC 1.00000000 0.84492557   Down
# 14    KIRC 1.00000000 1.28716875     Up
# 15    KIRP 0.79595894 0.72555996   Down
# 16    KIRP 0.40660162 0.36779100     Up
# 17    LIHC 0.24140621 2.30056367   Down
# 18    LIHC 1.00000000        Inf     Up
# 19    LUAD 0.43358678 0.57566138   Down
# 20    LUAD 1.00000000        Inf     Up
# 21    LUSC 0.56302028 0.60234805   Down
# 22    LUSC 0.12555228 0.13706568     Up
# 23    PRAD 0.56578365 1.35231279   Down
# 24    PRAD 1.00000000        Inf     Up
# 25    STAD 1.00000000 0.69399693   Down
# 26    STAD 0.24869425 0.18922361     Up
# 27    THCA 0.27052628 1.78519374   Down
# 28    THCA 1.00000000        Inf     Up



# Down
# nRDiffExp %>% group_by(Estatus) %>% arrange(Direction) %>% distinct(DISEASE, Symbol, .keep_all = T) %>%
#   group_by(Estatus, Direction) %>% count()
# 

# > fisher.test(matrix(c(42, 291, 182, 1316), 2, byrow = T))$p.value
# [1] 0.853275
# > fisher.test(matrix(c(42, 291, 182, 1316), 2, byrow = T))$estimate
# odds ratio 
# 1.043595


# Up
# nRDiffExp %>% group_by(Estatus) %>% arrange(desc(Direction)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>%
#   group_by(Estatus, Direction) %>% count()

# > fisher.test(matrix(c(14, 199, 68, 1430), 2, byrow = T))$p.value
# [1] 0.2273149
# > fisher.test(matrix(c(14, 199, 68, 1430), 2, byrow = T))$estimate
# odds ratio 
# 1.479066


########################
####down
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% arrange(Direction) %>% distinct(Symbol, .keep_all = T) %>% count(Direction)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% arrange(Direction) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(Direction)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


####up
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% arrange(desc(Direction)) %>% distinct(Symbol, .keep_all = T) %>% count(Direction)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot3 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% arrange(desc(Direction)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(Direction)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot4 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


pdf('/result/Section4/regulonScnaEnrich.pdf')

plot1
plot2
plot3
plot4

dev.off()



################################################Rho
########################
load(file = '/data/genesGisticPeakQvalue.RData')
nrGisticPeakQvalue <- rhoGisticPeakQvalue
nrGisticPeakQvalue <- nrGisticPeakQvalue %>% as.data.frame() %>% subset(Class == 'Small_GTPase_Rho') %>% 
  select(Approved.symbol, Disease, Direction) %>% distinct()


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% subset(Class == 'Small_GTPase_Rho') %>% 
  select(Approved.symbol) %>% distinct() %>% mutate(Direction = 'Neutral')


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

nRDiffExp <- nRDiffExp %>% left_join(nrGisticPeakQvalue, by = join_by(DISEASE == Disease, Symbol == Approved.symbol)) %>% 
  mutate(Direction = factor(Direction, levels = c('Deletion', 'Neutral', 'Amplification')))


#############
fisherT <- function(data) {
  
  
  downData <- subset(data, Estatus %in% c('Down', 'Neutral')) %>% mutate(Direction = ifelse(Direction == 'Deletion', 'Deletion', 'Other'))
  downData <- downData %>% arrange(Direction) %>% distinct(Symbol, .keep_all = T)
  
  downtable <- table(downData$Estatus, downData$Direction)
  if(nrow(downtable) == 2 & ncol(downtable) == 2){
    
    res <- fisher.test(downtable)
    downRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Down')
    
  }else{
    
    downRes <- data.frame(pValue = NA, OR = NA, status = 'Down')
  }  
  
  

  upData <- subset(data, Estatus %in% c('Up', 'Neutral')) %>% mutate(Direction = ifelse(Direction == 'Amplification', 'Amplification', 'Other'))
  upData <- upData %>% arrange(desc(Direction)) %>% distinct(Symbol, .keep_all = T)
  
  uptable <- table(upData$Estatus, upData$Direction)
  if(nrow(uptable) == 2 & ncol(uptable) == 2){
    
    res <- fisher.test(uptable)
    upRes <- data.frame(pValue = res$p.value, OR = res$estimate, status = 'Up') 
    
  }else{
    
    upRes <- data.frame(pValue = NA, OR = NA, status = 'Up')
  }  
  
  return(rbind(downRes, upRes))
  
}


fisherPvalue <- nRDiffExp %>% group_by(DISEASE) %>% do(fisherT(.))


#    DISEASE    pValue       OR status
# 1     BLCA 0.5294118 0.000000   Down
# 2     BLCA        NA       NA     Up
# 3     BRCA 1.0000000 0.000000   Down
# 4     BRCA        NA       NA     Up
# 5     COAD 1.0000000 1.000000   Down
# 6     COAD        NA       NA     Up
# 7     ESCA 1.0000000 0.000000   Down
# 8     ESCA 1.0000000      Inf     Up
# 9     HNSC 1.0000000 0.000000   Down
# 10    HNSC        NA       NA     Up
# 11    KICH        NA       NA   Down
# 12    KICH        NA       NA     Up
# 13    KIRC 1.0000000 0.000000   Down
# 14    KIRC        NA       NA     Up
# 15    KIRP 0.5147059 4.491444   Down
# 16    KIRP        NA       NA     Up
# 17    LIHC 1.0000000 0.000000   Down
# 18    LIHC        NA       NA     Up
# 19    LUAD 1.0000000 0.000000   Down
# 20    LUAD 1.0000000      Inf     Up
# 21    LUSC 0.6043956 1.911066   Down
# 22    LUSC        NA       NA     Up
# 23    PRAD 1.0000000 0.000000   Down
# 24    PRAD 1.0000000      Inf     Up
# 25    STAD 1.0000000 0.000000   Down
# 26    STAD 1.0000000      Inf     Up
# 27    THCA        NA       NA   Down
# 28    THCA        NA       NA     Up


# Down
# nRDiffExp %>% group_by(Estatus) %>% arrange(Direction) %>% distinct(DISEASE, Symbol, .keep_all = T) %>%
#   group_by(Estatus, Direction) %>% count()


# > fisher.test(matrix(c(5, 43, 36, 172), 2, byrow = T))$p.value
# [1] 0.2823046
# > fisher.test(matrix(c(5, 43, 36, 172), 2, byrow = T))$estimate
# odds ratio 
# 0.5566739


# Up
nRDiffExp %>% group_by(Estatus) %>% arrange(desc(Direction)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>%
  group_by(Estatus, Direction) %>% count()

# > fisher.test(matrix(c(0, 24, 6, 202), 2, byrow = T))$p.value
# [1] 1
# > fisher.test(matrix(c(0, 24, 6, 202), 2, byrow = T))$estimate
# odds ratio 
# 0


########################
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% arrange(Direction) %>% distinct(Symbol, .keep_all = T) %>% count(Direction)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot1 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% arrange(Direction) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(Direction)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot2 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 


####up
dataStat <- nRDiffExp %>% group_by(DISEASE, Estatus) %>% arrange(desc(Direction)) %>% distinct(Symbol, .keep_all = T) %>% count(Direction)

dataStat <- dataStat %>% group_by(DISEASE, Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(dataStat, by = c('DISEASE', 'Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot3 <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~ DISEASE, nrow = 1) + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



totalStat <- nRDiffExp %>% group_by(Estatus) %>% arrange(desc(Direction)) %>% distinct(DISEASE, Symbol, .keep_all = T) %>% count(Direction)

totalStat <- totalStat %>% group_by(Estatus) %>% summarise(total = sum(n)) %>% 
  left_join(totalStat, by = c('Estatus')) %>% mutate(perc = n/total, Direction = factor(Direction, levels = c('Amplification', 'Neutral', 'Deletion')))


plot4 <- ggplot(totalStat, aes(x = Estatus, y = perc, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'NR genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) 



pdf('/result/Section4/rhoScnaEnrich.pdf')

plot1
plot2
plot3
plot4

dev.off()