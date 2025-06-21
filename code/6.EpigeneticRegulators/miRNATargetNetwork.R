

library(dplyr)
#########################

miRnaTarByDisease <- readRDS(file = '/data/miRnaTarByDiseaseTargetScan.rds')
miRnaTarByDisease <- miRnaTarByDisease %>% mutate(miRNA = gsub('hsa-', '', miRNA))


network <- miRnaTarByDisease %>% group_by(miRNA, TargetGene) %>% count(name = 'width')

miRnaDegree <- miRnaTarByDisease %>% group_by(miRNA) %>% count(name = 'degree') %>% mutate(nodeType = 1, borderType = 1) %>% dplyr::rename(node = miRNA)
targetDegree <- miRnaTarByDisease %>% group_by(TargetGene) %>% count(name = 'degree') %>% mutate(nodeType = 2) %>% dplyr::rename(node = TargetGene)


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% distinct(Approved.symbol, .keep_all = TRUE) %>% select(Approved.symbol, Class) %>% 
  mutate(Class = recode(Class, DH_domain = 2, DOCKER_dom = 2, Rho_GDI = 3, RhoGAP_dom = 4, Small_GTPase_Rho = 5)) %>% dplyr::rename(borderType = Class)

targetDegree <- targetDegree %>% left_join(rhoFamilies, by = join_by(node == Approved.symbol))


nodeAttributes <- rbind.data.frame(miRnaDegree, targetDegree)

networkFilter <- subset(network, width > 3)


setwd('/result/Section4/')
write.table(network, file = 'miRnaTargerNetworkTargetScan.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(nodeAttributes, file = 'networkNodeAttributesTargetScan.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(networkFilter, file = 'miRnaTargerNetworkTargetScanFilter.txt', sep = '\t', col.names = T, row.names = F, quote = F)
