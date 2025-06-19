
library(dplyr)

#########################
# panCanSurPvalue <- readRDS(file = '/data/panCanRhoSurPvalue.rds')
panCanSurPvalue <-readRDS(file = '/data/panCanRhoSurUniMulLogPvalue.rds')

rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

rhoFamilies <- rhoFamilies %>% mutate(Class = recode(Class, DH_domain = 'GEF', DOCKER_dom = 'GEF', Rho_GDI = 'GDI', RhoGAP_dom = 'GAP', Small_GTPase_Rho = 'Rho'))

panCanSurPvalue <-merge.data.frame(panCanSurPvalue, rhoFamilies[, c('Approved.symbol', 'Class')], by.x = 'geneSymbol', by.y = 'Approved.symbol')


###########
selectSurPvalue <- subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05)
selectSurPvalue <- selectSurPvalue %>% mutate(hrDirect = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival'))

###########
setwd('/result/Section3')

# write.table(selectSurPvalue, sep = '\t', quote = F, row.names = F, col.names = T, file = 'rhoSurPvalue.txt')
write.table(selectSurPvalue, sep = '\t', quote = F, row.names = F, col.names = T, file = 'rhoSurUniMulLogPvalue.txt')


###########
disDysGCount <- selectSurPvalue %>% distinct(disease, geneSymbol, .keep_all = TRUE) %>% group_by(disease, hrDirect) %>% count(name = 'num')

write.table(disDysGCount, sep = '\t', quote = F, row.names = F, col.names = T, file = 'rhoSurAssoStatByDis.txt')


###########
totalSurAsso <- selectSurPvalue %>% group_by(Class, geneSymbol) %>% count(name = 'Total')
totalSurAsso <- selectSurPvalue %>% group_by(Class, geneSymbol, hrDirect) %>% count(name = 'Number') %>% 
  left_join(totalSurAsso, by = join_by(Class, geneSymbol)) %>% ungroup() %>%  group_by(Class) %>% arrange(desc(Total))


write.table(totalSurAsso, sep = '\t', quote = F, row.names = F, col.names = T, file = 'rhoSurAssoStat.txt')


#########################
panCanEnrichScoreSurPvalue <- readRDS(file = '/data/panCanRhoEnrichScoreSurPvalue.rds')

enrichScoreSurAsso <- panCanEnrichScoreSurPvalue %>% subset(coxPvalue < 0.05 & logrankPvalue < 0.05)

write.table(enrichScoreSurAsso, sep = '\t', quote = F, row.names = F, col.names = T, file = 'ssGseaScoreSurAsso.txt')


#########################
panCanSurPvalue <- readRDS(file = '/data/panCanRhoSurPvalue_GPL570.rds')
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

rhoFamilies <- rhoFamilies %>% mutate(Class = recode(Class, DH_domain = 'GEF', DOCKER_dom = 'GEF', Rho_GDI = 'GDI', RhoGAP_dom = 'GAP', Small_GTPase_Rho = 'Rho'))

panCanSurPvalue <-merge.data.frame(panCanSurPvalue, rhoFamilies[, c('Approved.symbol', 'Class')], by.x = 'geneSymbol', by.y = 'Approved.symbol') %>% 
  mutate(disease = stringr::str_to_title(disease))

###########
selectSurPvalue <- subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05)

setwd('/result/Section3')

write.table(selectSurPvalue, sep = '\t', quote = F, row.names = F, col.names = T, file = 'rhoSurPvalue_validate.txt')
