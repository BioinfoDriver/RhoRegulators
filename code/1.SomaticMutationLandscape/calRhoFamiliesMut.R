
#################
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')

mc3MutData$Hugo_Symbol[mc3MutData$Hugo_Symbol == 'HMHA1'] <- 'ARHGAP45'
mc3MutData$Hugo_Symbol[mc3MutData$Hugo_Symbol == 'C20orf95'] <- 'ARHGAP40' # missing


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

rhoFamilies$Class[rhoFamilies$Class %in% c('DH_domain', 'DOCKER_dom')] <- 'GEF'
rhoFamilies$Class[rhoFamilies$Class %in% c('RhoGAP_dom')] <- 'GAP'
rhoFamilies$Class[rhoFamilies$Class %in% c('Rho_GDI')] <- 'GDI'
rhoFamilies$Class[rhoFamilies$Class %in% c('Small_GTPase_Rho')] <- 'RHO'




rhoFamiliesMutData <- merge(mc3MutData, rhoFamilies, by.x = 'Hugo_Symbol', by.y = 'Approved.symbol')
rhoFamiliesMutData <- subset(rhoFamiliesMutData, 
                            Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',  'In_Frame_Ins', 
                                                          'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))

# Stat by patient
library('dplyr')
rhoFamiliesMutStatByDis <- rhoFamiliesMutData %>% group_by(DISEASE, SAMPLE_BARCODE) %>% 
  count(name = 'numOfMut') %>% group_by(DISEASE) %>% count(name = 'numOfPatWithMut')

rhoFamiliesMutStatByDisGroup <- rhoFamiliesMutData %>% group_by(DISEASE, Class, SAMPLE_BARCODE) %>% 
  count(name = 'numOfMut')  %>% group_by(DISEASE, Class)  %>% count(name = 'numOfPatWithMut')


disCount <- unique(rhoFamiliesMutData[, c('DISEASE', 'numOfPat')])

rhoFamiliesMutStatByDis <- merge(rhoFamiliesMutStatByDis, disCount, by = 'DISEASE')
rhoFamiliesMutStatByDisGroup <- merge(rhoFamiliesMutStatByDisGroup, disCount, by = 'DISEASE')

rhoFamiliesMutStatByDis <- rhoFamiliesMutStatByDis %>% mutate(altFre = numOfPatWithMut/numOfPat)
rhoFamiliesMutStatByDisGroup <- rhoFamiliesMutStatByDisGroup %>% mutate(altFre = numOfPatWithMut/numOfPat)



# Stat by gene
rhoFamiliesMutStatByDisGene <- rhoFamiliesMutData %>% group_by(DISEASE, SAMPLE_BARCODE, Hugo_Symbol) %>% count(name = 'numOfMut') %>% 
  group_by(DISEASE, Hugo_Symbol) %>% count(name = 'numOfPatWithMut') 


rhoFamiliesMutStatByDisGene <- merge(rhoFamiliesMutStatByDisGene, disCount, by = 'DISEASE')
rhoFamiliesMutStatByDisGene <- rhoFamiliesMutStatByDisGene %>% mutate(altFre = numOfPatWithMut/numOfPat)


rhoFamiliesMutStatByGene <- rhoFamiliesMutData %>% group_by(DISEASE, SAMPLE_BARCODE, Hugo_Symbol) %>% count(name = 'numOfMut') %>% 
  group_by(Hugo_Symbol) %>% count(name = 'numOfMut') 

rhoFamiliesMutStatByGene$altFre <- rhoFamiliesMutStatByGene$numOfMut/sum(disCount$numOfPat)


save(rhoFamiliesMutStatByDis, rhoFamiliesMutStatByDisGroup, rhoFamiliesMutStatByDisGene, rhoFamiliesMutStatByGene, file = '/data/rhoFamiliesMutStat.RData')



