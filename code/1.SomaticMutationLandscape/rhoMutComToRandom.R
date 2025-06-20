library('dplyr')

#########################
load(file = '/data/randomGenesList.RData')
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')
load(file = '/data/rhoFamiliesMutStat.RData')

# by disease
rgMutStatByDis <- lapply(randomGenesList, function(randomGenes){
  rgMutData <- merge(mc3MutData, randomGenes, by.x = 'Hugo_Symbol', by.y = 'gene_name')
  
  rgMutData <- subset(rgMutData, Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 
                                                               'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))
  
  rgMutStatByDis <- rgMutData %>% group_by(DISEASE, SAMPLE_BARCODE) %>% 
    count(name = 'numOfMut')  %>% group_by(DISEASE)  %>% count(name = 'numOfPatWithMut')
  
  return(rgMutStatByDis)
  
})

rgMutStatByDis <- Reduce(function(x, y) merge(x = x, y = y, by = 'DISEASE', all = TRUE), rgMutStatByDis)
colnames(rgMutStatByDis) <- c('DISEASE', paste0('rgRes', 1:10000))

comResByDis <- merge(rhoFamiliesMutStatByDis, rgMutStatByDis, by = 'DISEASE')

comResByDisPvalue <- sapply(seq(nrow(comResByDis)), function(i){
  
  res <- sum(comResByDis[i, 2] > comResByDis[i, 5:10004], na.rm = TRUE)
  return(res)
  
})
names(comResByDisPvalue) <- comResByDis$DISEASE
comResByDisPvalue <- 1 - comResByDisPvalue/10000


# by disease/group
rgMutStatListByDisGroup <- lapply(randomGenesListByGroup, function(rgGenesByGroup){
  
  rgMutStatByDisGroup <- lapply(rgGenesByGroup, function(rgGenes){
    
    rgMutData <- merge(mc3MutData, rgGenes, by.x = 'Hugo_Symbol', by.y = 'gene_name')
    
    rgMutData <- subset(rgMutData, Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 
                                                                 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))
    
    rgMutStatByDis <- rgMutData %>% group_by(DISEASE, SAMPLE_BARCODE) %>% 
      count(name = 'numOfMut')  %>% group_by(DISEASE)  %>% count(name = 'numOfPatWithMut')
    
    return(rgMutStatByDis)    
    
  })
  
  rgMutStatByDisGroup <- Reduce(function(x, y) merge(x = x, y = y, by = 'DISEASE', all = TRUE), rgMutStatByDisGroup)
  colnames(rgMutStatByDisGroup) <- c('DISEASE', paste0('rgRes', 1:10000))
  
  return(rgMutStatByDisGroup)
})


groupNames <- names(table(rhoFamiliesMutStatByDisGroup$Class))


comResByDisGroupPvalue <- lapply(seq(length(groupNames)), function(i){
  
  gName <- groupNames[i]
  nrMutStatByDisGroup <- subset(rhoFamiliesMutStatByDisGroup, Class == gName)
  
  
  comResByDisGroup <- merge(nrMutStatByDisGroup, rgMutStatListByDisGroup[[i]], by = 'DISEASE', all = TRUE)
  comResByDisGroup[is.na(comResByDisGroup)] <- 0
  
  comResByDisGroupPvalue <- sapply(seq(nrow(comResByDisGroup)), function(i){
    
    res <- sum(comResByDisGroup[i, 3] > comResByDisGroup[i, 6:10005])
    return(res)
    
  })
  names(comResByDisGroupPvalue) <- comResByDisGroup$DISEASE  
  comResByDisGroupPvalue <- 1 - comResByDisGroupPvalue/10000
  
  return(comResByDisGroupPvalue)
})

comResByDisGroupPvalue <- do.call(cbind, comResByDisGroupPvalue)
colnames(comResByDisGroupPvalue) <- groupNames


save(comResByDisPvalue, comResByDisGroupPvalue, file = '/data/rhoMutComToRandom.RData')


