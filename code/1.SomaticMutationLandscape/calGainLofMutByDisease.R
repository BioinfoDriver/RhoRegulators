library(dplyr)

#################
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')

mc3MutData$Hugo_Symbol[mc3MutData$Hugo_Symbol == 'HMHA1'] <- 'ARHGAP45'
mc3MutData$Hugo_Symbol[mc3MutData$Hugo_Symbol == 'C20orf95'] <- 'ARHGAP40' # missing


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% mutate(Class = ifelse(Approved.symbol %in% c('BCR', 'ABR'), 'RhoGAP_DH_dom', Class))%>% distinct()


rhoFamiliesMutData <- merge(mc3MutData, rhoFamilies, by.x = 'Hugo_Symbol', by.y = 'Approved.symbol')

rhoFamiliesMutData <- subset(rhoFamiliesMutData, 
                             Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
                                                           'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))


rhoFamiliesMutDataByGene <- split.data.frame(rhoFamiliesMutData, f = ~ DISEASE + Hugo_Symbol, drop = TRUE)


# Across the pan-cancer cohort, hotspot mutations were defined as missense or in-frame mutations
# at the same protein amino acid in > 2 patient samples.
# The fraction of hotspot mutations per gene was calculated as the total number of hotspot mutations
# over the total number of non-silent mutations found in that gene.

# The fraction of LoF mutations (defined as Frame_Shift_Ins, Frame_Shift_Del, Nonsense_Mutation, Nonstop_Mutation,
# Splice_Site, and Tanslation_Start_Site) per gene was calculated as the total number of LoF mutations over the total number of
# non-silent mutations in that gene.
#
# Genes with > 30% hotspot mutations, < 20% LoF mutations, and ≥ 5 unique hotspot mutation
# positions were identified as enriched with hotspot mutations, while genes containing > 30% LoF mutations, < 30% hotspot
# mutations, and ≥ 10 LoF mutations were identified as enriched with LoF mutations.

fracOfLoF <- lapply(rhoFamiliesMutDataByGene, function(geneMutData){
  
  lofMut <- subset(geneMutData, Variant_Classification %in% c('Frame_Shift_Ins', 'Frame_Shift_Del', 'Nonsense_Mutation', 'Nonstop_Mutation',
                                                              'Splice_Site', 'Translation_Start_Site'))
  
  numOfLofMut <- nrow(lofMut)
  numOfNonsMut <- nrow(geneMutData)
  lofMf <- numOfLofMut/numOfNonsMut
  
  res <- data.frame(numOfLofMut=numOfLofMut, 
                    lofMf = lofMf, 
                    numOfNonsMut = numOfNonsMut, 
                    Hugo_Symbol = geneMutData$Hugo_Symbol[1], 
                    DISEASE = geneMutData$DISEASE[1])
  return(res)
})

fracOfLoF <- do.call(rbind, fracOfLoF)



fracOfGainMut <- lapply(rhoFamiliesMutDataByGene, function(geneMutData){
  
  gainMut <- subset(geneMutData, Variant_Classification %in% c('In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation'))
  
  gainMut$mutPos <- gainMut$HGVSp
  gainMut$mutPos <- sapply(strsplit(gainMut$mutPos, split = '_'), function(x) x[1])
  gainMut$mutPos <- gsub("[A-Za-z]+$", '', gainMut$mutPos)
  
  
  mutPosi <- table(gainMut$mutPos)
  numOfgainMut <- sum(mutPosi[mutPosi >= 2])
  numOfHotMutPosi <- length(mutPosi[mutPosi >= 2])
  
  numOfNonsMut <- nrow(geneMutData)
  
  gainMf <- numOfgainMut/numOfNonsMut
  
  res <- data.frame(numOfgainMut = numOfgainMut, 
                    numOfHotMutPosi = numOfHotMutPosi, 
                    gainMf = gainMf, 
                    numOfNonsMut = numOfNonsMut,
                    Hugo_Symbol = geneMutData$Hugo_Symbol[1], 
                    DISEASE = geneMutData$DISEASE[1])
  
  return(res)
})

fracOfGainMut <- do.call(rbind, fracOfGainMut)


gainLoFMutStat <- merge.data.frame(fracOfGainMut, fracOfLoF, by = c('Hugo_Symbol', 'DISEASE', 'numOfNonsMut'))

gainLoFMutStat$Class <- rhoFamilies$Class[match(gainLoFMutStat$Hugo_Symbol,  rhoFamilies$Approved.symbol)]


saveRDS(gainLoFMutStat, file = '/data/gainLoFMutStatByDisease.rds')

