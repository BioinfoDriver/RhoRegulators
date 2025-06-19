
library(dplyr)
#########################
geneInfo <- read.csv(file = '/data/gene_with_protein_product.txt', sep = '\t', header = T, stringsAsFactors = F)
geneInfo <- geneInfo %>% select(symbol, hgnc_id, name, alias_symbol, prev_symbol, entrez_id, ensembl_gene_id, vega_id, location)



rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

rhoFamilies$Class[rhoFamilies$Class %in% c('DH_domain', 'DOCKER_dom')] <- 'GEF'
rhoFamilies$Class[rhoFamilies$Class %in% c('RhoGAP_dom')] <- 'GAP'
rhoFamilies$Class[rhoFamilies$Class %in% c('Rho_GDI')] <- 'GDI'
rhoFamilies$Class[rhoFamilies$Class %in% c('Small_GTPase_Rho')] <- 'RHO'

rhoFamilies <- rhoFamilies %>% select(-Locus.type) %>% left_join(geneInfo, 
                                                                 by = join_by(Approved.symbol == symbol, NCBI.Gene.ID == entrez_id))

write.table(rhoFamilies, file = '/result/Section1/rhoFamilies.txt', row.names = F, col.names = T, quote = F, sep = '\t')



#########################
gainLoFMutStat <- readRDS(file = '/data/gainLoFMutStat.rds')

gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('DH_domain', 'DOCKER_dom')] <- 'GEF'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('RhoGAP_dom')] <- 'GAP'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('Rho_GDI')] <- 'GDI'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('Small_GTPase_Rho')] <- 'RHO'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('RhoGAP_DH_dom')] <- 'GAP_GEF'


write.table(gainLoFMutStat, file = '/result/Section1/gainLoFMutStat.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)


#########################
load(file = '/data/rhoFamiliesMutSig2CV.RData')
load(file = '/data/rhoFamiliesMutStat.RData')


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies$Class[rhoFamilies$Class %in% c('DH_domain', 'DOCKER_dom')] <- 'GEF'
rhoFamilies$Class[rhoFamilies$Class %in% c('RhoGAP_dom')] <- 'GAP'
rhoFamilies$Class[rhoFamilies$Class %in% c('Rho_GDI')] <- 'GDI'
rhoFamilies$Class[rhoFamilies$Class %in% c('Small_GTPase_Rho')] <- 'RHO'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('RhoGAP_DH_dom')] <- 'GAP_GEF'


sigMut <- melt(mutSigGeneQvalue, id.vars = 'gene', value.name = 'qValue', variable.name = 'cancerType') %>% 
  mutate(cancerType = gsub('qvalue', '', cancerType)) %>% rename(Hugo_Symbol = gene, DISEASE = cancerType)

#########
rhoFamiliesMutStatByGene$DISEASE = 'PANCAN'
rhoFamiliesMutStatByGene$numOfPat = 8217
rhoFamiliesMutStatByGene <- rename(rhoFamiliesMutStatByGene, numOfPatWithMut = numOfMut)
rhoFamiliesMutStatByGene <- rhoFamiliesMutStatByGene[, colnames(rhoFamiliesMutStatByDisGene)]

nrMutSat <- rbind.data.frame(rhoFamiliesMutStatByGene, rhoFamiliesMutStatByDisGene)


#########
sigMut <- merge(sigMut, nrMutSat, by = c('Hugo_Symbol', 'DISEASE')) # %>% arrange(desc(altFre))

sigMut <- merge(sigMut, rhoFamilies[, c('Approved.symbol', 'Class')], by.x = 'Hugo_Symbol', by.y = 'Approved.symbol')

write.table(sigMut, file = '/result/Section1/rhoGTPasesMutFreq.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)

