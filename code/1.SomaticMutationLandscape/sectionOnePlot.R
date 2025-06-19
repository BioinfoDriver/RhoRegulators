
library('ggplot2')
library('reshape2')

#########################
gainLoFMutStat <- readRDS(file = '/data/gainLoFMutStat.rds')

# size
gainLoFMutStat$sizeLabel <- cut(gainLoFMutStat$numOfNonsMut, breaks = c(0, 36, 54, 86, 600), 
                                labels = c('n ≤ 36', '36 < n ≤ 54', '54 < n ≤ 86', 'n ≥ 86'))
# shape
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('DH_domain', 'DOCKER_dom')] <- 'GEF'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('RhoGAP_dom')] <- 'GAP'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('Rho_GDI')] <- 'GDI'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('Small_GTPase_Rho')] <- 'RHO'
gainLoFMutStat$Class[gainLoFMutStat$Class %in% c('RhoGAP_DH_dom')] <- 'GAP_GEF'


# LoF OR Hotspot
gainLoFMutStat$colLabel <- ifelse(gainLoFMutStat$gainMf > 0.3 & gainLoFMutStat$lofMf < 0.2 & gainLoFMutStat$numOfHotMutPosi > 5, 'HotspotMutG', 
                                  ifelse(gainLoFMutStat$gainMf < 0.2 & gainLoFMutStat$lofMf > 0.3 & gainLoFMutStat$numOfLofMut > 8, 'LofMutG', 'Other'))


gainLoFMutStat <- tibble::rownames_to_column(gainLoFMutStat, var = "GeneName")

colLabel <- c("#B7DBE3", "#DFE1E2", "#F5D9E6")
names(colLabel) <- c('LofMutG', 'OtherG', 'HotspotMutG')

sizeLabel <- c(1.5, 3, 4.5, 6)
names(sizeLabel) <- c('n ≤ 36', '36 < n ≤ 54', '54 < n ≤ 86', 'n ≥ 86')

shapeLabel <- c(14, 15, 16, 17, 18)
names(shapeLabel) <- c('GDI', 'GEF', 'GAP', 'GAP_GEF', 'RHO')


plot1 <- ggplot(gainLoFMutStat, aes(x = lofMf, y = gainMf, color = colLabel, size = sizeLabel)) +
  geom_point(aes(shape = Class)) + labs(x = "Fraction of LoF mutations", y = "Fracton of Hotspot mutations") + 
  scale_colour_manual(values = colLabel, guide = "legend", name = "Color") + 
  scale_size_manual(values = sizeLabel, guide = "legend", name = "Size") + 
  scale_shape_manual(values = shapeLabel, guide = "legend", name = "Shape") + 
  geom_text(data = subset(gainLoFMutStat, colLabel %in% c('LofMutG', 'HotspotMutG')),
            size = 3, check_overlap = TRUE, aes(label = GeneName), nudge_y = 0.02) + theme(panel.background = element_blank()) + 
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0, 0.2, 0.3, 0.5)) + 
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.3, 0.5, 0.7))


#########################
load(file = '/data/rhoFamiliesMutSig2CV.RData')
load(file = '/data/rhoFamiliesMutStat.RData')


rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies$Class[rhoFamilies$Class %in% c('DH_domain', 'DOCKER_dom')] <- 'GEF'
rhoFamilies$Class[rhoFamilies$Class %in% c('RhoGAP_dom')] <- 'GAP'
rhoFamilies$Class[rhoFamilies$Class %in% c('Rho_GDI')] <- 'GDI'
rhoFamilies$Class[rhoFamilies$Class %in% c('Small_GTPase_Rho')] <- 'RHO'


mutSigGeneQvalue <- melt(mutSigGeneQvalue, id.vars = 'gene', value.name = 'qValue', variable.name = 'cancerType')
mutSigGeneQvalue <- mutSigGeneQvalue %>% mutate(cancerType = gsub('qvalue', '', cancerType))


sigMut <- subset(mutSigGeneQvalue, qValue < 0.25)
sigMut <- rename(sigMut, Hugo_Symbol = gene, DISEASE = cancerType)

#########
rhoFamiliesMutStatByGene$DISEASE = 'PANCAN'
rhoFamiliesMutStatByGene$numOfPat = 8217
rhoFamiliesMutStatByGene <- rename(rhoFamiliesMutStatByGene, numOfPatWithMut = numOfMut)
rhoFamiliesMutStatByGene <- rhoFamiliesMutStatByGene[, colnames(rhoFamiliesMutStatByDisGene)]

nrMutSat <- rbind.data.frame(rhoFamiliesMutStatByGene, rhoFamiliesMutStatByDisGene)

#########
sigMut <- merge(sigMut, nrMutSat, by = c('Hugo_Symbol', 'DISEASE'))

sigMut <- merge(sigMut, rhoFamilies[, c('Approved.symbol', 'Class')], by.x = 'Hugo_Symbol', by.y = 'Approved.symbol')


# size
sigMut$sizeLable <- cut(sigMut$qValue, breaks = c(0, 1.0e-5, 1.0e-02, 5.0e-02, 0.25), right = FALSE, 
                        labels = c('q < 1.0e-5','q < 1.0e-02', 'q < 0.05', 'q < 0.25'))


sizeLable <- c(2, 4, 5, 6, 7)
names(sizeLable) <- rev(c('q < 1.0e-5','q < 1.0e-02', 'q < 0.05', 'q < 0.25'))

shapeLabel <- c(14, 15, 16, 18)
names(shapeLabel) <- c('GDI', 'GEF', 'GAP', 'RHO')


sigMut <- sigMut %>% group_by(Class) %>% arrange(qValue)
sigMut <- sigMut %>% mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = unique(Hugo_Symbol)), 
                            DISEASE = factor(DISEASE, 
                                             levels = c("PANCAN", "BLCA", "BRCA", "GBM", "HNSC", "LGG", "SKCM", "STAD", "UCEC", "UCS")))


plot2 <- ggplot(sigMut, aes(x = DISEASE, y = Hugo_Symbol, color = altFre, size = sizeLable)) +
  geom_point(aes(shape = Class)) +
  labs(x = "Cancer Types", y = "Significantly mutated genes") +
  scale_size_manual(values = sizeLable, guide = "legend", name = "MutSigCV q value") +
  scale_color_gradient(low = '#ECC38C', high = '#BD4146', name = 'Frequency') +
  scale_shape_manual(values = shapeLabel, guide = "legend", name = "Type") +
  theme(axis.text.x = element_text(angle = 60, hjust = 0.1, vjust = 0.5))

 
pdf(file = '/result/Section1/Figure1.1.pdf')
plot1
plot2
dev.off()


