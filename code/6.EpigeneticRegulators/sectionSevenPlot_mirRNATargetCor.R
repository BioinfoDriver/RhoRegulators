
library('dplyr')
library('ggplot2')
library('ggExtra')
library('Hmisc')
###############
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

load(file = '/data/panCanMiRnaExpData.RData')
# panCanTurMiRnaExp, panCanPairdTurMiRnaExp, panCanPairdNormMiRnaExp


load(file = '/data/tcgaExpSams.RData')
tagaExpData <- readRDS(file = '/data/tcgaFirehoseExpData.rds')
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
rhoFamilies <- rhoFamilies %>% distinct(Approved.symbol, .keep_all = TRUE) %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID))
rhoFamilies <- rhoFamilies %>% subset(NCBI.Gene.ID %in% rownames(tagaExpData$ACC$TPM))


tcgaRhoExp <- lapply(tagaExpData, function(expData){

  expData <- expData$TPM[rhoFamilies$NCBI.Gene.ID, ]
  rownames(expData) <- rhoFamilies$Approved.symbol

  return(expData)
})
tcgaRhoExp <- bind_cols(tcgaRhoExp)


panCanTurGeneExp <- tcgaRhoExp[, tumorExpSams]
panCanPairdTurGeneExp <- tcgaRhoExp[, paste(pairedSamsPats, '-01', sep = '')]
panCanPairdNormGeneExp <- tcgaRhoExp[, paste(pairedSamsPats, '-11', sep = '')]


###############
comSmas <- intersect(colnames(panCanTurMiRnaExp), colnames(panCanTurGeneExp))
turMiRnaExp <- panCanTurMiRnaExp[, comSmas]
turGeneExp <- panCanTurGeneExp[, comSmas]


tcgaPanCanSamples <- subset(tcgaPanCanSamples, SAMPLE_BARCODE %in% comSmas)


############
KIRCsamples <- subset(tcgaPanCanSamples, DISEASE == 'KIRC')
KIRCMiRnaExp <- turMiRnaExp['hsa-miR-142-5p', KIRCsamples$SAMPLE_BARCODE]
KIRCGeneExp <- turGeneExp['SYDE2', KIRCsamples$SAMPLE_BARCODE]


mirGeneExp <- as.data.frame(t(rbind.data.frame(KIRCMiRnaExp, KIRCGeneExp)))
colnames(mirGeneExp) <- c('miRNA', 'Target')


p2 <- ggplot(mirGeneExp, aes(x = miRNA, y = Target)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "Expression of hsa-miR-142-5p", y = "Expression of SYDE2") + 
  theme(legend.position = c(0.12,0.12), legend.background = element_blank(), legend.title = element_blank())


# > rcorr(as.matrix(mirGeneExp), type = "spearman")$r
# miRNA     Target
# miRNA   1.0000000 -0.3200288
# Target -0.3200288  1.0000000
# > rcorr(as.matrix(mirGeneExp), type = "spearman")$P
# miRNA       Target
# miRNA            NA 9.425465e-10
# Target 9.425465e-10           NA


############

comSmas <- Reduce(intersect, list(substr(colnames(panCanPairdTurMiRnaExp), 1, 12), substr(colnames(panCanPairdNormGeneExp), 1, 12),
                                  substr(colnames(panCanPairdTurGeneExp), 1, 12), substr(colnames(panCanPairdNormMiRnaExp), 1, 12)))

sams <- subset(tcgaPanCanSamples, DISEASE == 'KIRC' & PATIENT_BARCODE %in% comSmas)



tMirExp <- panCanPairdTurMiRnaExp["hsa-miR-142-5p", paste0(sams$PATIENT_BARCODE, '-01')]
colnames(tMirExp) <- substr(colnames(tMirExp), 1, 12)

tGeneExp <- panCanPairdTurGeneExp["SYDE2", paste0(sams$PATIENT_BARCODE, '-01')]
colnames(tGeneExp) <- substr(colnames(tGeneExp), 1, 12)

nMirExp <- panCanPairdNormMiRnaExp["hsa-miR-142-5p", paste0(sams$PATIENT_BARCODE, '-11')]
colnames(nMirExp) <- substr(colnames(nMirExp), 1, 12)

nGeneExp <- panCanPairdNormGeneExp["SYDE2", paste0(sams$PATIENT_BARCODE, '-11')]
colnames(nGeneExp) <- substr(colnames(nGeneExp), 1, 12)


plotdata <- bind_rows(as.data.frame(t(rbind.data.frame(tMirExp,tGeneExp))) %>% mutate(Type = 'Tumor'),
                  as.data.frame(t(rbind.data.frame(nMirExp,nGeneExp))) %>% mutate(Type = 'Normal'))
colnames(plotdata) <- c('miRNA', 'Target', 'Type')



p0 <- ggplot(plotdata,
             aes(miRNA, Target, colour = Type)) +
  geom_point() + theme_bw() + 
  scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set2')[1:2]) + 
  theme(legend.position = c(0.12,0.12), 
        legend.background = element_blank(), 
        legend.title = element_blank()) + labs(x = "Expression of hsa-miR-142-5p", y = "Expression of SYDE2")

p1 <- ggMarginal(p0, type = "boxplot", groupColour = TRUE, groupFill = TRUE)

# > wilcox.test(subset(plotdata, Type == 'Tumor')$miRNA, subset(plotdata, Type == 'Normal')$miRNA, paired = TRUE)$p.value
# [1] 5.112907e-11
# > wilcox.test(subset(plotdata, Type == 'Tumor')$Target, subset(plotdata, Type == 'Normal')$Target, paired = TRUE)$p.value
# [1] 1.125163e-09

############

pdf(file = '/result/Section4/SYDE2_mir142_KIRC_cor.pdf')
p1
p2
dev.off()



