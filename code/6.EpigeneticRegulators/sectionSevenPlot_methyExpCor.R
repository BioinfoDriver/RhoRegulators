
library('dplyr')
library('ggplot2')
library('ggExtra')
library('Hmisc')
library('RColorBrewer')
##########
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

load(file = '/data/panCanMethyData.RData')
# panCanTurMethy, panCanPairdTurMethy, panCanPairdNormMethy


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

##############

anno450k <- readRDS(file = '/data/methAnno450k.rds')

nrAnno450k <- subset(anno450k, RefGene_Name %in% 'PREX2')
nrpromoterAnno <- subset(nrAnno450k, RefGene_Group %in% c("5'UTR", "TSS1500", "TSS200"))


panCanTurMethy <- rownames_to_column(panCanTurMethy, var = 'Name')
TMethy <- panCanTurMethy %>% inner_join(nrpromoterAnno[, c('Name', 'RefGene_Name')], by = 'Name')

TMethy <- TMethy[!duplicated(TMethy$Name), ]
TMethy <- TMethy %>%  select(-RefGene_Name) %>% remove_rownames() %>% column_to_rownames(var = 'Name')  
TMethy <- colMeans(TMethy) %>% as.data.frame() %>% t()%>% as.data.frame()

rownames(TMethy) <- 'Methy'



panCanPairdNormMethy <- rownames_to_column(panCanPairdNormMethy, var = 'Name')
NMethy <- panCanPairdNormMethy %>% inner_join(nrpromoterAnno[, c('Name', 'RefGene_Name')], by = 'Name')


NMethy <- NMethy[!duplicated(NMethy$Name), ]
NMethy <- NMethy %>%  select(-RefGene_Name) %>% remove_rownames() %>% column_to_rownames(var = 'Name')  
NMethy <- colMeans(NMethy) %>% as.data.frame() %>% t()%>% as.data.frame()


rownames(NMethy) <- 'Methy'

TExp <- panCanTurGeneExp['PREX2', ]
NExp <- panCanPairdNormGeneExp['PREX2', ]

rownames(TExp) <- 'Exp'
rownames(NExp) <- 'Exp'


diseaseSamples <- subset(tcgaPanCanSamples, DISEASE == 'LUAD')

##############

comSmas <- intersect(intersect(colnames(TExp), colnames(TMethy)), diseaseSamples$SAMPLE_BARCODE)


methyGeneExp <- as.data.frame(t(rbind.data.frame(TMethy[, comSmas], TExp[, comSmas])))



p2 <- ggplot(methyGeneExp, aes(x = Methy, y = Exp)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "Methylation of PREX2", y = "Expression of PREX2") + 
  theme(legend.position = c(0.12,0.12), legend.background = element_blank(), legend.title = element_blank())

# > rcorr(as.matrix(methyGeneExp), type = "spearman")$r
# Methy        Exp
# Methy  1.0000000 -0.1836918
# Exp   -0.1836918  1.0000000
# > rcorr(as.matrix(methyGeneExp), type = "spearman")$P
# Methy          Exp
# Methy           NA 3.525536e-05
# Exp   3.525536e-05           NA


############
comSmas <- Reduce(intersect, list(substr(colnames(TMethy), 1, 12), substr(colnames(NMethy), 1, 12),
                                  substr(colnames(TExp), 1, 12), substr(colnames(NExp), 1, 12)))

sams <- subset(tcgaPanCanSamples, DISEASE == 'LUAD' & PATIENT_BARCODE %in% comSmas)



tGeneMethy <- TMethy[, paste0(sams$PATIENT_BARCODE, '-01')]
colnames(tGeneMethy) <- substr(colnames(tGeneMethy), 1, 12)

tGeneExp <- TExp[, paste0(sams$PATIENT_BARCODE, '-01')]
colnames(tGeneExp) <- substr(colnames(tGeneExp), 1, 12)

nGeneMethy <- NMethy[, paste0(sams$PATIENT_BARCODE, '-11')]
colnames(nGeneMethy) <- substr(colnames(nGeneMethy), 1, 12)

nGeneExp <- NExp[, paste0(sams$PATIENT_BARCODE, '-11')]
colnames(nGeneExp) <- substr(colnames(nGeneExp), 1, 12)


plotdata <- rbind(as.data.frame(t(rbind.data.frame(tGeneMethy,tGeneExp))) %>% mutate(Type = 'Tumor'),
                  as.data.frame(t(rbind.data.frame(nGeneMethy,nGeneExp))) %>% mutate(Type = 'Normal'))



p0 <- ggplot(plotdata,
             aes(Methy, Exp, colour = Type)) +
  geom_point() + theme_bw() + 
  scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set2')[1:2]) + 
  theme(legend.position = c(0.8, 0.80), 
        legend.background = element_blank(), 
        legend.title = element_blank()) + labs(x = "Methylation of PREX2", y = "Expression of PREX2")

p1 <- ggMarginal(p0, type = "boxplot", groupColour = TRUE, groupFill = TRUE)


# > wilcox.test(subset(plotdata, Type == 'Tumor')$Methy, subset(plotdata, Type == 'Normal')$Methy, paired = TRUE)$p.value
# [1] 5.340576e-05
# > wilcox.test(subset(plotdata, Type == 'Tumor')$Exp, subset(plotdata, Type == 'Normal')$Exp, paired = TRUE)$p.value
# [1] 0.000164032

############

pdf(file = '/result/Section4/LUAD_PREX2_methy_cor.pdf')
p1
p2
dev.off()





