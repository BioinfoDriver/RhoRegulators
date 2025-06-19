library(pheatmap)
library(tibble)
library(ggplot2)
library(ggpubr)
library(scatterpie)

#########################
load(file = '/data/tcgaRhoFamiliesDiffExp.RData') 
# deseq2DiffExp, limmaDiffExp

deseq2DiffExp <- deseq2DiffExp %>% mutate(rhoGroup = ifelse(Class %in% c('DH_domain', 'DOCKER_dom'), 'DH_DOCKER_dom', Class))

rhoGroupDysStat <- deseq2DiffExp %>% group_by(DISEASE, rhoGroup) %>% count(dysLabel)
rhoDysStat <- deseq2DiffExp %>% group_by(rhoGroup, Symbol) %>% count(dysLabel)

topDysGenes <- rhoDysStat %>% subset(dysLabel != 'Neutral') %>% group_by(rhoGroup, Symbol) %>% summarise(n = sum(n)) %>% 
  group_by(rhoGroup) %>% arrange(desc(n)) %>% slice_head(n = 6)

#########################
groupDysPiePlot <- rhoGroupDysStat %>% left_join(data.frame(DISEASE = unique(rhoGroupDysStat$DISEASE), x_locus = 1:14), by = 'DISEASE') %>% 
  left_join(data.frame(rhoGroup = unique(rhoGroupDysStat$rhoGroup), y_locus = 4:1), by = 'rhoGroup') %>% 
  mutate(dysLabel = factor(dysLabel, levels = c('Down', 'Neutral', 'Up')), value = n) %>% ungroup()


plot1 <- ggplot() + geom_scatterpie(aes(x=x_locus, y=y_locus), 
                           data=groupDysPiePlot, cols="dysLabel", long_format=TRUE) + coord_fixed() + 
  scale_fill_manual(values = c('Down' = 'blue', 'Neutral' = 'white', 'Up' = "red"), name = '') +
  scale_x_continuous(breaks = 1:14, labels = unique(groupDysPiePlot$DISEASE), ) + 
  scale_y_continuous(breaks = 1:4, labels = c('Rho', 'GDI', 'GAP' ,'GEF')) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 0.8), axis.text.y = element_text(angle = 60, vjust = -0.05))

#########################
rhoDysStatByDis <- deseq2DiffExp %>% group_by(DISEASE) %>% distinct(Symbol, .keep_all = T) %>% group_by(DISEASE) %>% count(dysLabel)

plot1.1 <- ggplot(rhoDysStatByDis, aes(x = DISEASE, y = n, fill = dysLabel)) +
  geom_bar(stat = "identity", position = "stack") + scale_y_continuous(breaks = seq(0, 175, 25), position = "left") +
  labs(y = 'Number of Dyregulated Genes') + scale_fill_manual(values = c("Down" = "blue", "Neutral" = 'white', "Up" = "red"), name = '') +
  theme(axis.title.x = element_blank(), legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black")) + 
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) 

ggsave(plot1.1, filename = '/result/Section3/rhoDysStatByDis.pdf')

#########################
logFcHeatmap <- deseq2DiffExp %>% subset(Symbol %in% topDysGenes$Symbol) %>% 
  reshape2::dcast(Symbol ~ DISEASE, value.var = "log2FoldChange") %>% column_to_rownames(var = "Symbol")


borderColors <- deseq2DiffExp %>% subset(Symbol %in% topDysGenes$Symbol) %>% 
  mutate(diffIndex = abs(log2FoldChange) > 1 & padj < 0.05) %>% reshape2::dcast(Symbol ~ DISEASE, value.var = "diffIndex") %>% 
  column_to_rownames(var = "Symbol")
borderColors <- t(apply(borderColors, 1, function(x) ifelse(x == TRUE, 'black', NA)))


rowAnno <- topDysGenes %>% select(rhoGroup, Symbol) %>% column_to_rownames(var = "Symbol")


plot2 <- pheatmap(logFcHeatmap[topDysGenes$Symbol, ], cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, 
         show_rownames = TRUE, number_color = "grey30", number_format = "%.2f", show_colnames = TRUE, 
         angle_col = 45, color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_number = 6, fontsize = 7, 
         cellheight = 12, cellwidth = 24, border = borderColors[topDysGenes$Symbol, ], annotation_row = rowAnno, na_col = "white")

#########################
topDysGenesBarplot <- topDysGenes %>% ungroup() %>% select(Symbol) %>% 
  inner_join(rhoDysStat, by = 'Symbol') %>% subset(dysLabel != 'Neutral') %>% 
  mutate(dysLabel = factor(dysLabel, levels = c('Up', 'Down')), Symbol = factor(Symbol, levels = rev(unique(Symbol))))


plot3 <- ggplot(topDysGenesBarplot, aes(x = Symbol, y = n, fill = dysLabel)) +
  geom_bar(stat = "identity", position = "stack") + scale_y_continuous(breaks = seq(0, 12, 3), position = "right") +
  labs(y = 'Number of cancers') + scale_fill_manual(values = c("Down" = "blue", "Up" = "red"), name = '') +
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black")) + coord_flip()

#########################
tcgaSsgseaRes <- readRDS(file = '/data/tcgaRhoFamiliesSsgsea2.0Res.rds')
load( file = '/data/tcgaExpSams.RData')
# tumorExpSams, normalExpSams, pairedSamsPats

tumorSamsScore <- tcgaSsgseaRes %>% subset(SampleIDs %in% tumorExpSams) %>% group_by(DISEASE, Groups) %>% mutate(meanScores = mean(Scores))

disOrder <- subset(tumorSamsScore, Groups == 'Small_GTPase_Rho') %>% arrange(meanScores) %>% ungroup()
disOrder <- unique(disOrder$DISEASE)


fillCol <- setNames(c('#734B27', '#364D99', '#168C42', '#B9A131', '#C7CAD9', '#C9A88D', '#EDB067', '#B0CC46', '#F5CED6', 
             '#DE6B71', '#98D6F0', '#92C9A4', '#D7EDF7', '#DC2425', '#F6DFC3', '#D92E88', '#4D2A80', '#B2202B', '#D07927', '#CFBFDC', '#0B477A', 
             '#791A1C', '#A64C93', '#967FB4', '#0FA095', '#1EA2DC', '#DEBC27', '#EEE33E', 
             '#1177A9', '#69779A', '#EEAAAE', '#CB96BF', '#ED8F27'), disOrder)


tumorSamsScore <- tumorSamsScore %>% mutate(DISEASE = factor(DISEASE, levels = disOrder))

plot4 <- ggplot(tumorSamsScore, aes(x = DISEASE, y = Scores, fill = DISEASE)) + 
  geom_violin() + scale_fill_manual(values = fillCol) +  
  # geom_jitter(width = 0.2, alpha = 0.6, size = 1) + 
  geom_point(aes(y = meanScores), shape = 95, size = 5, color = 'red') + 
  theme(axis.text.x = element_text(angle = 60, hjust = 0.8), legend.position = "none") + 
  xlab(NULL) + ylab('ssGsea Scores') + facet_wrap(~Groups, nrow = 3)

##########
disGroupMeanScore <- tcgaSsgseaRes %>% subset(SampleIDs %in% tumorExpSams) %>% group_by(DISEASE, Groups) %>% 
  summarise(meanScores = mean(Scores)) %>% reshape2::dcast(DISEASE~Groups, value.var = 'meanScores')


plot5 <- ggplot(disGroupMeanScore, aes(x = RhoGAP_dom, y = DH_DOCKER_dom, fill = DISEASE)) + 
  geom_point(shape = 21, size = 5, color = 'white') + scale_fill_manual(values = fillCol, guide = "legend", name = "Disease") + 
  geom_smooth(method = "lm", se = TRUE, color = '#ECC68C', fill = '#F5E4C8') + 
  geom_text(data = disGroupMeanScore, size = 3, check_overlap = TRUE, aes(label = DISEASE), nudge_y = 0.05) + 
  labs(x = 'GAP ssGsea Scores', y = 'GEF ssGsea Scores') + theme(legend.position = "none") 

# cor.test(~RhoGAP_dom + DH_DOCKER_dom, data = disGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.7022059
# pvalue = 9.69313e-06


plot6 <- ggplot(disGroupMeanScore, aes(x = DH_DOCKER_dom, y = Small_GTPase_Rho, fill = DISEASE)) + 
  geom_point(shape = 21, size = 5, color = 'white') + scale_fill_manual(values = fillCol, guide = "legend", name = "Disease") + 
  geom_smooth(method = "lm", se = TRUE, color = '#ECC68C', fill = '#F5E4C8') + 
  geom_text(data = disGroupMeanScore, size = 3, check_overlap = TRUE, aes(label = DISEASE), nudge_y = 0.05) + 
  labs(x = 'GEF ssGsea Scores', y = 'Rho ssGsea Scores') + theme(legend.position = "none") 

# cor.test(~DH_DOCKER_dom + Small_GTPase_Rho, data = disGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.3796791
# pvalue = 0.03003335

plot7 <- ggplot(disGroupMeanScore, aes(x = RhoGAP_dom, y = Small_GTPase_Rho, fill = DISEASE)) + 
  geom_point(shape = 21, size = 5, color = 'white') + scale_fill_manual(values = fillCol, guide = "legend", name = "Disease") + 
  geom_smooth(method = "lm", se = TRUE, color = '#ECC68C', fill = '#F5E4C8') + 
  geom_text(data = disGroupMeanScore, size = 3, check_overlap = TRUE, aes(label = DISEASE), nudge_y = 0.05) + 
  labs(x = 'GAP ssGsea Scores', y = 'Rho ssGsea Scores') + theme(legend.position = "none") 

# cor.test(~RhoGAP_dom + Small_GTPase_Rho, data = disGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.4141043
# pvalue = 0.01728471

########## unpaired
analyDis <- tcgaSsgseaRes %>% subset(SampleIDs %in% c(tumorExpSams, normalExpSams)) %>% 
  mutate(tumorTypes = ifelse(substr(SampleIDs, 14, 15) == '11', 'Normal', 'Tumor')) %>% group_by(DISEASE, tumorTypes) %>% 
  summarise(ncount = n_distinct(SampleIDs)) %>% subset(ncount >= 10) %>% count(DISEASE) %>% subset(n > 1)


unpairedSamSsgseaScore <- tcgaSsgseaRes %>% subset(SampleIDs %in% c(tumorExpSams, normalExpSams) & DISEASE %in% analyDis$DISEASE) %>% 
  mutate(tumorTypes = ifelse(substr(SampleIDs, 14, 15) == '11', 'Normal', 'Tumor'))


UnpairedTTest <- function(ssgseaScore){
  
  res <- t.test(Scores ~ tumorTypes, data = ssgseaScore)
  res <- data.frame(diffMeans = res$estimate[2]-res$estimate[1], pValues = res$p.value)
  
  return(res)
}

unpairedTTestRes <- unpairedSamSsgseaScore %>% group_by(Groups, DISEASE) %>% do(UnpairedTTest(.)) 
unpairedTTestRes$FDR <- p.adjust(unpairedTTestRes$pValues, method = 'BH')


########## paired
pairedSamSsgseaScore <- tcgaSsgseaRes %>% subset(PATIENT_BARCODE %in% pairedSamsPats & substr(SampleIDs, 14, 15) != '06') %>% 
  mutate(tumorTypes = ifelse(substr(SampleIDs, 14, 15) == '11', 'Normal', 'Tumor'))

PairedTTest <- function(ssgseaScore){
  
  ssgseaScore <- ssgseaScore %>% as.data.frame() %>% select(Scores, PATIENT_BARCODE, tumorTypes) %>% 
    reshape(direction = "wide", idvar = "PATIENT_BARCODE", timevar = "tumorTypes")

  res <- t.test(Pair(Scores.Tumor, Scores.Normal) ~ 1, data = ssgseaScore)
  res <- data.frame(diffMeans = res$estimate, pValues = res$p.value)
  
  return(res)
  
}

pairedTTestRes <- pairedSamSsgseaScore %>% group_by(Groups, DISEASE) %>% do(PairedTTest(.)) 
pairedTTestRes$FDR <- p.adjust(pairedTTestRes$pValues, method = 'BH')


DataTrans <- function(ssgseaScore){
  
  ssgseaScore <- ssgseaScore %>% as.data.frame() %>% select(Scores, PATIENT_BARCODE, tumorTypes) %>% 
    reshape(direction = "wide", idvar = "PATIENT_BARCODE", timevar = "tumorTypes")
  
  return(ssgseaScore)
  
}

pairedSamSsgseaScorePlot <- pairedSamSsgseaScore %>% group_by(Groups, DISEASE) %>% do(DataTrans(.))

colPalette <- c('#B2202B', '#69779A','#364D99', '#967FB4', '#1EA2DC', '#734B27', '#1177A9', '#EEAAAE', 
                '#C9A88D', '#DC2425', '#B9A131', '#C7CAD9', '#4D2A80', '#EEE33E')

plot8 <- pairedSamSsgseaScorePlot %>% ggpaired(cond1 = "Scores.Normal", cond2 = "Scores.Tumor", fill = "DISEASE", 
 color = "black", line.color = "gray80", line.size = 0.4, xlab = FALSE, palette = colPalette, 
 ylab = 'ssGsea Scores') + scale_x_discrete(labels = c('Normal', 'Tumor')) + geom_point(size = 1.2) +
  facet_wrap(facets = c('Groups', 'DISEASE'), nrow = 3) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 0.8, size = 6))

##########

pdf('/result/Section3/Figure3.1.pdf')
print(plot2)
dev.off()

pdf('/result/Section3/Figure3.2.pdf')
plot1
plot3
plot4
plot5
plot6
plot7
plot8
dev.off()


ggsave(plot = plot8, filename = '/result/Section3/rhoScoresdiff.pdf', width = 11)
