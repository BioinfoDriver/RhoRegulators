
library(tibble)
library(pheatmap)

##########################
load(file = '/data/rhoFamiliesMutStat.RData')
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')
load('/data/rhoMutComToRandom.RData')
nonsilentTBM <- readRDS(file = '/data/meanOfNonsilentMutBurden.rds')

##########################
heatmapData <- reshape2::dcast(rhoFamiliesMutStatByDisGroup, DISEASE ~ Class, value.var = "altFre")
heatmapData <- column_to_rownames(heatmapData, var = "DISEASE")

heatmapData[is.na(heatmapData)] <- 0
heatmapData <- as.data.frame(t(heatmapData))

borderColors <- apply(comResByDisGroupPvalue, 1, function(x) ifelse(x < 0.05, 'black', NA))

nonsilentTBM <- nonsilentTBM %>% column_to_rownames(var = "DISEASE") %>% mutate(TBM = round(n, 1), .keep = 'unused')

plot1 <- pheatmap(heatmapData, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, show_rownames = TRUE, 
                  number_color = "grey30", number_format = "%.2f", show_colnames = TRUE, angle_col = 45, 
                  fontsize_number = 5, fontsize = 7, 
                  color = colorRampPalette(c("white", "firebrick3"))(50), cellheight = 12, cellwidth = 14, border = borderColors, 
                  annotation_col = nonsilentTBM, annotation_colors = list(TBM = colorRampPalette(c("white", "#9EE092"))(100)))

##########################
nRMutStat <- rhoFamiliesMutStatByGene %>% mutate(DISEASE = 'ZPanCan', numOfPatWithMut = numOfMut, numOfPat = 8217)
nRMutStat <- rbind.data.frame(nRMutStat[, colnames(rhoFamiliesMutStatByDisGene)], rhoFamiliesMutStatByDisGene)

nRMutheatmap <- reshape2::dcast(nRMutStat, Hugo_Symbol~DISEASE, value.var = "altFre")
nRMutheatmap <- column_to_rownames(nRMutheatmap, var = "Hugo_Symbol")

nRMutheatmap[is.na(nRMutheatmap)] <- 0
# nRMutheatmap <- nRMutheatmap[apply(nRMutheatmap, 1, function(x) sum(x > 0.01)>15), ]
nRMutheatmap <- nRMutheatmap %>% arrange(desc(ZPanCan)) %>% slice_head(n = 25)


borderColors <- t(apply(nRMutheatmap, 1, function(x) ifelse(x > 0.05, 'black', NA)))

rhoFamilies <- rhoFamilies[match(rownames(nRMutheatmap), rhoFamilies$Approved.symbol), ] %>% arrange(Class)
nRMutheatmap <- nRMutheatmap[rhoFamilies$Approved.symbol, ]
borderColors <- borderColors[rhoFamilies$Approved.symbol, ]

plot2 <- pheatmap(nRMutheatmap, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, show_rownames = TRUE, 
                  number_color = "grey30", number_format = "%.2f", show_colnames = TRUE, angle_col = 45, 
                  fontsize_number = 5, fontsize = 7, 
                  color = colorRampPalette(c("white", "firebrick3"))(50), cellheight = 12, cellwidth = 12, border = borderColors)

##########################
nRMutPvalue <- as.data.frame(comResByDisPvalue) %>% rownames_to_column(var = 'DISEASE') %>% 
  mutate(freq = 1 - comResByDisPvalue, .keep = 'unused') %>% arrange(desc(freq))

fillCol <- setNames(c('#734B27', '#364D99', '#168C42', '#B9A131', '#C7CAD9', '#C9A88D', '#EDB067', '#B0CC46', '#F5CED6', 
                   '#DE6B71', '#98D6F0', '#92C9A4', '#D7EDF7', '#DC2425', '#F6DFC3', '#D92E88', '#4D2A80', '#B2202B', '#D07927', '#CFBFDC', '#0B477A', 
                   '#791A1C', '#A64C93', '#967FB4', '#0FA095', '#1EA2DC', '#DEBC27', '#EEE33E', 
                   '#1177A9', '#69779A', '#EEAAAE', '#CB96BF', '#ED8F27'), sort(nRMutPvalue$DISEASE))

plot3 <- ggplot(nRMutPvalue, aes(x = reorder(DISEASE, -freq), y = freq, fill = DISEASE)) +
  geom_bar(stat = "identity")  + labs(y = 'Normalized mutation load') + 
  scale_fill_manual(values = fillCol) + guides(fill = 'none')+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
  geom_hline(yintercept = c(0.05, 0.95), color = c('black', "red"), linetype = "dashed", linewidth = 1) 


pdf('/result/Section1/SFigure_MutBurden_1.pdf')
plot1
plot3
dev.off()


pdf('/result/Section1/SFigure_MutBurden_2.pdf')
plot2
dev.off()
