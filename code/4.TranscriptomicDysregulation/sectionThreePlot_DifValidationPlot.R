library(pheatmap)
library(tibble)
library(ggplot2)
library(ggpubr)
library(scatterpie)
library(dplyr)
#########################
load(file = '/data/difgene_GPL570.RData') 
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

deseq2DiffExp <- lapply(difexp1, rownames_to_column, var = 'Symbol') %>% bind_rows(.id = "DISEASE") %>% 
  mutate(DISEASE = str_to_title(DISEASE), dysLabel = ifelse(adj.P.Val < 0.05 & logFC > 0.585, 'Up', 
                                                            ifelse(adj.P.Val < 0.05 & logFC < -0.585, 'Down', 'Neutral'))) %>% 
  inner_join(rhoFamilies, by = c('Symbol' = 'Approved.symbol'))


deseq2DiffExp <- deseq2DiffExp %>% mutate(rhoGroup = ifelse(Class %in% c('DH_domain', 'DOCKER_dom'), 'DH_DOCKER_dom', Class))

rhoGroupDysStat <- deseq2DiffExp %>% group_by(DISEASE, rhoGroup) %>% count(dysLabel)
rhoDysStat <- deseq2DiffExp %>% group_by(rhoGroup, Symbol) %>% count(dysLabel)

topDysGenes <- rhoDysStat %>% subset(dysLabel != 'Neutral') %>% group_by(rhoGroup, Symbol) %>% summarise(n = sum(n)) %>% 
  group_by(rhoGroup) %>% arrange(desc(n)) %>% slice_head(n = 6)

#########################
groupDysPiePlot <- rhoGroupDysStat %>% left_join(data.frame(DISEASE = unique(rhoGroupDysStat$DISEASE), x_locus = 1:7), by = 'DISEASE') %>% 
  left_join(data.frame(rhoGroup = unique(rhoGroupDysStat$rhoGroup), y_locus = 4:1), by = 'rhoGroup') %>% 
  mutate(dysLabel = factor(dysLabel, levels = c('Down', 'Neutral', 'Up')), value = n) %>% ungroup()


plot1 <- ggplot() + geom_scatterpie(aes(x=x_locus, y=y_locus), 
                                    data=groupDysPiePlot, cols="dysLabel", long_format=TRUE) + coord_fixed() + 
  scale_fill_manual(values = c('Down' = 'blue', 'Neutral' = 'white', 'Up' = "red"), name = '') +
  scale_x_continuous(breaks = 1:7, labels = unique(groupDysPiePlot$DISEASE), ) + 
  scale_y_continuous(breaks = 1:4, labels = c('Rho', 'GDI', 'GAP' ,'GEF')) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 0.8), axis.text.y = element_text(angle = 60, vjust = -0.05))

#########################
rhoDysStatByDis <- deseq2DiffExp %>% group_by(DISEASE) %>% distinct(Symbol, .keep_all = T) %>% group_by(DISEASE) %>% count(dysLabel)

plot2 <- ggplot(rhoDysStatByDis, aes(x = DISEASE, y = n, fill = dysLabel)) +
  geom_bar(stat = "identity", position = "stack") + scale_y_continuous(breaks = seq(0, 175, 25), position = "left") +
  labs(y = 'Number of Dyregulated Genes') + scale_fill_manual(values = c("Down" = "blue", "Neutral" = 'white', "Up" = "red"), name = '') +
  theme(axis.title.x = element_blank(), legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black")) + 
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) 


#########################
logFcHeatmap <- deseq2DiffExp %>% subset(Symbol %in% topDysGenes$Symbol) %>% 
  reshape2::dcast(Symbol ~ DISEASE, value.var = "logFC") %>% column_to_rownames(var = "Symbol")


borderColors <- deseq2DiffExp %>% subset(Symbol %in% topDysGenes$Symbol) %>% 
  mutate(diffIndex = abs(logFC) > 0.585 & adj.P.Val < 0.05) %>% reshape2::dcast(Symbol ~ DISEASE, value.var = "diffIndex") %>% 
  column_to_rownames(var = "Symbol")
borderColors <- t(apply(borderColors, 1, function(x) ifelse(x == TRUE, 'black', NA)))


rowAnno <- topDysGenes %>% select(rhoGroup, Symbol) %>% column_to_rownames(var = "Symbol")


plot3 <- pheatmap(logFcHeatmap[topDysGenes$Symbol, ], cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, 
                  show_rownames = TRUE, number_color = "grey30", number_format = "%.2f", show_colnames = TRUE, 
                  angle_col = 45, color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_number = 6, fontsize = 7, 
                  cellheight = 12, cellwidth = 24, border = borderColors[topDysGenes$Symbol, ], annotation_row = rowAnno, na_col = "white")

#########################
topDysGenesBarplot <- topDysGenes %>% ungroup() %>% select(Symbol) %>% 
  inner_join(rhoDysStat, by = 'Symbol') %>% subset(dysLabel != 'Neutral') %>% 
  mutate(dysLabel = factor(dysLabel, levels = c('Up', 'Down')), Symbol = factor(Symbol, levels = rev(unique(Symbol))))


plot4 <- ggplot(topDysGenesBarplot, aes(x = Symbol, y = n, fill = dysLabel)) +
  geom_bar(stat = "identity", position = "stack") + scale_y_continuous(breaks = seq(0, 12, 3), position = "right") +
  labs(y = 'Number of cancers') + scale_fill_manual(values = c("Down" = "blue", "Up" = "red"), name = '') +
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black")) + coord_flip()

#########################
pdf('/result/Section3/validation/Figure3.1.pdf')
print(plot3)
dev.off()

pdf('/result/Section3/validation/Figure3.2.pdf')
plot1
plot2
plot4
dev.off()
