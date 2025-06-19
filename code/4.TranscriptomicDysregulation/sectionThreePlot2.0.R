
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
library(tibble)
library(ggbreak)

#########################
panCanSurPvalue <- readRDS(file = '/data/panCanRhoSurPvalue.rds')
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

rhoFamilies <- rhoFamilies %>% mutate(Class = ifelse(Class %in% c('DH_domain', 'DOCKER_dom'), 'DH_DOCKER_dom', Class))
panCanSurPvalue <-merge.data.frame(panCanSurPvalue, rhoFamilies[, c('Approved.symbol', 'Class')], by.x = 'geneSymbol', by.y = 'Approved.symbol')
  

###########

selectSurPvalue <- subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05)
selectSurPvalue <- selectSurPvalue %>% group_by(Class, geneSymbol) %>% count() %>% group_by(Class) %>% arrange(desc(n)) %>% slice_head(n = 5) %>%
  ungroup() %>% left_join(panCanSurPvalue, by='geneSymbol') %>% subset(coxPvalue < 0.05) %>% 
  mutate(geneSymbol = factor(geneSymbol, levels = rev(unique(geneSymbol))))



selectSurPvalue <- selectSurPvalue %>% mutate(hrDirect = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival'), 
                                              logrankPcateg = ifelse(logrankPvalue < 0.05, 'logrank p < 0.05', 'logrank p ≥ 0.05')) 
selectSurPvalue <- selectSurPvalue %>% mutate(coxPcateg = cut(coxPvalue, breaks = c(0, 0.00001, 0.001, 0.01, 0.05), 
                                                              labels = c('cox p < 1e-05', '1e-05 < cox p ≤ 0.001', '0.001 < cox p ≤ 0.01', '0.01 < cox p < 0.05')))


fillLabel <- c("pink", "lightblue")
names(fillLabel) <- c('High level-> worse survival', 'High level-> better survival')
sizeLable <- c(3, 4, 5, 6)
names(sizeLable) <- c('0.01 < cox p < 0.05', '0.001 < cox p ≤ 0.01', '1e-05 < cox p ≤ 0.001', 'cox p < 1e-05')
colLabel <- c("black", "white")
names(colLabel) <- c('logrank p < 0.05', 'logrank p ≥ 0.05')



plot1 <- ggplot(selectSurPvalue, aes(x = disease, y = geneSymbol, color = logrankPcateg, size = coxPcateg, fill = hrDirect)) +
  geom_point(shape = 21, stroke = 1.5) + labs(x = NULL, y = NULL) + 
  scale_colour_manual(values = colLabel, guide = "legend", name = "Color") + 
  scale_size_manual(values = sizeLable, guide = "legend", name = "Size") + 
  scale_fill_manual(values = fillLabel, guide = "legend", name = "Fill") + 
  theme(legend.position = "top", legend.justification = "left", axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


###########
selectSurPvalue <- subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05)
selectSurPvalue <- selectSurPvalue %>% mutate(hrDirect = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival'))

disDysGCount <- selectSurPvalue %>% distinct(disease, geneSymbol, .keep_all = TRUE) %>% group_by(disease, hrDirect) %>% count(name = 'num')

showGenes <- selectSurPvalue %>% group_by(Class, geneSymbol) %>% count() %>% group_by(Class) %>% arrange(desc(n)) %>% slice_head(n = 5)
showGenes <- showGenes$geneSymbol
dysGcount <- selectSurPvalue %>% group_by(geneSymbol, hrDirect) %>% count(name = 'num') %>% subset(geneSymbol %in% showGenes) %>% 
  mutate(geneSymbol = factor(geneSymbol, levels = rev(showGenes)))


plot2 <- ggplot(disDysGCount, aes(x = disease, y = num, fill = hrDirect)) +
  geom_bar(stat = "identity")  + labs(y = 'Number of dysregulated genes') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "top", legend.justification = "left") + 
  scale_y_break(breaks=c(50, 70), ticklabels = c(70, 80, 90), scales = 0.2, space = 0.2)


plot3 <- ggplot(dysGcount, aes(x = geneSymbol, y = num, fill = hrDirect)) +
  geom_bar(stat = "identity")  + labs(y = 'Number of dysregulated genes') + 
  scale_y_continuous(position = "right") + 
  theme(legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black"), 
        axis.title.y = element_blank()) + coord_flip()


#########################
panCanEnrichScoreSurPvalue <- readRDS(file = '/data/panCanRhoEnrichScoreSurPvalue.rds')

enrichScoreSurPlot <- panCanEnrichScoreSurPvalue %>% subset(coxPvalue < 0.5)


enrichScoreSurPlot <- enrichScoreSurPlot %>% mutate(hrDirect = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival'), 
                                              logrankPcateg = ifelse(logrankPvalue < 0.05, 'logrank p < 0.05', 'logrank p ≥ 0.05')) 
enrichScoreSurPlot <- enrichScoreSurPlot %>% mutate(coxPcateg = cut(coxPvalue, breaks = c(0, 0.001, 0.01, 0.05, 0.5), right = FALSE, 
                                                              labels = c('cox p < 0.001', '0.001 ≤ cox p < 0.01', '0.01 ≤ cox p < 0.05', 'cox p ≥ 0.05')))


fillLabel <- c("pink", "lightblue")
names(fillLabel) <- c('High level-> worse survival', 'High level-> better survival')
sizeLable <- c(2, 4, 5, 6)
names(sizeLable) <- c('cox p ≥ 0.05', '0.01 ≤ cox p < 0.05', '0.001 ≤ cox p < 0.01', 'cox p < 0.001')
colLabel <- c("black", "white")
names(colLabel) <- c('logrank p < 0.05', 'logrank p ≥ 0.05')



plot4 <- ggplot(enrichScoreSurPlot, aes(x = disease, y = geneSymbol, color = logrankPcateg, size = coxPcateg, fill = hrDirect)) +
  geom_point(shape = 21, stroke = 1.5) + labs(x = NULL, y = NULL) + 
  scale_colour_manual(values = colLabel, guide = "legend", name = "Color") + 
  scale_size_manual(values = sizeLable, guide = "legend", name = "Size") + 
  scale_fill_manual(values = fillLabel, guide = "legend", name = "Fill") + 
  theme(legend.position = "top", legend.justification = "left", axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


###########

panCanEnrichGroupSurPvalue <- readRDS(file = '/data/panCanRhoEnrichGroupSurPvalue.rds')

enrichGroupSurPlot <- panCanEnrichGroupSurPvalue %>% subset(coxPvalue < 1.0)


enrichGroupSurPlot <- enrichGroupSurPlot %>% mutate(hrDirect = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival'), 
                                                    logrankPcateg = ifelse(logrankPvalue < 0.05, 'logrank p < 0.05', 'logrank p ≥ 0.05')) 
enrichGroupSurPlot <- enrichGroupSurPlot %>% mutate(coxPcateg = cut(coxPvalue, breaks = c(0, 0.001, 0.01, 0.05, 1.0), right = FALSE, 
                                                                    labels = c('cox p < 0.001', '0.001 ≤ cox p < 0.01', '0.01 ≤ cox p < 0.05', 'cox p ≥ 0.05')))


fillLabel <- c("pink", "lightblue")
names(fillLabel) <- c('High level-> worse survival', 'High level-> better survival')
sizeLable <- c(2, 4, 5, 6)
names(sizeLable) <- c('cox p ≥ 0.05', '0.01 ≤ cox p < 0.05', '0.001 ≤ cox p < 0.01', 'cox p < 0.001')
colLabel <- c("black", "white")
names(colLabel) <- c('logrank p < 0.05', 'logrank p ≥ 0.05')



plot5 <- ggplot(enrichGroupSurPlot, aes(x = disease, y = geneSymbol, color = logrankPcateg, size = coxPcateg, fill = hrDirect)) +
  geom_point(shape = 21, stroke = 1.5) + labs(x = NULL, y = NULL) + 
  scale_colour_manual(values = colLabel, guide = "legend", name = "Color") + 
  scale_size_manual(values = sizeLable, guide = "legend", name = "Size") + 
  scale_fill_manual(values = fillLabel, guide = "legend", name = "Fill") + 
  theme(legend.position = "top", legend.justification = "left", axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


pdf('/result/Section4/Figure4.1.pdf')
plot1
plot2
plot3
plot4
plot5
dev.off()

