
library('ggplot2')
library('reshape2')

#########################
gainLoFMutStat <- readRDS(file = '/data/gainLoFMutStatByDisease.rds')

cols <- setNames(c('#734B27', '#364D99', '#168C42', '#B9A131', '#C7CAD9', '#C9A88D', '#EDB067', '#B0CC46', '#F5CED6', 
                      '#DE6B71', '#98D6F0', '#92C9A4', '#D7EDF7', '#DC2425', '#F6DFC3', '#D92E88', '#4D2A80', '#B2202B', '#D07927', '#CFBFDC', '#0B477A', 
                      '#791A1C', '#A64C93', '#967FB4', '#0FA095', '#1EA2DC', '#DEBC27', '#EEE33E', 
                      '#1177A9', '#69779A', '#EEAAAE', '#CB96BF', '#ED8F27'), sort(unique(gainLoFMutStat$DISEASE)))

##########PIK3R1
mutStat <- subset(gainLoFMutStat, Hugo_Symbol %in% c('PIK3R1'))

# size
mutStat$sizeLabel <- cut(mutStat$numOfNonsMut, breaks = c(0, 5, 10, 30, 100), right = F,
                                labels = c('n < 5', '5 ≤ n < 10', '10 ≤ n < 30', 'n ≥ 30'))

sizeLabel <- c(1.5, 3, 4.5, 6)
names(sizeLabel) <- c('n < 5', '5 ≤ n < 10', '10 ≤ n < 30', 'n ≥ 30')


plot1 <- ggplot(mutStat, aes(x = lofMf, y = gainMf, color = DISEASE, size = sizeLabel)) +
  geom_point() + labs(x = "Fraction of LoF mutations", y = "Fracton of Hotspot mutations") + 
  scale_size_manual(values = sizeLabel, guide = "legend", name = "Size") + scale_color_manual(values = cols) +
  geom_text(size = 3, check_overlap = TRUE, aes(label = DISEASE), nudge_y = 0.02) + guides(color = "none") + 
  scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.1, 0.3, 0.5, 0.7, 1.0)) + 
  scale_y_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) + 
  theme(legend.position = "bottom", legend.justification = "left") + ggtitle("PIK3R1")


##########RHOA
mutStat <- subset(gainLoFMutStat, Hugo_Symbol %in% c('RHOA'))

# size
mutStat$sizeLabel <- cut(mutStat$numOfNonsMut, breaks = c(0, 5, 10, 20, 100), right = F,
                         labels = c('n < 5', '5 ≤ n < 10', '10 ≤ n < 20', 'n ≥ 20'))

sizeLabel <- c(1.5, 3, 4.5, 6)
names(sizeLabel) <- c('n < 5', '5 ≤ n < 10', '10 ≤ n < 20', 'n ≥ 20')


plot2 <- ggplot(mutStat, aes(x = lofMf, y = gainMf, color = DISEASE, size = sizeLabel)) +
  geom_point() + labs(x = "Fraction of LoF mutations", y = "Fracton of Hotspot mutations") + 
  scale_size_manual(values = sizeLabel, guide = "legend", name = "Size") + scale_color_manual(values = cols) +
  geom_text(size = 3, check_overlap = TRUE, aes(label = DISEASE), nudge_y = 0.02) + guides(color = "none") + 
  scale_x_continuous(breaks = c(0, 0.1, 0.3, 0.5, 0.7, 1.0)) + 
  scale_y_continuous( breaks = c(0, 0.1, 0.3, 0.5, 0.7, 1.0)) + 
  theme(legend.position = "bottom", legend.justification = "left") + ggtitle("RHOA")


##########RHOB
mutStat <- subset(gainLoFMutStat, Hugo_Symbol %in% c('RHOB'))

# size
mutStat$sizeLabel <- cut(mutStat$numOfNonsMut, breaks = c(0, 5, 10, 20, 100), right = F,
                         labels = c('n < 5', '5 ≤ n < 10', '10 ≤ n < 20', 'n ≥ 20'))

sizeLabel <- c(1.5, 3, 4.5, 6)
names(sizeLabel) <- c('n < 5', '5 ≤ n < 10', '10 ≤ n < 20', 'n ≥ 20')


plot3 <- ggplot(mutStat, aes(x = lofMf, y = gainMf, color = DISEASE, size = sizeLabel)) +
  geom_point() + labs(x = "Fraction of LoF mutations", y = "Fracton of Hotspot mutations") + 
  scale_size_manual(values = sizeLabel, guide = "legend", name = "Size") + scale_color_manual(values = cols) +
  geom_text(size = 3, check_overlap = TRUE, aes(label = DISEASE), nudge_y = 0.02) + guides(color = "none") + 
  scale_x_continuous(breaks = c(0, 0.1, 0.3, 0.5, 0.7, 1.0)) + 
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.1, 0.3, 0.5, 0.7)) + 
  theme(legend.position = "bottom", legend.justification = "left") + ggtitle("RHOB")

##########RAC1
mutStat <- subset(gainLoFMutStat, Hugo_Symbol %in% c('RAC1'))

# size
mutStat$sizeLabel <- cut(mutStat$numOfNonsMut, breaks = c(0, 5, 10, 20, 100), right = F,
                         labels = c('n < 5', '5 ≤ n < 10', '10 ≤ n < 20', 'n ≥ 20'))

sizeLabel <- c(1.5, 3, 4.5, 6)
names(sizeLabel) <- c('n < 5', '5 ≤ n < 10', '10 ≤ n < 20', 'n ≥ 20')
 

plot4 <- ggplot(mutStat, aes(x = lofMf, y = gainMf, color = DISEASE, size = sizeLabel)) +
  geom_point() + labs(x = "Fraction of LoF mutations", y = "Fracton of Hotspot mutations") + 
  scale_size_manual(values = sizeLabel, guide = "legend", name = "Size") + scale_color_manual(values = cols) +
  geom_text(size = 3, check_overlap = TRUE, aes(label = DISEASE), nudge_y = 0.02) + guides(color = "none") + 
  scale_x_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.3, 0.4)) + 
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.1, 0.3, 0.5, 0.7)) + 
  theme(legend.position = "bottom", legend.justification = "left") + ggtitle("RAC1")



##########ARHGAP35
mutStat <- subset(gainLoFMutStat, Hugo_Symbol %in% c('ARHGAP35'))

# size
mutStat$sizeLabel <- cut(mutStat$numOfNonsMut, breaks = c(0, 5, 10, 30, 100), right = F,
                         labels = c('n < 5', '5 ≤ n < 10', '10 ≤ n < 30', 'n ≥ 30'))

sizeLabel <- c(1.5, 3, 4.5, 6)
names(sizeLabel) <- c('n < 5', '5 ≤ n < 10', '10 ≤ n < 30', 'n ≥ 30')


plot5 <- ggplot(mutStat, aes(x = lofMf, y = gainMf, size = sizeLabel, color = DISEASE)) +
  geom_point() + labs(x = "Fraction of LoF mutations", y = "Fracton of Hotspot mutations") + 
  scale_size_manual(values = sizeLabel, guide = "legend", name = "Size") + scale_color_manual(values = cols) +
  geom_text(size = 3, check_overlap = TRUE, aes(label = DISEASE), nudge_y = 0.003) + guides(color = "none") + 
  scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.1, 0.3, 0.5, 0.7, 1.0)) + 
  scale_y_continuous(limits = c(0, 0.04), breaks = c(0, 0.01, 0.02, 0.03, 0.04)) + 
  theme(legend.position = "bottom", legend.justification = "left") + ggtitle("ARHGAP35")


#############
pdf('/result/Section1/SFigure_HotLofGenes.pdf')
plot1
plot2
plot3
plot4
plot5
dev.off()


