library('dplyr')
library(ggplot2)

####################
load(file = '/data/tcgaExpSams.RData')
# tumorExpSams, normalExpSams, pairedSamsPats
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

####################
# paired samples
disSamPairs <- subset.data.frame(tcgaPanCanSamples, PATIENT_BARCODE %in% pairedSamsPats) 
disSamPairsStat <- disSamPairs %>% group_by(DISEASE) %>% count()

fillCol <- setNames(c('#734B27', '#364D99', '#168C42', '#B9A131', '#C7CAD9', '#C9A88D', '#EDB067', '#B0CC46', '#F5CED6', 
                      '#DE6B71', '#98D6F0', '#92C9A4', '#D7EDF7', '#DC2425', '#F6DFC3', '#D92E88', '#4D2A80', '#B2202B', 
                      '#D07927', '#CFBFDC', '#0B477A', '#791A1C', '#A64C93', '#967FB4', '#0FA095', '#1EA2DC', '#DEBC27', 
                      '#EEE33E', '#1177A9', '#69779A', '#EEAAAE', '#CB96BF', '#ED8F27'), sort(unique(tcgaPanCanSamples$DISEASE)))

plot1 <- ggplot(disSamPairsStat, aes(x = reorder(DISEASE, -n), y = n, fill = DISEASE)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = n), vjust = -0.5) + labs(y = 'Number of Samples') + 
  theme(axis.title.x = element_blank(), panel.background = element_rect(fill = "white", color = "black")) + 
  scale_fill_manual(values = fillCol, guide = 'none')

####################
tumorSamsStat <- subset.data.frame(tcgaPanCanSamples, SAMPLE_BARCODE %in% tumorExpSams) %>% group_by(DISEASE) %>% count()

plot2 <- ggplot(tumorSamsStat, aes(x = reorder(DISEASE, -n), y = n, fill = DISEASE)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = n), vjust = -0.5, size = 3) + labs(y = 'Number of Samples') + 
  theme(axis.title.x = element_blank(), panel.background = element_rect(fill = "white", color = "black")) + 
  scale_fill_manual(values = fillCol, guide = 'none')

####
pdf(file = '/result/Section3/tcgaExpSamStat.pdf')
plot1
plot2
dev.off()
