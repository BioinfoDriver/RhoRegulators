library(ggpubr)
library(ggplot2)
library(dbplyr)
library(reshape2)
library(ggcorrplot)
#####################

valiDataSsgseaRes <- readRDS(file = '/data/rhoFam.GPL570CaseCtrl.rds')

##################### paired
pairedSams <- valiDataSsgseaRes %>% subset(!is.na(patientID)) %>% group_by(cancer_type, GEO_number, patientID) %>% 
  count(patientID) %>% subset.data.frame(n == 6)

pairedSamsSsgseaRes <- valiDataSsgseaRes %>% inner_join(pairedSams, join_by(cancer_type, GEO_number, patientID))

pairedSams <- pairedSamsSsgseaRes %>% group_by(cancer_type, GEO_number, patientID) %>% 
  summarise(ncount = n_distinct(tissue_type1)) %>% subset(ncount > 1)

pairedSamsSsgseaRes <- pairedSamsSsgseaRes %>% inner_join(pairedSams, join_by(cancer_type, GEO_number, patientID))

# pairedSamsSsgseaRes %>% group_by(cancer_type, GEO_number) %>% 
#   summarise(ncount = n_distinct(patientID)) %>% group_by(cancer_type) %>% summarise(sum(ncount))
# 1 breast                 26
# 2 colon                  96
# 3 kidney                180
# 4 liver                 138
# 5 lung                  154
# 6 stomach               129
# 7 thyroid                44
# totally 757


PairedTTest <- function(ssgseaScore){
  
  ssgseaScore <- ssgseaScore %>% as.data.frame() %>% select(Scores, patientID, tissue_type1) %>% 
    reshape(direction = "wide", idvar = "patientID", timevar = "tissue_type1")
  
  res <- t.test(Pair(Scores.tumor, Scores.normal) ~ 1, data = ssgseaScore)
  
  res <- data.frame(diffMeans = res$estimate, pValues = res$p.value)
  
  return(res)
  
}

pairedTTestRes <- pairedSamsSsgseaRes %>% group_by(Groups, cancer_type) %>% do(PairedTTest(.)) 
pairedTTestRes$FDR <- p.adjust(pairedTTestRes$pValues, method = 'BH')


##################### paired plot
DataTrans <- function(ssgseaScore){
  
  ssgseaScore <- ssgseaScore %>% as.data.frame() %>% select(Scores, patientID, tissue_type1) %>% 
    reshape(direction = "wide", idvar = "patientID", timevar = "tissue_type1")
  
  return(ssgseaScore)
  
}

pairedSamSsgseaScorePlot <- pairedSamsSsgseaRes %>% group_by(Groups, cancer_type) %>% do(DataTrans(.))

colPalette <- c('#69779A', '#364D99', '#1177A9', '#DC2425', '#C9A88D', '#4D2A80', '#EEE33E')


plot1 <- pairedSamSsgseaScorePlot %>% ggpaired(cond1 = "Scores.normal", cond2 = "Scores.tumor", fill = "cancer_type", 
                                               color = "black", line.color = "gray95", line.size = 0.4, xlab = FALSE, 
                                               palette = colPalette, ylab = 'ssGsea Scores') + 
  scale_x_discrete(labels = c('Normal', 'Tumor'))+ geom_point(size = 1.2) +  
  facet_wrap(facets = c('Groups', 'cancer_type'), nrow = 3) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 0.8, size = 6))



########## unpaired
UnpairedTTest <- function(ssgseaScore){
  
  res <- t.test(Scores ~ tissue_type1, data = ssgseaScore)
  res <- data.frame(diffMeans = res$estimate[2]-res$estimate[1], pValues = res$p.value)
  
  return(res)
}

UnpairedTTestRes <- valiDataSsgseaRes %>% group_by(Groups, cancer_type) %>% do(UnpairedTTest(.)) 
UnpairedTTestRes$FDR <- p.adjust(UnpairedTTestRes$pValues, method = 'BH')


#####################
disGroupMeanScore <- valiDataSsgseaRes %>% subset(tissue_type1 == 'tumor') %>% group_by(cancer_type, Groups) %>% 
  summarise(meanScores = mean(Scores)) %>% reshape2::dcast(cancer_type~Groups, value.var = 'meanScores')


# cor.test(~RhoGAP_dom + DH_DOCKER_dom, data = disGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.9642857
# pvalue = 0.002777778

# cor.test(~DH_DOCKER_dom + Small_GTPase_Rho, data = disGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.8928571
# pvalue = 0.01230159

# cor.test(~RhoGAP_dom + Small_GTPase_Rho, data = disGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.9285714
# pvalue = 0.006746032


norGroupMeanScore <- valiDataSsgseaRes %>% subset(tissue_type1 == 'normal') %>% group_by(cancer_type, Groups) %>% 
  summarise(meanScores = mean(Scores)) %>% reshape2::dcast(cancer_type~Groups, value.var = 'meanScores')

# cor.test(~RhoGAP_dom + DH_DOCKER_dom, data = norGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.9285714
# pvalue = 0.006746032

# cor.test(~DH_DOCKER_dom + Small_GTPase_Rho, data = norGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.8571429
# pvalue = 0.02380952

# cor.test(~RhoGAP_dom + Small_GTPase_Rho, data = norGroupMeanScore, method = 'spearman', alternative = 'two.sided')
# rho = 0.7857143
# pvalue = 0.04801587


##################### ssgsea score correlation by disease
tumorSamsScore <- valiDataSsgseaRes %>% subset(tissue_type1 == 'tumor') %>% split.data.frame(f = ~ cancer_type)

turRhoSsGseaScoreCor <- sapply(tumorSamsScore, function(gseaScore){
  
  gseaScore <- gseaScore %>% dcast(SampleIDs~Groups, value.var = 'Scores') %>% column_to_rownames(var = 'SampleIDs')
  corRes <- cor(x = gseaScore, method = "spearman")
  
  corRes <- setNames(object = c(corRes['DH_DOCKER_dom', 'RhoGAP_dom'], corRes['DH_DOCKER_dom', 'Small_GTPase_Rho'], 
                                corRes['RhoGAP_dom', 'Small_GTPase_Rho']), nm = c('GEF_GAP', 'GEF_Rho', 'GAP_Rho'))
  return(corRes)
})

########
turRhoSsGseaScoreCorPvalue <- sapply(tumorSamsScore, function(gseaScore){
  
  gseaScore <- gseaScore %>% dcast(SampleIDs~Groups, value.var = 'Scores') %>% column_to_rownames(var = 'SampleIDs')
  corPvalue <- corrplot::cor.mtest(mat = gseaScore, method = 'spearman', conf.level = 0.95)$p
  
  corPvalue <- setNames(object = c(corPvalue['DH_DOCKER_dom', 'RhoGAP_dom'], corPvalue['DH_DOCKER_dom', 'Small_GTPase_Rho'], 
                                   corPvalue['RhoGAP_dom', 'Small_GTPase_Rho']), nm = c('GEF_GAP', 'GEF_Rho', 'GAP_Rho'))  
  return(corPvalue)
})

########
plot2 <- ggcorrplot(corr = t(turRhoSsGseaScoreCor), method = "circle", type = "full", hc.order = FALSE, tl.cex = 8, 
                    lab = TRUE, lab_size = 3, p.mat = t(turRhoSsGseaScoreCorPvalue), sig.level = 0.05, colors = c("#6D9EC1", "white", "#E46726")) + 
  theme(legend.position = "bottom", legend.direction = "horizontal")



#####################
normalSamsScore <- valiDataSsgseaRes %>% subset(tissue_type1 == 'normal') %>% split.data.frame(f = ~ cancer_type)


norRhoSsGseaScoreCor <- sapply(normalSamsScore, function(gseaScore){
  
  gseaScore <- gseaScore %>% dcast(SampleIDs~Groups, value.var = 'Scores') %>% column_to_rownames(var = 'SampleIDs')
  corRes <- cor(x = gseaScore, method = "spearman")
  
  corRes <- setNames(object = c(corRes['DH_DOCKER_dom', 'RhoGAP_dom'], corRes['DH_DOCKER_dom', 'Small_GTPase_Rho'], 
                                corRes['RhoGAP_dom', 'Small_GTPase_Rho']), nm = c('GEF_GAP', 'GEF_Rho', 'GAP_Rho'))
  return(corRes)
})

########
norRhoSsGseaScoreCorPvalue <- sapply(normalSamsScore, function(gseaScore){
  
  gseaScore <- gseaScore %>% dcast(SampleIDs~Groups, value.var = 'Scores') %>% column_to_rownames(var = 'SampleIDs')
  corPvalue <- corrplot::cor.mtest(mat = gseaScore, method = 'spearman', conf.level = 0.95)$p
  
  corPvalue <- setNames(object = c(corPvalue['DH_DOCKER_dom', 'RhoGAP_dom'], corPvalue['DH_DOCKER_dom', 'Small_GTPase_Rho'], 
                                   corPvalue['RhoGAP_dom', 'Small_GTPase_Rho']), nm = c('GEF_GAP', 'GEF_Rho', 'GAP_Rho'))  
  return(corPvalue)
})

########
plot3 <- ggcorrplot(corr = t(norRhoSsGseaScoreCor), method = "circle", type = "full", hc.order = FALSE, tl.cex = 8, 
                    lab = TRUE, lab_size = 3, p.mat = t(norRhoSsGseaScoreCorPvalue), sig.level = 0.05, colors = c("#6D9EC1", "white", "#E46726")) + 
  theme(legend.position = "bottom", legend.direction = "horizontal")


#####################
pdf('/result/Section3/validation/rhoSsGseaScoreCor.pdf')
plot1
plot2
plot3

dev.off()


#####################
turRhoSsGseaScoreCor <- turRhoSsGseaScoreCor %>% melt() %>% dplyr::rename(GroupCor = Var1, Disease = Var2, TurCor = value)
norRhoSsGseaScoreCor <- norRhoSsGseaScoreCor %>% melt() %>% dplyr::rename(GroupCor = Var1, Disease = Var2, NorCor = value)
rhoSsGseaScoreCor <- merge(turRhoSsGseaScoreCor, norRhoSsGseaScoreCor, by = c('GroupCor', 'Disease'))


WilcoxTest <- function(corData){
  pairedRes <- wilcox.test(Pair(TurCor, NorCor) ~ 1, data = corData)
  unpairedRes <- wilcox.test(x = corData$TurCor, y = corData$NorCor)
  
  pValue <- data.frame(pairedPvalue = pairedRes$p.value, unpairedPvalue = unpairedRes$p.value)
  return(pValue)
}

# rhoSsGseaScoreCor %>% group_by(GroupCor) %>% do(WilcoxTest(.))

# GroupCor pairedPvalue unpairedPvalue
# <fct>           <dbl>          <dbl>
#   1 GEF_GAP         0.688          0.805
# 2 GEF_Rho         0.219          0.165
# 3 GAP_Rho         0.109          0.128


colPalette <- c('#92C9A4', "#6D9EC1", "#E46726")

plot6 <- rhoSsGseaScoreCor %>% ggpaired(cond1 = "NorCor", cond2 = "TurCor", line.color = 'gray95', fill = 'GroupCor',
                                        line.size = 0.4, xlab = FALSE, palette = colPalette, 
                                        ylab = 'Correlation coefficient') + geom_point(size = 1) +
  scale_x_discrete(labels = c('Normal', 'Tumor')) + facet_wrap(facets = c('GroupCor'), ncol = 3) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 0.8, size = 8))



ggsave(filename = '/result/Section3/validation/rhoSsGseaScoreCorCompare.pdf', plot = plot6)


