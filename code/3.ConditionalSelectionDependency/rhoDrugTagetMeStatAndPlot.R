library(dplyr)
######################
selectRes <- readRDS(file = '/data/rhoSelectScoreRes.RData')
durgTargetGenes <- read.csv(file = '/data/oncokb_biomarker_drug_associations.tsv', header = T, sep = '\t', stringsAsFactors = F)
geneInfo <- read.csv(file= '/data/gene_with_protein_product.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')


rhoFamilies <- rhoFamilies %>% distinct(Approved.symbol, .keep_all = TRUE) %>% 
  mutate(Class = ifelse(Class %in% c('DH_domain', 'DOCKER_dom'), 'DH_DOCKER_domain', Class))


durgTargetGenes <- durgTargetGenes %>% distinct(Gene, .keep_all = TRUE)

sigselectRes <- selectRes %>% subset(FDR & APC_good) %>% mutate(SFE_1 = gsub("_.*$", "", SFE_1), SFE_2 = gsub("_.*$", "", SFE_2)) %>% 
  subset(SFE_1 %in% durgTargetGenes$Gene | SFE_2 %in% durgTargetGenes$Gene) %>% subset(!(SFE_1 %in% durgTargetGenes$Gene & SFE_2 %in% durgTargetGenes$Gene))


sigselectRes <- sigselectRes %>% inner_join(geneInfo[, c('symbol', 'location')], by = join_by(SFE_1 == symbol)) %>% 
  rename(SFE_1_location = location) %>% inner_join(geneInfo[, c('symbol', 'location')], by = join_by(SFE_2 == symbol)) %>% 
  rename(SFE_2_location = location) 


# > dim(sigselectRes)
# [1] 827  28

######################
# # Filtering
# sigselectRes <- sigselectRes %>% subset(SFE_1_location != SFE_2_location)
# # > dim(sigselectRes)
# # [1] 776   28
# 
# # Filtering
# sigselectRes <- sigselectRes %>% mutate(SFE_1_location_1 = gsub("\\..*$", "", SFE_1_location), SFE_2_location_1 = gsub("\\..*$", "", SFE_2_location))
# 
# sigselectRes <- sigselectRes %>% subset(SFE_1_location_1 != SFE_2_location_1)
# # > dim(sigselectRes)
# # [1] 665  30
# 
# 
# # Filtering
# sigselectRes <- sigselectRes %>% mutate(SFE_1_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_1_location_1), 
#                                         SFE_2_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_2_location_1))
# 
# sigselectRes <- sigselectRes %>% subset(SFE_1_location_2 != SFE_2_location_2)
# # > dim(sigselectRes)
# # [1] 339  32

######################
# Filtering
sigselectRes <- sigselectRes %>% subset(!((SFE_1_location == SFE_2_location) & 
                                            (!int_type %in% c('AMP - DEL', 'AMP - MUT', 'MUT - MUT'))))
# > dim(sigselectRes)
# [1] 776  30




# Filtering
sigselectRes <- sigselectRes %>% mutate(SFE_1_location_1 = gsub("\\..*$", "", SFE_1_location), SFE_2_location_1 = gsub("\\..*$", "", SFE_2_location))

sigselectRes <- sigselectRes %>% subset(!((SFE_1_location_1 == SFE_2_location_1) & 
                                            (!int_type %in% c('AMP - DEL', 'AMP - MUT', 'MUT - MUT'))))
# > dim(sigselectRes)
# [1] 665  30


# Filtering
sigselectRes <- sigselectRes %>% mutate(SFE_1_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_1_location_1),
                                        SFE_2_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_2_location_1),
                                        SFE_1_location_2_subloca = as.numeric(stringr::str_split(SFE_1_location_1, pattern = 'p|q', simplify = T)[, 2]), 
                                        SFE_2_location_2_subloca = as.numeric(stringr::str_split(SFE_2_location_1, pattern = 'p|q', simplify = T)[, 2]))

sigselectRes <- sigselectRes %>% subset(!((SFE_1_location_2 == SFE_2_location_2) & 
                                            (abs(SFE_2_location_2_subloca - SFE_1_location_2_subloca) < 6) & 
                                            (!int_type %in% c('AMP - DEL', 'AMP - MUT', 'MUT - MUT'))))
# > dim(sigselectRes)
# [1] 446  34

# 
# write.table(x = sigselectRes, row.names = F, col.names = T, sep = '\t', quote = F,
#             file = '/result/Section2/rhoDrugTagetMe.txt')

######################
library(ComplexHeatmap)


showSigselectRes <- sigselectRes %>% subset(freq_1 > 0.05 | freq_2 > 0.05 | SFE_1 %in% c('ARHGEF1', 'ARHGDIB'))
meGenePairs <- stringr::str_split(showSigselectRes$name, pattern = ' - ', simplify = TRUE) %>% as.data.frame() %>% rename(Gene1 = V1, Gene2 = V2)
showSigselectRes <- cbind.data.frame(meGenePairs, showSigselectRes)

showSigselectRes <- showSigselectRes %>% mutate(Gene_1 = ifelse(SFE_1 %in% durgTargetGenes$Gene, Gene1, Gene2), 
                                                Gene_2 = ifelse(SFE_1 %in% durgTargetGenes$Gene, Gene2, Gene1)) %>% 
  mutate(APC = ifelse(direction == 'CO', -APC, APC))



showData <- showSigselectRes %>% tidyr::pivot_wider(id_cols = Gene_1, names_from = Gene_2, values_from = APC) %>% 
  column_to_rownames(var = 'Gene_1')


colNames <- stringr::str_split(string = colnames(showData), pattern = '_', simplify = T) %>% as.data.frame() %>% 
  rename(Approved.symbol = V1, altStatus = V2) %>% 
  left_join(rhoFamilies[, c('Approved.symbol', 'Class')], by = join_by(Approved.symbol)) %>% 
  mutate(altName = paste(Approved.symbol, altStatus, sep = '_')) %>% arrange(Class)

showData <- showData[, colNames$altName]



colors<- circlize::colorRamp2(breaks = c(-0.2, 0, 0.01), colors = c("#8BB77B", "white","#80539A"))

p1 <- Heatmap(showData, cluster_rows = FALSE, cluster_columns = FALSE, na_col = "white", 
                              show_heatmap_legend = F, col = colors,)

lgd <- list(Legend(title = "SELECT score", 
                   col_fun = colors, 
                   at = c(-0.2, -0.1, 0, 0.005, 0.01), 
                   direction = "horizontal"))

pdf(file = '/result/Section2/rhoDrugTagetMe.pdf')
draw(p1, annotation_legend_list = lgd, annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()


