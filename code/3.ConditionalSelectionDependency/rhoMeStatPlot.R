
library(dplyr)
######################
selectRes <- readRDS(file = '/data/rhoSelectScoreRes.RData')
geneInfo <- read.csv(file= '/data/gene_with_protein_product.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')


rhoFamilies <- rhoFamilies %>% distinct(Approved.symbol, .keep_all = TRUE) %>% 
  mutate(Class = ifelse(Class %in% c('DH_domain', 'DOCKER_dom'), 'DH_DOCKER_domain', Class))

rhoFamilies <- geneInfo %>% select(symbol, location) %>% inner_join(rhoFamilies, by = join_by(symbol == Approved.symbol))

sigselectRes <- selectRes %>% subset(FDR & APC_good) %>% mutate(SFE_1 = gsub("_.*$", "", SFE_1), SFE_2 = gsub("_.*$", "", SFE_2)) %>% 
  subset(SFE_1 %in% rhoFamilies$symbol & SFE_2 %in% rhoFamilies$symbol)

sigselectRes <- sigselectRes %>% left_join(rhoFamilies[, c('symbol', 'location', 'Class')], by = join_by(SFE_1 == symbol)) %>% 
  rename(SFE_1_location = location, SFE_1_Class = Class) %>% left_join(rhoFamilies[, c('symbol', 'location', 'Class')], by = join_by(SFE_2 == symbol)) %>% 
  rename(SFE_2_location = location, SFE_2_Class = Class) 

# > dim(sigselectRes)
# [1] 347  30


##############
# # Filtering
# sigselectRes <- sigselectRes %>% subset(SFE_1_location != SFE_2_location)
# # > dim(sigselectRes)
# # [1] 327  30
# 
# # Filtering
# sigselectRes <- sigselectRes %>% mutate(SFE_1_location_1 = gsub("\\..*$", "", SFE_1_location), SFE_2_location_1 = gsub("\\..*$", "", SFE_2_location))
# 
# sigselectRes <- sigselectRes %>% subset(SFE_1_location_1 != SFE_2_location_1)
# # > dim(sigselectRes)
# # [1] 270  32
# 
# 
# # Filtering
# sigselectRes <- sigselectRes %>% mutate(SFE_1_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_1_location_1),
#                                         SFE_2_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_2_location_1))
# 
# sigselectRes <- sigselectRes %>% subset(SFE_1_location_2 != SFE_2_location_2)
# # > dim(sigselectRes)
# # [1] 121  34
# 
# 
# # write.table(x = sigselectRes, row.names = F, col.names = T, sep = '\t', quote = F,
# #             file = '/result/Section2/rhoMe.txt')
# 

##############

# Filtering
sigselectRes <- sigselectRes %>% subset(!((SFE_1_location == SFE_2_location) & (int_type != 'MUT - MUT')))
# > dim(sigselectRes)
# [1] 327  30

# Filtering
sigselectRes <- sigselectRes %>% mutate(SFE_1_location_1 = gsub("\\..*$", "", SFE_1_location), SFE_2_location_1 = gsub("\\..*$", "", SFE_2_location))

sigselectRes <- sigselectRes %>% subset(!((SFE_1_location_1 == SFE_2_location_1) & (int_type != 'MUT - MUT')))
# > dim(sigselectRes)
# [1] 270  32


# Filtering
sigselectRes <- sigselectRes %>% mutate(SFE_1_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_1_location_1),
                                        SFE_2_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_2_location_1),
                                        SFE_1_location_2_subloca = as.numeric(stringr::str_split(SFE_1_location_1, pattern = 'p|q', simplify = T)[, 2]), 
                                        SFE_2_location_2_subloca = as.numeric(stringr::str_split(SFE_2_location_1, pattern = 'p|q', simplify = T)[, 2]))

sigselectRes <- sigselectRes %>% subset(!((SFE_1_location_2 == SFE_2_location_2) & 
 (abs(SFE_2_location_2_subloca - SFE_1_location_2_subloca) < 6) & (int_type != 'MUT - MUT')))
# > dim(sigselectRes)
# [1] 170  36


# write.table(x = sigselectRes, row.names = F, col.names = T, sep = '\t', quote = F,
#             file = '/result/Section2/rhoMe.txt')


######################
library(ComplexHeatmap)


# showSigselectRes <- sigselectRes %>% subset(freq_1 > 0.03 | freq_2 > 0.03)

selectGeneA <- sigselectRes %>% group_by(SFE_1_Class) %>% arrange(desc(freq_1)) %>% slice_head(n = 6)
selectGeneB <- sigselectRes %>% group_by(SFE_2_Class) %>% arrange(desc(freq_2)) %>% slice_head(n = 7)
showSigselectRes <- sigselectRes %>% subset(name %in% selectGeneA$name | name %in% selectGeneB$name)


meGenePairs <- stringr::str_split(showSigselectRes$name, pattern = ' - ', simplify = TRUE) %>% as.data.frame() %>% rename(Gene1 = V1, Gene2 = V2)
showSigselectRes <- cbind.data.frame(meGenePairs, showSigselectRes)

showSigselectRes <- showSigselectRes %>% mutate(APC = ifelse(direction == 'CO', -APC, APC))

showData <- showSigselectRes %>% pivot_wider(id_cols = Gene1, names_from = Gene2, values_from = APC) %>% 
  column_to_rownames(var = 'Gene1')


rowNames <- stringr::str_split(string = rownames(showData), pattern = '_', simplify = T) %>% as.data.frame() %>% 
  rename(symbol = V1, altStatus = V2) %>% left_join(rhoFamilies[, c('symbol', 'Class')], by = join_by(symbol)) %>% 
  mutate(altName = paste(symbol, altStatus, sep = '_')) %>% arrange(Class)

colNames <- stringr::str_split(string = colnames(showData), pattern = '_', simplify = T) %>% as.data.frame() %>% 
  rename(symbol = V1, altStatus = V2) %>% left_join(rhoFamilies[, c('symbol', 'Class')], by = join_by(symbol)) %>% 
  mutate(altName = paste(symbol, altStatus, sep = '_')) %>% arrange(Class)

showData <- showData[rowNames$altName, colNames$altName]

colors<- circlize::colorRamp2(breaks = c(-0.2, 0, 0.006), colors = c("#8BB77B", "white","#80539A"))

p1 <- Heatmap(showData, cluster_rows = FALSE, cluster_columns = FALSE, na_col = "white", 
              show_heatmap_legend = F, col = colors, border = 'gray')

lgd <- list(Legend(title = "SELECT score", 
                   col_fun = colors, 
                   at = c(-0.2, -0.1, 0, 0.003, 0.006), 
                   direction = "horizontal"))

pdf(file = '/result/Section2/rhoMe_2.pdf')
draw(p1, annotation_legend_list = lgd, annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()



