library('dplyr')

load(file = '/data/tcgaPairedSamsDiffExp.RData')
# pairedSamsDeseq2DiffExp, pairedSamsLimmaDiffExp
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

#########################
# deseq2DiffExp <- subset(pairedSamsDeseq2DiffExp, abs(log2FoldChange) > 1.0 & padj < 0.05) %>% split.data.frame(f = ~ DISEASE)
# limmaDiffExp <- subset(pairedSamsLimmaDiffExp, abs(logFC) > 1.0 & adj.P.Val < 0.05) %>% split.data.frame(f = ~ DISEASE)

# sapply(names(deseq2DiffExp), function(disease){
#   
#   comDiffNum <- length(intersect(deseq2DiffExp[[disease]]$Symbol, limmaDiffExp[[disease]]$Symbol))
#   deseqUniDiffNum <- length(setdiff(deseq2DiffExp[[disease]]$Symbol, limmaDiffExp[[disease]]$Symbol))
#   limmaUniDiffNum <- length(setdiff(limmaDiffExp[[disease]]$Symbol, deseq2DiffExp[[disease]]$Symbol))
#   
#   return(c('comDiffNum' = comDiffNum, 'deseqUniDiffNum' = deseqUniDiffNum, 'limmaUniDiffNum' = limmaUniDiffNum))
# })
# 
#                 BLCA BRCA COAD ESCA HNSC KICH KIRC KIRP LIHC LUAD LUSC PRAD STAD THCA
# comDiffNum      2139 2131 2048 1710 1439 2961 2676 2249 1042 2344 3843  802 1407 1250
# deseqUniDiffNum 1948 1804 2213 2151 2030 2371 1952 1976 2195 1689 2160 1389 1782 1643
# limmaUniDiffNum  186   16  138 4404  909 1376   25   52 1176  121  216   12  738  125


#########################
# missing-----did not contain in firehose expression matrix
# 343578
# ARHGAP40
deseq2DiffExp <- pairedSamsDeseq2DiffExp %>% inner_join(rhoFamilies[, c('NCBI.Gene.ID', 'Class')], 
                                                        by = join_by(geneID == NCBI.Gene.ID), relationship = "many-to-many")
limmaDiffExp <- pairedSamsLimmaDiffExp %>% inner_join(rhoFamilies[, c('NCBI.Gene.ID', 'Class')], 
                                                      by = join_by(geneID == NCBI.Gene.ID), relationship = "many-to-many")


deseq2DiffExp <- deseq2DiffExp %>% mutate(dysLabel = ifelse(log2FoldChange > 1.0 & padj < 0.05, 'Up', 
                                                            ifelse(log2FoldChange < -1.0 & padj < 0.05, 'Down', 'Neutral')))

limmaDiffExp <- limmaDiffExp %>% mutate(dysLabel = ifelse(logFC > 1.0 & adj.P.Val < 0.05, 'Up', 
                                                            ifelse(logFC < -1.0 & adj.P.Val < 0.05, 'Down', 'Neutral')))


save(deseq2DiffExp, limmaDiffExp, file = '/data/tcgaRhoFamiliesDiffExp.RData') 





