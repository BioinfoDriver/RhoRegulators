
#########################
load(file = '/data/tcgaRhoFamiliesDiffExp.RData') 
# deseq2DiffExp, limmaDiffExp

##########
deseq2DiffExp <- deseq2DiffExp %>% mutate(Class = recode(Class, 
  'Rho_GDI' = 'GDI', 'RhoGAP_dom' = 'GAP', 'Small_GTPase_Rho' = 'Rho', 'DH_domain' = 'GEF', 'DOCKER_dom' = 'GEF'))

write.table(x = deseq2DiffExp, file = '/result/Section3/rhoFamiliesDiffExp.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

##########
rhoGroupDysStat <- deseq2DiffExp %>% group_by(DISEASE, Class) %>% count(dysLabel)

write.table(x = rhoGroupDysStat, file = '/result/Section3/rhoGroupDisDysStat.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

##########

rhoDysStat <- deseq2DiffExp %>% group_by(Class, Symbol) %>% count(dysLabel)

write.table(x = rhoDysStat, file = '/result/Section3/rhoDysStat.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


######################### validation
load(file = '/data/difgene_GPL570.RData') 
rhoFamilies <- readRDS(file = '/data/rhoFamilies.rds')

deseq2DiffExp <- lapply(difexp1, rownames_to_column, var = 'Symbol') %>% bind_rows(.id = "DISEASE") %>% 
  mutate(DISEASE = stringr::str_to_title(DISEASE), dysLabel = ifelse(adj.P.Val < 0.05 & logFC > 0.585, 'Up', 
                                                            ifelse(adj.P.Val < 0.05 & logFC < -0.585, 'Down', 'Neutral'))) %>% 
  inner_join(rhoFamilies, by = c('Symbol' = 'Approved.symbol'))


deseq2DiffExp <- deseq2DiffExp %>% 
  mutate(Class = recode(Class, 'Rho_GDI' = 'GDI', 'RhoGAP_dom' = 'GAP', 
                        'Small_GTPase_Rho' = 'Rho', 'DH_domain' = 'GEF', 'DOCKER_dom' = 'GEF'))

write.table(x = deseq2DiffExp, file = '/result/Section3/validation/rhoFamiliesDiffExp.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)



##########
rhoGroupDysStat <- deseq2DiffExp %>% group_by(DISEASE, Class) %>% count(dysLabel)

write.table(x = rhoGroupDysStat, file = '/result/Section3/validation/rhoGroupDisDysStat.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

##########

rhoDysStat <- deseq2DiffExp %>% group_by(Class, Symbol) %>% count(dysLabel)

write.table(x = rhoDysStat, file = '/result/Section3/validation/rhoDysStat.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)




#########################
tcgaSsgseaRes <- readRDS(file = '/data/tcgaRhoFamiliesSsgsea2.0Res.rds')
load( file = '/data/tcgaExpSams.RData')

pairedSamSsgseaScore <- tcgaSsgseaRes %>% subset(PATIENT_BARCODE %in% pairedSamsPats & substr(SampleIDs, 14, 15) != '06') %>% 
  mutate(tumorTypes = ifelse(substr(SampleIDs, 14, 15) == '11', 'Tumor', 'Normal'))

PairedTTest <- function(ssgseaScore){
  
  ssgseaScore <- ssgseaScore %>% as.data.frame() %>% select(Scores, PATIENT_BARCODE, tumorTypes) %>% 
    reshape(direction = "wide", idvar = "PATIENT_BARCODE", timevar = "tumorTypes")
  
  res <- t.test(Pair(Scores.Tumor, Scores.Normal) ~ 1, data = ssgseaScore)
  res <- data.frame(diffMeans = res$estimate, pValues = res$p.value)
  
  return(res)
  
}

pairedTTestRes <- pairedSamSsgseaScore %>% group_by(Groups, DISEASE) %>% do(PairedTTest(.)) 
pairedTTestRes$FDR <- p.adjust(pairedTTestRes$pValues, method = 'BH')


write.table(x = pairedTTestRes, file = '/result/Section3/rhoFamiliesSsgseaDiff.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

