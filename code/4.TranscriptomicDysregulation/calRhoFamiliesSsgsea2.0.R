
library('ssGSEA2')

# setwd('/data/ssGSEAdata')
# 
# files <- list.files(path = ".", pattern = '.gct', recursive = F)
# 
# 
# ssGseaRes <- lapply(files, function(filename){
#   
#   dir.create(path = strsplit(filename, split = '_')[[1]][1])
#   
#   res = run_ssGSEA2(input.ds = filename,
#                     output.prefix = strsplit(filename, split = '_')[[1]][1],
#                     gene.set.databases = "rhoFamilies.gmt",
#                     output.directory = paste0('./', strsplit(filename, split = '_')[[1]][1]),
#                     statistic = "area.under.RES",
#                     output.score.type = "NES")  
#   
#   return(res)
# })
# 
# 
# 
# saveRDS(ssGseaRes, file = '/data/tcgaRhoFamiliesSsgsea2.0InterRes.rds')
# 

setwd('/data/ssGSEAdata')
files <- list.files(path = ".")
files <- files[-grep(pattern = "\\.(gct|gmt|log)$", files)]


tcgaSsgseaRes <- lapply(files, function(file){
  
  setwd(dir = paste0('/data/ssGSEAdata/', file))
  
  scores <- parse_gctx(fname = paste0(file, '-scores.gct'))@mat
  scores <- reshape::melt(data = scores) %>% rename(Groups = X1, SampleIDs = X2, Scores = value)
  
  pvalues <- parse_gctx(fname = paste0(file, '-pvalues.gct'))@mat
  pvalues <- reshape::melt(data = pvalues) %>% rename(Groups = X1, SampleIDs = X2, Pvalues = value)  
  
  fdr.pvalues <- parse_gctx(fname = paste0(file, '-fdr-pvalues.gct'))@mat
  fdr.pvalues <- reshape::melt(data = fdr.pvalues) %>% rename(Groups = X1, SampleIDs = X2, FDR = value)  
  
  res <- Reduce(function(x, y) inner_join(x, y, by = join_by(Groups, SampleIDs)), list(scores, pvalues, fdr.pvalues))
  
  return(res)
})

tcgaSsgseaRes <- do.call(rbind, tcgaSsgseaRes)


tcgaSsgseaRes <- tcgaSsgseaRes %>% mutate(PATIENT_BARCODE = substr(SampleIDs, 1, 12))
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

tcgaSsgseaRes <- tcgaSsgseaRes %>% left_join(tcgaPanCanSamples[, c('PATIENT_BARCODE', 'DISEASE')], by = 'PATIENT_BARCODE')


saveRDS(tcgaSsgseaRes, file = '/data/tcgaRhoFamiliesSsgsea2.0Res.rds')
