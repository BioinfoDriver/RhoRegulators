
setwd('/data/OncoKB/')

cancerGeneList <- read.csv(file = 'cancerGeneList.tsv', header = T, sep = '\t', stringsAsFactors = F)
load(file = '/data/rhoFamiliesMutSig2CV.RData')

mutSigGeneQvalue <- mutSigGeneQvalue %>% melt() %>% subset(value < 0.25) %>% arrange(desc(value))
mutSigGene <- unique(mutSigGeneQvalue$gene)

# subset(cancerGeneList, Hugo.Symbol %in% mutSigGene)
# Hugo.Symbol Entrez.Gene.ID  GRCh37.Isoform GRCh37.RefSeq  GRCh38.Isoform GRCh38.RefSeq Is.Oncogene Is.Tumor.Suppressor.Gene
# 75       PIK3R1           5295 ENST00000274335   NM_181523.2 ENST00000521381   NM_181523.2          No                      Yes
# 265        RAC1           5879 ENST00000356142   NM_018890.3 ENST00000356142   NM_018890.3         Yes                       No
# 271        RHOA            387 ENST00000418115   NM_001664.2 ENST00000418115   NM_001664.2         Yes                       No
# 344      PIK3R2           5296 ENST00000222254   NM_005027.3 ENST00000222254   NM_005027.3          No                      Yes
# 381    ARHGAP35           2909 ENST00000404338   NM_004491.4 ENST00000404338   NM_004491.4         Yes                      Yes
# 589        SOS1           6654 ENST00000402219   NM_005633.3 ENST00000402219   NM_005633.3         Yes                       No
# X..of.occurrence.within.resources..Column.J.P. OncoKB.Annotated MSK.IMPACT MSK.HEME FOUNDATION.ONE FOUNDATION.ONE.HEME Vogelstein
# 75                                               7              Yes        Yes      Yes            Yes                 Yes        Yes
# 265                                              5              Yes        Yes      Yes            Yes                  No         No
# 271                                              5              Yes        Yes      Yes             No                 Yes         No
# 344                                              4              Yes        Yes      Yes             No                 Yes         No
# 381                                              3              Yes        Yes       No             No                  No         No
# 589                                              3              Yes        Yes      Yes             No                  No         No
# COSMIC.CGC..v99.                                           Gene.Aliases
# 75               Yes                                        GRB1, p85-ALPHA
# 265              Yes                                 Rac-1, TC-25, p21-Rac1
# 271              Yes                             ARH12, ARHA, RHOH12, Rho12
# 344               No                                                   P85B
# 381              Yes GRF-1, GRLF1, KIAA1722, P190A, p190ARhoGAP, p190RhoGAP
# 589               No                                                  GINGF