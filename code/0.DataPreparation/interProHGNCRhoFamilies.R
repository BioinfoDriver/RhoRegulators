library('dplyr')

######################
hgncRhoFamilies <- readRDS(file = '/data/HGNCRhoFamilies.rds')
interProRhoFamilies <- readRDS(file = '/data/InterProRhoFamilies.rds')
geneInfo <- read.csv(file= '/data/gene_with_protein_product.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)


# setdiff(interProRhoFamilies$ApprovedSymbol, geneInfo$symbol)
# "BARGIN"
# BARGIN probably is SH3BP1
# We describe a brain-specific RhoGAP splice variant, BARGIN (BGIN), 
# which comprises a combination of BAR, GAP, and partial CIN phosphatase domains spliced from adjacent SH3BP1 and CIN gene loci.


interProRhoFamilies <- interProRhoFamilies %>% inner_join(geneInfo[, c('symbol', 'entrez_id', 'locus_type')], by = c('ApprovedSymbol' = 'symbol'))


# subset(hgncRhoFamilies, Approved.symbol %in% setdiff(hgncRhoFamilies$Approved.symbol, interProRhoFamilies$ApprovedSymbol),
#        select = c('Approved.symbol', 'Locus.type', 'classesOfRho'))
# 
# Approved.symbol                Locus.type                        classesOfRho
#       ARHGAP16P                pseudogene         RhoGTPaseActivatingProteins
#          ALS2CL gene with protein product RhoGuanineNucleotideExchangeFactors
#       ARHGEF34P                pseudogene RhoGuanineNucleotideExchangeFactors
#        ARHGEF35 gene with protein product RhoGuanineNucleotideExchangeFactors


# subset(interProRhoFamilies, ApprovedSymbol %in% setdiff(interProRhoFamilies$ApprovedSymbol, hgncRhoFamilies$Approved.symbol),
#        select = c('ApprovedSymbol', 'Class'))
# 
# ApprovedSymbol            Class
#          MYO9A       RhoGAP_dom
#         PIK3R2       RhoGAP_dom
#         PIK3R1       RhoGAP_dom
#         INPP5B       RhoGAP_dom
#           OCRL       RhoGAP_dom
#          MYO9B       RhoGAP_dom
#         RALBP1       RhoGAP_dom
#          SYDE2       RhoGAP_dom
#          SYDE1       RhoGAP_dom
#        DEPDC1B       RhoGAP_dom
#          ARAP3       RhoGAP_dom
#          ARAP2       RhoGAP_dom
#          ARAP1       RhoGAP_dom
#        RACGAP1       RhoGAP_dom
#          RHOT2 Small_GTPase_Rho
#          RHOT1 Small_GTPase_Rho


######################
hgncRhoFamilies$Class[hgncRhoFamilies$Group.name == 'RhoGDI protein family'] <- 'Rho_GDI'
hgncRhoFamilies$Class[hgncRhoFamilies$Group.name == 'Rho family GTPases'] <- 'Small_GTPase_Rho'
hgncRhoFamilies$Class[hgncRhoFamilies$Group.name == 'Dbl family Rho GEFs'] <- 'DH_domain'
hgncRhoFamilies$Class[hgncRhoFamilies$Group.name == 'DOCK family Rho GEFs'] <- 'DOCKER_dom'
hgncRhoFamilies$Class[hgncRhoFamilies$Group.name == 'Rho GTPase activating proteins'] <- 'RhoGAP_dom'
hgncRhoFamilies$Class[hgncRhoFamilies$Group.name == 'AH/BAR family Rho GTPase activating proteins'] <- 'RhoGAP_dom'


######################
rhoFamilies <- interProRhoFamilies[, c('ApprovedSymbol', 'entrez_id', 'locus_type', 'Class')] %>% 
  rename(Approved.symbol = ApprovedSymbol, NCBI.Gene.ID = entrez_id, Locus.type = locus_type) %>% 
  rbind.data.frame(hgncRhoFamilies[, c('Approved.symbol', 'NCBI.Gene.ID', 'Locus.type', 'Class')]) %>% 
  distinct(Approved.symbol, NCBI.Gene.ID, Locus.type, Class)


# table(rhoFamilies$Class, rhoFamilies$Locus.type)
#                  gene with protein product pseudogene
# DH_domain                               72          1
# DOCKER_dom                              11          0
# Rho_GDI                                  3          0
# RhoGAP_dom                              65          1
# Small_GTPase_Rho                        22          0


# pseudogene, Mitochondrial Rho proteins(atypical Rho GTPases)
rhoFamilies <- subset(rhoFamilies, Locus.type != 'pseudogene' & Approved.symbol != 'RHOT2' & Approved.symbol != 'RHOT1')

saveRDS(rhoFamilies, file = '/data/rhoFamilies.rds')

