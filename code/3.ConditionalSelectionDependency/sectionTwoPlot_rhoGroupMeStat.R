
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
################
rhoMe <- read.table(file = '/result/Section2/rhoMe.txt', sep = '\t', header = T, stringsAsFactors = F)

rhoGroupMeStat <- rhoMe %>% mutate(SFE_1_Class = str_replace_all(SFE_1_Class, 
                              c("RhoGAP_dom" = "GAP", "Small_GTPase_Rho" = "Rho", "Rho_GDI" = "GDI", "DH_DOCKER_domain" = 'GEF')), 
                 SFE_2_Class = str_replace_all(SFE_2_Class, 
                                               c("RhoGAP_dom" = "GAP", "Small_GTPase_Rho" = "Rho", "Rho_GDI" = "GDI", "DH_DOCKER_domain" = 'GEF'))) %>% 
  mutate(MeType = paste(SFE_1_Class, SFE_2_Class, sep = '_')) %>% group_by(MeType) %>% count(direction) %>% 
  dcast(MeType~direction, value.var = 'n', fill = 0)

################
withinGroup <- rhoGroupMeStat %>% subset(MeType %in% c('Rho_Rho', 'GAP_GAP', 'GEF_GEF', 'GDI_GDI')) %>% 
  summarize(sum_CO = sum(CO), sum_ME = sum(ME))

betweenGroup <- rhoGroupMeStat %>% subset(!MeType %in% c('Rho_Rho', 'GAP_GAP', 'GEF_GEF', 'GDI_GDI')) %>% 
  summarize(sum_CO = sum(CO), sum_ME = sum(ME))


################
rhoMeSum <- rbind.data.frame(withinGroup, betweenGroup) %>% mutate(Group = c('WithinGroup', 'betweenGroup'))

fisher.test(rhoMeSum[, c('sum_ME', 'sum_CO')])$p.value
# 0.03224476
fisher.test(rhoMeSum[, c('sum_ME', 'sum_CO')])$estimate
# 2.156768
chisq.test(rhoMeSum[, c('sum_CO', 'sum_ME')])$p.value
# 0.03591509


plot1 <- ggplot(data = melt(rhoMeSum, id = 'Group'), aes(x = Group, y = value, fill = variable)) +
  geom_bar(stat = "identity") + geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 4) +
  labs(x = NULL, y = "Number of ME/CO variations", fill = "Patterns") +
  scale_fill_brewer(palette = "Set1") + theme(legend.position = "bottom")

ggsave(plot1, filename = '/result/Section2/rhoGroupMeStat.pdf')


