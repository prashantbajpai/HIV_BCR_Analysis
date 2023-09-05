library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(writexl)
library(patchwork)
source('/Users/prashant/R_themes/theme_publication.R')
source('/Users/prashant/R_themes/theme_publication_rnaseq.R')

mapped_mabs <- read_xlsx("mabs_mapping_50_perc_cutoff_with_indel.xlsx", sheet = 1)
mapped_mabs$vgene_ins[is.na(mapped_mabs$vgene_ins)] <- "NA"
mapped_mabs$vgene_del[is.na(mapped_mabs$vgene_del)] <- "NA"
mapped_mabs$isIndel <- ifelse(grepl("do not cause frameshift", mapped_mabs$vgene_ins)|grepl("do not cause frameshift", mapped_mabs$vgene_del), "Has Indel", "No Indel")

savedat <- split(mapped_mabs, f = mapped_mabs$mab_id)
write_xlsx(savedat, "suppliment/mapped_mabs_seq_info.xlsx")

summary_mabs_heavy <- mapped_mabs %>%
  group_by(sample_id, mab_id) %>%
  summarise(meanSHM = mean(mu_freq), counts = n())
summary_mabs$isIndel <- "No Indel"  

for (i in 1:nrow(summary_mabs)) {
  test <- mapped_mabs
  test <- test[test$sample_id == summary_mabs$sample_id[i],]
  test <- test[test$mab_id == summary_mabs$mab_id[i],]
  if(any(grepl("Has Indel", test$isIndel))){summary_mabs$isIndel[i] <- "Has Indel"}
}

cols <- c("#9d9dff", "#c7e9c0", "#ff9d9d", "#1414ff", "#238b45", "#ff1414")

ggplot(summary_mabs, aes(x = mab_id, y = meanSHM, color = sample_id)) + geom_point(aes(size = counts, shape = isIndel), stroke = 1) + theme_Publication() + scale_colour_manual(values = cols, name = "") + scale_shape_manual(name = "",values = c(2, 1)) + scale_size_continuous(breaks = c(50, seq(0, 10000, 2000)), labels = c(50, seq(0, 10000, 2000)), range = c(1, 6)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank()) 
  
ggsave("plot_ver1/mabs_mapped.png", width = 15, height = 7, units = "in")
ggsave("plot_ver1/mabs_mapped.pdf", width = 15, height = 7, units = "in")

##80% cutoff

mapped_mabs_80 <- readRDS("rds/bnabs_80perc_match_alignment.rds")












