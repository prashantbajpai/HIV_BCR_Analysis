library(alakazam)
library(shazam)
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(writexl)
library(ggrepel)
library(patchwork)
source('/Users/pbajpai/Documents/Rthemes/theme_publication.R')

########
all_files <- list.files("public_mabs_data/", pattern = "set")
setwd("public_mabs_data/")
hiv_mabs1 <- lapply(all_files, function(x){
  dat <- read_xls(x)
  dat <- dat[!is.na(dat$`J-GENE and allele`),]
  dat$`CDR3-IMGT length` <- as.numeric(dat$`CDR3-IMGT length`)
  return(dat)
})
setwd("..")
hiv_mabs1 <- bind_rows(hiv_mabs1)

geneusage_transformation = function(x){
  out = gsub('\\*.*', '', x)
}
clean.ambgs.calls = function(df){
  df$v_gene = gsub('Homsap ', '', df$`V-GENE and allele`)
  df$v_gene = gsub('\\*.*', '', df$v_gene)
  df$v_family = gsub('-.*', '', df$v_gene)
  df$j_gene = gsub('Homsap ', '', df$`J-GENE and allele`)
  df$j_gene = gsub('\\*.*', '', df$j_gene)
  df$j_family = gsub('-.*', '', df$j_gene)
  return(df)
}

hiv_mabs1 <- clean.ambgs.calls(hiv_mabs1)
hiv_mabs1 = hiv_mabs1[c("Sequence ID", "AA JUNCTION", "v_gene", "j_gene")]
hiv_mabs1$`AA JUNCTION` <- gsub(" .*", "",hiv_mabs1$`AA JUNCTION`)
hiv_mabs1$cdrl3_aa = stringr::str_sub(hiv_mabs1$`AA JUNCTION`, 2, -2)
hiv_mabs1$cdrl3_aa_length = nchar(hiv_mabs1$cdrl3_aa)

hiv_mabs1$clone_id = paste(hiv_mabs1$v_gene, hiv_mabs1$j_gene,
                           hiv_mabs1$cdrl3_aa_length, sep = ',')

hiv_mabs1$seq_id = gsub('_.*', '', hiv_mabs1$`Sequence ID`)
hiv_mabs1$seq_id = gsub(' .*', '', hiv_mabs1$seq_id)
hiv_mabs1 = hiv_mabs1[c("v_gene", "j_gene", "cdrl3_aa_length", 
                        "cdrl3_aa", "clone_id", "seq_id")]
#####
hiv_mabs2 = read_xlsx('ab_metadata/ab_germlines-updated_PB.xlsx', sheet = 1)
hiv_mabs2 = hiv_mabs2[!is.na(hiv_mabs2$`Light CDR3 sequence`),]
hiv_mabs2 = hiv_mabs2[!is.na(hiv_mabs2$`Light V`),]
hiv_mabs2 = hiv_mabs2[!is.na(hiv_mabs2$`Light J`),]
hiv_mabs2$v_gene = paste('IG', hiv_mabs2$`Light V`, sep = '')
hiv_mabs2$j_gene = paste('IG', hiv_mabs2$`Light J`, sep = '')
hiv_mabs2$cdrl3_aa = hiv_mabs2$`Light CDR3 sequence`
hiv_mabs2$cdrl3_aa_length = nchar(hiv_mabs2$cdrl3_aa)
hiv_mabs2$check = hiv_mabs2$cdrl3_aa_length - as.numeric(hiv_mabs2$`Light CDR3 length`)
hiv_mabs2 %>%
  mutate(across(c("v_gene", "j_gene"), 
                geneusage_transformation)) -> hiv_mabs2
hiv_mabs2$clone_id = paste(hiv_mabs2$v_gene, hiv_mabs2$j_gene,
                           hiv_mabs2$cdrl3_aa_length, sep = ',')
hiv_mabs2$seq_id = hiv_mabs2$Antibody
hiv_mabs2 = hiv_mabs2[c("v_gene", "j_gene", "cdrl3_aa_length", 
                        "cdrl3_aa", "clone_id", "seq_id")]
#####

hiv_mabs = rbind(hiv_mabs1, hiv_mabs2)
rm(hiv_mabs1, hiv_mabs2)

#hiv_mabs = hiv_mabs[!duplicated(hiv_mabs$cdrl3_aa),]

#filter the mabs mapped in heavy chain
mapped_mabs <- read_xlsx("mabs_mapping_50_perc_cutoff_with_indel.xlsx", sheet = 1)
mapped_mabs <- unique(mapped_mabs$mab_id)

hiv_mabs <- hiv_mabs[hiv_mabs$seq_id %in% mapped_mabs,]
hiv_mabs_list = split(hiv_mabs, f = hiv_mabs$clone_id)

#####
#light chain data
abstar_files_lc = list.files('abstar_output_LC/', pattern = '.txt')
setwd('abstar_output_LC/')
filnames = gsub('.txt', '', abstar_files_lc)
abstar_dat_lc = lapply(seq_along(abstar_files_lc), function(i){
  dat = read.csv(abstar_files_lc[i])
  dat$group_id = filnames[i]
  return(dat)
})
setwd('..')
abstar_dat_lc = bind_rows(abstar_dat_lc)
abstar_dat_lc <- abstar_dat_lc[abstar_dat_lc$chain == "lambda",]

#Kappa chain
abstar_files_kc = list.files('abstar_output_KC/', pattern = '.txt')
setwd('abstar_output_KC/')
filnames = gsub('.txt', '', abstar_files_kc)
abstar_dat_kc = lapply(seq_along(abstar_files_kc), function(i){
  dat = read.csv(abstar_files_kc[i])
  dat$group_id = filnames[i]
  return(dat)
})
setwd('..')
abstar_dat_kc = bind_rows(abstar_dat_kc)
abstar_dat_kc <- abstar_dat_kc[abstar_dat_kc$chain == "kappa",]

abstar_dat_lightC <- rbind(abstar_dat_lc, abstar_dat_kc)
abstar_dat_lightC$cdrl3_aa <- translateDNA(abstar_dat_lightC$cdr3_nt, trim = T)
abstar_dat_lightC$cdrl3_aa_length <- nchar(abstar_dat_lightC$cdrl3_aa)

abstar_dat_lightC$clone_id = paste(abstar_dat_lightC$v_gene, abstar_dat_lightC$j_gene, abstar_dat_lightC$cdrl3_aa_length, sep = ',')

abstar_dat_lightC <- abstar_dat_lightC[c("v_gene", "j_gene", "cdrl3_aa_length", 
                            "cdrl3_aa", "clone_id", "seq_id")]

abstar_dat_lightC_list = split(abstar_dat_lightC, f = abstar_dat_lightC$clone_id)

common_clones = intersect(names(abstar_dat_lightC_list), names(hiv_mabs_list))
####
#mapping
#ref = hiv_mabs_list[[c("IGKV1-33,IGKJ1,5")]]
#ref = ref[!duplicated(ref$cdrl3_aa),]

#test = abstar_dat_lightC_list[["IGKV1-33,IGKJ1,5"]]

custom_fun = function(ref, mab){
  mab_align = lapply(ref, function(x){
    mymax = round(nchar(mab)*0.5)
    matched_data = matchPattern(mab, x, 
                                max.mismatch = mymax, min.mismatch = 0)
    matched_num = nmismatch(mab, matched_data)
    out = NA
    if(length(matched_num)==1){
      out = matched_num
    }
    if(length(matched_num)>1){
      out = sort(unique(matched_num))[1]
      print(paste('wtf', mab, sep = ' '))
    }
    return(out)
  })
}

cl <- makeCluster(detectCores())
clusterExport(cl, c('abstar_dat_lightC_list', 'custom_fun', 'hiv_mabs_list'))
clusterEvalQ(cl, library(Biostrings))

check = parLapply(cl,common_clones, function(x){
  seq_data = abstar_dat_lightC_list[[x]]
  mabs = hiv_mabs_list[[x]]
  mabs = mabs[!duplicated(mabs$cdrl3_aa),]
  test = seq_data$cdrl3_aa
  names(test) = seq_data$seq_id
  out = lapply(mabs$cdrl3_aa, function(mab){
    res = custom_fun(test, mab)
  })
  names(out) = mabs$cdrl3_aa
  return(out)
})
names(check) = common_clones

dat = lapply(check, function(list1){
  #print(i)
  #list1 = check[[i]]
  out = lapply(list1, function(x){
    x = bind_rows(x, .id = 'cdrl3_aa')
  })
})
dat = lapply(dat, function(x){
  out = bind_rows(x, .id = 'cdrl3_aa')
})

mabs = lapply(seq_along(check), function(i){
  print(i)
  x = check[[i]]
  x = unlist(x)
  x = as.data.frame(x)
  x$id = rownames(x)
  rownames(x) = NULL
  x$subject = names(check)[i]
  colnames(x)[1] = c('matches')
  x = x[!is.na(x$matches),]
  return(x)
})

mabs = bind_rows(mabs)
colnames(mabs)[colnames(mabs) == "subject"] <- "clone_id"

mapped_mabs_light <- merge(mabs, hiv_mabs, by = "clone_id")
colnames(mapped_mabs_light)[colnames(mapped_mabs_light) == "seq_id"] <- "mAb_id"

#merge sample info
sample_data <- rbind(abstar_dat_kc, abstar_dat_lc)
sample_data$mu_freq <- 100 - sample_data$var_identity_aa
sample_data <- sample_data[c("seq_id", "group_id", "mu_freq")]

trim_seq_id = function(x){
  x = stringi::stri_reverse(x)
  x = sub("^[^:]*:", "", x)
  x = stringi::stri_reverse(x)
  return(x)
}

sample_data$sequence_id_new <- sapply(sample_data$seq_id, trim_seq_id)

#add trimmed sequence info to mapped mabs for merging
mapped_mabs_light$id <- sub("^[^.]*.", "", mapped_mabs_light$id)
mapped_mabs_light$sequence_id_new <- sapply(mapped_mabs_light$id, trim_seq_id)

mapped_mabs_light <- merge(sample_data, mapped_mabs_light, by = "sequence_id_new")
#drop sequences from 330_2018
test <- grep("330_2018", mapped_mabs_light$group_id)
set.seed(10)
#keep 20000 rows
test <- sample(test, 99978 - 20000)
mapped_mabs_light_sampled <- mapped_mabs_light[-test,]

summary_mabs <- mapped_mabs_light_sampled %>%
  group_by(group_id, mAb_id) %>%
  summarise(meanSHM = mean(mu_freq), counts = n())

cols <- c("#9d9dff", "#c7e9c0", "#ff9d9d", "#1414ff", "#238b45", "#ff1414")

ggplot(summary_mabs, aes(x = mAb_id, y = meanSHM, color = group_id)) + geom_point(aes(size = counts), stroke = 1) + theme_Publication() + scale_colour_manual(values = cols, name = "") + scale_shape_manual(name = "",values = c(2, 1)) + scale_size_continuous(breaks = c(50, seq(0, 8000, 2000)), labels = c(50, seq(0, 8000, 2000)), range = c(1, 6)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank()) 

ggsave("plot_ver1/mabs_mapped_light.png", width = 15, height = 7, units = "in")
ggsave("plot_ver1/mabs_mapped_light.pdf", width = 15, height = 7, units = "in")



###correlation
p1 <- summary_mabs_heavy[c("mab_id", "counts")]
colnames(p1)[2] <- "counts_heavy"
p2 <- summary_mabs[c("mAb_id", "counts")]
colnames(p2) <- c("mab_id","counts_light")
plotdat <- merge(p1, p2, by = "mab_id")

plotdat <- plotdat %>%
  group_by(mab_id) %>%
  summarise(counts_heavy = sum(counts_heavy), counts_light = sum(counts_light))

ggplot(plotdat, aes(x = counts_heavy, y = counts_light)) + geom_point() + stat_cor() + geom_smooth(method = "lm")

##
#Venn
p1 <- summary_mabs_heavy[summary_mabs_heavy$sample_id == "330_2018",]
p2 <- summary_mabs[summary_mabs$group_id == "330_2018",]

plotdat <- list("mAbs mapped to \nheavy chain" = unique(p1$mab_id), 
                "mAbs mapped to \nlight chain" = unique(p2$mAb_id))

ggVennDiagram(plotdat, label_alpha = 0, edge_size = 0, set_size = 5) + scale_fill_gradient(low="grey", high = "#e41a1c") + theme(legend.position = "none")

ggsave("plot_ver1/venn_mabs_mapped_heavy_light_330_2018.pdf", width = 6, height = 5, units = "in")

p1 <- summary_mabs_heavy[summary_mabs_heavy$sample_id == "329_2018",]
p2 <- summary_mabs[summary_mabs$group_id == "329_2018",]

plotdat <- list("mAbs mapped to \nheavy chain" = unique(p1$mab_id), 
                "mAbs mapped to \nlight chain" = unique(p2$mAb_id))

ggVennDiagram(plotdat, label_alpha = 0, edge_size = 0, set_size = 5) + scale_fill_gradient(low="grey", high = "#e41a1c") + theme(legend.position = "none")

ggsave("plot_ver1/venn_mabs_mapped_heavy_light_330_2018.pdf", width = 6, height = 5, units = "in")