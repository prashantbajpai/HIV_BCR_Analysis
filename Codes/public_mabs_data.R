library(alakazam)
library(shazam)
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(writexl)
library(ggrepel)
library(patchwork)
source('/Users/prashant/R_themes/theme_publication.R')

########
hiv_mabs1 = readChangeoDb('ab_metadata/hiv_public_mabs_db-pass.tsv')
dist_ham <- distToNearest(hiv_mabs1, sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=4)

output <- findThreshold(dist_ham$dist_nearest, method="density")
#########
hiv_mabs1 = readChangeoDb('ab_metadata/hiv_public_mabs_db-pass_clone-pass_germ-pass.tsv')
hiv_mabs1 = hiv_mabs1[c("sequence_id", "junction", "germline_v_call", 
                        "germline_d_call", "germline_j_call")]
hiv_mabs1$cdrh3_aa = translateDNA(hiv_mabs1$junction, trim = T)
hiv_mabs1$cdrh3_aa_length = nchar(hiv_mabs1$cdrh3_aa)

geneusage_transformation = function(x){
  out = gsub('\\*.*', '', x)
}

hiv_mabs1 %>%
  mutate(across(c("germline_v_call", "germline_d_call", "germline_j_call"), 
                geneusage_transformation)) -> hiv_mabs1

hiv_mabs1$clone_id = paste(hiv_mabs1$germline_v_call, hiv_mabs1$germline_j_call,
                           hiv_mabs1$cdrh3_aa_length, sep = ',')

hiv_mabs1$seq_id = gsub('_.*', '', hiv_mabs1$sequence_id)
hiv_mabs1$seq_id = gsub(' .*', '', hiv_mabs1$seq_id)
hiv_mabs1 = hiv_mabs1[c("germline_v_call", "germline_j_call", "cdrh3_aa_length", 
                        "cdrh3_aa", "cdrh3_aa_length", "clone_id", "seq_id")]
#####
hiv_mabs2 = read_xlsx('ab_metadata/ab_germlines-updated_PB.xlsx', sheet = 1)
hiv_mabs2 = hiv_mabs2[!is.na(hiv_mabs2$`Heavy CDR3 sequence`),]
hiv_mabs2 = hiv_mabs2[!is.na(hiv_mabs2$`Heavy V`),]
hiv_mabs2 = hiv_mabs2[!is.na(hiv_mabs2$`Heavy J`),]
hiv_mabs2$germline_v_call = paste('IG', hiv_mabs2$`Heavy V`, sep = '')
hiv_mabs2$germline_j_call = paste('IG', hiv_mabs2$`Heavy J`, sep = '')
hiv_mabs2$cdrh3_aa = hiv_mabs2$`Heavy CDR3 sequence`
hiv_mabs2$cdrh3_aa_length = nchar(hiv_mabs2$cdrh3_aa)
hiv_mabs2$check = hiv_mabs2$cdrh3_aa_length - as.numeric(hiv_mabs2$`Heavy CDR3 length`)
hiv_mabs2 %>%
  mutate(across(c("germline_v_call", "germline_j_call"), 
                geneusage_transformation)) -> hiv_mabs2
hiv_mabs2$clone_id = paste(hiv_mabs2$germline_v_call, hiv_mabs2$germline_j_call,
                           hiv_mabs2$cdrh3_aa_length, sep = ',')
hiv_mabs2$seq_id = hiv_mabs2$Antibody
hiv_mabs2 = hiv_mabs2[c("germline_v_call", "germline_j_call", "cdrh3_aa_length", 
                        "cdrh3_aa", "cdrh3_aa_length", "clone_id", "seq_id")]
#####
hiv_mabs3 = read_xlsx('ab_metadata/ab_germlines-updated_PB.xlsx', sheet = 2)
hiv_mabs3$sequence_id = gsub('.* ', '', hiv_mabs3$sequence_id)
hiv_mabs3 = hiv_mabs3[!is.na(hiv_mabs3$cdrh3_aa),]
hiv_mabs3 = hiv_mabs3[!is.na(hiv_mabs3$v_call),]
hiv_mabs3 = hiv_mabs3[!is.na(hiv_mabs3$j_call),]
hiv_mabs3$germline_v_call = paste('IGHV', hiv_mabs3$v_call, sep = '')
hiv_mabs3$germline_j_call = paste('IGHV', hiv_mabs3$j_call, sep = '')
hiv_mabs3$cdrh3_aa_length = as.numeric(gsub(' aa', '', hiv_mabs3$cdrh3_aa_length))
hiv_mabs3$cdrh3_aa_length2 = nchar(hiv_mabs3$cdrh3_aa)
hiv_mabs3$check = hiv_mabs3$cdrh3_aa_length - hiv_mabs3$cdrh3_aa_length2
hiv_mabs3 %>%
  mutate(across(c("germline_v_call", "germline_j_call"), 
                geneusage_transformation)) -> hiv_mabs3
hiv_mabs3$clone_id = paste(hiv_mabs3$germline_v_call, hiv_mabs3$germline_j_call,
                           hiv_mabs3$cdrh3_aa_length, sep = ',')
hiv_mabs3$seq_id = hiv_mabs3$sequence_id
hiv_mabs3 = hiv_mabs3[c("germline_v_call", "germline_j_call", "cdrh3_aa_length", 
                        "cdrh3_aa", "cdrh3_aa_length", "clone_id", "seq_id")]
#####
hiv_mabs = rbind(hiv_mabs1, hiv_mabs2, hiv_mabs3)
rm(hiv_mabs1, hiv_mabs2, hiv_mabs3)

hiv_mabs <- hiv_mabs[!duplicated(hiv_mabs$seq_id),]
hiv_mabs$cdrh3_aa_length <- NULL
############

#Light chain info of mabs
all_files <- list.files("public_mabs_data/", pattern = "set")
setwd("public_mabs_data/")
hiv_mabs1_light <- lapply(all_files, function(x){
  dat <- read_xls(x)
  dat <- dat[!is.na(dat$`J-GENE and allele`),]
  dat$`CDR3-IMGT length` <- as.numeric(dat$`CDR3-IMGT length`)
  return(dat)
})
setwd("..")
hiv_mabs1_light <- bind_rows(hiv_mabs1_light)

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

hiv_mabs1_light <- clean.ambgs.calls(hiv_mabs1_light)
hiv_mabs1_light = hiv_mabs1_light[c("Sequence ID", "AA JUNCTION", "v_gene", "j_gene")]
hiv_mabs1_light$`AA JUNCTION` <- gsub(" .*", "",hiv_mabs1_light$`AA JUNCTION`)
hiv_mabs1_light$cdrl3_aa = stringr::str_sub(hiv_mabs1_light$`AA JUNCTION`, 2, -2)
hiv_mabs1_light$cdrl3_aa_length = nchar(hiv_mabs1_light$cdrl3_aa)

hiv_mabs1_light$clone_id = paste(hiv_mabs1_light$v_gene, hiv_mabs1_light$j_gene,
                                 hiv_mabs1_light$cdrl3_aa_length, sep = ',')

hiv_mabs1_light$seq_id = gsub('_.*', '', hiv_mabs1_light$`Sequence ID`)
hiv_mabs1_light$seq_id = gsub(' .*', '', hiv_mabs1_light$seq_id)
hiv_mabs1_light = hiv_mabs1_light[c("v_gene", "j_gene", "cdrl3_aa_length", 
                                    "cdrl3_aa", "clone_id", "seq_id")]
#####
hiv_mabs2_light = read_xlsx('ab_metadata/ab_germlines-updated_PB.xlsx', sheet = 1)
hiv_mabs2_light = hiv_mabs2_light[!is.na(hiv_mabs2_light$`Light CDR3 sequence`),]
hiv_mabs2_light = hiv_mabs2_light[!is.na(hiv_mabs2_light$`Light V`),]
hiv_mabs2_light = hiv_mabs2_light[!is.na(hiv_mabs2_light$`Light J`),]
hiv_mabs2_light$v_gene = paste('IG', hiv_mabs2_light$`Light V`, sep = '')
hiv_mabs2_light$j_gene = paste('IG', hiv_mabs2_light$`Light J`, sep = '')
hiv_mabs2_light$cdrl3_aa = hiv_mabs2_light$`Light CDR3 sequence`
hiv_mabs2_light$cdrl3_aa_length = nchar(hiv_mabs2_light$cdrl3_aa)
hiv_mabs2_light$check = hiv_mabs2_light$cdrl3_aa_length - as.numeric(hiv_mabs2_light$`Light CDR3 length`)
hiv_mabs2_light %>%
  mutate(across(c("v_gene", "j_gene"), 
                geneusage_transformation)) -> hiv_mabs2_light
hiv_mabs2_light$clone_id = paste(hiv_mabs2_light$v_gene, hiv_mabs2_light$j_gene,
                                 hiv_mabs2_light$cdrl3_aa_length, sep = ',')
hiv_mabs2_light$seq_id = hiv_mabs2_light$Antibody
hiv_mabs2_light = hiv_mabs2_light[c("v_gene", "j_gene", "cdrl3_aa_length", 
                                    "cdrl3_aa", "clone_id", "seq_id")]
#####

hiv_mabs_light = rbind(hiv_mabs1_light, hiv_mabs2_light)

hiv_mabs_light_heavy_c <- merge(hiv_mabs, hiv_mabs_light, by = "seq_id")

write_xlsx(hiv_mabs_light_heavy_c, "ab_metadata/hiv_mabs_combined_heavy_light.xlsx")
