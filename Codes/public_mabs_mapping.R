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

hiv_mabs = hiv_mabs[!duplicated(hiv_mabs$cdrh3_aa),]
hiv_mabs_list = split(hiv_mabs, f = hiv_mabs$clone_id)
############

all_data = readChangeoDb('imgt_combined_HC/combined_clone-pass_germ-pass.tsv')
all_data = all_data[!duplicated(all_data$sequence_id),]

#all_data = all_data[sample(seq_len(nrow(all_data)), 10000),]
#all_data = all_data[!duplicated(all_data$sequence_id),]
all_data$group_id = gsub('_sample.*', '', all_data$group_id)
all_data$cdrh3_aa = translateDNA(all_data$junction, trim = T)
all_data$cdrh3_aa_length = nchar(all_data$cdrh3_aa)

all_data %>%
  mutate(across(c("germline_v_call", "germline_d_call", "germline_j_call"), 
                geneusage_transformation)) -> all_data

all_data$clone_id = paste(all_data$germline_v_call, all_data$germline_j_call,
                          all_data$cdrh3_aa_length, sep = ',')

all_data_list = split(all_data, f = all_data$clone_id)

common_clones = intersect(names(all_data_list), hiv_mabs$clone_id)
test = all_data[[1]]

ref = hiv_mabs1_list[[c("IGHV1-18,IGHJ4,16")]]
ref = ref[!duplicated(ref$cdrh3_aa),]

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
clusterExport(cl, c('all_data_list', 'custom_fun', 'hiv_mabs_list'))
clusterEvalQ(cl, library(Biostrings))

# check = parLapply(cl,common_clones, function(x){
#   seq_data = all_data_list[[x]]
#   mabs = hiv_mabs_list[[x]]
#   mabs = mabs[!duplicated(mabs$cdrh3_aa),]
#   test = seq_data$cdrh3_aa
#   names(test) = seq_data$sequence_id
#   out = lapply(mabs$cdrh3_aa, function(mab){
#     res = custom_fun(test, mab)
#   })
#   names(out) = mabs$cdrh3_aa
#   return(out)
# })
# names(check) = common_clones
# 
# dat = lapply(check, function(list1){
#   #print(i)
#   #list1 = check[[i]]
#   out = lapply(list1, function(x){
#     x = bind_rows(x, .id = 'cdrh3_aa')
#   })
# })
# dat = lapply(dat, function(x){
#   out = bind_rows(x, .id = 'cdrh3_aa')
# })
# dat = bind_rows(dat, .id = 'clone_id')

###########################################################################
hiv_mabs = read_xlsx('ab_metadata/hiv_mabs_combined.xlsx', sheet = 1)
hiv_mabs_list = split(hiv_mabs, f = hiv_mabs$seq_id)

all_data = readRDS('rds/combined_data_with_clones_22_feb.rds')

mapped_data = lapply(seq_along(hiv_mabs_list), function(i){
  print(i)
  x = hiv_mabs_list[[i]]
  outdat = all_data[all_data$clone_id %in% x$clone_id,]
  outdat = outdat[c("sequence_id", "cdrh3_aa", "cdrh3_aa_length", "clone_id")]
  outdat$mapped_mab = x$seq_id
  outdat$mapped_mab_aa = x$cdrh3_aa
  return(outdat)
})

mapped_data_df = bind_rows(mapped_data)
mapped_data_df = mapped_data_df[!duplicated(mapped_data_df$cdrh3_aa),]

########
#add annotation and indel to the data
setwd('imgt_output_HC/')
summary_files = list.files(pattern = 'Summary', recursive = T)

summary_dat = lapply(summary_files, function(x){
  dat = read.csv(x, sep = '\t')
  dat = dat[dat$V.DOMAIN.Functionality %in% 
              c('productive', 'productive (see comment)'),]
  dat = dat[c("Sequence.ID", "V.REGION.identity..", "V.REGION.insertions",
              "V.REGION.deletions")]
  dat$mu_freq = 100 - dat$V.REGION.identity..
  colnames(dat) = c('sequence_id', 'vgene_identity', 'vgene_ins', 'vgene_del',
                    'mu_freq')
  group_id = gsub('/.*', '', x)
  dat$group_id = group_id
  return(dat)
})
setwd('..')
summary_dat_df = bind_rows(summary_dat)
summary_dat_df = summary_dat_df[!duplicated(summary_dat_df$sequence_id),]

trim_seq_id = function(x){
  x = stringi::stri_reverse(x)
  x = sub("^[^:]*:", "", x)
  x = stringi::stri_reverse(x)
  return(x)
}

summary_dat_df$sequence_id_new = sapply(summary_dat_df$sequence_id, trim_seq_id)
write.csv(summary_dat_df, 'all_sequences_indel_info.csv')

all_data_fil = all_data[c('sequence_id', 'germline_v_call', 'germline_j_call',
                          'germline_d_call', 'group_id', 'sequence')]
all_data_fil$sequence_id_new = sapply(all_data_fil$sequence_id, trim_seq_id)

mapped_data_df$sequence_id_new = sapply(mapped_data_df$sequence_id, trim_seq_id)
mapped_data_df = merge(mapped_data_df, summary_dat_df, by = 'sequence_id_new')
mapped_data_df = merge(mapped_data_df, all_data_fil, by = 'sequence_id_new')
mapped_data_df = mapped_data_df[c('sequence_id_new', 'sequence','germline_v_call', 
                                  'germline_d_call', 'germline_j_call',
                                  'cdrh3_aa', 'cdrh3_aa_length', 'mu_freq',
                                  'clone_id', 'mapped_mab', 'mapped_mab_aa',
                                  'vgene_ins', 'vgene_del', 'group_id.y')]

write_xlsx(mapped_data_df, 'ab_metadata/mapped_mabs_without_identity_cutoff.xlsx')

###############
#give clonotype data for 330_2018
dat_330_2018 = readChangeoDb('imgt_combined_HC/combined_clone-pass_germ-pass.tsv')
dat_330_2018$group_id = gsub('_sample.*', '', dat_330_2018$group_id)
dat_330_2018 = dat_330_2018[dat_330_2018$group_id == '330_2018',]
dat_330_2018 = dat_330_2018[!duplicated(dat_330_2018$sequence),]

clones_dat = countClones(dat_330_2018)

dat_330_2018_clones = dat_330_2018[!duplicated(dat_330_2018$clone_id),]
dat_330_2018_clones = merge(dat_330_2018_clones, clones_dat, by = 'clone_id')
dat_330_2018_clones$sequence_id_new = sapply(dat_330_2018_clones$sequence_id,
                                             trim_seq_id)
dat_330_2018_clones = merge(dat_330_2018_clones, summary_dat_df, 
                            by = 'sequence_id_new')
dat_330_2018_clones$cdrh3_aa = translateDNA(dat_330_2018_clones$junction, trim = T)
dat_330_2018_clones$cdrh3_aa_length = nchar(dat_330_2018_clones$cdrh3_aa)

dat_330_2018_clones = dat_330_2018_clones[c('sequence_id_new', 'germline_v_call', 
                                  'germline_d_call', 'germline_j_call',
                                  'cdrh3_aa', 'cdrh3_aa_length', 'mu_freq',
                                  'clone_id', 'vgene_ins', 'vgene_del', 'group_id.x',
                                  'seq_count')]

write_xlsx(dat_330_2018_clones, '330_2018_clonotypes_sample.xlsx')


test = dat_330_2018
test$sequence_id = sapply(test$sequence_id, trim_seq_id)
test$cdrh3_aa = translateDNA(test$junction, trim = T)
test1 = test[test$cdrh3_aa == 'ARDVLAFAFQDYYHGMDV',]
test1 = test[test$clone_id == '88705',]

########################
##load mapping data from server
mapped_dat = readRDS('rds/hiv_mabs_mapping_0.5.rds')
names(mapped_dat) = hiv_mabs$seq_id

mabs = lapply(seq_along(mapped_dat), function(i){
  print(i)
  x = mapped_dat[[i]]
  x = unlist(x)
  x = as.data.frame(x)
  x$id = rownames(x)
  rownames(x) = NULL
  x$subject = names(mapped_dat)[i]
  colnames(x)[1] = c('matches')
  x = x[!is.na(x$matches),]
  return(x)
})

mabs = bind_rows(mabs)

mabs = readRDS('rds/hiv_mabs_mapping_0.5_df.rds')
colnames(mabs)[2] = 'sequence_id'

all_data_fil = all_data[c('sequence_id', 'germline_v_call', 'germline_j_call',
                          'germline_d_call', 'group_id', 'sequence', 
                          'cdrh3_aa', 'cdrh3_aa_length')]

all_data_mapped = merge(all_data_fil, mabs, by = 'sequence_id')
colnames(all_data_mapped)[colnames(all_data_mapped) == 'subject'] = 'seq_id'
all_data_mapped = merge(all_data_mapped, hiv_mabs, by = 'seq_id')

colnames(all_data_mapped) = c('seq_id', 'sequence_id', 'germline_v_call', 
                              'germline_j_call', 'germline_d_call',
                              'group_id', 'sequence', 'cdrh3_aa', 'cdrh3_aa_length',
                              'num_mismatch', 'mab_vcall', 'mab_jcall', 
                              'mab_cdrh3_aa_length', 'mab_cdrh3_aa', 'clone_id')

all_data_mapped$mab_vj = paste(all_data_mapped$mab_vcall, all_data_mapped$mab_jcall,
                               sep = ',')

all_data_mapped$seq_vj = paste(all_data_mapped$germline_v_call, 
                               all_data_mapped$germline_j_call, sep = ',')

all_data_mapped$seq_cloneid = paste(all_data_mapped$germline_v_call, 
                                    all_data_mapped$germline_j_call,
                                    all_data_mapped$cdrh3_aa_length,
                                    sep = ',')

all_data_mapped$filter1 = ifelse(all_data_mapped$seq_vj == all_data_mapped$mab_vj,
                                 'Yes', 'No')
all_data_mapped$filter2 = ifelse(all_data_mapped$clone_id == all_data_mapped$seq_cloneid,
                                 'Yes', 'No')
  
all_data_mapped_vj_filter = all_data_mapped[all_data_mapped$filter1 == 'Yes',]  
  
all_data_mapped_clone_filter = all_data_mapped[all_data_mapped$filter2 == 'Yes',]  
all_data_mapped_clone_filter$sequence_id_new = 
  sapply(all_data_mapped_clone_filter$sequence_id, trim_seq_id)

#####
#map indel data
indel_dat = summary_dat_df
indel_dat = merge(all_data_mapped_clone_filter, indel_dat, by = 'sequence_id_new')

indel_dat = indel_dat[c('sequence_id_new', 'sequence_id.x', 'sequence', 
                        'germline_v_call', 'germline_d_call', 'germline_j_call',
                        'group_id.x', 'cdrh3_aa', 'cdrh3_aa_length', 
                        'vgene_ins', 'vgene_del', 'mu_freq', 'seq_id', 
                        'mab_vcall', 'mab_jcall', 'mab_cdrh3_aa', 
                        'mab_cdrh3_aa_length')]

colnames(indel_dat) = c('sequence_id_new', 'sequence_id', 'sequence', 
                        'v_call', 'd_call', 'j_call',
                        'sample_id', 'cdrh3_aa', 'cdrh3_aa_length', 
                        'vgene_ins', 'vgene_del', 'mu_freq', 'mab_id', 
                        'mab_vcall', 'mab_jcall', 'mab_cdrh3_aa', 
                        'mab_cdrh3_aa_length')

mabs_summary = indel_dat[c('mab_id', 'sample_id')]
mabs_summary = split(mabs_summary, f = mabs_summary$mab_id)
mabs_out = lapply(seq_along(mabs_summary), function(i){
  x = mabs_summary[[i]]
  out = data.frame(table(x$sample_id))
  out$comb = paste(out$Var1, ' (', out$Freq, ')', sep = '')
  #out$comb = paste(out$comb, collapse = ', ')
  out2 = data.frame('mab_id' = names(mabs_summary)[i], 
                    'summary' = paste(out$comb, collapse = ', '))
})

mabs_out = bind_rows(mabs_out)

write_xlsx(list(mapped_mabs = indel_dat, mabs_summary = mabs_out), 
           'mabs_mapping_50_perc_cutoff_with_indel.xlsx')



