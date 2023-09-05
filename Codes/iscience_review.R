all_dat_hc = readRDS('rds/all_dat_cdr3_shm_HC.rds')
all_dat_hc$mu_freq = round(all_dat_hc$mu_freq*100, 2)

all_dat_hc$group_id[all_dat_hc$group_id == '330_2018_sample3'] = '330_2018'

clean.ambgs.calls = function(df){
  df$v_gene = gsub('Homsap ', '', df$v_call)
  df$v_gene = gsub('\\*.*', '', df$v_gene)
  df$v_family = gsub('-.*', '', df$v_gene)
  df_clones_split = split(df, df$clone_id)
  df_clones_gene_corrected = lapply(df_clones_split, function(x){
    v_genes = x$v_gene[order(table(x$v_gene), decreasing = T)]
    if(length(v_genes)>1){print(paste('Found Ambiguity, Selecting', v_genes[1], 
                                      'from', paste(v_genes, collapse = ',')))}
    x$v_gene = v_genes[1]
    v_families = x$v_family[order(table(x$v_family), decreasing = T)]
    if(length(v_families)>1){print(paste('Found Ambiguity, Selecting', v_families[1], 
                                         'from', paste(v_families, collapse = ',')))}
    x$v_family = v_families[1]
    return(x)
  })
  df_clones_gene_corrected = bind_rows(df_clones_gene_corrected)
  return(df_clones_gene_corrected)
}

all_dat_hc_cleaned <- clean.ambgs.calls(all_dat_hc)
all_dat_hc_cleaned$j_gene <- gsub('Homsap ', '', all_dat_hc_cleaned$j_call)
all_dat_hc_cleaned$j_gene = gsub('\\*.*', '', all_dat_hc_cleaned$j_gene)

test <- all_dat_hc_cleaned[all_dat_hc_cleaned$v_gene == "IGHV4-59",]
test <- test[test$j_gene == "IGHJ6",]
test <- test[test$cdr3_aa_length == 18,]

abstar_files = list.files('abstar_output_HC/', pattern = '.txt')
setwd('abstar_output_HC/')
filnames = gsub('.txt', '', abstar_files)
abstar_dat_hc = lapply(seq_along(abstar_files), function(i){
  dat = read.csv(abstar_files[i])
  dat$group_id = filnames[i]
  return(dat)
})
setwd('..')
abstar_dat_hc = bind_rows(abstar_dat_hc)
colnames(abstar_dat_hc)[1] = 'sequence_id'

#table of sequences
dat <- as.data.frame(table(all_dat_hc$group_id))
write_xlsx(dat, "~/Downloads/iscience reviews/Table S1.xlsx")




