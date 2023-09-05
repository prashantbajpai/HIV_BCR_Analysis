library(alakazam)
library(shazam)
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(writexl)
library(ggrepel)
library(patchwork)
#source('/Users/prashant/R_themes/theme_publication.R')

all_dat_hc = readRDS('rds/all_dat_cdr3_shm_HC.rds')
all_dat_hc$mu_freq = round(all_dat_hc$mu_freq*100, 2)

all_dat_hc$group_id[all_dat_hc$group_id == '330_2018_sample3'] = '330_2018'
#all_dat_hc$sequence_id = gsub('.{11}$', '', all_dat_hc$sequence_id)
all_dat_hc$v_region_insertion = 'NA'
all_dat_hc$v_region_deletion = 'NA'
#add indels info
setwd('imgt_output_HC/')
indel_summary = list.files(pattern = 'Summary', recursive = T)
indel_dat = lapply(indel_summary, function(x){
  dat = read.delim(x)
  dat = dat[c("Sequence.ID", "V.REGION.insertions", "V.REGION.deletions")]
  return(dat)
})
indel_dat = bind_rows(indel_dat)
colnames(indel_dat)[1] = 'sequence_id'
library(foreach)
library(doParallel)
indel_dat_added = foreach(i = seq_along(indel_dat$sequence_id)) %dopar% {
  print(i)
  all_dat_hc$v_region_insertion[grep(indel_dat$sequence_id[i], 
                                     all_dat_hc$sequence_id)] = 
    indel_dat$V.REGION.insertions[i]
  all_dat_hc$v_region_deletion[grep(indel_dat$sequence_id[i], 
                                    all_dat_hc$sequence_id)] = 
    indel_dat$V.REGION.deletions[i]
  return(NULL)
}
setwd('..')
saveRDS(all_dat_hc, 'rds/all_dat_cdr3_shm_HC_indels.rds')
saveRDS(indel_dat_added, 'rds/indel_dat_added.rds')

##################################
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
abstar_dat_hc = abstar_dat_hc[c('sequence_id', 'isotype')]

merged_dat = merge(abstar_dat_hc, all_dat_hc, by = 'sequence_id')
merged_dat$clone_id = as.character(merged_dat$clone_id)

cdr3_plotdat = merged_dat[!is.na(merged_dat$cdr3_aa_length),]
cdr3_means = cdr3_plotdat %>%
  group_by(group_id) %>%
  summarise(mean = mean(cdr3_aa_length), .groups = 'drop')

cdr3_means$mean = round(cdr3_means$mean, 2)
write_xlsx(cdr3_means, "suppliment/cdr3_mean_heavychain.xlsx")

p1 = ggplot(cdr3_plotdat, aes(x=cdr3_aa_length, fill = group_id)) + 
  geom_density(bw = 0.8, alpha = 0.8) +
  facet_wrap(group_id~., ncol = 1, strip.position = 'left') +
  geom_vline(data = cdr3_means, aes(xintercept = mean), 
             linetype=2) +
  geom_text_repel(data = cdr3_means, aes(x = mean, y = 0.12, label = mean),
                  size = 8) +
  xlab('CDR3 Length') + ylab('Density') + ylab('') +
  coord_cartesian(xlim = c(1, 28)) +
  theme_Publication() +
  scale_fill_manual(name = '', values = c("#9d9dff", "#c7e9c0", "#ff9d9d", "#1414ff", "#238b45", "#ff1414")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 14),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(face = 'bold', size = 16),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 16))

p_final = p1 + plot_annotation(title = 'CDR3 Length Distribution') &
  theme(plot.title = element_text(face = 'bold', size = 20, hjust = 0.5), 
        legend.position = 'none')

ggsave('plot_ver1/cdr3_length_HC.png', p_final, width = 9.5, height = 11, units = 'in')
ggsave('plot_ver1/cdr3_length_HC.pdf', p_final, width = 9.5, height = 11, units = 'in')

#####
#SHM

shm <- merged_dat %>% group_by(group_id) %>% summarise(median_shm = median(mu_freq))
write_xlsx(shm, "suppliment/shm_all_groups.xlsx")

shm = merged_dat

p2 = ggplot(shm, aes(x = group_id, y = mu_freq, fill = group_id)) +
  geom_violin(bw = 0.8) + 
  stat_summary(fun = median, geom = 'point', size = 4, show.legend = F) +
  theme_Publication() + 
  scale_fill_manual(name = '', values = c("#9d9dff", "#c7e9c0", "#ff9d9d", "#1414ff", "#238b45", "#ff1414")) +
  coord_cartesian(ylim = c(0,18)) +
  ylab('SHM (%)') + xlab('') +
  theme(axis.text.x = element_text(size = 16, angle = -60, hjust = 0, face = 'bold', vjust = 0.5),
        legend.key.size = unit(1,'cm'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        legend.position = "none")

ggsave('plot_ver1/shm_hc.png', p2, width = 10, height = 6, units = 'in')
ggsave('plot_ver1/shm_hc.pdf', p2, width = 10, height = 6, units = 'in')

#merge IgAs, IgGs, remove unknown

shm$isotype = gsub('[0-9]', '', shm$isotype)
shm = shm[!shm$isotype == 'unknown',]

p3 = ggplot(shm, aes(x = group_id, y = mu_freq, fill = group_id)) +
  geom_violin(bw = 0.8) + 
  stat_summary(fun = median, geom = 'point', size = 3, show.legend = F) +
  facet_wrap(isotype~., scales = 'free') +
  theme_Publication() + 
  scale_fill_manual(name = '', values = c("#9d9dff", "#c7e9c0", "#ff9d9d", "#1414ff", "#238b45", "#ff1414")) +
  coord_cartesian(ylim = c(0,30)) +
  ylab('SHM (%)') + xlab('') +
  theme(axis.text.x = element_text(size = 16, angle = -60, hjust = 0, face = 'bold', vjust = 0.5),
        legend.key.size = unit(1,'cm'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        legend.position = 'none',
        strip.text = element_text(size = 18))

ggsave('plot_ver1/shm_by_isotype_hc.png', p3, width = 14, height = 10, units = 'in')
ggsave('plot_ver1/shm_by_isotype_hc.pdf', p3, width = 14, height = 10, units = 'in')

shm <- shm %>% group_by(group_id, isotype) %>% summarise(median_shm = median(mu_freq), .groups = "drop")
shm <- split(shm, f = shm$isotype)
write_xlsx(shm, "suppliment/shm_all_groups_byisotype.xlsx")
##gene usage
vcall = countGenes(merged_dat, gene = 'v_call', mode = 'gene', groups = 'group_id')
vcall$seq_freq = round(vcall$seq_freq*100, 2)
vcall$call = 'v_gene'
dcall = countGenes(merged_dat, gene = 'd_call', mode = 'gene', groups = 'group_id')
dcall$seq_freq = round(dcall$seq_freq*100, 2)
dcall$call = 'd_gene'
jcall = countGenes(merged_dat, gene = 'j_call', mode = 'gene', groups = 'group_id')
jcall$seq_freq = round(jcall$seq_freq*100, 2)
jcall$call = 'j_gene'
genecall_dat = rbind(vcall, dcall, jcall)

#save data
write_xlsx(list("vcall" = vcall, "dcall" = dcall, "jcall" = jcall), "suppliment/genecall_heavychain.xlsx")

geneusage_plots = function(df, title,legendpos){
  df = df[order(df$seq_freq, decreasing = T), ]
  df$gene = factor(df$gene, levels = unique(df$gene))
  #remove genes with low freq (sum of P+S < 1%)
  df %>%
    group_by(gene) %>%
    summarise(freqsum = sum(seq_freq), .groups = 'drop') -> freqsum
  freqsum = freqsum[freqsum$freqsum>5, ]
  df = df[df$gene %in% freqsum$gene,]
  p1 = ggplot(df, aes(x = gene, y = seq_freq, fill = group_id)) +
    geom_bar(stat = 'identity', position = position_dodge(preserve = 'single'), 
             color = 'black') +
    labs(title = title) + ylab('Gene Usage (%)') + xlab('') +
    theme_Publication() + 
    scale_fill_manual(name = '', values = c("#9d9dff", "#c7e9c0", "#ff9d9d", "#1414ff", "#238b45", "#ff1414")) +
    theme(axis.text.x = element_text(angle = 90, size = 20, vjust = 0.5, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 20),
          legend.key.size = unit(1,'cm'),
          legend.text = element_text(size = 20),
          legend.position = legendpos,
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 22)) +
    guides(fill=guide_legend(nrow=1,byrow=TRUE))
  return(p1)
}

genecall_plasma_vgene = genecall_dat[genecall_dat$call == 'v_gene', ]
genecall_plasma_dgene = genecall_dat[genecall_dat$call == 'd_gene', ]
genecall_plasma_jgene = genecall_dat[genecall_dat$call == 'j_gene', ]

p1 = geneusage_plots(genecall_plasma_vgene, 'V Gene Usage', 'none')
p2 = geneusage_plots(genecall_plasma_dgene, 'D Gene Usage', 'none')
p3 = geneusage_plots(genecall_plasma_jgene, 'J Gene Usage', 'none')

p_final = p1 + {p2 + p3 + plot_layout(ncol = 2, widths = c(1, 0.4))} + 
  plot_layout(nrow = 2, guides = 'collect') + 
  plot_annotation(title = 'Gene Usage Heavy Chain') &
  theme(plot.title = element_text(face = 'bold', size = 20, hjust = 0.5), 
        legend.position = 'none') 

ggsave('plot_ver1/gene_usage_HC.png', p_final, width = 20, height = 16, units = 'in')
ggsave('plot_ver1/gene_usage_HC.pdf', p_final, width = 20, height = 16, units = 'in')

####
#Upset plot
clones <- merged_dat
clones <- split(clones, f = clones$group_id)
leastcount <- unname(sort(sapply(clones, function(x)nrow(x)))[1])
clones_sampled <- lapply(clones, function(x){
  set.seed(10)
  x <- x[sample(nrow(x), leastcount),]
})
clones <- bind_rows(clones_sampled)
#remove singlets
clones = countClones(merged_dat, groups = c('group_id'))
#filter singlet clones
clones = clones[clones$seq_count>1,]
clones = split(clones, f = clones$group_id)

merged_dat_singlet_rem = split(merged_dat, f = merged_dat$group_id)
merged_dat_singlet_rem = lapply(seq_along(merged_dat_singlet_rem), function(i){
  clones_to_keep = clones[[i]]$clone_id
  dat = merged_dat_singlet_rem[[i]]
  dat = dat[dat$clone_id %in% clones_to_keep,]
  return(dat)
})
merged_dat_singlet_rem = bind_rows(merged_dat_singlet_rem)

# #sample similar number of sequences from plasmablasts
# merged_dat_singlet_rem_split = split(merged_dat_singlet_rem, 
#                                      f = merged_dat_singlet_rem$group_id)
# nsample = sapply(merged_dat_singlet_rem_split, function(x)dim(x)[1])
# nsample = min(nsample)
# merged_dat_singlet_rem_split = lapply(merged_dat_singlet_rem_split, function(dat){
#   set.seed(100)
#   dat = dat[sample(seq_len(nrow(dat)), nsample),]
#   return(dat)
# })

#all_dat_sampled = bind_rows(merged_dat_singlet_rem_split)

clones = countClones(merged_dat_singlet_rem, groups = c('group_id'))
clones %>%
  group_by(clone_id) %>%
  summarise(seq_count = mean(seq_count), seq_freq = mean(seq_freq)) ->
  clones

all_dat_sampled_fil = merged_dat_singlet_rem[c('group_id', 'v_call', 'd_call', 
                                               'j_call','cdr3', 'cdr3_aa_length', 
                                               'mu_freq', 'clone_id')]

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

all_dat_sampled_fil = clean.ambgs.calls(all_dat_sampled_fil)

all_dat_sampled_fil %>%
  group_by(group_id, clone_id, v_family) %>%
  summarise(cdr3_aa_length = mean(cdr3_aa_length), mu_freq = mean(mu_freq),
            .groups = 'drop') -> dat_sum

library(tidyr)
library(ComplexHeatmap)
library(ggbeeswarm)
library(ComplexUpset)

listset = split(dat_sum, f = dat_sum$group_id)
listset = lapply(listset, function(x){x = as.character(x$clone_id)})
upsetmat = list_to_matrix(listset)
upsetmat = as.data.frame(upsetmat)
upsetmat$clone_id = rownames(upsetmat)

dat_sum %>%
  group_by(clone_id, cdr3_aa_length, v_family) %>%
  summarise(mu_freq = mean(mu_freq)) -> dat_sum

plotdat = Reduce(function(x, y) merge(x, y, by = 'clone_id'), 
                 list(clones, dat_sum, upsetmat))

plotdat$threhold = cut(plotdat$mu_freq, breaks = c(0, 5, 10, 15, Inf))
plotdat$threhold[which(is.na(plotdat$threhold))] = '(0,5]'
#make upset plots
combination = grep('329|330', colnames(plotdat), value = T)

#set_size(8,3)
set.seed(100)
p1 = upset(plotdat, combination, width_ratio = 0.2, 
           name = 'Sample ID',
           themes = upset_modify_themes(
             list(
               'intersections_matrix'=theme(text=element_text(size=18, face = 'bold'))
             )),
           sort_intersections=FALSE,
           intersections = list(
             "329_2015", "329_2016", "329_2018", 
             "330_2015", "330_2016", "330_2018",
             c("329_2015", "329_2016", "329_2018"),
             c("330_2015", "330_2016", "330_2018"),
             c("329_2015", "329_2016"), 
             c("329_2016", "329_2018"),
             c("330_2015", "330_2016"),
             c("330_2016", "330_2018"),
             c("329_2015", "330_2015"),
             c("329_2016", "330_2016"),
             c("329_2018", "330_2018")),
           base_annotations = list(
             'Intersection size' = (intersection_size(text = list(size = 9))
                                    + ylab('Number of \n Clonotypes')
                                    + theme(axis.title.y = element_text(face = 'bold',
                                                                        size = 18),
                                            axis.text.y = element_text(size = 16))
             )),
           set_sizes = upset_set_size(aes(width = 0.1)) + ylab('Number of Clonotypes') 
           + geom_text(aes(label=..count..), hjust=1.1, stat='count', 
                       fontface = 'bold', size = 7)
           + expand_limits(y=65000)
           + theme(axis.text.x = element_blank(), 
                   axis.title = element_text(size = 18, face = 'bold')),
           sort_sets = F,
           annotations = list(
             'V-Gene Distribution'=(
               ggplot(mapping = aes(fill = v_family))
               + geom_bar(stat='count', position=position_fill(reverse = T), na.rm=TRUE)
               + ylab('V-Family Distribution')
               + scale_y_continuous(labels=scales::percent_format())
               + scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb',
                                              '#e78ac3', '#a6d854', '#ffd92f',
                                              '#e5c494'), name = 'V-Family')
               + scale_color_manual(values=c('show'='black', 'hide'='transparent'), 
                                    guide=FALSE)
               + theme(axis.title.y = element_text(face = 'bold', size = 18),
                       legend.text = element_text(size = 16),
                       legend.key.size = unit(1, 'cm'),
                       legend.title = element_text(size = 16),
                       axis.text.y = element_text(size = 16))
             ),
             'CDR3 Length'=(
               ggplot(mapping = aes(y=cdr3_aa_length, size = seq_count, 
                                    color = threhold))
               + scale_color_manual(values = c('#e6f5d0','#4d9221', 
                                               '#de77ae', '#8e0152'),
                                    name = 'Mean V-Region Mutation',
                                    labels = c('<5', '5-10', '10-15', '>15'))
               + geom_quasirandom(na.rm=TRUE, alpha = 0.8)
               + scale_size_continuous(name = 'Read Count', breaks = 
                                         c(1, 10, 100, 500, 1000, 5000), 
                                       range = c(2, 11)) 
               + guides(colour = guide_legend(override.aes = list(size=10)))
               + theme(legend.text = element_text(size = 16),
                       axis.title.y = element_text(face = 'bold', size = 18),
                       axis.text.y = element_text(size = 16),
                       legend.title = element_text(size = 16))
             )
           )
)


ggsave('plot_ver1/upset_plot_hc.pdf', p1, width = 30, height = 20, units = 'in')
ggsave('plot_ver1/upset_plot_hc.png', p1, width = 30, height = 20, units = 'in')

####
seqlogodat = plotdat
seqlogodat %>%
  rowwise() %>%
  mutate('rowsum' = sum(c_across(7:12)),
         'rowsum_329' = sum(c_across(7:9)),
         'rowsum_330' = sum(c_across(10:12))) -> seqlogodat

clones_2016 = seqlogodat[seqlogodat$`329_2016` == 1 & 
                          seqlogodat$`330_2016` == 1 &
                          seqlogodat$rowsum == 2 & 
                          seqlogodat$rowsum_329 == 1 &
                          seqlogodat$rowsum_330 == 1,]

clone.seqs = function(clones_df, complete_dat = merged_dat_singlet_rem){
  common_clones = unique(clones_df$clone_id)
  common_clones_dat = complete_dat[complete_dat$clone_id %in% common_clones,]
  common_clones_dat = common_clones_dat[order(common_clones_dat$v_call, 
                                              decreasing = F),]
  listorder = unique(common_clones_dat$clone_id)
  common_clones_dat$cdr3_aa = translateDNA(common_clones_dat$cdr3, trim = T)
  common_clones_dat = split(common_clones_dat, f = common_clones_dat$clone_id)
  common_clones_dat = common_clones_dat[listorder]
  return(common_clones_dat)
}

clones_2016 = clone.seqs(clones_2016)

library(ggseqlogo)
library(ggpubr)

custom.seqlogo = function(myclones){
  all_plots = lapply(seq_along(myclones), function(i){
    x = myclones[[i]]
    x$v_gene = gsub('Homsap ', '', x$v_call)
    x$v_gene = gsub('\\*.*', '', x$v_gene)
    cdr3_clones = split(x, f = x$group_id)
    mytitle = paste(unique(x$v_gene), ' (', unique(x$clone_id), ')', sep = '')
    print(names(cdr3_clones))
    p1 = ggseqlogo(cdr3_clones[[names(cdr3_clones)[1]]]$cdr3_aa, seq_type = 'aa', method = 'prob') + 
      labs(title = paste(names(cdr3_clones)[1], ' (SHM: ', 
                         round(mean(cdr3_clones[[names(cdr3_clones)[1]]]$mu_freq), 2), ')', sep = '')) + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.text = element_text(size = 16),
            axis.title.y = element_text(size = 18))
    p2 = ggseqlogo(cdr3_clones[[names(cdr3_clones)[2]]]$cdr3_aa, seq_type = 'aa', method = 'prob') +
      labs(title = paste(names(cdr3_clones)[2], ' (SHM: ', 
                         round(mean(cdr3_clones[[names(cdr3_clones)[2]]]$mu_freq), 2), ')', sep = '')) + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.text = element_text(size = 16),
            axis.title.y = element_text(size = 18))
    # p3 = ggseqlogo(cdr3_clones[['Memory']]$cdr3_aa, seq_type = 'aa', method = 'prob') +
    #   labs(title = paste('Memory B-Cell', ' (SHM: ', 
    #                      round(mean(cdr3_clones[['Memory']]$mu_freq), 2), ')', sep = '')) + 
    #   theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
    #         axis.text = element_text(size = 16),
    #         axis.title.y = element_text(size = 18))
    p_final = p1 + p2 + plot_layout(nrow = 2, guides = 'collect') + 
      plot_annotation(title = mytitle) & 
      theme(legend.position = 'bottom', 
            plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))
    return(p_final)
  })
  
}

all_plots = custom.seqlogo(clones_2016)
p1 = ggarrange(plotlist = all_plots[1:4], nrow = 2, ncol = 2)
p2 = ggarrange(plotlist = all_plots[5:8], nrow = 2, ncol = 2)
p3 = ggarrange(plotlist = all_plots[9:12], nrow = 2, ncol = 2)

ggsave('plot_ver1/seqlogo/shared_clones_2016_1.png', p1, width = 30, height = 25, 
       units = 'in')
ggsave('plot_ver1/seqlogo/shared_clones_2016_2.png', p2, width = 30, height = 25, 
       units = 'in')
ggsave('plot_ver1/seqlogo/shared_clones_2016_3.png', p3, width = 30, height = 25, 
       units = 'in')

bnabs = read_xlsx('HIV_bnabs_sanders_virology.xlsx', sheet = 1)

aligndat = merged_dat
aligndat$cdr3_aa = translateDNA(aligndat$cdr3, trim = F)

library(Biostrings)

align_matches = function(bnabs, cdr3_data){
  library(Biostrings)
  mismatch_data = data.frame(matrix(data = NA, nrow = nrow(bnabs), 
                                    ncol = nrow(cdr3_data)))
  rownames(mismatch_data) = bnabs$Antibody
  colnames(mismatch_data) = cdr3_data$sequence_id
  for(i in 1:nrow(bnabs)){
    for(j in seq_along(cdr3_data$cdr3_aa)){
      mymax = round(nchar(bnabs$HCDR3a[i])*0.2)
      matched_data = matchPattern(bnabs$HCDR3a[i], cdr3_data$cdr3_aa[j], 
                                  max.mismatch = mymax, min.mismatch = 0)
      matched_num = nmismatch(bnabs$HCDR3a[i], matched_data)
      if(length(matched_num)==1){
        mismatch_data[i,j] = matched_num
      }
      if(length(matched_num)>1){
        mismatch_data[i,j] = sort(unique(matched_num))[1]
        #print('wtf')
      }
    }
  }
  return(mismatch_data)
}

mymatches = align_matches(bnabs = bnabs, cdr3_data = test)

saveRDS(mymatches, 'rds/bnabs_80perc_match_alignment.rds')

mymatches$rowsums = apply(mymatches, 1, function(x)sum(is.na(x)))
mymatches$rowsums = mymatches$rowsums - 934304

test = mymatches[mymatches$rowsums !=0, ]
test$rowsums = NULL
test = unlist(test)
test1 = na.omit(test)

test = merged_dat[merged_dat$sequence_id == 
                    'M03837:141:000000000-CKNLB:1:1119:23600:2142:TAGCTT;0.974669',]
coi = merged_dat[merged_dat$clone_id == '147082',]

##
mymatches_pct64 = align_matches(bnabs = bnabs[16,], cdr3_data = aligndat)

sum(!is.na(unlist(mymatches_pct64)))

###########
#Mapped mAbs exploration May 26 2023
mymatches <- readRDS("rds/bnabs_80perc_match_alignment.rds")


#######upset with sampling
####
#Upset plot
clones <- merged_dat
clones <- split(clones, f = clones$group_id)
leastcount <- unname(sort(sapply(clones, function(x)nrow(x)))[1])
clones_sampled <- lapply(clones, function(x){
  set.seed(10)
  x <- x[sample(nrow(x), leastcount),]
})
clones <- bind_rows(clones_sampled)
#remove singlets
clones_counts = countClones(clones, groups = c('group_id'))
#filter singlet clones
clones_counts = clones_counts[clones_counts$seq_count>1,]

#add indel info to the clones
clones_indel <- split(clones_counts, f = clones_counts$clone_id)
test <- clones_indel[[1]]
test <- clones[clones$clone_id == test$clone_id,]
test <- test[c("sequence_id", "clone_id")]

test$indel = 'No Indel'

indel_ins <- indel_dat[indel_dat$V.REGION.insertions != "",]
indel_ins <- indel_ins[grep("do not cause frameshift", indel_ins$V.REGION.insertions, ignore.case = T),]
indel_del <- indel_dat[indel_dat$V.REGION.deletions != "",]
indel_del <- indel_del[grep("do not cause frameshift", indel_del$V.REGION.insertions, ignore.case = T),]
indel_fil <- rbind(indel_ins, indel_del)
indel_fil <- indel_fil[!duplicated(indel_fil$sequence_id),]
indel_fil$indel <- "Has Indel"
indel_fil <- indel_fil[c("sequence_id", "indel")]

library(foreach)
library(doParallel)

clones_indel <- lapply(seq_along(clones_indel), function(i){
  print(i)
  test <- clones_indel[[i]]
  test <- clones[clones$clone_id == test$clone_id,]
  test <- test[c("sequence_id", "clone_id")]
  test$indel = 'No Indel'
  indel_dat_added = foreach(j = seq_along(indel_fil$sequence_id)) %dopar% {
    #print(i)
    test$indel[grep(indel_fil$sequence_id[j], 
                    test$sequence_id)] = 
      indel_fil$indel[j]
    return(NULL)
  }
  return(test)
})

clones <- clones_counts
clones = split(clones, f = clones$group_id)

merged_dat_singlet_rem = split(merged_dat, f = merged_dat$group_id)
merged_dat_singlet_rem = lapply(seq_along(merged_dat_singlet_rem), function(i){
  clones_to_keep = clones[[i]]$clone_id
  dat = merged_dat_singlet_rem[[i]]
  dat = dat[dat$clone_id %in% clones_to_keep,]
  return(dat)
})
merged_dat_singlet_rem = bind_rows(merged_dat_singlet_rem)

# #sample similar number of sequences from plasmablasts
# merged_dat_singlet_rem_split = split(merged_dat_singlet_rem, 
#                                      f = merged_dat_singlet_rem$group_id)
# nsample = sapply(merged_dat_singlet_rem_split, function(x)dim(x)[1])
# nsample = min(nsample)
# merged_dat_singlet_rem_split = lapply(merged_dat_singlet_rem_split, function(dat){
#   set.seed(100)
#   dat = dat[sample(seq_len(nrow(dat)), nsample),]
#   return(dat)
# })

#all_dat_sampled = bind_rows(merged_dat_singlet_rem_split)

clones = countClones(merged_dat_singlet_rem, groups = c('group_id'))
clones %>%
  group_by(clone_id) %>%
  summarise(seq_count = mean(seq_count), seq_freq = mean(seq_freq)) ->
  clones

all_dat_sampled_fil = merged_dat_singlet_rem[c('group_id', 'v_call', 'd_call', 
                                               'j_call','cdr3', 'cdr3_aa_length', 
                                               'mu_freq', 'clone_id')]

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

all_dat_sampled_fil = clean.ambgs.calls(all_dat_sampled_fil)

all_dat_sampled_fil %>%
  group_by(group_id, clone_id, v_family) %>%
  summarise(cdr3_aa_length = mean(cdr3_aa_length), mu_freq = mean(mu_freq),
            .groups = 'drop') -> dat_sum

library(tidyr)
library(ComplexHeatmap)
library(ggbeeswarm)
library(ComplexUpset)

listset = split(dat_sum, f = dat_sum$group_id)
listset = lapply(listset, function(x){x = as.character(x$clone_id)})
upsetmat = list_to_matrix(listset)
upsetmat = as.data.frame(upsetmat)
upsetmat$clone_id = rownames(upsetmat)

dat_sum %>%
  group_by(clone_id, cdr3_aa_length, v_family) %>%
  summarise(mu_freq = mean(mu_freq)) -> dat_sum

plotdat = Reduce(function(x, y) merge(x, y, by = 'clone_id'), 
                 list(clones, dat_sum, upsetmat))

plotdat$threhold = cut(plotdat$mu_freq, breaks = c(0, 5, 10, 15, Inf))
plotdat$threhold[which(is.na(plotdat$threhold))] = '(0,5]'
plotdat$label <- ""
for (i in seq_along(plotdat$label)) {
  if(plotdat$mu_freq[i] >= 10){
    plotdat$label[i] <- plotdat$clone_id[i]
  }
}
#make upset plots
combination = grep('329|330', colnames(plotdat), value = T)

#set_size(8,3)
set.seed(100)
p1 = upset(plotdat, combination, width_ratio = 0.2, 
           name = 'Sample ID',
           themes = upset_modify_themes(
             list(
               'intersections_matrix'=theme(text=element_text(size=18, face = 'bold'))
             )),
           sort_intersections=FALSE,
           intersections = list(
             "329_2015", "329_2016", "329_2018", 
             "330_2015", "330_2016", "330_2018",
             c("329_2015", "329_2016", "329_2018"),
             c("330_2015", "330_2016", "330_2018"),
             c("329_2015", "329_2016"), 
             c("329_2016", "329_2018"),
             c("330_2015", "330_2016"),
             c("330_2016", "330_2018"),
             c("329_2015", "330_2015"),
             c("329_2016", "330_2016"),
             c("329_2018", "330_2018")),
           base_annotations = list(
             'Intersection size' = (intersection_size(text = list(size = 9))
                                    + ylab('Number of \n Clonotypes')
                                    + theme(axis.title.y = element_text(face = 'bold',
                                                                        size = 18),
                                            axis.text.y = element_text(size = 16))
             )),
           set_sizes = upset_set_size(aes(width = 0.1)) + ylab('Number of Clonotypes') 
           + geom_text(aes(label=..count..), hjust=1.1, stat='count', 
                       fontface = 'bold', size = 7)
           + expand_limits(y=65000)
           + theme(axis.text.x = element_blank(), 
                   axis.title = element_text(size = 18, face = 'bold')),
           sort_sets = F,
           annotations = list(
             'V-Gene Distribution'=(
               ggplot(mapping = aes(fill = v_family))
               + geom_bar(stat='count', position=position_fill(reverse = T), na.rm=TRUE)
               + ylab('V-Family Distribution')
               + scale_y_continuous(labels=scales::percent_format())
               + scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb',
                                              '#e78ac3', '#a6d854', '#ffd92f',
                                              '#e5c494'), name = 'V-Family')
               + scale_color_manual(values=c('show'='black', 'hide'='transparent'), 
                                    guide=FALSE)
               + theme(axis.title.y = element_text(face = 'bold', size = 18),
                       legend.text = element_text(size = 16),
                       legend.key.size = unit(1, 'cm'),
                       legend.title = element_text(size = 16),
                       axis.text.y = element_text(size = 16))
             ),
             'CDR3 Length'=(
               ggplot(mapping = aes(y=cdr3_aa_length, size = seq_count, 
                                    color = threhold, label = label))
               + scale_color_manual(values = c('#e6f5d0','#4d9221', 
                                               '#de77ae', '#8e0152'),
                                    name = 'Mean V-Region Mutation',
                                    labels = c('<5', '5-10', '10-15', '>15'))
               + geom_quasirandom(na.rm=TRUE, alpha = 0.8)
               + geom_text_repel(box.padding = 0.5, max.overlaps = Inf, color = "black", size = 3) 
               + scale_size_continuous(name = 'Read Count', breaks = 
                                         c(1, 10, 100, 500, 1000, 5000), 
                                       range = c(2, 11)) 
               + guides(colour = guide_legend(override.aes = list(size=10)))
               + theme(legend.text = element_text(size = 16),
                       axis.title.y = element_text(face = 'bold', size = 18),
                       axis.text.y = element_text(size = 16),
                       legend.title = element_text(size = 16))
             )
           )
)


ggsave('plot_ver1/upset_plot_hc_sampled.pdf', p1, width = 30, height = 20, units = 'in')
ggsave('plot_ver1/upset_plot_hc_sampled_1.png', p1, width = 30, height = 20, units = 'in')
####
###seqlogo of clones of interest

####
clones_ofint = plotdat[plotdat$clone_id %in% c("221871", "116705", "299920", "252697", "10611", "254244", "373110", "39960", "386635"),]

clone.seqs = function(clones_df, complete_dat = merged_dat_singlet_rem){
  common_clones = unique(clones_df$clone_id)
  common_clones_dat = complete_dat[complete_dat$clone_id %in% common_clones,]
  common_clones_dat = common_clones_dat[order(common_clones_dat$v_call, 
                                              decreasing = F),]
  listorder = unique(common_clones_dat$clone_id)
  common_clones_dat$cdr3_aa = translateDNA(common_clones_dat$cdr3, trim = T)
  common_clones_dat = split(common_clones_dat, f = common_clones_dat$clone_id)
  common_clones_dat = common_clones_dat[listorder]
  return(common_clones_dat)
}

clones_ofint = clone.seqs(clones_ofint)

#check for indels in these clones
clones <- merged_dat
clones <- split(clones, f = clones$group_id)
leastcount <- unname(sort(sapply(clones, function(x)nrow(x)))[1])
clones_sampled <- lapply(clones, function(x){
  set.seed(10)
  x <- x[sample(nrow(x), leastcount),]
})
clones <- bind_rows(clones_sampled)

trim_seq_id = function(x){
  x = stringi::stri_reverse(x)
  x = sub("^[^:]*:", "", x)
  x = stringi::stri_reverse(x)
  return(x)
}

clones$sequence_id_new <- sapply(clones$sequence_id, trim_seq_id)
clones <- clones[clones$clone_id %in% c("221871", "116705", "299920", "252697", "10611", "254244", "373110", "39960", "386635"),]
indel_fil$sequence_id_new <- sapply(indel_fil$sequence_id, trim_seq_id)
test <- merge(clones, indel_fil, by = "sequence_id_new")
#299920 has indel

library(ggseqlogo)
library(ggpubr)

custom.seqlogo = function(myclones){
  all_plots = lapply(seq_along(myclones), function(i){
    x = myclones[[i]]
    x$v_gene = gsub('Homsap ', '', x$v_call)
    x$v_gene = gsub('\\*.*', '', x$v_gene)
    cdr3_clones = split(x, f = x$group_id)
    mytitle = paste(unique(x$v_gene), ' (', unique(x$clone_id), ')', sep = '')
    print(names(cdr3_clones))
    p1 = ggseqlogo(cdr3_clones[[names(cdr3_clones)[1]]]$cdr3_aa, seq_type = 'aa', method = 'prob') + 
      labs(title = paste(names(cdr3_clones)[1], ' (SHM: ', 
                         round(mean(cdr3_clones[[names(cdr3_clones)[1]]]$mu_freq), 2), ')', sep = '')) + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.text = element_text(size = 16),
            axis.title.y = element_text(size = 18))
    p2 = ggseqlogo(cdr3_clones[[names(cdr3_clones)[2]]]$cdr3_aa, seq_type = 'aa', method = 'prob') +
      labs(title = paste(names(cdr3_clones)[2], ' (SHM: ', 
                         round(mean(cdr3_clones[[names(cdr3_clones)[2]]]$mu_freq), 2), ')', sep = '')) + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.text = element_text(size = 16),
            axis.title.y = element_text(size = 18))
    # p3 = ggseqlogo(cdr3_clones[['Memory']]$cdr3_aa, seq_type = 'aa', method = 'prob') +
    #   labs(title = paste('Memory B-Cell', ' (SHM: ', 
    #                      round(mean(cdr3_clones[['Memory']]$mu_freq), 2), ')', sep = '')) + 
    #   theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
    #         axis.text = element_text(size = 16),
    #         axis.title.y = element_text(size = 18))
    p_final = p1 + p2 + plot_layout(nrow = 2, guides = 'collect') + 
      plot_annotation(title = mytitle) & 
      theme(legend.position = 'none', 
            plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))
    return(p_final)
  })
  
}

all_plots = custom.seqlogo(clones_ofint)
names(all_plots) <- names(clones_ofint)

lapply(seq_along(all_plots), function(i){
  x <- all_plots[[i]]
  ggsave(paste("plot_ver1/ggseqlogo_",names(all_plots)[i], ".pdf", sep = ""), x, width = 18, height = 9, units = "in")
  return(NULL)
})









