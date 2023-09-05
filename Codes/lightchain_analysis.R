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

#CDR3 length
cdr3_plotdat = abstar_dat_lightC[!is.na(abstar_dat_lightC$cdr3_length),]
cdr3_plotdat <- cdr3_plotdat[cdr3_plotdat$cdr3_length <= 25,]
cdr3_means = cdr3_plotdat %>%
  group_by(group_id, chain) %>%
  summarise(mean = mean(cdr3_length), .groups = 'drop')

cdr3_means$mean = round(cdr3_means$mean, 2)
write_xlsx(cdr3_means, "suppliment/cdr3_mean_lightchain.xlsx")

p1 = ggplot(cdr3_plotdat, aes(x=cdr3_length, fill = group_id)) + 
  geom_density(bw = 0.8, alpha = 0.8) +
  facet_grid(group_id~chain) +
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

ggsave('plot_ver1/cdr3_length_LC.png', p1, width = 9.5, height = 9.5, units = 'in')
ggsave('plot_ver1/cdr3_length_LC.pdf', p1, width = 9.5, height = 9.5, units = 'in')

##gene usage
vcall_lightC = countGenes(abstar_dat_lightC, gene = 'v_gene', mode = 'gene', groups = 'group_id')
vcall_lightC$seq_freq = round(vcall_lightC$seq_freq*100, 2)
vcall_lightC$call = 'v_gene'

jcall_lightC = countGenes(abstar_dat_lightC, gene = 'j_gene', mode = 'gene', groups = 'group_id')
jcall_lightC$seq_freq = round(jcall_lightC$seq_freq*100, 2)
jcall_lightC$call = 'j_gene'
genecall_lightC_dat = rbind(vcall_lightC, jcall_lightC)

write_xlsx(list("vcall" = vcall_lightC, "jcall" = jcall_lightC), "suppliment/genecall_lightchain.xlsx")
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

genecall_plasma_vgene_lightC = genecall_lightC_dat[genecall_lightC_dat$call == 'v_gene', ]
#genecall_plasma_dgene_lightC = genecall_lightC_dat[genecall_lightC_dat$call == 'd_gene', ]
genecall_plasma_jgene_lightC = genecall_lightC_dat[genecall_lightC_dat$call == 'j_gene', ]

p1 = geneusage_plots(genecall_plasma_vgene_lightC, 'V Gene Usage', 'none')
p2 = geneusage_plots(genecall_plasma_jgene_lightC, 'J Gene Usage', 'none')

p_final = p1 + p2 + plot_layout(ncol = 2, widths = c(1, 0.4)) + 
  plot_annotation(title = 'Gene Usage Light Chain') &
  theme(plot.title = element_text(face = 'bold', size = 20, hjust = 0.5), 
        legend.position = 'none') 

ggsave('plot_ver1/gene_usage_lightC.png', p_final, width = 18, height = 8, units = 'in')
ggsave('plot_ver1/gene_usage_lightC.pdf', p_final, width = 18, height = 8, units = 'in')
