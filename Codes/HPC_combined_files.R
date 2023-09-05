library(alakazam)

all_files  = list.files(pattern = '.tsv')

filnames = gsub('(.*)_db.*', '\\1', all_files)

all_dat = lapply(seq_along(all_files), function(i){
  dat = readChangeoDb(all_files[i])
  dat$group_id = filnames[i]
  dat = dat[dat$productive == TRUE, ]
  return(dat)
})

all_dat = bind_rows(all_dat)

saveRDS(all_dat, 'combined_productive.rds')

write.table(all_dat, 'combined_data/combined_productive.tsv', sep = '\t', 
            row.names = F)

##calculate cloning threshold

dist_ham <- distToNearest(all_dat, sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=120)

output <- findThreshold(dist_ham$dist_nearest, method="density")
threshold <- output@threshold
#0.07398438

p1 <- ggplot(subset(dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=threshold, color="firebrick", linetype=2)
