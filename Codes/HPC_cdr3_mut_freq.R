library(alakazam)
library(shazam)

path = 'imgt_output_HC/makedb/combined_data/clones/germlines/combined_clone-pass_germ-pass.tsv'

all_dat = readChangeoDb('combined_clone_germ-pass.tsv')


test = split(all_dat, f = all_dat$group_id)

clones = sapply(test, function(x){
  x = unique(x$clone_id)
})

#pick 330_2018_sample3 which has highest overlapping clones with other subjects
# approx 78600

all_dat = all_dat[all_dat$group_id %in% c('329_2015', '329_2016', '329_2018', 
                                          '330_2015', '330_2016', 
                                          '330_2018_sample3'), ]


all_dat = aminoAcidProperties(all_dat, seq="cdr3", nt=TRUE, trim=TRUE, 
                              label="cdr3", property = 'length')

#remove NAs
all_dat = all_dat[!is.na(all_dat$cdr3_aa_length), ]

all_dat_shm = observedMutations(all_dat, 
                                sequenceColumn="sequence_alignment",
                                germlineColumn="germline_alignment_d_mask",
                                regionDefinition=IMGT_V, frequency=TRUE,
                                combine=TRUE, nproc=120)


