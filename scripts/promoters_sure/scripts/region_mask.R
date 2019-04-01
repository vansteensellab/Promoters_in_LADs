#!/usr/bin/env Rscript
library(rtracklayer)

opt = commandArgs(trailingOnly = TRUE)
chrom_sizes = read.table(opt[2], row.names=1)
gencode_gr = import.gff(opt[1])
gencode_gr = gencode_gr[gencode_gr$type=='gene']

follow_vec = follow(gencode_gr, ignore.strand=T)
precede_vec = precede(gencode_gr, ignore.strand=T)

regions = data.frame(row.names = gencode_gr$ID,
                     strand = strand(gencode_gr))
left = ifelse(is.na(follow_vec), -Inf,
              (start(gencode_gr) - end(gencode_gr)[follow_vec])/2)
right = ifelse(is.na(precede_vec), Inf,
               (start(gencode_gr)[precede_vec] - end(gencode_gr))/2)

regions$start = ifelse(regions$strand=='+', -left, -right)
regions$end = ifelse(regions$strand=='+', right, left)

write.table(regions, opt[3], sep='\t', col.names=F, quote=F)
