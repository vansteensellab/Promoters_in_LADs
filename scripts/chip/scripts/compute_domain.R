#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='Process counts into normalized mean and geometric means')

parser$add_argument('--idxstats', help='bam files of samples for calculating coverage')
parser$add_argument('--counts', help='count file with counts in region')


argv <- parser$parse_args()

save(argv, file='test.Rdata')


stats = read.table(argv$idxstats, stringsAsFactors=F, header=T)
stats = stats[-grep('[*]|chrM', stats$seqnames), ]
total = colSums(stats[,3:ncol(stats)])

counts = read.table(argv$counts, stringsAsFactors=F, header=T)

exp_counts = counts[,-c(1:2, grep('input', colnames(counts))), drop=F]
exp_total = total[-c(grep('input', names(total)))]

norm_list = lapply(1:ncol(exp_counts), function(i){
    t_vec = c(exp_total[i], total['input'])
    norm = t(rbind(exp_counts[,i], counts$input) / t_vec * min(t_vec))
    return(norm[,1] / norm[,2])
})

gm_mean = function(x){
  prod(x)^(1/length(x))
}

norm_table = do.call(cbind, norm_list)

mean_table = data.frame(domain=counts$name,
                        mean=rowMeans(norm_table),
                        gm_mean = apply(norm_table, 1, gm_mean))

write.table(mean_table, quote=F, sep='\t', row.names=F)
