#!/usr/bin/env Rscript

library(argparse)
suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(rtracklayer)))

parser <- ArgumentParser(description='Process counts into normalized mean and geometric means')

parser$add_argument('--counts', help='count file with counts per region')
parser$add_argument('--bed', help='bed file with regions of interest')
parser$add_argument('--window', help='window around intergration')
parser$add_argument('--out', help='output')

argv <- parser$parse_args()

save(argv, file='test.Rdata')

count_vec = unlist(strsplit(argv$counts, ','))
# count_vec = c('/home/t.v.schaik/mydata/proj/tests/results/ts171030_PipelineImprovements/results_alldata/counts/K562_r1_LMNB1-gatc.counts.txt.gz',
#               '/home/t.v.schaik/mydata/proj/tests/results/ts171030_PipelineImprovements/results_alldata/counts/K562_r1_Dam16-gatc.counts.txt.gz')



count_list = lapply(count_vec, function(x){
                        fread(paste('zcat', x), stringsAsFactors=F,
                              col.names=c('seqnames', 'start', 'end', 'count'))})
chr_count = lapply(count_list, function(dt){dt[,.(total=sum(count)),by=.(seqnames)]})

map_gr = import.bed(argv$bed)
# map_gr = import.bed('sites/10_A.bed')

gatc_gr = GRanges(count_list[[1]]$seqnames, IRanges(count_list[[1]]$start,
                                                    count_list[[1]]$end))

o = findOverlaps(map_gr, gatc_gr, maxgap=as.numeric(argv$window)/2)
# o = findOverlaps(map_gr, gatc_gr, maxgap=20000)

gatc_list = lapply(count_list, function(x){
    dt = data.table(name=map_gr$name[queryHits(o)],
                    seqnames=x$seqnames[subjectHits(o)],
                    x[subjectHits(o), 'count'])
    setkey(dt, name)
    sum_dt = dt[,.(seqnames=seqnames[1], sum=sum(count)), by=.(name)]
})
total_list = lapply(1:length(gatc_list), function(i){
    dt = chr_count[[i]][gatc_list[[i]], on=.(seqnames)]
    dt$total
})
gatc_dt = data.table(name=gatc_list[[1]]$name,
                     do.call(cbind, lapply(gatc_list, function(x){x$sum})),
                     do.call(cbind, total_list))


gm_mean = function(x){
      prod(x)^(1/length(x))
}

R = length(count_vec)/2
norm_list = lapply(1:R, function(i){
    i_vec = 0:3 * R + 1 + i
    ct = gatc_dt[,.SD,.SDcols=c(1,i_vec)]
    colnames(ct)[2:5] = c('exp', 'ctrl', 't_exp', 't_ctrl')
    norm_c=ct[,.(exp=exp/t_exp * min(t_exp, t_ctrl) + 1,
                 ctrl=ctrl/t_ctrl * min(t_exp, t_ctrl) + 1), by=name]
    norm_c[,exp/ctrl]
})

norm_table = data.table(name=gatc_dt$name,
                        do.call(cbind, norm_list))

mean_table = norm_table[,.(mean=log2(mean(unlist(.SD))),
                           gm_mean = log2(gm_mean(.SD))), by=name]

write.table(mean_table, file=argv$out, quote=F, sep='\t', row.names=F)
