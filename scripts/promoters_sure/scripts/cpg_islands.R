#!/usr/bin/env Rscript
library(argparse)
library(data.table)
parser <- ArgumentParser(description='Process counts into normalized mean and geometric means')

parser$add_argument('-i', help='input')
parser$add_argument('-o', help='count file with counts in region')

argv <- parser$parse_args()
save(argv, file='test.Rdata')

input_df = read.table(argv$i, stringsAsFactors=F, sep='\t')

prom_col = c('seqnames', 'start', 'end', 'transcript_id',
            'score', 'strand', 'cpg_seqnames', 'cpg_start',
            'cpg_end', 'cpg_name', 'length', 'cpgNum', 'cgNum',
            'perCpg', 'perGc', 'obsExp')

enh_col = prom_col[!prom_col%in%c('score', 'strand')]

if (ncol(input_df)==length(prom_col)){
    colnames(input_df) = prom_col
} else {
    colnames(input_df) = enh_col
}

col_vec = grep('length', colnames(input_df)):ncol(input_df)
input_df[,col_vec][input_df[,col_vec]=='.'] = '0'

input_df[,col_vec] = lapply(col_vec, function(i){as.numeric(input_df[,i])})

input_dt = data.table(input_df)
setkey(input_dt, 'transcript_id')
result = input_dt[,.SD[which.max(length),],
                       by=transcript_id][,c('transcript_id', 'cpg_name', 'length',
                                           'cpgNum', 'cgNum', 'perCpg')]

fwrite(result, argv$o, sep='\t', quote=F, row.names=F)
