#!/usr/bin/env Rscript
library(argparse)
library(data.table)

parser <- ArgumentParser(description='Process counts into normalized mean and geometric means')

parser$add_argument('--exp', nargs='+', help='bam files of experiment')
parser$add_argument('--exp-label', nargs='+', dest='labels',
                    help='labels for experiments')
parser$add_argument('--input', help='bam file with input for calculating coverage')
parser$add_argument('--out', help='output file')

argv <- parser$parse_args()
save(argv, file='test.Rdata')

fread_idx <- function(file_name){
    fread(paste('samtools idxstats', file_name),
          col.names=c('seqnames', 'length', 'N_mapped', 'N_unmapped'),
          stringsAsFactors=F, key='seqnames')
}

exp_idx = lapply(argv$exp, fread_idx)

input_idx = fread_idx(argv$input)

seqnames = unique(c(unlist(lapply(exp_idx, function(x){x$seqnames}),
                  input_idx$seqnames)))

exp_mapped = do.call(cbind, lapply(exp_idx, function(x){x[seqnames, 'N_mapped']}))
if (all(!is.null(argv$labels))){
    colnames(exp_mapped) = argv$labels
} else{
    colnames(exp_mapped) = paste0('rep', 1:length(exp_idx))
}

result = data.frame(seqnames, input_idx[seqnames, 'length'],
                    exp_mapped, input_idx[,'N_mapped'])
colnames(result)[ncol(result)] = 'input'
write.table(result, argv$out, quote=F, sep='\t', row.names=F)
