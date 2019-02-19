#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='Process counts into normalized mean and geometric means')

parser$add_argument('--exp', help='bam files of samples for calculating coverage')
parser$add_argument('--counts', help='count file with counts in region')
parser$add_argument('--input', help='bam file with input for calculating coverage')
parser$add_argument('--bed', help='bed file with regions of interest')

argv <- parser$parse_args()

save(argv, file='test.Rdata')
exp_vec = unlist(strsplit(argv$exp, ','))
idxstats <- function(bam_file){
    line_vec = system(paste('samtools idxstats', bam_file),
                      intern=T)
    total_matrix = do.call(rbind, strsplit(line_vec, '\t'))
    total_df = data.frame(total_matrix, stringsAsFactors=F)
    colnames(total_df) = c('seqnames', 'seq_length', 'n_mapped', 'n_unmapped')
    total_df$seq_length = as.numeric(total_df$seq_length)
    total_df$n_mapped = as.numeric(total_df$n_mapped)
    total_df$n_unmapped = as.numeric(total_df$n_unmapped)
    rownames(total_df) = total_df$seqnames
    return(total_df)
}

exp_stats = lapply(exp_vec, idxstats)
names(exp_stats) = gsub('.*/(.*)/(.*)[.]bam', '\\1_\\2', exp_vec)
names(exp_stats) = gsub('-', '.', names(exp_stats))
input_stats = idxstats(argv$input)

count_table = read.table(argv$counts, header=T, comment.char='',
                         stringsAsFactors=F)

colnames(count_table)[1] = "seqnames"

gm_mean = function(x){
  prod(x)^(1/length(x))
}

norm_list = lapply(names(exp_stats), function(s){
    m_ip = sum(input_stats[, 'n_mapped'])
    m_s = sum(exp_stats[[s]][, 'n_mapped'])
    m_min = min(m_ip, m_s)
    ct = count_table[,c('Input', s)]
    norm_c = t(t(ct) / c(m_ip, m_s) * m_min + 1)
    norm_c[,s] / norm_c[, 'Input']
})

norm_table = do.call(cbind, norm_list)

gm_mean_vec = apply(norm_table, 1, gm_mean)
mean_vec = rowMeans(norm_table)

bed_table = read.table(argv$bed, row.names=4, stringsAsFactors=F)
match_vec = match(paste(bed_table[,1], bed_table[,2]),
                  paste(count_table[,1], count_table[,2]))

mean_table = data.frame(transcript_id=rownames(bed_table),
                        mean=mean_vec[match_vec],
                        gm_mean=gm_mean_vec[match_vec])

write.table(mean_table, quote=F, sep='\t', row.names=F)
