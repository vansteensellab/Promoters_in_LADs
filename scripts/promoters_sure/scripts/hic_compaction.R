#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(GENOVA)))
suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(sqldf)))


parser <- ArgumentParser(description=paste('Calculate compaction scores from',
                                           'HiC data in a window around a site',
                                           'of interest by fitting a slope on',
                                           'the median HiC signal at each bin distance',
                                           'measure in the window up to a max distance.'))

parser$add_argument('--signal', '-s',
                    help='HiC signal file (GENOVA)')
parser$add_argument('--indices', '-i',
                    help='HiC indeces file (GENOVA)')
parser$add_argument('--window', '-w', type='integer',
                    help='size of the window used')
parser$add_argument('--max-dist', '-d', type='integer',
                    help='maximum distance between bins used to calculate the slope')
parser$add_argument('--bed', '-b',
                    help='bed file with sites of interest')
parser$add_argument('--out', '-o',
                    help='output file')

argv <- parser$parse_args()
save(argv, file='test.Rdata')


shift = argv$window / 2
command = paste("awk -vOFS='\t' '{",
                "    start=$2-", shift, ";",
                "    end=$3+", shift, ";",
                "    print $1, start<0 ? 0 : start, end, $4, $5, $6",
                "}'", argv$bed, '| bedtools intersect -wb -a - -b ',
                argv$indices)

bed = fread(command, stringsAsFactors=F, select=c(1:4, 10),
            col.names=c('seqnames', 'start', 'end', 'names', 'bin'))
setkey(bed, names)


sig = fread(argv$sig, stringsAsFactors=F,
            col.names=c('bin1', 'bin2', 'signal'))

resolution = max(bed[, end-start])
bin_dist = argv$max_dist/resolution
sig = sig[bin2-bin1 < eval(bin_dist) & bin2!=bin1, ]

setkey(sig, bin1)
combined = bed[, sig[bin1%in%bin,][bin2%in%bin, ], by=names]

combined$dist = (combined$bin2 - combined$bin1) * resolution
median_dt = combined[,list(median=median(signal)), by=list(names, dist)]

run_lm <- function(dt){
    fit = lm(log10(median) ~ log10(dist), dt)
    return(list(r2=summary(fit)$r.squared, slope=fit$coefficients[2]))
}


lm_fit = median_dt[,run_lm(.SD),by=names]

write.table(lm_fit, file=argv$out, quote=F, row.names=F, sep='\t')
