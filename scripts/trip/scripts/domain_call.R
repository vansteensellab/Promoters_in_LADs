#!/usr/bin/env Rscript

library(argparse)
suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(rtracklayer)))

parser <- ArgumentParser(description=paste('call sites in domain, calculate',
                                           'distance to border, domain width'))

parser$add_argument('--regions', help='bed file with regions of interest')
parser$add_argument('--domains', help='bed files with domains ("," seperated)')
parser$add_argument('--out', help='output')

argv <- parser$parse_args()
save(argv, file='test.Rdata')

region_gr = import.bed(argv$regions)


domain_vec = strsplit(argv$domains, ',')
domain_list = lapply(domain_vec, import.bed)

if (length(domain_list) > 1){
    domain_gr = Reduce(intersect, domain_list)
} else {
    domain_gr = domain_list[[1]]
}


o = findOverlaps(region_gr, domain_gr)
d = distanceToNearest(region_gr, domain_gr)

dt = data.table(name = region_gr$name,
                distToNearest = mcols(d)$distance,
                domain_index = NaN)

dt$domain_index[queryHits(o)] = subjectHits(o)
dt[, distToStart := start(eval(region_gr)) -
                    start(eval(domain_gr))[domain_index]]

dt[, distToEnd := end(eval(domain_gr))[domain_index] -
                      end(eval(region_gr))]

dt[,border_distance := ifelse(is.na(distToStart), distToNearest,
                              ifelse(distToEnd<distToStart, distToEnd,
                                     distToStart))]


dt[,domain_width := width(eval(domain_gr))[domain_index]]

fwrite(dt[,.(name, domain_index, border_distance, domain_width)],
       file=argv$out, sep='\t', quote=F, na="NA")
