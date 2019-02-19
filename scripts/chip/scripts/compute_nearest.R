#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(rtracklayer)))
args = commandArgs(trailingOnly=T)
#
# args = c('/DATA/usr/c.leemans/projects/trip/cl20180409_TRIP_K562_evsr/sites/hPGK_B.bed',
#          '/DATA/usr/c.leemans/projects/chip_snake/hg38/peaks/K562/Schmidl2015/H3K4ME1_chipmentation_shared.txt')
#          '/DATA/usr/c.leemans/projects/chip_snake/hg38/domains/K562/Snyder2012/MYC_ENCSR000EGJ_shared.txt')
df_a = try(read.table(args[1]))
df_b = try(read.table(args[2]))

if (!inherits(df_a, "try-error") &
    !inherits(df_b, "try-error")){
    gr_a = GRanges(df_a[,1], IRanges(df_a[,2], df_a[,3]))
    gr_b = GRanges(df_b[,1], IRanges(df_b[,2], df_b[,3]))
    if (ncol(df_a) > 6 && all(df_b[,6] %in% c(-1, 1, '+', '-', '-1', '1'))){
        strand(gr_a) = df_a[,6]
    }
    if (ncol(df_b) > 6 && all(df_b[,6] %in% c(-1, 1, '+', '-', '-1', '1'))){
        strand(gr_b) = df_b[,6]
    }

    nearest = distanceToNearest(gr_a, gr_b)

    df_a$distance = Inf
    df_a$distance[queryHits(nearest)] = mcols(nearest)$distance

    write.table(df_a, row.names=F, col.names=F, quote=F, sep='\t')
}
