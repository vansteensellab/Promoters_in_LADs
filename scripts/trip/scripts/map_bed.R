#!/usr/bin/env Rscript
library(rtracklayer)
library(plyr)
options(scipen=999)
args = commandArgs(trailingOnly=T)
#
# args = c('/DATA/usr/c.leemans/projects/trip/cl20180409_TRIP_K562_evsr/mapping/11_A_r1.2.table',
#          '/DATA/usr/c.leemans/projects/trip/cl20180409_TRIP_K562_evsr/mapping/11_A_r2.2.table')

in_list = lapply(args, read.table, stringsAsFactors=F, header=T)

in_table = do.call(rbind, in_list)

barcode_table = ddply(in_table, .(barcode), function(data){
                            seqname = data$seqname[1]
                            ori = data$ori[1]
                            start = data$start_pos[1]
                            total = sum(data$total_mapped)
                            if (all(data$seqname==seqname) && all(data$ori==ori) &&
                                    all(data$start_pos==start)){
                                return(c(seqname, start-1, start, ori, total))
                            }})
barcode_table = barcode_table[with(barcode_table, order(-as.numeric(V5))),]

bc_table = ddply(barcode_table, .(V1, V2, V3, V4), function(x){x[1,]})

bc_table = bc_table[with(bc_table, order(V1, as.numeric(V2))),]

bed_table = data.frame(bc_table[,2:4],
                       bc_table[,1],
                       '0',
                       bc_table[,5])

write.table(bed_table, col.names=F, row.names=F, quote=F, sep='\t')
