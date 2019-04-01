#!/usr/bin/env Rscript

library(argparse)
suppressMessages(suppressWarnings(require(data.table)))

selection = args[1]
link = args[2]
dist = args[3]
data <- fread("file:///dev/stdin", fill=T)


dt = data[data[, list(I=.I[which.min(V17)]), by=c('V1', 'V6', 'V7', 'V11')]$I,]

dt = dt[V17<(dist), ]

write.table(dt[,c(1,3,4,6,7)], file=selection, sep='\t', quote=F, row.names=F,
            col.names=F)
write.table(dt[,c(4,11)], file=link, sep='\t', quote=F, row.names=F,
            col.names=F)
