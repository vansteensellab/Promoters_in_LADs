#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(rtracklayer)))
args = commandArgs(trailingOnly=T)
#
# args = c('/DATA/usr/c.leemans/projects/chip_snake/hg38/domains/K562/Snyder2012/MYC_ENCSR000EGJ_r1_SRR502421_analysis.bed',
#          '/DATA/usr/c.leemans/projects/chip_snake/hg38/domains/K562/Snyder2012/MYC_ENCSR000EGJ_r2_SRR502422_analysis.bed')

if (length(args) > 1){
    in_list = lapply(args, function(fn){
                         df = tryCatch(read.table(fn, stringsAsFactors=F),
                                       error=function(e) NULL)
                         if (is.null(df)){return(NULL)}
                         if (ncol(df) > 3){
                             rownames(df) = df[,4]
                         } else {
                             rownames(df) = range(1,nrow(df))
                         }
                         return(df)
                     })
    if (any(unlist(lapply(in_list, is.null)))){
        df = NULL
    } else {
        in_gr = lapply(in_list, function(df){
                           gr = GRanges(seqnames = df[,1], IRanges(df[,2], df[,3]))
                           names(gr) = rownames(df)
                           return(gr)
                       })

        shared_gr = in_gr[[1]]
        mcols(shared_gr) = DataFrame(matrix('', nrow=length(shared_gr),
                                            ncol=length(in_gr)))
        mcols(shared_gr)[,1] = names(in_gr[[1]])
        for (i in 2:length(in_gr)){
            o = findOverlaps(shared_gr, in_gr[[i]])
            shared_gr = shared_gr[queryHits(o)]
            other_gr = in_gr[[i]][subjectHits(o)]

            mcols(shared_gr)[,i] = names(other_gr)

            start(shared_gr) = unlist(apply(cbind(start(shared_gr), start(other_gr)), 1, max))
            end(shared_gr) = unlist(apply(cbind(end(shared_gr), end(other_gr)), 1, min))
        }
        if (length(shared_gr) == 0){
            df = data.frame()
        } else{
            df = data.frame(data.frame(shared_gr)[,1:3],
                            name=paste0("peak_intersect", 1:length(shared_gr)),
                            score="0", strand=".")
            for (i in 1:length(in_gr)){
                name_vec = mcols(shared_gr)[,i]
                n = ncol(in_list[[i]])
                if (n > 6){
                    df = data.frame(df, in_list[[i]][name_vec, c(4,7:n)])
                } else{
                    df = data.frame(df, in_list[[i]][name_vec, 4])
                }
            }

        }
    }
} else {
    df = tryCatch(read.table(args[1], stringsAsFactors=F),
                  error=function(e) NULL)
}


write.table(df, col.names=F, row.names=F, quote=F, sep='\t')
