#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(rtracklayer)))



parser <- ArgumentParser(description=paste('create gff for selected transcripts'))

parser$add_argument('--gff', '-g',
                    help='gff file')
parser$add_argument('--transcript', '-t',
                    help='file with transcripts')
parser$add_argument('--column', '-c', type='integer',
                    help='column number of the transcript-id')
parser$add_argument('--out', '-o',
                    help='output file')

argv <- parser$parse_args()
save(argv, file='test.Rdata')

transcript_vec = read.table(argv$transcript, stringsAsFactors=F)[,argv$column]

gff_gr = import.gff(argv$gff)


gff_selection = gff_gr[gff_gr$transcript_id %in% transcript_vec]

export(gff_selection, argv$out, format='bed15')
