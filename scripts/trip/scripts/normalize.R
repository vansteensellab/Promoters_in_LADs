#!/usr/bin/R
args = commandArgs(trailingOnly=T)
cDNA_input = read.table(args[1], col.names=c('barcode','count', 'other_list'),
                        stringsAsFactors=F, fill=T)
gDNA_input = read.table(args[2], col.names=c('barcode','count', 'other_list'),
                        stringsAsFactors=F, fill=T)

# cDNA_input = read.table("../cl20180102_TRIP_K562_evsr/cDNA/6_B_r1.starcode.count",
#                         col.names=c('barcode','count', 'other_list'),
#                         stringsAsFactors=F, fill=T)
#
# gDNA_input = read.table("../cl20180102_TRIP_K562_evsr/gDNA/6_B_r1.starcode.count",
#                         col.names=c('barcode','count', 'other_list'),
#                         stringsAsFactors=F, fill=T)

if (length(args) == 4){
  spike_input = read.table(args[3], col.names=c('barcode','count', 'other_list'),
                           stringsAsFactors=F, fill=T)
  spike_sum = sum(spike_input$count)
  output = args[4]
} else if (length(args) == 3){
  output = args[3]
}


gDNA_data = gDNA_input[gDNA_input$count > 0,]
match_vec = match(cDNA_input$barcode, gDNA_data$barcode)
cDNA_data = cDNA_input[!is.na(match_vec), ]

count_table = cbind.data.frame(cDNA_data$count + 1,
                               gDNA_data$count[match_vec[!is.na(match_vec)]])
col_sum = colSums(count_table)
cpm_table = t(t(count_table) / col_sum * min(col_sum))
normalized_by_gDNA = cpm_table[,1] / cpm_table[,2]
output_data = cbind.data.frame(cDNA_data$barcode, count_table, cpm_table,
                               normalized_by_gDNA)
colnames(output_data) = c('barcode', 'cDNA_count', 'gDNA_count', 'cDNA_cpm', 'gDNA_cpm', 'normalized_by_gDNA')
if (length(args)==4){
  output_data$normalized_by_spike = output_data$normalized_by_gDNA / spike_sum * 1000000
}

write.table(output_data, file=output, row.names=F, quote=F, sep='\t')
