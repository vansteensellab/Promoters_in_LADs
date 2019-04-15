args <- commandArgs(trailingOnly = TRUE)

source(paste0(args[5], "hiddenDomains.R"))


##print(paste("treat.bin.file:", args[1]))
##print(paste("control.bin.file:", args[2]))
##print(paste("out.file.name:", args[3]))
##print(paste("chr.file:", args[4]))
##print(paste("path:", args[5]))

chr.data <- read.delim(file=args[4], sep="\t", stringsAsFactors=FALSE, header=FALSE)


hiddenDomains(treat.bin.file=args[1], control.bin.file=args[2], out.file.name=args[3], chr.names=chr.data[,1])