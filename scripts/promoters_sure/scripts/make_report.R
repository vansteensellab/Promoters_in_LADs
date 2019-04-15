#!/usr/bin/env Rscript
library(knitr)
opt = commandArgs(trailingOnly = TRUE)
root_dir = getwd()
opts_knit$set(root.dir = root_dir)
knit2html(opt[[1]], output=opt[[2]])
