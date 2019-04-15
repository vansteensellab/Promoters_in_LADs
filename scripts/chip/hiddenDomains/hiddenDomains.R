##========================================================================
##
## Author: Joshua Starmer <josh.starmer@gmail.com>, 2014
## 
## Copyright (C) 2014, Joshua Starmer
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
## 
##========================================================================

library("HiddenMarkov")
library("depmixS4")

failed.to.converge <- FALSE

hiddenDomains <- function(treat.bin.file=NULL, control.bin.file=NULL, normalize=TRUE, chr.names=NULL, skip.chrs=NULL, min.prob=0.6, max.read.count=200, min.read.count=-10, out.file.name="results.txt", debug=FALSE, treat.data=NULL, control.data=NULL) {
  
  if (is.null(chr.names)) {
    cat("You must provide a non-null value to the \'chr.names\' argument.\n")
    cat("You can either use \'mouse\', \'human\', or specify a vector of strings.\n")    
    return("Exiting")
  }
  
  if (is.null(treat.data)) {
    print(paste("Reading treatment file:", treat.bin.file))
    treat.data <- read.delim(file=treat.bin.file, sep="\t",
                             stringsAsFactors=TRUE, header=TRUE)
  } else {
    print("Using treatment data provided...") 
  }
  
  if (!is.null(control.bin.file) || !is.null(control.data)) {
    if (is.null(control.data)) {
      print(paste("Reading control file:", control.bin.file))
      control.data <- read.delim(file=control.bin.file, sep="\t",
                                 stringsAsFactors=TRUE, header=TRUE)
    } else {
      print("Using control data provided...") 
    }
    print("Merging treatment and control data")
    merged.data <- merge(treat.data, control.data, 
                         by.x=c("chr", "pos"),
                         by.y=c("chr", "pos"), all.x=TRUE)
    ## add pseudo counts to control bins that have no reads
    merged.data[is.na(merged.data$count.y),]$count.y <- 1
    if (normalize) {      
      ## normalize reads...
      print("Normalizing reads")
      norm.factor <- 1000000
      
      data.read.count <- sum(merged.data$count.x)
      control.data.read.count <- sum(merged.data$count.y)
      
      smaller.count <- data.read.count
      if(smaller.count > control.data.read.count) {
        smaller.count <- control.data.read.count 
      }
      
      if (smaller.count < 10000000) {
        norm.factor <- 1000000
      } else if (smaller.count < 100000000) {
        norm.factor <- 10000000 
      } else if (smaller.count < 1000000000) {
        norm.factor <- 100000000
      }
      
      merged.data$count.x <- merged.data$count.x / 
        (data.read.count / norm.factor)
      merged.data$count.y <- merged.data$count.y / 
        (control.data.read.count / norm.factor)
      if (debug) {
        treat.scale.factor <- (data.read.count / norm.factor)
        control.scale.factor <- (control.data.read.count / norm.factor)
        cat("treat count total: ")
        cat(data.read.count, "\n")
        cat("treat.scale.factor: ")
        cat(treat.scale.factor, "\n")
        cat("control count total: ")
        cat(control.data.read.count, "\n")
        cat("control.scale.factor: ")
        cat(control.scale.factor, "\n")
      }
    }
    
    data <- data.frame(id=merged.data$id.x,                       
                       chr=merged.data$chr,
                       pos=merged.data$pos,
                       count=(merged.data$count.x - merged.data$count.y))
    
    ##when I set the min value to 0, depmix throws an error, so don't do this
    ##data[data$count < 0,]$count = 0
  } else {
    data <- treat.data 
  }
  
  if (length(chr.names) == 1) {
    if (chr.names == "mouse") {
      chr.names = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                    "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                    "chrX", "chrY")
    } else if (chr.names == "human") {
      chr.names = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                    "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                    "chr20", "chr21", "chr22", "chrX", "chrY")
    }
  }
  
  #####################
  ##
  ## Estimate HMM parameters
  ##
  #####################
  all.parameters <- data.frame(chr=character(),
                               pr1=numeric(),
                               pr2=numeric(),
                               s1tos1=numeric(),
                               s1tos2=numeric(),
                               s2tos1=numeric(),
                               s2tos2=numeric(),
                               s1mean=numeric(),
                               s1sd=numeric(),
                               s2mean=numeric(),
                               s2sd=numeric())
  
  failed.chrs <- c()
  for (chr.name in chr.names) {
    if (chr.name %in% skip.chrs) {
      cat("Skipping", chr.name, " because you told me to.\n")
      next
    }
    set.seed(1)
    cat("\n")
    print(paste("Estimating parameters for:", chr.name))
    chr.data <- data[data$chr == chr.name,]
    
    if ((nrow(chr.data) == 0) || (max(chr.data$count) < 1)) {
      cat("WARNING!!!  There was not enough data to call domains on chromosome:", chr.name, "\n")
      cat("Skipping", chr.name, "\n")
      next 
    }
    
    if (max(chr.data$count) > max.read.count) {
      ## I have to do this because depmix will die if a bin has too many
      ## reads in it.
      chr.data[chr.data$count > max.read.count,]$count <- max.read.count
    }
    if (min(chr.data$count) < min.read.count) {
      chr.data[chr.data$count < min.read.count,]$count <- min.read.count
    }
    
    ## NOTE: each chromosome has so much data that the priors made from
    ## one chromosome will be washed out by the data, so you might as well
    ## start from scratch each time.  It's easier, and the results are
    ## the same.
    
    #################################################
    ##
    ## This is where we estiomate parameters for the HMM
    ##
    #################################################
    failed.to.converge <<- FALSE
    
    #####################
    ##
    ## Trying library: depmix
    ##
    #####################
    print("Trying the depmix library")
    parameters <- tryDepmixEstimate(chr.data, min.prob, debug)
    
    if(!failed.to.converge) {
      #print("these are the depmix parameters...")
      #print(parameters)
      if (!is.null(parameters) && length(parameters)==10) {
        all.parameters <- rbind(all.parameters, 
                                data.frame(chr=chr.name, parameters))
      }
      next      
    } else {
      #print("depmix failed... trying hiddenMarkov") 
    }
    
    #####################
    ##
    ## Trying library: HiddenMarkov
    ##
    #####################
    print("Trying the HiddenMarkov library")
    failed.to.converge <<- FALSE
    parameters <- tryHiddenEstimate(chr.data, min.prob)
    
    if(!failed.to.converge) {
      #print("these are the hiddenMarkov parameters...")
      #print(parameters)
      if (!is.null(parameters) && length(parameters)==10) {
        all.parameters <- rbind(all.parameters, 
                                data.frame(chr=chr.name, parameters))
      }
      next
    } else {
      failed.chrs <- c(failed.chrs, chr.name)
    }
  }
  
  #####################
  ##
  ## now determine the enriched/depleted states for each chromosome.
  ##
  #####################
  print("Determing enriched/depleted states for each chromosome")
  ## 1) average the parameter estimates.
  
  delta <- c(1, 0)
  #delta <- c(0, 1)
  
  Pi <- matrix(c(mean(all.parameters$s1tos1),
                 mean(all.parameters$s1tos2),
                 mean(all.parameters$s2tos1),
                 mean(all.parameters$s2tos2)), 
               nrow=2, byrow=TRUE)
  pm <- list(mean=c(mean(all.parameters$s1mean), 
                    mean(all.parameters$s2mean)),
             sd=c(mean(all.parameters$s1sd),
                  mean(all.parameters$s2sd)))
  
  avg.parameters = list(delta=delta, Pi=Pi, pm=pm)
  
  if (debug) {
    print("parameters:")
    print(all.parameters)
    print(avg.parameters)
  }
  
  final.results <- data.frame(id=character(),
                              chr=character(),
                              pos=numeric(),
                              count=numeric(),
                              state=numeric(),
                              enriched.pr=numeric())
  
  failed.chr.results <- c()
  for (chr.name in chr.names) {
    if (chr.name %in% skip.chrs) {
      next
    }
    set.seed(1)
    cat("\n")
    print(paste("Working on:", chr.name))
    chr.data <- data[data$chr == chr.name,]
    
    if ((nrow(chr.data) == 0) || (max(chr.data$count) < 1)) {
      next 
    }
    
    if (max(chr.data$count) > max.read.count) {
      ## I have to do this because depmix will die if a bin has too many
      ## reads in it.
      chr.data[chr.data$count > max.read.count,]$count <- max.read.count
    }
    if (min(chr.data$count) < min.read.count) {
      chr.data[chr.data$count < min.read.count,]$count <- min.read.count
    }
    
    if (!(chr.name %in% failed.chrs)) {
      Pi <- matrix(c(all.parameters[all.parameters$chr==chr.name,]$s1tos1,
                     all.parameters[all.parameters$chr==chr.name,]$s1tos2,
                     all.parameters[all.parameters$chr==chr.name,]$s2tos1,
                     all.parameters[all.parameters$chr==chr.name,]$s2tos2), 
                   nrow=2, byrow=TRUE)
      pm <- list(mean=c(all.parameters[all.parameters$chr==chr.name,]$s1mean, 
                        all.parameters[all.parameters$chr==chr.name,]$s2mean),
                 sd=c(all.parameters[all.parameters$chr==chr.name,]$s1sd,
                      all.parameters[all.parameters$chr==chr.name,]$s2sd))
      
      parameters = list(delta=delta, Pi=Pi, pm=pm)
    } else {
      cat("Using avgerage parameters for ", chr.name, "\n")
      parameters = avg.parameters 
    }
    
    x.with.fit <- dthmm(chr.data$count, parameters$Pi,
                        parameters$delta, "norm", parameters$pm, discrete=TRUE)
    states <- Viterbi(x.with.fit)
    
    x.results <- forwardback(x.with.fit$x, parameters$Pi, parameters$delta, "norm", parameters$pm)
    
    hmm.results <- data.frame(chr.data)
    
    enriched.pr <- exp((x.results$logalpha[,2] + x.results$logbeta[,2]) - x.results$LL)
    hmm.results <- cbind(hmm.results, 
                         state=states, 
                         enriched.pr=round(enriched.pr, digits=4))
    
    ## round the counts...
    hmm.results$count <- round(hmm.results$count, digits=4)
    
    ## remove bins with NaN for probabilities...
    hmm.results <- hmm.results[!is.nan(hmm.results$enriched.pr), ]
    
    if (debug) {
      print(head(hmm.results))
    }
    
    
    
    if (nrow(hmm.results) == 0) {
      cat("\n\nWARNING!!!\n")
      cat("Failed to find any domains on ", chr.name, "\n")
      failed.chr.results <- c(failed.chr.results, chr.name)
    } else {
      
      ## we want bins with a postive number of reads
      hmm.results <- hmm.results[hmm.results$count > 0, ]
      
#       state1.mean <- mean(hmm.results[(hmm.results$state==1 & 
#                                          hmm.results$enriched.pr == 0), "count"])
#       state1.count <- nrow(hmm.results[(hmm.results$state==1 & 
#                                           hmm.results$enriched.pr == 0), "count"])
#       state2.mean <- mean(hmm.results[(hmm.results$state==2 &
#                                          hmm.results$enriched.pr == 1), "count"])
#       state2.count <- nrow(hmm.results[(hmm.results$state==2 & 
#                                           hmm.results$enriched.pr == 0), "count"])
      
#       if (debug) {
#         print(paste("state1.mean:", state1.mean))
#         print(paste("state2.mean:", state2.mean))
#       }
#       
#       if (stat1.count > 10) {
#         if (!is.nan(state1.mean) & (
#           is.nan(state2.mean) | 
#           (state1.mean > state2.mean))) {
#           
#           hmm.results$enriched.pr = 1 - hmm.results$enriched.pr
#         }
#       }
      if (parameters$pm$sd[1] > parameters$pm$sd[2]) {
        if(debug) {
          print("flipping the posterior probabilities...")
#          print(paste("parameters$pm$sd[1]:", parameters$pm$sd[1]))
#          print(paste("parameters$pm$sd[2]:", parameters$pm$sd[2]))
        }
        hmm.results$enriched.pr = 1 - hmm.results$enriched.pr
      }
      ## we only want bins that have a good posterior probability
      hmm.results <- hmm.results[hmm.results$enriched.pr > min.prob, ]
      
      final.results <- rbind(final.results, hmm.results)
    }
  }
  
  
  write.table(final.results, file=out.file.name, sep="\t", quote=FALSE, row.names=FALSE)
  cat("\n")
  cat("Hooray!  All done!  Check", out.file.name, "for results\n")
  cat("\n")
  
  if (length(failed.chr.results) > 0) {
    cat("\n")
    cat("Warning, there were no results for the following chromosomes:\n")
    print(failed.chr.results)
    cat("\nConsider changing max.read.count to a smaller value, like 10\n")
    cat("If you want to be more rational about setting that parameter,\n")
    cat("eyeball the normalized bigWig track for the data in the UCSC\n")
    cat("Genome Browser to get an idea of what cutoff is appropriate for peaks.\n")
    cat("\n")
  }
  
}

tryDepmixEstimate <- function(chr.data, min.prob, debug) {
  chr.data.mean <- mean(chr.data$count)
  chr.data.sd <- sd(chr.data$count)
  if (debug) {
    print(paste("Mean:", chr.data.mean))
    print(paste("SD:", chr.data.sd))
  }
  
  return(tryCatch({
    
    mod <- depmix(response = count ~ 1, data=chr.data, nstates=2,
                  instart=c(0, 1))#,
                  # respstart=c(chr.data.mean*3, chr.data.sd, 
                  #          chr.data.mean, chr.data.sd))
    
    fm <- fit(mod, verbose=TRUE, emcontrol=em.control(random.start=TRUE))
    #fm <- fit(mod, verbose=TRUE)
    #summary(fm)
    
    parameter.df <- data.frame(pr1=getpars(fm)[1],
                               pr2=getpars(fm)[2],
                               s1tos1=getpars(fm)[3],
                               s1tos2=getpars(fm)[4],
                               s2tos1=getpars(fm)[5],
                               s2tos2=getpars(fm)[6],
                               s1mean=getpars(fm)[7],
                               s1sd=getpars(fm)[8],
                               s2mean=getpars(fm)[9],
                               s2sd=getpars(fm)[10])
    
    if (getpars(fm)[7] < getpars(fm)[9]) {
      #print("State 2 was the enriched state...")
      
    } else {
      #print("State 1 was the enriched state...") 
      
      parameter.df <- data.frame(pr1=getpars(fm)[2],
                                 pr2=getpars(fm)[1],
                                 s1tos1=getpars(fm)[6],
                                 s1tos2=getpars(fm)[5],
                                 s2tos1=getpars(fm)[4],
                                 s2tos2=getpars(fm)[3],
                                 s1mean=getpars(fm)[9],
                                 s1sd=getpars(fm)[10],
                                 s2mean=getpars(fm)[7],
                                 s2sd=getpars(fm)[8])
      
    }
    
    parameter.df
    
  }, warning = function(war) {
    # do nothing
  }, error = function(err) {
    ##print(err, "\n")
    failed.to.converge <<- TRUE
  }, finally = {
    # do nothing
  }) ## end of tryCatch(
  ) ## end of return(
  
}

tryHiddenEstimate <- function(chr.data, min.prob) {
  chr.data.mean <- mean(chr.data$count)
  chr.data.sd <- sd(chr.data$count)
  
  return(tryCatch({
    
    Pi <- matrix(c(0.99, 0.01, 0.03, 0.97), byrow=TRUE, nrow=2)
    delta <- c(1,0)
    x <- dthmm(chr.data$count, Pi, delta, "norm", 
               list(mean=c(chr.data.mean, 3*chr.data.mean), 
                    sd=c(chr.data.sd, 3*chr.data.sd)), 
               discrete = TRUE)
    y <- BaumWelch(x, bwcontrol(prt=FALSE))
    #print(summary(y))
    
    pr1 <- y$delta[1]
    pr2 <- y$delta[2]
    s1tos1 <- y$Pi[1,1]
    s1tos2 <- y$Pi[1,2]
    s2tos1 <- y$Pi[2,1]
    s2tos2 <- y$Pi[2,2]
    s1mean <- y$pm$mean[1]
    s1sd   <- y$pm$sd[1]
    s2mean <- y$pm$mean[2]
    s2sd   <- y$pm$sd[2]
    
    parameters <- c()
    if (y$pm$mean[1] < y$pm$mean[2]) {
      #print("State 2 was the enriched state...")
      parameters <- data.frame(pr1=pr1,
                               pr2=pr2,
                               s1tos1=s1tos1,
                               s1tos2=s1tos2,
                               s2tos1=s2tos1,
                               s2tos2=s2tos2,
                               s1mean=s1mean,
                               s1sd=s1sd,
                               s2mean=s2mean,
                               s2sd=s2sd)
    } else {
      #print("State 1 was the enriched state...") 
      parameters <- data.frame(pr1=pr2,
                               pr2=pr1,
                               s1tos1=s2tos2,
                               s1tos2=s2tos1, 
                               s2tos1=s1tos2,
                               s2tos2=s1tos1,
                               s1mean=s2mean,
                               s1sd=s2sd,
                               s2mean=s1mean,
                               s2sd=s1sd)
    }
    
    parameters    
  }, warning = function(war) {
    # do nothing
  }, error = function(err) {
    failed.to.converge <<- TRUE 
  }, finally = {
    # do nothing
  }) ## end of tryCatch(
  ) ## end return(
}
