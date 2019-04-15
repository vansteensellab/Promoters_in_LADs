
basedir = "/DATA/usr/c.leemans/projects/Promoters_in_LADs"

library(ggpubr)
library(ggplot2)
library(rtracklayer)
library(plyr)
library(reshape2)
library(ggbeeswarm)
library(knitr)
library(ppcor)
library(data.table)
library(gdata)
library(cowplot)
library(scales)
library(grid)
## get a table with matching sets
## table = complete table to take matching sets from
## class_col = column name of class of interest
## class = name of class to match the set on
## order_on = column name to order on

matchSet <- function(table, class_col, class, order_on, bs=10){
  o_vec = order(table[,order_on])
  o_table = table[o_vec, ]
  set_A = which(o_table[,class_col]==class)
  n = length(o_vec)
  bin_n = ceiling((n - set_A[1] - 1) / bs)
  seq_vec = seq(n-bin_n*bs, n, bs)
  set_B = c()
  for(i in 1:(length(seq_vec)-1)){
      sub_table = o_table[(seq_vec[i] + 1):seq_vec[i + 1], ]
      sub_A = which(sub_table[,class_col]==class)
      if (length(sub_A) < bs/2){
          sub_B = sample(which(sub_table[,class_col]!=class), length(sub_A))
      } else {
          sub_B = which(sub_table[,class_col]!=class)
      }
      set_B = c(set_B, sub_B + seq_vec[i])
  }
  ## can also return o_table[c(setA, setB), ]
  ## but this way order is perserved.
  i_vec = o_vec[c(set_A, set_B)]
  return(table[i_vec[order(i_vec)], ])
}


COLi<-"#00BBFF11" #dot color for iLAD promoters
COL_lad = COL_nad = c("#FF0000", "#0077FF")
names(COL_lad)<-c('LAD', 'iLAD')

names(COL_nad)<-c('NAD', 'iNAD')

#color vector for plotting:
COL_class<-c("#A020F0", "#FFA500", "#006400", "#7e7e7e", "#0077FF")
names(COL_class)<-c("repressed", "escaper", "inactive", 'boundary', 'iLAD')

COL<-c("#A020F0", "#FFA500", "#006400")
names(COL)<-c("repressed", "escaper", "inactive")

P_pro_tssr = read.table(paste0(basedir, '/data/promoter_expression/expression/',
                               'gencode.v27_proseq_tssr.tsv'),
                        stringsAsFactors=F, header=T)

for (col in colnames(P_pro_tssr)){
    P_pro_tssr[which(is.na(P_pro_tssr[,col])),col] = 0
}


P_pro_body = read.table(paste0(basedir, '/data/lad_promoters/expression/',
                               'gencode.v27_proseq_body.tsv'),
                        stringsAsFactors=F, header=T)

for (col in colnames(P_pro_body)){
    P_pro_body[which(is.na(P_pro_body[,col])),col] = 0
}

P_exp = read.table(paste0(basedir, '/data/lad_promoters/expression/',
                          'gencode.v27_stranded_expression.txt.gz'),
                       stringsAsFactors=F, header=T, row.names=1)

for (col in colnames(P_exp)){
    P_exp[which(is.na(P_exp[,col])),col] = 0
}

P_tss = read.table(paste0(basedir, '/data/lad_promoters/selection/',
                          'gencode.v27_fantom_selection.txt'),
                   stringsAsFactors=F, row.names=3,
                   col.names=c('seqnames', 'tss', 'transcript_id', 'strand',
                               'gene_id', 'max_fantom', 'tissues_expressed'))

P_domain = read.table(paste0(basedir, '/data/lad_promoters/domains/',
                             'gencode.v27_domains.txt'),
                      stringsAsFactors=F, header=T, row.names=1)

P = data.frame(P_tss, P_exp[rownames(P_tss), ], P_domain[rownames(P_tss), ])
P$K562_GROcap = ifelse(P$strand=='+',
                       P$GROcap_plus,
                       P$GROcap_min)
P$K562.B1_sense = ifelse(P$strand=='+',
                         P$K562.B1_plus,
                         P$K562.B1_min)
P$K562.B2_sense = ifelse(P$strand=='+',
                         P$K562.B2_plus,
                         P$K562.B2_min)
P$HT1080.B1_sense = ifelse(P$strand=='+',
                           P$HT1080.B1_plus,
                           P$HT1080.B1_min)
P$HT1080.B2_sense = ifelse(P$strand=='+',
                           P$HT1080.B2_plus,
                           P$HT1080.B2_min)

pro_tssr_match = match(rownames(P), P_pro_tssr$transcript_id)
P$K562_PROseq_tssr = ifelse(P$strand=='+',
                            rowMeans(P_pro_tssr[pro_tssr_match,c("rep1_plus",
                                                                 "rep2_plus")]),
                            rowMeans(P_pro_tssr[pro_tssr_match,c("rep1_min",
                                                                 "rep2_min")]))
P$K562_PROseq_tssr[is.na(P$K562_PROseq_tssr)] = 0
pro_body_match = match(rownames(P), P_pro_body$transcript_id)
P$K562_PROseq_body = ifelse(P$strand=='+',
                           rowMeans(P_pro_body[pro_body_match,c("rep1_plus",
                                                            "rep2_plus")]),
                          rowMeans(P_pro_body[pro_body_match,c("rep1_min",
                                                           "rep2_min")]))
P$K562_PROseq_body[is.na(P$K562_PROseq_body)] = 0

P$K562_SuRE = rowMeans(P[,c('K562.B1_sense', 'K562.B2_sense')])
P$HT1080_SuRE = rowMeans(P[,c('HT1080.B1_sense', 'HT1080.B2_sense')])

sd_jit = min(P$K562_GROcap[P$K562_GROcap>0])
jit = rnorm(nrow(P), sd = sd_jit / 20)

P$K562_GROcap_jitter = log10(P$K562_GROcap + jit + sd_jit / 2)

pseudo_log10 <- function(val_vec){
    Pseud=min(val_vec[val_vec > 0], na.rm=TRUE)/2
    val_vec = val_vec + Pseud
    return(log10(val_vec))
}

for (col in c('K562_SuRE', 'HT1080_SuRE', 'K562_GROcap', 'K562_PROseq_tssr',
              'K562_PROseq_body')){
    P[,col] = pseudo_log10(P[,col])
}

create_RM <-function(data, x, y, ad){
    #then calculate running mean for iLAD promoters:
    #sort by SuRE and then random for ties
    o = order(data[,x],sample(c(1:nrow(data))))

    x_sorted = data[o,x]
    y_sorted = data[o,y]
    ad_sorted = data[o,ad]

    n<-60 #number of windows
    w<-501 #window width (number of datapoints); if n*w > nrow(P) then windows overlap
    s<-round(seq(from=w/2+0.0001, to=nrow(data)-w/2, length.out=n))
    RM<-data.frame(x.low=rep(NA,n), x.mean=rep(NA,n), x.hi=rep(NA,n), y.AD=rep(NA,n), y.iAD=rep(NA,n))
    RM$x.low=x_sorted[s-floor(w/2)]
    for(i in 1:n){RM$x.mean[i]=mean(x_sorted[(s[i]-floor(w/2)):(s[i]+floor(w/2))], na.rm=TRUE)}
    RM$x.hi=x_sorted[s+floor(w/2)]
    for(i in 1:n)
      {t<-data.frame(AD=ad_sorted[(s[i]-floor(w/2)):(s[i]+floor(w/2))],
                     y=y_sorted[(s[i]-floor(w/2)):(s[i]+floor(w/2))])
       RM$y.AD[i]<-mean(t$y[t$AD==1], na.rm=TRUE)
       RM$y.iAD[i]<-mean(t$y[t$AD==0], na.rm=TRUE)
      }
    #add first datapoint (SuRE equals pseudocount)
    RM1<-RM[0,] #empty df
    RM1[1,]<-c(rep(min(x_sorted),3), mean(y_sorted[x_sorted==min(x_sorted) & ad_sorted==1]), mean(y_sorted[x_sorted==min(x_sorted) & ad_sorted==0]))
    RM<-rbind(RM1, RM)
    rm(RM1)
    return(RM)
}

RM_LAD = create_RM(P, 'K562_SuRE', 'K562_GROcap', ad='K562_LAD')


P$LRS_LAD<- P$K562_GROcap - approx(x=RM_LAD$x.mean, y=RM_LAD$y.iAD,
                                   xout=P$K562_SuRE, rule=2)$y

classify <- function(sure, exp, lrs, ad, exp_cut){
    INACT<- sure < -0.3 & ad & exp < exp_cut #inactive
    NREP<- sure > 0 & lrs > -0.5 & ad & exp > exp_cut #not repressed
    REP<- sure > 0.3 & lrs < -1 & ad  & exp < exp_cut #repressed
    Pcnts<-c(length(which(REP)), length(which(NREP)), length(which(INACT)))
    names(Pcnts)<-c("repressed", "escaper", "inactive")
    BND <- ad & !INACT & !NREP & !REP
    class = rep(NA, length(sure))
    class[ad==0] = 'iLAD'
    class[INACT]<-"inactive"
    class[NREP]<-"escaper"
    class[REP]<-"repressed"
    class[BND] <- "boundary"
    return(factor(class, levels=c('iLAD', 'escaper', 'repressed', 'inactive', 'boundary')))
}

P$class_LAD = classify(P$K562_SuRE, P$K562_GROcap, P$LRS_LAD, P$K562_LAD, -2)

lad_names = c(LAD=paste0('LAD; n=', table(P$K562_LAD)['1']),
              iLAD=paste0('iLAD; n=', table(P$K562_LAD)['0']))
P$K562_LAD_n = factor(ifelse(P$K562_LAD==1, lad_names['LAD'], lad_names['iLAD']))
COL_lad_n = COL_lad
names(COL_lad_n) = lad_names


mask_regions = fread(paste0(basedir, '/data/lad_promoters/window/',
                            'gencode.v27_mask.txt'),
                     header=F, col.names=c('gene_id', 'strand', 'start', 'end'),
                     stringsAsFactors=F)
mask_regions[,transcript_id := eval(rownames(P))[match(gene_id, eval(P$gene_id))]]
setkey(mask_regions, 'transcript_id')

P$dist_to_nearest = abs(mask_regions[rownames(P), 'start'])
lad_class_vec = c('inactive', 'repressed', 'escaper')

p_bed = import.bed(paste0(basedir, '/data/lad_promoters/selection/',
                          'gencode.v27_fantom_selection.bed'))

width = width(p_bed[match(rownames(P), p_bed$name)])

P_big = P[width>20000,]
P_far = P_big[which(P_big$class_LAD %in% lad_class_vec |
                    P_big$dist_to_nearest > 20000), ]

p_matched = matchSet(P_far[P_far$class_LAD%in%c('iLAD', 'escaper'), ], 'class_LAD',
                     'escaper', 'K562_GROcap', 20)

p_matched_SuRE = matchSet(P_far[P_far$class_LAD%in%c('iLAD', 'escaper'), ], 'class_LAD',
                          'escaper', 'K562_SuRE')

p_sort = data.frame(p_bed[match(rownames(P), p_bed$name)])
escaper_data = data.frame(p_sort[, c('seqnames', 'start', 'end', 'strand')],
                          P[, c('gene_id', 'class_LAD', 'K562_SuRE',
                                'K562_GROcap', 'LRS_LAD')])


p_subset = rbind(p_matched,
                 P_big[P_big$class_LAD%in%c('repressed', 'inactive'), ])

save(p_subset, 'cl20181127_promoter_matched_sets.RData')

figS2D = ggplot(p_matched, aes(x=K562_GROcap, color=class_LAD)) +
            geom_density(adjust=0.8) +
            xlab('log10(GRO-cap)') +
            scale_color_manual(values=COL_class) +
            ggtitle('matched set GRO-cap') +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(size = .5, color = "grey"),
                  legend.background = element_rect(fill="transparent"),
                  text = element_text(size=14),
                  legend.justification =c(1,1),
                  legend.position=c(0.05,0.95),
                  legend.title=element_blank())

figS2F = ggplot(p_matched_SuRE, aes(x=K562_SuRE, color=class_LAD)) +
            geom_density(adjust=0.8) +
            xlab('log10(SuRE)') +
            scale_color_manual(values=COL_class) +
            ggtitle('matched set SuRE') +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(size = .5, color = "grey"),
                  legend.background = element_rect(fill="transparent"),
                  text = element_text(size=14),
                  legend.justification =c(1,1),
                  legend.position=c(0.05,0.95),
                  legend.title=element_blank())


## DAM-ID window based on bigwig track containing counts centred on each GATC-fragment.
damid_vec = list.files(paste0(basedir, '/data/lad_promoters/window/'),
                       pattern='centre.r.*22000_100.txt', full.names=T)


runsum_bins <- function(dt, k, step){
    runsum_list = lapply(seq(1,ncol(dt),step), function(i){
        if (i - k/2 < 1){
            start = 1
            end = k
        } else if (i + k/2 > ncol(dt)){
            start = ncol(dt) - k + 1
            end = ncol(dt)
        } else {
            start = i - k/2
            end = i + k/2 - 1
        }
        rowSums(dt[, start:end])
    })
    do.call(cbind, runsum_list)
}




k=20
step=2
damid_list = lapply(damid_vec, function(file_name, P, k, step){
    input_dt = fread(paste('zcat', file_name), skip=1,
                     stringsAsFactors=F, sep='\t',
                     key='V4')
    runsum_dt = runsum_bins(input_dt[rownames(P),7:ncol(input_dt)], k, step)
    nuc_step = 100 * step
    pos_vec = ((-22000 / nuc_step) : (22000 / nuc_step - 1)) * nuc_step
    colnames(runsum_dt) = pos_vec
    return(runsum_dt)
}, p_subset, k, step)

names(damid_list) = gsub('.*v27_(.*).centre.(.*?)_.*', 'K562_\\2_\\1', damid_vec)
names(damid_list) = gsub('[.]', '_', names(damid_list))


damid_stats = read.table(paste0(basedir, '/data/damid_statistics.txt'),
                         header=T, row.names=1)

lmnb1_vec = grep('LMNB1', names(damid_list), value=T)
lmnb1_vec = lmnb1_vec[order(gsub(".*_r([0-9])_.*", '\\1', lmnb1_vec))]
dam_vec = grep('Dam', names(damid_list), value=T)
dam_vec = dam_vec[order(gsub(".*_r([0-9])_.*", '\\1', dam_vec))]


damid_norm_list = lapply(1:2, function(i){
    lmnb1_name = lmnb1_vec[i]
    dam_name = dam_vec[i]
    lmnb1 = damid_list[[lmnb1_name]]
    dam = damid_list[[dam_name]]
    n_lmnb1 = damid_stats[lmnb1_name,'used']
    n_dam = damid_stats[dam_name,'used']
    n_min = min(n_lmnb1, n_dam)

    lmnb1_norm = lmnb1 / n_lmnb1 * n_min + 1
    dam_norm = dam / n_dam * n_min + 1
    lmnb1_norm / dam_norm
})

damid_dt = data.table(class = p_subset$class,
                      log2(do.call('+', damid_norm_list) / 2))

mean_dt = damid_dt[, lapply(.SD, median),by=class]
mean_melt = melt(mean_dt,variable.name='pos', id.vars='class',
                 value.name='log2')
mean_melt[,pos:=as.numeric(as.character(pos))]

w = k * 100

window=data.frame(pos=c(-20000, -20000 + w), log2=c(0,0))
title = paste('Lamin-B1 Dam-ID k=', k, 'step=',step)
fig2A = ggplot(mean_melt, aes(x=pos, y=log2, color=class)) +
        geom_line(size=1.2) +
        geom_line(data=window, color='black', size=1.5) +
        ylab('Lamin-B1/Dam') +
        scale_color_manual(values=COL_class) +
        theme_bw() +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5, face='bold'),
              legend.justification =c(0,1),
              legend.position=c(0.05,0.95),
              legend.background = element_rect(fill="transparent"),
              legend.title=element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x=element_blank()) +
        geom_vline(xintercept=0, linetype='dotdash') +
        coord_cartesian(xlim=c(-20000,20000))


group_mean <- function(signal_data, P, class_name='class_LAD',
                       start=-22000, end=22000, step=200){
    pos_vec = (start / step) : (end / step - 1) + .5
    data = data.frame(class=P[, class_name], signal_data[rownames(P), ])
    mean_data = ddply(data, .(class), function(x){
         apply(x[,2:ncol(x)],2,function(y){
             if (!all(is.na(y)) & !all(y[!is.na(y)]==y[!is.na(y)][1])){
                 t = t.test(y)
                 r = as.vector(c(t$estimate, t$conf.int))
             } else{
                 r = rep(NaN, 3)
             }
             return(r)
         })
    })
    class_vec = unique(mean_data$class)
    i_vec = 1:length(class_vec) * 3
    result = data.frame(class=rep(class_vec, ncol(signal_data)),
                        mean=do.call(c, mean_data[i_vec-2, -1]),
                        min=do.call(c, mean_data[i_vec-1, -1]),
                        max=do.call(c, mean_data[i_vec, -1]),
                        pos=rep(pos_vec * step, each=length(class_vec)))
    return(result)
}

group_median <- function(signal_data, P, class_name='class_LAD',
                       start=-22000, end=22000, step=200){
    pos_vec = (start / step) : (end / step - 1) + .5
    data = data.frame(class=P[, class_name], signal_data[rownames(P), ])
    median_data = ddply(data, .(class), function(x){
         m = t(apply(x[,2:ncol(x)],2,function(y){
             return(quantile(y, probs=c(0.1,0.25,0.5,0.75,0.9), na.rm=T))
            }))
         m = data.frame(m)
         colnames(m) = c('m.1', 'm.25', 'median', 'm.75', 'm.9')
         m$pos = pos_vec * step
         return(m)
    })
    return(median_data)
}


file_vec = list.files('../../chip_snake/test_coverage/window/', full.names=T,
                      pattern='gencode.*22000_100.txt')
# file_vec = grep('H3K9|H3K36', file_vec, value=T)
base_vec = gsub('.*/(.*).txt.gz', '\\1', file_vec)


file_matrix = do.call(rbind, strsplit(base_vec, '_'))
file_df = data.frame(file_matrix, stringsAsFactors=F)

colnames(file_df) = c('region_source', 'origin', 'target', 'rep', 'ID',
                      'refpoint', 'a', 'b', 'binsize')
file_df$file_name = file_vec

salzberg_file = file_df$origin == 'Salzberg2017' & file_df$refpoint == 'TSS'
schmidl_H3K36_file = file_df$origin == 'Schmidl2015' & file_df$target == 'H3K36me3' &
                file_df$ID == 'chipmentation'
pol2as2_file = file_df$target == 'POL2AS2'

k=20
step=2

mean_table = ddply(file_df[salzberg_file | schmidl_H3K36_file | pol2as2_file, ],
                   .(region_source, origin, target, ID, refpoint),
      function(x, P){
          binsize = as.numeric(x$binsize[1])
          start = -as.numeric(x$a[1])
          end = as.numeric(x$b[1])


          runsum_list = mclapply(x$file_name, function(file_name, P){
                          input_dt = fread(paste('zcat', file_name), skip=1,
                                           stringsAsFactors=F, sep='\t',
                                           key='V4')
                          runsum_dt = runsum_bins(input_dt[rownames(P),7:ncol(input_dt)],
                                                  k, step)

                          nuc_step = 100 * step
                          pos_vec = ((-22000 / nuc_step) : (22000 / nuc_step - 1)) * nuc_step
                          colnames(runsum_dt) = pos_vec
                          return(runsum_dt)
                      }, P, mc.cores=3)

         ctrl = grep('control', x$file_name)


         idx_stats = read.table(paste0('../../chip_snake/hg38/idxstats/K562/',
                                       x$origin[1], '/', x$target[1], '_',
                                       x$ID[1], '.txt'), header=T,
                                stringsAsFactors=F)

         n_mapped = colSums(idx_stats[-grep('chrM|[*]', idx_stats$seqnames), -(1:2) ])
         i_input = grep('input', names(n_mapped))
         n_exp = n_mapped[-i_input]
         n_input = n_mapped[i_input]

         ctrl_dt = runsum_list[[ctrl]]
         norm_list = lapply(1:length(n_exp), function(i){
             exp_df = runsum_list[-ctrl][[i]]
             n_min = min(n_exp[i], n_input)
             exp_norm = exp_df / n_exp[i] * n_min * 100 + 1
             ctrl_norm = ctrl_dt / n_input * n_min * 100 + 1
             log2(exp_norm / ctrl_norm)
         })
         chip_dt = data.table(class = P$class_LAD,
                               do.call('+', norm_list) / length(norm_list))
         mean_dt = chip_dt[, lapply(.SD, mean),by=class]

         mean_melt = melt(mean_dt,variable.name='pos', id.vars='class',
                          value.name='log2')
         mean_melt[,pos:=as.numeric(as.character(pos))]
         return(mean_melt)
     }, p_subset)
mean_table$class = as.character(mean_table$class)
u = unique(mean_table[,c('origin', 'target', 'refpoint', 'ID')])


pdf('cl20181109_H3K36me2_POL2AS2_TES.pdf')
for (target in c('H3K36me3', 'POL2AS2')) {
    df = mean_table[which(mean_table$target==target &
                          mean_table$refpoint=='TES'), ]
    print(ggplot(df, aes(x=pos, y=log2)) +
                geom_line(aes(color=class)) +
                scale_color_manual(values=COL_class) +
                scale_fill_manual(values=COL_class) +
                theme_bw() +
                ggtitle(paste(target, 'ChIP TES')) +
                theme(plot.title = element_text(hjust = 0.5, face='bold'),
                      legend.justification =c(0,1),
                      legend.position=c(0.05,0.95),
                      legend.background = element_rect(fill="transparent"),
                      legend.title=element_blank(),
                      panel.grid.minor = element_blank(),
                            axis.title.x=element_blank()) +
                ylab('log2(ChIP/Input)') +
                geom_vline(xintercept=0, linetype='dotdash') +
                geom_vline(xintercept=3000, linetype='dotdash') +
                coord_cartesian(xlim=c(-20000,20000)))
}
dev.off()

salzberg = u[which(u$origin=='Salzberg2017'), ]


i = which(salzberg$target=='H3K9me2' & salzberg$refpoint=='TSS')
df = mean_table[which(mean_table$origin == salzberg[i,'origin'] &
                      mean_table$target == salzberg[i,'target'] &
                      mean_table$refpoint == salzberg[i,'refpoint'] &
                      mean_table$ID == salzberg[i,'ID']),]

figS2B = ggplot(df, aes(x=pos, y=log2)) +
            geom_line(aes(color=class)) +
            scale_color_manual(values=COL_class) +
            scale_fill_manual(values=COL_class) +
            theme_bw() +
            ggtitle('H3K9me2 ChIP') +
            theme(plot.title = element_text(hjust = 0.5, face='bold'),
                  legend.justification =c(0,1),
                  legend.position=c(0.05,0.95),
                  legend.background = element_rect(fill="transparent"),
                  legend.title=element_blank(),
                  panel.grid.minor = element_blank(),
                        axis.title.x=element_blank()) +
            ylab('log2(ChIP/Input)') +
            geom_vline(xintercept=0, linetype='dotdash') +
            coord_cartesian(xlim=c(-20000,20000))


i = which(salzberg$target=='H3K9me3' & salzberg$refpoint=='TSS')
df = mean_table[which(mean_table$origin == salzberg[i,'origin'] &
                      mean_table$target == salzberg[i,'target'] &
                      mean_table$refpoint == salzberg[i,'refpoint'] &
                      mean_table$ID == salzberg[i,'ID']),]
ggplot(df, aes(x=pos, y=log2)) +
          geom_line(aes(color=class)) +
          scale_color_manual(values=COL_class) +
          scale_fill_manual(values=COL_class) +
          theme_bw() +
          ggtitle('H3K9me3 ChIP') +
          theme(plot.title = element_text(hjust = 0.5, face='bold'),
                legend.justification =c(0,1),
                legend.position=c(0.05,0.95),
                legend.background = element_rect(fill="transparent"),
                legend.title=element_blank(),
                panel.grid.minor = element_blank(),
                      axis.title.x=element_blank()) +
          ylab('log2(ChIP/Input)') +
          geom_vline(xintercept=0, linetype='dotdash') +
          coord_cartesian(xlim=c(-20000,20000))

ggplot(df, aes(x=pos, y=log2(control))) +
    geom_line(aes(color=class)) +
    scale_color_manual(values=COL_class) +
    scale_fill_manual(values=COL_class) +
    theme_bw() +
    ggtitle('H3K9me3 ChIP') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
        legend.justification =c(0,1),
        legend.position=c(0.05,0.95),
        legend.background = element_rect(fill="transparent"),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
              axis.title.x=element_blank()) +
    ylab('log2(ChIP/Input)') +
    geom_vline(xintercept=0, linetype='dotdash') +
    coord_cartesian(xlim=c(-20000,20000))
figS2C = ggplot(df, aes(x=pos, y=log2)) +
            geom_line(aes(color=class)) +
            scale_color_manual(values=COL_class) +
            scale_fill_manual(values=COL_class) +
            theme_bw() +
            ggtitle('H3K9me3 ChIP') +
            theme(plot.title = element_text(hjust = 0.5, face='bold'),
                  legend.justification =c(0,1),
                  legend.position=c(0.05,0.95),
                  legend.background = element_rect(fill="transparent"),
                  legend.title=element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.title.x=element_blank()) +
            ylab('log2(ChIP/Input)') +
            geom_vline(xintercept=0, linetype='dotdash') +
            coord_cartesian(xlim=c(-20000,20000))



schmidl_H3K36 = u[which(u$origin=='Schmidl2015' & u$target=='H3K36me3'), ]


i = which(schmidl_H3K36$ID=='chipmentation' & schmidl_H3K36$refpoint=='TSS')
df = mean_table[which(mean_table$origin == schmidl_H3K36[i,'origin'] &
                      mean_table$target == schmidl_H3K36[i,'target'] &
                      mean_table$refpoint == schmidl_H3K36[i,'refpoint'] &
                      mean_table$ID == schmidl_H3K36[i,'ID']),]

fig2E = ggplot(df, aes(x=pos, y=log2)) +
            geom_line(aes(color=class)) +
            scale_color_manual(values=COL_class) +
            scale_fill_manual(values=COL_class) +
            theme_bw() +
            ggtitle('H3K36me3 ChIP') +
            theme(plot.title = element_text(hjust = 0.5, face='bold'),
                  legend.justification =c(0,1),
                  legend.position=c(0.05,0.95),
                  legend.background = element_rect(fill="transparent"),
                  legend.title=element_blank(),
                  panel.grid.minor = element_blank(),
                        axis.title.x=element_blank()) +
            ylab('log2(ChIP/Input)') +
            geom_vline(xintercept=0, linetype='dotdash') +
            coord_cartesian(xlim=c(-10000,20000))


pos_vec = ((-22000 / 100) : (22000 / 100 - 1) + .5) * 100


file_vec = list.files('../cl20180517_hg38_snakemake/window/',
                      pattern='TTseq.*_100.txt',
                      full.names=T)
file_df = cbind(grep('plus', file_vec, value=T), grep('min', file_vec, value=T))
mean_list = apply(file_df, 1, function(file_vec, P){
                input_plus = fread(paste('zcat', file_vec[1]), skip=1,
                                 stringsAsFactors=F, sep='\t',
                                 key='V4')
                input_min = fread(paste('zcat', file_vec[2]), skip=1,
                                 stringsAsFactors=F, sep='\t',
                                 key='V4')
                input_sense = rbind(input_plus[input_plus$V6=='+', ],
                                    input_min[input_min$V6=='-', ])
                input_sense$class = P[input_sense$V4, 'class_LAD']
                setkey(input_sense, 'V4')
                means = ddply(input_sense[rownames(P), ], .(class), function(x){
                               means= Rle(colMeans(x[,7:(ncol(x)-1)], na.rm=T))
                               as.numeric(runmean(means, k=3, endrule='constant'))
                           })
                colnames(means)[-1] = pos_vec
                mean_melt = melt(means,variable.name='pos', id.vars='class')
                mean_melt$pos = as.numeric(as.character(mean_melt$pos))
                return(mean_melt)
            }, p_subset)

means = data.frame(mean_list[[1]][, 1:2],
                 mean = rowMeans(do.call(cbind, lapply(mean_list,
                                                       function(x){x[,3]}))))
fig2C = ggplot(means, aes(x=pos, y=mean, color=class)) +
         scale_color_manual(values=COL_class) +
         theme_bw() +
         ggtitle('TT-seq') +
         ylab('coverage') +
         theme(plot.title = element_text(hjust = 0.5, face='bold'),
               legend.justification =c(0,1),
               legend.position=c(0.05,0.95),
               legend.background = element_rect(fill="transparent"),
               legend.title=element_blank(),
               panel.grid.minor = element_blank(),
               axis.title.x=element_blank()) +
         geom_vline(xintercept=0, linetype='dotdash') +
         geom_line() +
         coord_cartesian(xlim=c(-10000,20000))

dnase_file = list.files('../cl20180517_hg38_snakemake/window/',
                      pattern='DNAse.*_100.txt',
                      full.names=T)
input_dnase = fread(paste('zcat', dnase_file), skip=1,
                    stringsAsFactors=F, sep='\t',
                    key='V4')
input_dnase[,class:=P[V4,'class_LAD']]

dnase_means = ddply(input_dnase[rownames(p_subset), ], .(class), function(x){
               means= Rle(colMeans(x[,7:(ncol(x)-1)], na.rm=T))
               as.numeric(runmean(means, k=3, endrule='constant'))
           })
colnames(dnase_means)[-1] = pos_vec
dnase_melt = melt(dnase_means,variable.name='pos', id.vars='class')
dnase_melt$pos = as.numeric(as.character(dnase_melt$pos))

figS2A = ggplot(dnase_melt, aes(x=pos, y=value, color=class)) +
         scale_color_manual(values=COL_class) +
         theme_bw() +
         ggtitle('DNAse-seq') +
         ylab('coverage') +
         theme(plot.title = element_text(hjust = 0.5, face='bold'),
               legend.justification =c(0,1),
               legend.position=c(0.05,0.95),
               legend.background = element_rect(fill="transparent"),
               legend.title=element_blank(),
               panel.grid.minor = element_blank(),
               axis.title.x=element_blank()) +
         geom_vline(xintercept=0, linetype='dotdash') +
         geom_line() +
         coord_cartesian(xlim=c(-20000,20000))


tssr_file_vec = list.files('../../chip_snake/cl20181108_lad_repression/means/', full.names=T,
                            pattern='gencode.v27_tssr')
base_vec = gsub('.*/(.*).txt', '\\1', tssr_file_vec)

file_matrix = do.call(rbind, strsplit(base_vec, '_'))
tssr_file_df = data.frame(file_matrix, stringsAsFactors=F)

colnames(tssr_file_df) = c('region_source', 'region', 'chip_source', 'target', 'ID')
tssr_file_df$file_name = tssr_file_vec
tssr_file_df = tssr_file_df[file.info(tssr_file_df$file_name)$size > 1, ]



body_file_vec = list.files('../../chip_snake/cl20181108_lad_repression/means/', full.names=T,
                            pattern='gencode.v27_body')
base_vec = gsub('.*/(.*).txt', '\\1', body_file_vec)

file_matrix = do.call(rbind, strsplit(base_vec, '_'))
body_file_df = data.frame(file_matrix, stringsAsFactors=F)

colnames(body_file_df) = c('region_source', 'region', 'chip_source', 'target', 'ID')
body_file_df$file_name = body_file_vec
body_file_df = body_file_df[file.info(body_file_df$file_name)$size > 1, ]


tssr_pol2 = read.table(tssr_file_df[tssr_file_df$ID=='ENCSR388QZF', 'file_name'],
                       header=T, row.names=1)
body_pol2 = read.table(body_file_df[body_file_df$ID=='ENCSR388QZF', 'file_name'],
                       header=T, row.names=1)

pol2_df = data.frame(class=p_subset$class_LAD,
                     tssr=tssr_pol2[rownames(p_subset), 'mean'],
                     body=body_pol2[rownames(p_subset), 'mean'])

pol2_melt = melt(pol2_df, variable.name='region')


pol2_esc_ilad = pol2_melt[pol2_melt$class%in%c('escaper', 'iLAD'), ]
quant_summary = ddply(pol2_esc_ilad, .(class, region), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$value,probs=c(0.25, 0.75))
                         df = data.frame(region, x$region[1], rbind(q, q))

                         colnames(df) = c('x', 'region', 'down', 'up')
                         return(df)
                     })

median_summary = ddply(pol2_esc_ilad, .(class, region), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$value)
                          df = data.frame(region, x$region[1], rbind(q, q))

                          colnames(df) = c('x', 'region','med')
                          return(df)
                      })
p = ggplot(pol2_esc_ilad, aes(x=class, y=log2(value), color=class)) +
    geom_quasirandom() +
    theme_bw() +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    scale_color_manual(values=COL_class) +
    facet_wrap(~region) +
    ggtitle('POL2 ChIP occupancy') +
    ylab('log2(ChIP/input)') +
    guides(color=F) +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())

for (class in unique(median_summary$class)){
    for (region in unique(median_summary$region)){
        quant_df = quant_summary[quant_summary$class==class &
                                 quant_summary$region==region, ]
        median_df = median_summary[median_summary$class==class &
                                   median_summary$region==region, ]
        p = p + geom_line(data=quant_df, aes(x=x, y=log2(down)), size=1, color='black') +
                 geom_line(data=median_df, aes(x=x, y=log2(med)), size=1, color='black') +
                 geom_line(data=quant_df, aes(x=x, y=log2(up)), size=1, color='black')
    }
}

fig2D = p


p_proseq = p_subset[p_subset$class_LAD %in% c('escaper', 'iLAD'),
                    c('class_LAD', 'K562_PROseq_tssr', 'K562_PROseq_body')]

colnames(p_proseq) = c('class', 'tssr', 'body')
p_pro_melt = melt(p_proseq, variable.name='region')

quant_summary = ddply(p_pro_melt, .(class, region), function(x){
                         prom_adj = which(levels(x$class)==x$class[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$value,probs=c(0.25, 0.75))
                         df = data.frame(region, x$region[1], rbind(q, q))

                         colnames(df) = c('x', 'region', 'down', 'up')
                         return(df)
                     })

median_summary = ddply(p_pro_melt, .(class, region), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$value)
                          df = data.frame(region, x$region[1], rbind(q, q))

                          colnames(df) = c('x', 'region','med')
                          return(df)
                      })
figS2G = ggplot(p_pro_melt, aes(x=class, y=value, color=class)) +
    geom_quasirandom() +
    theme_bw() +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    scale_color_manual(values=COL_class) +
    facet_wrap(~region) +
    ggtitle('PROseq') +
    ylab('log10(PROseq)') +
    guides(color=F) +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())

for (class in unique(median_summary$class)){
    for (region in unique(median_summary$region)){
        quant_df = quant_summary[quant_summary$class==class &
                                 quant_summary$region==region, ]
        median_df = median_summary[median_summary$class==class &
                                   median_summary$region==region, ]
        figS2G = figS2G + geom_line(data=quant_df, aes(x=x, y=down), size=1, color='black') +
                     geom_line(data=median_df, aes(x=x, y=med), size=1, color='black') +
                     geom_line(data=quant_df, aes(x=x, y=up), size=1, color='black')
    }
}




body_H3K36_file = body_file_df[which(body_file_df$chip_source=='Schmidl2015' &
                                     body_file_df$target=='H3K36me3' &
                                     body_file_df$ID=='chipmentation'), 'file_name']
body_H3K36 = read.table(body_H3K36_file, header=T, row.names=1)

wilcox_H3K36 = wilcox.test(body_H3K36[rownames(p_subset)[p_subset$class=='escaper'], 'mean'],
                           body_H3K36[rownames(p_subset)[p_subset$class=='iLAD'], 'mean'])


fig2E = fig2E + geom_text(data=data.frame(pos=10000,log2=0),color='black',
                  label=sprintf('p = %0.2g', wilcox_H3K36$p.value))



rep1 = fread('../raw_data/K562_rna_rep1_ENCFF004LGY.tsv')
rep1$gene_base = gsub('[.][0-9]+', '', rep1$gene_id)
setkey(rep1, 'gene_base')
rep2 = fread('../raw_data/K562_rna_rep2_ENCFF222NCB.tsv')
rep2$gene_base = gsub('[.][0-9]+', '', rep2$gene_id)
setkey(rep2, 'gene_base')
subset_id = gsub('[.][0-9]+', '', p_subset$gene_id)
pme_fpkm = cbind(rep1[subset_id, pme_FPKM], rep2[subset_id, pme_FPKM])
p_subset$RNA_seq = rowMeans(pme_fpkm)
p_esc_ilad = p_subset[p_subset$class_LAD%in%c('iLAD', 'escaper'), ]
p_esc_ilad$class_LAD = factor(p_esc_ilad$class_LAD)
quant_summary = ddply(data.table(p_esc_ilad), .(class_LAD), function(x){
                        print(x$class_LAD[1])
                          prom_adj = which(levels(x$class_LAD)==x$class_LAD[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$RNA_seq,probs=c(0.2, 0.75), na.rm=T)
                         df = data.frame(region, rbind(q, q))
                         colnames(df) = c('x', 'down', 'up')
                         return(df)
                     })
median_summary = ddply(data.table(p_esc_ilad), .(class_LAD), function(x){
                          prom_adj = which(levels(x$class_LAD)==x$class_LAD[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$RNA_seq, na.rm=T)
                          df = data.frame(region, rbind(q, q))
                          colnames(df) = c('x','med')
                          return(df)
                      })
mean_summary = aggregate(RNA_seq ~ class_LAD, mean, data=p_esc_ilad)

fc = data.frame(fc=paste0(round(median_summary[1,'med']/median_summary[3,'med']), '-fold'),
                class_LAD='iLAD',
                RNA_seq=10^3)

line_df = data.frame(class_LAD=c(1,1,2,2),
                     RNA_seq=c(600, 800, 800, 600))

p = ggplot(p_subset[p_subset$class_LAD%in%c('escaper', 'iLAD'), ],
           aes(x=class_LAD, y=log10(RNA_seq), color=class_LAD)) +
        geom_quasirandom() +
        stat_compare_means(method='wilcox.test',
                           comparisons=list(c('escaper', 'iLAD'))) +
        # geom_text(data=fc, aes(label=fc), color='black', nudge_x=0.5) +
        # geom_line(data=line_df, color='black') +
        scale_color_manual(values=COL_class) +
        theme_bw() +
        ggtitle('polyA RNA-seq') +
        guides(color=F) +
        ylab('log10(FPKM)') +
        theme(plot.title = element_text(hjust = 0.5, face='bold'),
              panel.grid.minor = element_blank(),
              axis.title.x=element_blank())
for (class in unique(median_summary$class)){
    quant_df = quant_summary[quant_summary$class==class, ]
    median_df = median_summary[median_summary$class==class, ]
    p = p + geom_line(data=quant_df, aes(x=x, y=log10(down)), size=1, color='black') +
             geom_line(data=median_df, aes(x=x, y=log10(med)), size=1, color='black') +
             geom_line(data=quant_df, aes(x=x, y=log10(up)), size=1, color='black')
}
fig2B = p

subset_id = gsub('[.][0-9]+', '', p_matched_SuRE$gene_id)
pme_fpkm = cbind(rep1[subset_id, pme_FPKM], rep2[subset_id, pme_FPKM])
p_matched_SuRE$RNA_seq = rowMeans(pme_fpkm)
p_matched_SuRE$class_LAD = factor(p_matched_SuRE$class_LAD)
quant_summary = ddply(data.table(p_matched_SuRE), .(class_LAD), function(x){
                        print(x$class_LAD[1])
                          prom_adj = which(levels(x$class_LAD)==x$class_LAD[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$RNA_seq,probs=c(0.2, 0.75), na.rm=T)
                         df = data.frame(region, rbind(q, q))
                         colnames(df) = c('x', 'down', 'up')
                         return(df)
                     })
median_summary = ddply(data.table(p_matched_SuRE), .(class_LAD), function(x){
                          prom_adj = which(levels(x$class_LAD)==x$class_LAD[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$RNA_seq, na.rm=T)
                          df = data.frame(region, rbind(q, q))
                          colnames(df) = c('x','med')
                          return(df)
                      })
mean_summary = aggregate(RNA_seq ~ class_LAD, mean, data=p_esc_ilad)

fc = data.frame(fc=paste0(round(median_summary[1,'med']/median_summary[3,'med']), '-fold'),
                class_LAD='iLAD',
                RNA_seq=10^3)

line_df = data.frame(class_LAD=c(1,1,2,2),
                     RNA_seq=c(600, 800, 800, 600))

p = ggplot(p_matched_SuRE, aes(x=class_LAD, y=log10(RNA_seq), color=class_LAD)) +
        geom_quasirandom() +
        stat_compare_means(method='wilcox.test',
                           comparisons=list(c('escaper', 'iLAD'))) +
        # geom_text(data=fc, aes(label=fc), color='black', nudge_x=0.5) +
        # geom_line(data=line_df, color='black') +
        scale_color_manual(values=COL_class) +
        theme_bw() +
        ggtitle('polyA RNA-seq') +
        guides(color=F) +
        ylab('log10(FPKM)') +
        theme(plot.title = element_text(hjust = 0.5, face='bold'),
              panel.grid.minor = element_blank(),
              axis.title.x=element_blank())
for (class in unique(median_summary$class)){
    quant_df = quant_summary[quant_summary$class==class, ]
    median_df = median_summary[median_summary$class==class, ]
    p = p + geom_line(data=quant_df, aes(x=x, y=log10(down)), size=1, color='black') +
             geom_line(data=median_df, aes(x=x, y=log10(med)), size=1, color='black') +
             geom_line(data=quant_df, aes(x=x, y=log10(up)), size=1, color='black')
}
figS2E = p













myc_df = read.table(paste0('../../chip_snake/test_coverage/means/',
                           'gencode.v27_1kb_Snyder2012_MYC_ENCSR000EGJ.txt'),
                    header=T, row.names=1, stringsAsFactors=F)

myc_data = data.frame(class=p_esc_ilad$class_LAD, myc_df[rownames(p_esc_ilad), ])


quant_summary = ddply(data.table(myc_data), .(class), function(x){
                        print(x$class_LAD[1])
                          prom_adj = which(levels(x$class)==x$class[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$mean,probs=c(0.2, 0.75), na.rm=T)
                         df = data.frame(region, rbind(q, q))
                         colnames(df) = c('x', 'down', 'up')
                         return(df)
                     })
median_summary = ddply(data.table(myc_data), .(class), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$mean, na.rm=T)
                          df = data.frame(region, rbind(q, q))
                          colnames(df) = c('x','med')
                          return(df)
                      })

fig2F = ggplot(myc_data, aes(x=class, y=log2(mean), color=class)) +
    geom_quasirandom() +
    scale_color_manual(values=COL_class) +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    theme_bw() +
    ggtitle('MYC ChIP') +
    guides(color=F) +
    ylab('log2(ChIP/Input)') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())
for (class in unique(median_summary$class)){
    quant_df = quant_summary[quant_summary$class==class, ]
    median_df = median_summary[median_summary$class==class, ]
    fig2F = fig2F + geom_line(data=quant_df, aes(x=x, y=log2(down)), size=1, color='black') +
        geom_line(data=median_df, aes(x=x, y=log2(med)), size=1, color='black') +
        geom_line(data=quant_df, aes(x=x, y=log2(up)), size=1, color='black')
}
fig2F





brd4_df = read.table(paste0('../../chip_snake/test_coverage/means/',
                           'gencode.v27_1kb_Tyler2017_BRD4_GSE88747.txt'),
                    header=T, row.names=1, stringsAsFactors=F)

brd4_data = data.frame(class=p_esc_ilad$class_LAD, brd4_df[rownames(p_esc_ilad), ])


quant_summary = ddply(data.table(brd4_data), .(class), function(x){
                        print(x$class_LAD[1])
                          prom_adj = which(levels(x$class)==x$class[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$mean,probs=c(0.2, 0.75), na.rm=T)
                         df = data.frame(region, rbind(q, q))
                         colnames(df) = c('x', 'down', 'up')
                         return(df)
                     })
median_summary = ddply(data.table(brd4_data), .(class), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$mean, na.rm=T)
                          df = data.frame(region, rbind(q, q))
                          colnames(df) = c('x','med')
                          return(df)
                      })

figS2K = ggplot(brd4_data, aes(x=class, y=log2(mean), color=class)) +
    geom_quasirandom() +
    scale_color_manual(values=COL_class) +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    theme_bw() +
    ggtitle('BRD4 ChIP') +
    guides(color=F) +
    ylab('log2(ChIP/Input)') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())
for (class in unique(median_summary$class)){
    quant_df = quant_summary[quant_summary$class==class, ]
    median_df = median_summary[median_summary$class==class, ]
    figS2K = figS2K + geom_line(data=quant_df, aes(x=x, y=log2(down)), size=1, color='black') +
        geom_line(data=median_df, aes(x=x, y=log2(med)), size=1, color='black') +
        geom_line(data=quant_df, aes(x=x, y=log2(up)), size=1, color='black')
}
figS2K



nelfe_df = read.table(paste0('../../chip_snake/test_coverage/means/',
                           'gencode.v27_1kb_Struhl2011_NELFE_ENCSR000DOF.txt'),
                    header=T, row.names=1, stringsAsFactors=F)

nelfe_data = data.frame(class=p_esc_ilad$class_LAD, nelfe_df[rownames(p_esc_ilad), ])


quant_summary = ddply(data.table(nelfe_data), .(class), function(x){
                        print(x$class_LAD[1])
                          prom_adj = which(levels(x$class)==x$class[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$mean,probs=c(0.2, 0.75), na.rm=T)
                         df = data.frame(region, rbind(q, q))
                         colnames(df) = c('x', 'down', 'up')
                         return(df)
                     })
median_summary = ddply(data.table(nelfe_data), .(class), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$mean, na.rm=T)
                          df = data.frame(region, rbind(q, q))
                          colnames(df) = c('x','med')
                          return(df)
                      })

figS2I = ggplot(nelfe_data, aes(x=class, y=log2(mean), color=class)) +
    geom_quasirandom() +
    scale_color_manual(values=COL_class) +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    theme_bw() +
    ggtitle('NELFE ChIP') +
    guides(color=F) +
    ylab('log2(ChIP/Input)') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())
for (class in unique(median_summary$class)){
    quant_df = quant_summary[quant_summary$class==class, ]
    median_df = median_summary[median_summary$class==class, ]
    figS2I = figS2I + geom_line(data=quant_df, aes(x=x, y=log2(down)), size=1, color='black') +
        geom_line(data=median_df, aes(x=x, y=log2(med)), size=1, color='black') +
        geom_line(data=quant_df, aes(x=x, y=log2(up)), size=1, color='black')
}
figS2I


mllt1_df = read.table(paste0('../../chip_snake/test_coverage/means/',
                           'gencode.v27_1kb_Snyder2016_MLLT1_ENCSR675LRO.txt'),
                    header=T, row.names=1, stringsAsFactors=F)

mllt1_data = data.frame(class=p_esc_ilad$class_LAD, mllt1_df[rownames(p_esc_ilad), ])


quant_summary = ddply(mllt1_data, .(class), function(x){
                        print(x$class[1])
                         prom_adj = which(levels(x$class)==x$class[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         print(region)
                         q = quantile(x$mean,probs=c(0.2, 0.75), na.rm=T)
                         df = data.frame(region, rbind(q, q))
                         print(df)
                         colnames(df) = c('x', 'down', 'up')
                         return(df)
                     })
median_summary = ddply(mllt1_data, .(class), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$mean, na.rm=T)
                          df = data.frame(region, rbind(q, q))
                          colnames(df) = c('x','med')
                          return(df)
                      })

figS2L = ggplot(mllt1_data, aes(x=class, y=log2(mean), color=class)) +
    geom_quasirandom() +
    scale_color_manual(values=COL_class) +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    theme_bw() +
    ggtitle('MLLT1 ChIP') +
    guides(color=F) +
    ylab('log2(ChIP/Input)') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())
for (class in unique(median_summary$class)){
    quant_df = quant_summary[quant_summary$class==class, ]
    median_df = median_summary[median_summary$class==class, ]
    figS2L = figS2L + geom_line(data=quant_df, aes(x=x, y=log2(down)), size=1, color='black') +
        geom_line(data=median_df, aes(x=x, y=log2(med)), size=1, color='black') +
        geom_line(data=quant_df, aes(x=x, y=log2(up)), size=1, color='black')
}
figS2L




larp_df = read.table(paste0('../../chip_snake/test_coverage/means/',
                           'gencode.v27_1kb_Snyder2016_LARP7_ENCSR288MOZ.txt'),
                    header=T, row.names=1, stringsAsFactors=F)

larp_data = data.frame(class=p_esc_ilad$class_LAD, larp_df[rownames(p_esc_ilad), ])


quant_summary = ddply(data.table(larp_data), .(class), function(x){
                        print(x$class_LAD[1])
                          prom_adj = which(levels(x$class)==x$class[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$mean,probs=c(0.2, 0.75), na.rm=T)
                         df = data.frame(region, rbind(q, q))
                         colnames(df) = c('x', 'down', 'up')
                         return(df)
                     })
median_summary = ddply(data.table(larp_data), .(class), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$mean, na.rm=T)
                          df = data.frame(region, rbind(q, q))
                          colnames(df) = c('x','med')
                          return(df)
                      })

figS2J = ggplot(larp_data, aes(x=class, y=log2(mean), color=class)) +
    geom_quasirandom() +
    scale_color_manual(values=COL_class) +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    theme_bw() +
    ggtitle('LARP7 ChIP') +
    guides(color=F) +
    ylab('log2(ChIP/Input)') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())
for (class in unique(median_summary$class)){
    quant_df = quant_summary[quant_summary$class==class, ]
    median_df = median_summary[median_summary$class==class, ]
    figS2J = figS2J + geom_line(data=quant_df, aes(x=x, y=log2(down)), size=1, color='black') +
        geom_line(data=median_df, aes(x=x, y=log2(med)), size=1, color='black') +
        geom_line(data=quant_df, aes(x=x, y=log2(up)), size=1, color='black')
}
figS2J


gtf2f_df = read.table(paste0('../../chip_snake/test_coverage/means/',
                            'gencode.v27_1kb_Snyder2017_GTF2F1_ENCSR189VXS.txt'),
                    header=T, row.names=1, stringsAsFactors=F)

gtf2f_data = data.frame(class=p_esc_ilad$class_LAD, gtf2f_df[rownames(p_esc_ilad), ])


quant_summary = ddply(data.table(gtf2f_data), .(class), function(x){
                        print(x$class_LAD[1])
                          prom_adj = which(levels(x$class)==x$class[1])
                         region = c(prom_adj - 0.1, prom_adj + 0.1)
                         q = quantile(x$mean,probs=c(0.2, 0.75), na.rm=T)
                         df = data.frame(region, rbind(q, q))
                         colnames(df) = c('x', 'down', 'up')
                         return(df)
                     })
median_summary = ddply(data.table(gtf2f_data), .(class), function(x){
                          prom_adj = which(levels(x$class)==x$class[1])
                          region = c(prom_adj  - 0.2, prom_adj  + 0.2)
                          q = median(x$mean, na.rm=T)
                          df = data.frame(region, rbind(q, q))
                          colnames(df) = c('x','med')
                          return(df)
                      })

figS2H = ggplot(gtf2f_data, aes(x=class, y=log2(mean), color=class)) +
    geom_quasirandom() +
    scale_color_manual(values=COL_class) +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    theme_bw() +
    ggtitle('GTF2F1 ChIP') +
    guides(color=F) +
    ylab('log2(ChIP/Input)') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())
for (class in unique(median_summary$class)){
    quant_df = quant_summary[quant_summary$class==class, ]
    median_df = median_summary[median_summary$class==class, ]
    figS2H = figS2H + geom_line(data=quant_df, aes(x=x, y=log2(down)), size=1, color='black') +
        geom_line(data=median_df, aes(x=x, y=log2(med)), size=1, color='black') +
        geom_line(data=quant_df, aes(x=x, y=log2(up)), size=1, color='black')
}
figS2H


## one worry is that maybe the greylist approach I use is to stringent for these
## datasets. I see a lot of greylisted areas in the tracks, this could lead to
## underrepresentation of the actual fold-differences

encode_df = read.table('../encode_elongation_test.txt',
                       col.names=c('transcript_id', 'GTF2F1', 'MLLT1'),
                        stringsAsFactors=F)

t_match = match(rownames(p_esc_ilad), encode_df$transcript_id)
encode_data = data.frame(class=p_esc_ilad$class_LAD,
                         encode_df[t_match, 2:3])

pdf('elongation_test.pdf', useDingbats=F)

ggplot(encode_data, aes(x=class, y=GTF2F1, color=class)) +
    geom_quasirandom() +
    scale_color_manual(values=COL_class) +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    theme_bw() +
    ggtitle('GTF2F1 ChIP') +
    guides(color=F) +
    ylab('log2(ChIP/Input)') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())

ggplot(encode_data, aes(x=class, y=MLLT1, color=class)) +
    geom_quasirandom() +
    scale_color_manual(values=COL_class) +
    stat_compare_means(method='wilcox.test',
                       comparisons=list(c('escaper', 'iLAD'))) +
    theme_bw() +
    ggtitle('MLLT1 ChIP') +
    guides(color=F) +
    ylab('log2(ChIP/Input)') +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())

dev.off()


pdf('Figure_2_LMNB1_and_elongation.pdf', useDingbats=F, width=8, height=6)
plot_grid(fig2A, fig2C, fig2B, fig2D, fig2E, fig2F, nrow=2,
          rel_widths=c(1,1,.6), align='vh', axis='l',
          labels=c('A', 'C', 'B', 'D', 'E', 'F'))

dev.off()


pdf('Figure_S2_H3K9_dip_and_dnase.pdf', useDingbats=F, width=10, height=8)
plot_grid(figS2A, figS2B, figS2C, figS2D, figS2E, figS2F, figS2G, figS2H,
          figS2I, figS2J, figS2K, figS2L, nrow=3, align='v', axis='l', labels='AUTO')

dev.off()
