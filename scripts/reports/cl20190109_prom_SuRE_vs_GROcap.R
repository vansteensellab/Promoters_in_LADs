
library(ggpubr)
library(ggplot2)
library(rtracklayer)
library(plyr)
library(reshape2)
library(ggbeeswarm)
library(knitr)
library(ppcor)
library(cowplot)
library(data.table)
library(ggExtra)
library(gridExtra)
## get a table with matching sets
## table = complete table to take matching sets from
## class_col = column name of class of interest
## class = name of class to match the set on
## order_on = column name to order on
matchSet <- function(table, class_col, class, order_on){
  o_vec = order(table[,order_on])
  o_table = table[o_vec, ]
  setA = which(o_table[,class_col]==class)
  setB = c(setA + 1, setA -1)
  ## check if setB is all within the possible indexes
  setB = setB[setB %in% 1:length(o_vec)]
  ## can also return o_table[unique(c(setA, setB)), ]
  ## but this way order is perserved.
  i_vec = o_vec[unique(c(setA, setB))]
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


P_pro_body = read.table(paste0(basedir, '/data/lad_promoters/expression/',
                               'gencode.v27_proseq_body.tsv'),
                        stringsAsFactors=F, header=T)

P_exp = read.table(paste0(basedir, '/data/lad_promoters/expression/',
                          'gencode.v27_stranded_expression.txt.gz'),
                       stringsAsFactors=F, header=T, row.names=1)

for (col in colnames(P_exp)){
    P_exp[which(is.na(P_exp[,col])),col] = 0
}
for (col in colnames(P_proseq)){
    P_proseq[which(is.na(P_proseq[,col])),col] = 0
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
P$K562_CAGE_FANTOM = ifelse(P$strand=='+',
                            rowMeans(P[,c("CAGE.FANTOM.rep1_plus",
                                          "CAGE.FANTOM.rep2_plus",
                                          "CAGE.FANTOM.rep3_plus")]),
                            rowMeans(P[,c("CAGE.FANTOM.rep1_min",
                                          "CAGE.FANTOM.rep2_min",
                                          "CAGE.FANTOM.rep3_min")]))
P$K562_PROseq = ifelse(P$strand=='+',
                       rowMeans(P[,c("PROseq.rep1_plus",
                                     "PROseq.rep2_plus")]),
                       rowMeans(P[,c("PROseq.rep1_min",
                                     "PROseq.rep2_min")]))

proseq_match = match(rownames(P), P_proseq$transcript_id)
P$K562_PROseq_200 = ifelse(P$strand=='+',
                           rowMeans(P_proseq[proseq_match,c("rep1_plus",
                                                            "rep2_plus")]),
                          rowMeans(P_proseq[proseq_match,c("rep1_min",
                                                           "rep2_min")]))



P$K562_GROcap = ifelse(P$strand=='+',
                       P$GROcap_plus,
                       P$GROcap_min)
P$K562.B1_sense = ifelse(P$strand=='+',
                         P$K562.B1_plus,
                         P$K562.B1_min)
P$K562.B2_sense = ifelse(P$strand=='+',
                         P$K562.B2_plus,
                         P$K562.B2_min)

P$K562_SuRE = rowMeans(P[,c('K562.B1_sense', 'K562.B2_sense')])

sd_jit = min(P$K562_GROcap[P$K562_GROcap>0])
jit = rnorm(nrow(P), sd = sd_jit / 20)

P$K562_GROcap_jitter = log10(P$K562_GROcap + jit + sd_jit / 2)

sd_jit = min(P$K562_CAGE_FANTOM[P$K562_CAGE_FANTOM>0])
jit = rnorm(nrow(P), sd = sd_jit / 20)
P$K562_CAGE_jitter = log10(P$K562_CAGE_FANTOM + jit + sd_jit / 2)

sd_jit = min(P$K562_PROseq[P$K562_PROseq>0])
jit = rnorm(nrow(P), sd = sd_jit / 20)
P$K562_PROseq_jitter = log10(P$K562_PROseq + jit + sd_jit / 2)
P$K562_PROseq_200_jitter = log10(P$K562_PROseq_200 + jit + sd_jit / 2)


pseudo_log10 <- function(val_vec){
    Pseud=min(val_vec[val_vec > 0], na.rm=TRUE)/2
    val_vec = val_vec + Pseud
    return(log10(val_vec))
}

for (col in c('K562_PROseq', 'K562_CAGE_FANTOM', 'K562_SuRE',
              'K562_PROseq_200', 'K562_GROcap')){
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

RM_CAGE = create_RM(P, 'K562_SuRE', 'K562_CAGE_FANTOM', ad='K562_LAD')
RM_PROseq = create_RM(P, 'K562_SuRE', 'K562_PROseq', ad='K562_LAD')

P$LRS_LAD<- P$K562_GROcap - approx(x=RM_LAD$x.mean, y=RM_LAD$y.iAD,
                                   xout=P$K562_SuRE, rule=2)$y

classify <- function(sure, exp, lrs, ad, exp_cut){
    INACT<- sure< -0.3 & ad & exp< exp_cut #inactive
    NREP<- sure> 0 & lrs > -0.5 & ad & exp> exp_cut #not repressed
    REP<- sure> 0.3 & lrs < -1 & ad  & exp< exp_cut #repressed
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

write.table(P[,c()])

lad_names = c(LAD=paste0('LAD; n=', table(P$K562_LAD)['1']),
              iLAD=paste0('iLAD; n=', table(P$K562_LAD)['0']))
P$K562_LAD_n = factor(ifelse(P$K562_LAD==1, lad_names['LAD'], lad_names['iLAD']))
COL_lad_n = COL_lad
names(COL_lad_n) = lad_names

print('correlation iLADs:')

p_cor = cor(P$K562_SuRE[P$K562_LAD==0], P$K562_GROcap[P$K562_LAD==0], method='pearson')
print(paste('pearson:', round(p_cor, 2)))
sp_cor = cor(P$K562_SuRE[P$K562_LAD==0], P$K562_GROcap[P$K562_LAD==0], method='spearman')
print(paste('spearman:', round(sp_cor, 2)))

print('correlation LADs:')

p_cor = cor(P$K562_SuRE[P$K562_LAD==1], P$K562_GROcap[P$K562_LAD==1], method='pearson')
print(paste('pearson:', round(p_cor, 2)))
sp_cor = cor(P$K562_SuRE[P$K562_LAD==1], P$K562_GROcap[P$K562_LAD==1], method='spearman')
print(paste('spearman:', round(sp_cor, 2)))


P$GROcap_0 = P$K562_GROcap < -3
P$SuRE_0 = P$K562_SuRE > 0
m = matchSet(P[P$SuRE_0 | P$K562_LAD==0, ], 'K562_LAD', 1, 'K562_SuRE')

no_gro = table(m[,c('GROcap_0', 'K562_LAD')])
print(no_gro)
perc = t(t(no_gro) / colSums(no_gro) * 100)
round(perc)

fisher.test(no_gro)$p.value

p_df = data.frame(lad=c('LAD', 'iLAD'), inactive=perc[2,c('1','0')])

figS1A = ggplot(p_df, aes(x=lad, y=inactive, fill=lad)) +
    geom_bar(stat='identity') +
    ylab('% of genes inactive') +
    theme_bw() +
    ylim(0, 100) +
    guides(fill=FALSE) +
    theme(axis.title.x=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = .5, color = "grey"),
          text = element_text(size=14)) +
    geom_text(aes(label = paste(round(inactive), "%")),
              vjust = 5) +
    scale_fill_manual(values=COL_lad)

m$lad = ifelse(m$K562_LAD==1, 'LAD', 'iLAD')
figS1B = ggplot(m, aes(x=K562_SuRE, color=lad)) +
    theme_bw() +
    geom_density(adjust=0.8) +
    xlab('log10(SuRE)') +
    scale_color_manual(values=COL_lad) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = .5, color = "grey"),
          legend.background = element_rect(fill="transparent"),
          text = element_text(size=14),
          legend.justification =c(1,1),
          legend.position=c(0.05,0.95),
          legend.title=element_blank())
pdf('cl20181017_GROcap_percentage_inactive.pdf', useDingbats=F, width=6, height=3)
plot_grid(figS1A, figS1B, labels='AUTO',rel_widths=c(1,1.5), align='h', nrow=1)
dev.off()

RM_melt = melt(RM_LAD, measure.vars=c('y.iAD', 'y.AD'))
RM_melt$variable = ifelse(RM_melt$variable=='y.AD',
                          lad_names['LAD'], lad_names['iLAD'])
ilad_melt = RM_melt[grep('iLAD', RM_melt$variable),]

figA = ggplot(P, aes(x=K562_SuRE, y=K562_GROcap_jitter, color=K562_LAD_n)) +
                  geom_point(data=P[P$K562_LAD==0, ], size=0.25, alpha=0.1) +
                  theme_bw() +
                  geom_line(data=RM_LAD, aes(x=x.mean, y=y.iAD), color='black',  size=1) +
                  geom_line(data=RM_LAD, aes(x=x.mean, y=y.iAD), color=COL_lad['iLAD'], size=0.5) +
                  labs(y='log10(GROcap)', x='log10(SuRE)') +
                  theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_line(size = .5, color = "grey"),
                        axis.line = element_line(size=.7, color = "black"),
						legend.background = element_rect(fill="transparent"),
                        text = element_text(size=14),
						legend.justification =c(0,1),
		 			    legend.position=c(0.05,0.95),
						legend.title=element_blank()) +
                  scale_color_manual(values=COL_lad_n) +
                  coord_equal(ratio=1) +
          		  guides(colour = guide_legend(override.aes = list(size=1,linetype=0, alpha=1)))


figB = ggplot(P, aes(x=K562_SuRE, y=K562_GROcap_jitter, color=K562_LAD_n)) +
                  geom_point(data=P[P$K562_LAD==0, ], size=0.25, alpha=0.1) +
                  geom_point(data=P[P$K562_LAD==1, ], size=0.8, alpha=0.5) +
                  theme_bw() +
                  geom_line(data=RM_LAD, aes(x=x.mean, y=y.AD), color='black',  size=1) +
                  geom_line(data=RM_LAD, aes(x=x.mean, y=y.iAD), color='black',  size=1) +
                  geom_line(data=RM_melt, aes(x=x.mean, y=value, color=variable), size=0.5) +
                  labs(y='log10(GROcap)', x='log10(SuRE)') +
                  theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_line(size = .5, color = "grey"),
                        axis.line = element_line(size=.7, color = "black"),
						legend.background = element_rect(fill="transparent"),
                        text = element_text(size=14),
						legend.justification =c(0,1),
		 			    legend.position=c(0.05,0.95),
						legend.title=element_blank()) +
                  scale_color_manual(values=COL_lad_n) +
                  coord_equal(ratio=1) +
          		  guides(colour = guide_legend(override.aes = list(size=1,linetype=0, alpha=1)))

class_names = paste0(levels(P$class_LAD), '; n=',table(P$class_LAD))
names(class_names) = levels(P$class_LAD)
P$class_LAD_n = P$class_LAD
levels(P$class_LAD_n) = class_names
COL_class_LAD_n = COL_class[names(class_names)]
names(COL_class_LAD_n) = class_names


y_line = approx(x=RM_LAD$x.mean, y=RM_LAD$y.ilad, xout=0.3, rule=2)$y - 1
p_classes = P[which(P$class_LAD %in% c('inactive', 'escaper', 'repressed')),]


figC = ggplot(P[P$K562_LAD==0, ], aes(x=K562_SuRE, y=K562_GROcap_jitter,
                                          color=class_LAD_n)) +
    geom_line(data=RM_LAD[RM_LAD$y.iAD > -1.5,],
              aes(x=x.mean, y=y.iAD - 0.5), color='black',
              linetype='dotdash', size=0.5, show.legend=F) +
    geom_line(data=RM_LAD[RM_LAD$x.mean > 0.3 & RM_LAD$y.iAD < -1,],
           aes(x=x.mean, y=y.iAD - 1), color='black',
           linetype='dotdash', size=0.5, show.legend=F) +
    geom_segment(x=0.3, xend=0.3, y=y_line, yend=min(P$K562_GROcap_jitter),
              linetype='dotdash', color='black', size=0.5) +
    geom_segment(x=-0.3, xend=-0.3, y=-2, yend=min(P$K562_GROcap_jitter),
              linetype='dotdash', color='black', size=0.5) +
    geom_hline(yintercept=-2, linetype='dotdash', size=0.5) +
    geom_point(size=0.1, alpha=0.1) +
    geom_point(data=p_classes, size=0.8) +
    theme_bw() +
    geom_line(data=RM_LAD, aes(x=x.mean, y=y.iAD),
             color='black', size=1) +
    geom_line(data=RM_LAD, aes(x=x.mean, y=y.iAD),
             color=COL_lad['iLAD'], size=0.5) +
    labs(y='log10(GROcap)', x='log10(SuRE)') +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = .5, color = "grey"),
          axis.line = element_line(size=.7, color = "black"),
          legend.background = element_rect(fill="transparent"),
          text = element_text(size=14),
          legend.justification =c(0,1),
          legend.position=c(0.05,0.95),
          legend.title=element_blank()) +
    scale_color_manual(values=COL_class_LAD_n) +
    coord_equal(ratio=1)

RM_CAGE_melt = melt(RM_CAGE, measure.vars=c('y.iAD', 'y.AD'))
RM_CAGE_melt$variable = ifelse(RM_CAGE_melt$variable=='y.AD',
                               lad_names['LAD'], lad_names['iLAD'])
RM_proseq_melt = melt(RM_PROseq, measure.vars=c('y.iAD', 'y.AD'))
RM_proseq_melt$variable = ifelse(RM_proseq_melt$variable=='y.AD',
                                 lad_names['LAD'], lad_names['iLAD'])

p_lad_cage = ggplot(P, aes(x=K562_SuRE, y=K562_CAGE_jitter, color=K562_LAD_n)) +
                  geom_point(data=P[P$K562_LAD==0, ], size=0.25, alpha=0.1) +
                  geom_point(data=P[P$K562_LAD==1, ], size=0.8, alpha=0.5) +
                  theme_bw() +
                  geom_line(data=RM_CAGE, aes(x=x.mean, y=y.AD), color='black',  size=1) +
                  geom_line(data=RM_CAGE, aes(x=x.mean, y=y.iAD), color='black',  size=1) +
                  geom_line(data=RM_CAGE_melt, aes(x=x.mean, y=value, color=variable), size=0.5) +
                  labs(y='log10(CAGE)', x='log10(SuRE)') +
                  theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_line(size = .5, color = "grey"),
                        axis.line = element_line(size=.7, color = "black"),
						legend.background = element_rect(fill="transparent"),
                        text = element_text(size=14),
						legend.justification =c(0,1),
		 			    legend.position=c(0.05,0.95),
						legend.title=element_blank()) +
                  scale_color_manual(values=COL_lad_n) +
                  coord_equal(ratio=1) +
          		  guides(colour = guide_legend(override.aes = list(size=1,linetype=0, alpha=1)))

p_lad_proseq = ggplot(P, aes(x=K562_SuRE, y=K562_PROseq_jitter, color=K562_LAD_n)) +
                  geom_point(data=P[P$K562_LAD==0, ], size=0.25, alpha=0.1) +
                  geom_point(data=P[P$K562_LAD==1, ], size=0.8, alpha=0.5) +
                  theme_bw() +
                  geom_line(data=RM_PROseq, aes(x=x.mean, y=y.AD), color='black',  size=1) +
                  geom_line(data=RM_PROseq, aes(x=x.mean, y=y.iAD), color='black',  size=1) +
                  geom_line(data=RM_proseq_melt, aes(x=x.mean, y=value, color=variable), size=0.5) +
                  labs(y='log10(CAGE)', x='log10(SuRE)') +
                  theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_line(size = .5, color = "grey"),
                        axis.line = element_line(size=.7, color = "black"),
						legend.background = element_rect(fill="transparent"),
                        text = element_text(size=14),
						legend.justification =c(0,1),
		 			    legend.position=c(0.05,0.95),
						legend.title=element_blank()) +
                  scale_color_manual(values=COL_lad_n) +
                  coord_equal(ratio=1) +
          		  guides(colour = guide_legend(override.aes = list(size=1,linetype=0, alpha=1)))


p_cage = ggplot(P[P$K562_LAD==0, ], aes(x=K562_SuRE, y=K562_CAGE_jitter,
                                          color=class_LAD_n, fill=class_LAD_n)) +
    geom_point(size=0.1, alpha=0.1) +
    geom_point(data=p_classes, size=0.5, alpha=1) +
    theme_bw() +
    labs(y='log10(CAGE)', x='log10(SuRE)') +
    ggtitle('GROcap classification on CAGE data') +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = .5, color = "grey"),
          axis.line = element_line(size=.7, color = "black"),
          legend.background = element_rect(fill="transparent"),
          text = element_text(size=14),
          legend.justification =c(0,1),
          legend.position=c(0.05,0.95),
          legend.title=element_blank()) +
    scale_color_manual(values=COL_class_LAD_n) +
    coord_equal(ratio=1)

p_cageclass = ggplot(P[P$K562_LAD==0, ], aes(x=K562_SuRE, y=K562_CAGE_jitter,
                                          color=class_LAD_n, fill=class_LAD_n)) +
    geom_point(data=p_classes, size=0.5, alpha=1) +
    geom_point(data=P[P$K562_LAD==0, ], size=0.1, alpha=0.1) +
    theme_bw() +
    labs(y='log10(CAGE)', x='log10(SuRE)') +
    ggtitle('GROcap classification on CAGE data') +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = .5, color = "grey"),
          axis.line = element_line(size=.7, color = "black"),
          legend.background = element_rect(fill="transparent"),
          text = element_text(size=14),
          legend.justification =c(0,1),
          legend.position=c(0.05,0.95),
          legend.title=element_blank()) +
    scale_color_manual(values=COL_class_LAD_n) +
    coord_equal(ratio=1)


p_proseq = ggplot(p_classes, aes(x=K562_SuRE, y=K562_PROseq_jitter,
                                          color=class_LAD_n, fill=class_LAD_n)) +
    geom_point(data=P[P$K562_LAD==0, ], size=0.1, alpha=0.1) +
    geom_point(data=p_classes, size=0.5, alpha=1) +
    theme_bw() +
    labs(y='log10(PROseq)', x='log10(SuRE)') +
    ggtitle('GROcap classification on PROseq data') +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = .5, color = "grey"),
          axis.line = element_line(size=.7, color = "black"),
          legend.background = element_rect(fill="transparent"),
          text = element_text(size=14),
          legend.justification =c(0,1),
          legend.position=c(0.05,0.95),
          legend.title=element_blank()) +
    scale_color_manual(values=COL_class_LAD_n) +
    coord_equal(ratio=1)



p_proclass = ggplot(p_classes, aes(x=K562_SuRE, y=K562_PROseq_jitter,
                                          color=class_LAD_n, fill=class_LAD_n)) +
    geom_point(data=p_classes, size=0.5, alpha=1) +
    geom_point(data=P[P$K562_LAD==0, ], size=0.1, alpha=0.1) +
    theme_bw() +
    labs(y='log10(PROseq)', x='log10(SuRE)') +
    ggtitle('GROcap classification on PROseq data') +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = .5, color = "grey"),
          axis.line = element_line(size=.7, color = "black"),
          legend.background = element_rect(fill="transparent"),
          text = element_text(size=14),
          legend.justification =c(0,1),
          legend.position=c(0.05,0.95),
          legend.title=element_blank()) +
    scale_color_manual(values=COL_class_LAD_n) +
    coord_equal(ratio=1)

pdf('cl20190104_CAGE_PROseq.pdf', width=8, height=5, useDingbats=F)

grid.arrange(p_lad_cage, p_lad_proseq, nrow=1)

p1 = ggMarginal(p_cage, margins="y", type="density", yparams(alpha=0.1),
                groupFill=T, groupColour=T, size=2)
grid.arrange(p1)
p2 = ggMarginal(p_cageclass, margins="y", type="density", yparams(alpha=0.1),
                groupFill=T, groupColour=T, size=2)
grid.arrange(p2)
p3 = ggMarginal(p_proseq, margins="y", type="density", yparams(alpha=0.1),
                groupFill=T, groupColour=T, size=2)
grid.arrange(p3)
p4 = ggMarginal(p_proclass, margins="y", type="density", yparams(alpha=0.1),
                groupFill=T, groupColour=T, size=2)
grid.arrange(p4)
dev.off()



figD = ggplot(P[P$class_LAD!='boundary', ], aes(x=class_LAD, y=tissues_expressed, color=class_LAD)) +
    geom_quasirandom() +
    theme_bw() +
    scale_color_manual(values=COL_class) +
    xlab('promoter class') +
    ylab('number of tissues expressed') +
    guides(color=F) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 color='black', geom = "crossbar", width = 0.4) +
    theme(plot.title = element_text(hjust = 0.5, face='bold'),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank())

pdf('cl20181017_SuRE_vs_GROcap_prom.pdf', useDingbats=F, width=10, height=7)
plot_grid(figA, figB, figC, figD, labels='AUTO', nrow=2)
dev.off()
