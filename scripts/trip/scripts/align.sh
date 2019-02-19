#!/bin/bash
INPUT=$1
LOG=$2
THREADS=$3
OPTIONS=$4
INDEX=$5
OUTPUT=$6

gunzip -c $INPUT | \
    awk -vLOG=$LOG '
        BEGIN{
            hit = 0
            short = 0
        }
        {
            step=NR%4
            if (step==0 && length(a[2])>6){
                hit += 1
                for (i in a){
                    print a[i]
                }
                print $0
            } else if (step!=0){
                a[step]=$0;
            } else {
                short += 1
            }
        } END {
            print "filtering before mapping with bowtie2:" > LOG;
            printf "%i\ttreads; of these:\n", hit+short > LOG;
            printf "  %i (%2.2f%%) were long enough (> 6bp)\n", hit, hit/(hit+short)*100 > LOG;
            printf "  %i (%2.2f%%) were too short (<= 6bp)\n\n", short, short/(hit+short)*100 > LOG;
            print "stats from bowtie2:" > LOG;
        }' | bowtie2 -p $THREADS $OPTIONS -x $INDEX -U - | samtools view -Sb - > $OUTPUT ;

samtools flagstat $OUTPUT > $LOG
