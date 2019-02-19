#!/bin/bash

zcat $1 | \
    awk -vOFS='\t' -vW=$2 '{
        start=$2 - W
        end=$3 + W
        if (end > start){
            print $1, $2 - W, $3 + W
        }
    }' | bedtools intersect -c -a $3 -b -
