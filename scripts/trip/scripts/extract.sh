#!/bin/sh

bwtool extract bed $1 /dev/stdin /dev/stdout | awk -F'[,_\t|/]' 'BEGIN{
    print "barcode\tsignal\tcount\ttotal"
  }{
    i=$4"\t"($8 + $9 + $10 + $11)/4
    count[i] += $5
    total[i] = $6
  }END{
    for (i in count){
      print i"\t"count[i]"\t"total[i]
    }
  }'