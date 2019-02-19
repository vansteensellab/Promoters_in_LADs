#!/bin/bash

bedtools closest -t all -D a -a $1 -b - | awk 'BEGIN {
print "barcode\tCpG_name\tdistance\tcount\ttotal"
}{
  count[$4"\t"$10"\t"$11] += $5
  total[$4"\t"$10"\t"$11] = $6
}END{
  for (i in count){
    print i"\t"count[i]"\t"total[i]
  }
}'
