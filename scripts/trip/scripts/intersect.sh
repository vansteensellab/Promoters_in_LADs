bedtools intersect -a $1 -b - -wa -wb | awk -v s="$2" 'BEGIN {
print "barcode\t"s"\tcount\ttotal"
}{
  count[$4"\t"$10] += $5
  total[$4"\t"$10] = $6
}END{
  for (i in count){
    print i"\t"count[i]"\t"total[i]
  }
}'
