#!/bin/bash
usage="$(basename "$0") [-h] -l <...> -e <...>

script takes a tab seperated file with gene-id's in the first column and peak
names in the second, and count file from fantom. It then reports the sum of
reads normalized over library size for each column (cell-types/tissues)
[(count1 / total1) + (count2/total2) + etc.], and the number of columns in
which count > 0 (number of tissues expressed).

where:
    -h  show this help text
    -l  two column link file between identifiers
    -e  file with expression (fantom data)"



while getopts ':h:l:e:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    e) exp="$OPTARG"
       ;;
    l) link="$OPTARG"
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
done
shift $((OPTIND - 1))


awk -vOFS='\t' '{
        if (NR==FNR){
            arr[$2][$1] = 1
        } else {
            if ($1=="01STAT:MAPPED"){
                for (i=2;i<NF;i++){
                    total[i] = $i
                }
            } else if ($1 in arr) {
                sum=0
                number=0
                for (i=2;i<NF;i++){
                    sum+=$i/total[i]
                    if ($i > 0){
                        number++
                    }
                }
                for (name in arr[$1]){
                    print name, sum, number
                }
            }
        }
    }' <(zcat $link) \
    <(zcat $exp) | \
    gzip -c
