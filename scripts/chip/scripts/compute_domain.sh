#!/bin/bash
usage="$(basename "$0") [-h] -t <...> -u <INT> -d <INT>

script to return a bed-file with regions around reference-points.

where:
    -h  show this help text
    -b  bam files
    -d  bed file with domains
    -l  labels
    -o  output file
    -t  threads
    "

threads=1

while getopts ':h:d:b:o:l:?t:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    d) domain="$OPTARG"
       ;;
    b) bam_list="$OPTARG"
       ;;
    l) labels="$OPTARG"
       ;;
    o) out="$OPTARG"
       ;;
    t) threads="$OPTARG"
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

name_arr=$(zcat $domain | awk '{ arr[$4]==1 }END{ for (a in arr) { print a }}')


echo "name total_length "$labels | sed 's/ /\t/g'  > $out
for name in $name_arr
do
    tmp_in=$(tempfile)".bed"
    tmp_npz=$(tempfile)".npz"
    tmp_out=$(tempfile)
    zcat $domain | grep -P "[ \t]"$name"$" | \
        awk -vOFS='\t' '{
                if (NR==1){
                    chr=$1
                    start=$2
                    end=$3
                } else {
                    if (chr==$1 && end==$2){
                        end = $3
                    } else {
                        print chr, start, end
                        chr=$1
                        start=$2
                        end=$3
                    }
                }
            }END{
                print chr, start, end
            }'  > $tmp_in
    multiBamSummary BED-file -p $threads \
                                 --BED $tmp_in \
                                 -b $bam_list \
                                 -out $tmp_npz \
                                 --outRawCounts $tmp_out
    tail -n+2 $tmp_out | \
        awk -vOFS='\t' -vname=$name \
            '{
                len+=$3 - $2
                for (i=4;i<=NF;i++){
                    cnt[i] = cnt[i] + $i
                }
             } END {
                 for (i=1;i<=length(cnt);i++){
                     len = len "\t" cnt[i + 3]
                 }
                 print name, len
             }' >> $out
    rm $tmp_in $tmp_npz $tmp_out
done
