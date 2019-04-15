#!/bin/bash
usage="$(basename "$0") [-h] -t <...> -u <INT> -d <INT>

script to return a bed-file with regions around reference-points.

where:
    -h  show this help text
    -t  text file with: chromosome name in 1st column, TSS position in
        2nd column, identifier in 3rd column, and strand in 4th column
        (other columns are ignored)
    -g  genome in fasta format
    -u  upstream region from TSS
    -d  downstream region from TSS"



while getopts ':h:t:g:u:d:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    t) tss="$OPTARG"
       ;;
    g) genome="$OPTARG"
       ;;
    u) up="$OPTARG"
       ;;
    d) down="$OPTARG"
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




awk -v OFS='\t' -v up=$up -v down=$down \
        '{
            if ($4 == ""){
                print $1, $2 - up, $2 + down, $3, ".", "."
            }
            else{
                if ($4 == "+"){
                    start = $2 - up
                    end = $2 + down
                } else{
                    start = $2 - down
                    end = $2 + up
                }
                    print $1, start, end, $3, ".", $4
            }}' $tss | \
    bedtools nuc -C -pattern CG -fi $genome -bed /dev/stdin | \
    awk -vOFS='\t' \
        '{
            if (NR==1){
                print "seqnames", "start", "end", "name",
                      "CpG", "C+G", "exp_CpG", "CpG_OE"
            } else {
                CGdinuc=$16
                CGsum=$10 + $11
                l=$15
                CGE=((CGsum/2)*(CGsum/2))/l
                if (CGE==0){
                    print $1, $2, $3, $4, CGdinuc, CGsum, CGE, "NA"
                } else {
                    print $1, $2, $3, $4, CGdinuc, CGsum, CGE, CGdinuc/CGE
                }
            }
        }'
