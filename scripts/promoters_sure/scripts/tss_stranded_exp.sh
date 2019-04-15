#!/bin/bash
usage="$(basename "$0") [-h] -t <...> -u <INT> -d <INT> -p <bigwig> -m <bigwig>

script to calculate stranded expression around TSS's.

where:
    -h  show this help text
    -t  text file with: chromosome name in 1st column, TSS position in
        2nd column, identifier in 3rd column, and strand in 4th column
        (other columns are ignored)
    -u  upstream region to include in calculation
    -d  downstream region to include in calculation
    -p  bigWig with signal for positive strand
    -m  bigWig with signal for negative strand"



while getopts ':h:t:u:d:p:m:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    t) tss="$OPTARG"
       ;;
    u) up="$OPTARG"
       ;;
    d) down="$OPTARG"
       ;;
    p) plus="$OPTARG"
       ;;
    m) min="$OPTARG"
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


function bed {
    awk -vOFS='\t' -v u=$1 -v d=$2 '{
                    if ($4 == "+"){
                        start = $2 - u - 1 #-1 for conversion to 0-based bed
                        end = $2 + d
                    } else {
                        start = $2 - d - 1
                        end = $2 + u
                    }
                    print $1, start, end, $3, "0", $4
                }' $3
}


paste <(bed $up $down $tss | \
        bwtool summary -skip-median -with-sum -keep-bed /dev/stdin \
                              $plus /dev/stdout ) \
      <(bed $up $down $tss | \
        bwtool summary -skip-median -with-sum -keep-bed /dev/stdin \
                              $min /dev/stdout ) | \
    awk -vOFS='\t' 'BEGIN{
                        print "name", "plus", "min", "sense", "antisense"
                    }{
                        if ($6 ~ /[+-]/){
                            print $4 , $12 , $24 , $6=="+"?$12:$24, $6=="+"?$24:$12
                        } else {
                            print $4 , $12 , $24
                        }

                    }' | \
    gzip -c
