#!/bin/bash
usage="$(basename "$0") [-h] -t <...> -u <INT> -d <INT>

script to return a bed-file with regions around reference-points.

where:
    -h  show this help text
    -t  text file with: chromosome name in 1st column, TSS position in
        2nd column, identifier in 3rd column, and strand in 4th column
        (other columns are ignored)
    -u  upstream region to include in calculation
    -d  downstream region to include in calculation "



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

awk -vOFS='\t' -v u=$up -v d=$down '{
                if ($4 == "+"){
                    start = $2 - u - 1 #-1 for conversion to 0-based bed
                    end = $2 + d
                } else {
                    start = $2 - d - 1
                    end = $2 + u
                }
                print $1, start, end, $3, "0", $4==""?".":$4
            }' $tss
