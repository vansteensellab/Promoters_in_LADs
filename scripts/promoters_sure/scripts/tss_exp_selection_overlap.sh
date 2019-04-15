#!/bin/bash
usage="$(basename "$0") [-h] -t <bed> -e <bed> -s <...> -l <...>

script to select unique TSS positions of transcripts and select only those
overlapping with expression peaks (e.g. CAGE).
Intersect with exact TSS position and peak bed is obtained, the first transcript
at each position is used in case of overlapping TSS's.

where:
    -h  show this help text
    -t  7-column bed file with TSS's of transcripts with column 4 containing
        transcript-ID and column 7 containing gene-ID
    -e  bed file with expression peaks
    -d  max distance between elements
    -s  output file for 7 column bed file containing the selected TSS's
    -l  link table to connect id's of expression peaks with transcript id's"



while getopts ':h:t:e:s:d:l:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    t) tss="$OPTARG"
       ;;
    e) exp="$OPTARG"
       ;;
    s) selection="$OPTARG"
       ;;
    l) link="$OPTARG"
       ;;
    d) dist="$OPTARG"
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


## For multiple transcripts coming from the same gene, we want to select transcription start
## sites at least 500bp apart.

## select unique transcript start sites which overlap with a cage peak.
## CAGE peaks have at least 1 transcript in one of the tissues.
## (multiple transcripts of same gene can start at same position we don't want those).


bedtools closest -s -d -wb -a $tss -b <(bedtools sort -i $exp) | \
    awk -v OFS='\t' -v sel=$selection -v link=$link -v dist=$dist '{
        if ($8 > 0 && $17 < dist){
            name=$1 FS $3 FS $6 FS $7 FS $11
            if (!(name in seen)){
                location=$1 FS $3 FS $6 FS $7
                if (!(location in loc)){
                    print $1, $3, $4, $6, $7 | "gzip > " sel
                }
                loc[location] = 1
                print $4, $11 | "gzip > " link
            }
            seen[name]=1
        }}'
