#!/bin/bash
usage="$(basename "$0") [-h] [-s <...>] [-t <...>]

script to convert a bed-file into gff3 format

where:
    -h|--help        show this help text
    -s|--source      source; 2nd field of gff file. (e.g. Cufflinks,Maker) [default: .]
    -t|--type        data type; 3rd field of gff file. (e.g. TF_binding_site) [default: match]
    -b|--bed         bed file
    -n|--narrowPeak  bed is in narrowPeak format"

SOURCE="."
TYPE="match"


if [[ ${1^^} == "-H" || ${1^^} == "--HELP" ]]; then
    echo "$usage"; exit 1
else
  while [[ $# -gt 0 ]]; do
    opt="$1"
    shift;
    case "$opt" in
      "-s"|"--source"     ) SOURCE="$1"; shift;;
      "-t"|"--type"       ) TYPE="$1"; shift;;
      "-b"|"--bed"        ) BED="$1"; shift;;
      "-n"|"--narrowPeak" ) NARROW=true;;
      *                   ) echo "ERROR: Invalid option: \""$opt"\"" >&2
                            exit 1;;
    esac
  done
fi


if [[ "$BED" == "" ]]; then
  echo "ERROR: Bed file is required (-b)." >&2
  exit 1
fi

awk -vOFS='\t' -vtype=$TYPE -vsource=$SOURCE -vnarrow=NARROW \
  '{
     id=""
     strand="."
     score=0
     phase="."
     if ($4 != ""){
       id="ID="$4
     }
     if ($5 != ""){
       score=$5
     }
     if ($6 != ""){
       strand=$6
     }
     if (narrow == "true"){
         id=id";signalValue="$7";pValue="$8";qValue="$9";peak="$10
     }
     print $1, source, type, $2 + 1, $3, score, strand, phase, id
   }' $BED
