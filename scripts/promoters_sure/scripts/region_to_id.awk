#!/usr/bin/awk -f
{
    if (NR==1){
        gsub($1"\t"$2"\t"$3, "ID", $0)
        print($0)
    }
    if (NR==FNR){
        arr[$1"\t"$2"\t"$3] = $0
    } else {
        line=arr[$1"\t"$2"\t"$3]
        gsub($1"\t"$2"\t"$3, $4, line)
        print line
    }
}
