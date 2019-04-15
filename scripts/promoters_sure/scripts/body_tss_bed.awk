#!/usr/bin/awk -f
BEGIN{
    OFS="\t"
}
{
    if ($3 - $2 + body_end > body_start){
        if ($6 == "+"){
            print $1, $2 + tssr_start, $2 + tssr_end, $4, $5, $6 > tss
            print $1, $2 - 500, $2 + 500, $4, $5, $6 > kb
            print $1, $2 + body_start, $3 + body_end, $4, $5, $6 > body

        } else {
            print $1, $3 - tssr_end, $3 - tssr_start, $4, $5, $6 > tss
            print $1, $3 - 500, $3 + 500, $4, $5, $6 > kb
            print $1, $2 - body_end, $3 - body_start, $4, $5, $6 > body
        }
    }
}
