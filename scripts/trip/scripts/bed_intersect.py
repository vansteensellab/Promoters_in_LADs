import os
# just a simple wrapper to perform an awk call
# in snakemake I will have to deal with escaping both {} and "
# Arguments:
#    - bed file with mapping locations
#    - bed file with data track for comparison
#    - output file
command = ("nice -18 awk '{if ($1!=\"*\") {print $0}}' %s | "
           "bedtools intersect -a - -b "
           "%s -wa -wb | awk 'BEGIN {"
           "print \"barcode\tchrom_state\tcount\ttotal\""
           "}{"
           "  count[$4\"\\t\"$10] += $5;"
           "  total[$4\"\\t\"$10] = $6"
           "}END{ "
           "  for (i in count){ "
           "    print i\"\\t\"count[i]\"\\t\"total[i]"
           "  }"
           "}' > %s")
os.system(command % (snakemake.input[0], snakemake.input[1],
                     snakemake.output[0]))
