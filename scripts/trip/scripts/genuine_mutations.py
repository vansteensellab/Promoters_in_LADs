import subprocess
import re

raw_mut = snakemake.input[0]
starcode = snakemake.input[1]
mapping = snakemake.input[2]

not_genuine = snakemake.output[0]
genuine_mapped = snakemake.output[1]
genuine_unmapped = snakemake.output[2]


starcode_set = set()
with open(starcode) as in_starcode:
    for line in in_starcode.readlines():
        starcode_set.add(line.split('\t')[0])

mapping_set = set()
with open(mapping) as in_mapping:
    line = in_mapping.readline()
    for line in in_mapping.readlines():
        mapping_set.add(line.split('\t')[0])

out_not_genuine = open(not_genuine, 'w')
out_mapped = open(genuine_mapped, 'w')
out_unmapped = open(genuine_unmapped, 'w')
with open(raw_mut) as in_raw:
    for line in in_raw.readlines():
        line_split = line.strip().split('\t')
        bc = line_split[0]
        out_line = '\t'.join(line_split[0:3])
        if bc in starcode_set:
            if bc in mapping_set:
                print(out_line, file=out_mapped)
            else:
                print(out_line, file=out_unmapped)
        else:
            print(out_line, file=out_not_genuine)
out_not_genuine.close()
out_mapped.close()
out_unmapped.close()
