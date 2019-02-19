import gzip
import subprocess
import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

barcode_file = str(snakemake.input)
output_file = snakemake.output[0]

target = snakemake.params.target
spacer_list = snakemake.params.spacer_list
gap_list = snakemake.params.gap_list


with gzip.open(barcode_file) as f_in:
    with open(output_file, 'w') as f_out:
        for byte in f_in.readlines():
            line = str(byte, 'utf-8')
            line_split = line.strip().split('\t')
            barcode = line_split[2]
            dna_str = line_split[1]

            match_list = [re.search(spacer, dna_str) for spacer in spacer_list]
            indel_list = [match.start() - gap
                          for gap,match in zip(gap_list, match_list)
                          if match is not None]
            if len(indel_list) > 0:
                score = indel_list[0]
            else:
                score = 'NA'

            target_match = re.search(target, dna_str)
            if target_match is not None:
                status = 'wt'
                score = 0 if score is 'NA' else score
            elif score is not 'NA':
                status = 'ins' if score > 0 else 'del' if score != 0 else 'wt_point_mut'
            else:
                status = 'not_clear'
            print('\t'.join((barcode, status, str(score), dna_str)), file=f_out)
