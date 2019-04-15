import scipy.sparse
import pysam
import itertools
import pickle

# chrom_sizes = 'GRCm38_68.chrom.sizes'
# snp_sites = 'CAST_129S1_isec/sites.txt'
# bam_in = '../../damid/results/mapped/cl18_Dam_129_r1_comb_filtered.bam'


# bed_ambiguous = open('../../damid/cl18_Dam_129_r1_ambiguous.bed', 'w')
# bed_nosnp = open('../../damid/cl18_Dam_129_r1_nosnp.bed', 'w')
# bed_cast = open('../../damid/cl18_Dam_129_r1_CAST.bed', 'w')
# bed_129s1 = open('../../damid/cl18_Dam_129_r1_129S1.bed', 'w')

def split_map(b_file, bed_ambiguous, bed_nosnp, bed_cast, bed_129s1, chrom_dict):
    for line in b_file:
        if not line.is_unmapped:
            chrom = line.reference_name
            if chrom in chrom_dict:
                line_snps = chrom_dict[chrom][:,line.get_reference_positions()]
                cx = line_snps.tocoo()
                allele_count = [0,0]
                score = 0
                if len(cx.row) > 0:
                    seq = line.query_alignment_sequence
                    for allele, pos, base in itertools.zip_longest(cx.row, cx.col, cx.data):
                        if base.find(seq[pos-1]) > -1:
                            allele_count[allele] += 1
                    if sum(allele_count) > 0:
                        score = max(allele_count) / sum(allele_count)
                        if allele_count[0] > allele_count[1]:
                            bed_out = bed_cast
                        elif allele_count[1] > allele_count[0]:
                            bed_out = bed_129s1
                        else:
                            bed_out = bed_ambiguous
                    else:
                        score = 0
                        bed_out = bed_ambiguous
                else:
                    bed_out = bed_nosnp
                strand = '-' if line.is_reverse else '+'
                out_line = "%s\t%i\t%i\t%s\t%g\t%s\t%i\t%i" % (chrom,
                                                               line.reference_start,
                                                               line.reference_end,
                                                               line.query_name,
                                                               round(score,3),
                                                               strand,
                                                               allele_count[0],
                                                               allele_count[1])
                print(out_line, file=bed_out)


if __name__ == '__main__':
    snakein = snakemake.input
    snakeout = snakemake.output
    bed_ambiguous = open(snakeout.bed_ambiguous, 'w')
    bed_nosnp = open(snakeout.bed_nosnp, 'w')
    bed_cast = open(snakeout.bed_cast, 'w')
    bed_129s1 = open(snakeout.bed_129S1, 'w')
    b_file = pysam.AlignmentFile(snakein.bam, 'rb')
    chrom_dict = pickle.load(open(snakein.chrom_dict, 'rb'))
    split_map(b_file, bed_ambiguous, bed_nosnp, bed_cast, bed_129s1, chrom_dict)
    bed_ambiguous.close()
    bed_nosnp.close()
    bed_cast.close()
    bed_129s1.close()
