import Levenshtein
import subprocess

if __name__=='__main__':
    snakein = snakemake.input
    snakeout = snakemake.output
    vcf_dict = snakemake.params.vcf
    command = ("samtools faidx {fasta} {chr}:{start}-{end} | "
               "bcftools consensus {vcf}")
    with open(snakein.table,'r') as fin:
        with open(snakeout[0], 'w') as fout:
            print('barcode\treverse_read\tforward_read', file=fout)
            line = next(fin)
            for line in fin:
                line_split = line.split('\t')
                line_split[-1] = line_split[-1].strip()
                barcode = line_split[0]
                if line_split[1] != '' and line_split[8] != '':
                    chr_dict = {'rev':line_split[1], 'fwd':line_split[8]}
                    strand_dict = {'rev':line_split[2], 'fwd':line_split[9]}
                    pos_dict = {'rev':int(line_split[3]), 'fwd':int(line_split[10])}
                    seq_dict = {'rev':line_split[15], 'fwd':line_split[16]}
                    calls = {'rev':'unknown', 'fwd':'unknown'}
                elif line_split[1] != '':
                    chr_dict = {'rev':line_split[1]}
                    strand_dict = {'rev':line_split[2]}
                    pos_dict = {'rev':int(line_split[3])}
                    seq_dict = {'rev':line_split[8]}
                    calls = {'rev':'unknown', 'fwd':'missing'}
                elif line_split[8] != '':
                    chr_dict = {'fwd':line_split[8]}
                    strand_dict = {'fwd':line_split[9]}
                    pos_dict = {'fwd':int(line_split[10])}
                    seq_dict = {'fwd':line_split[16]}
                    calls = {'rev':'missing', 'fwd':'unknown'}
                for r in chr_dict.keys():
                    start = pos_dict[r] if strand_dict[r]=='+' else pos_dict[r] - len(seq_dict[r])
                    end = pos_dict[r] + len(seq_dict[r]) if strand_dict[r]=='+' else pos_dict[r]
                    seq_dist = {}
                    for species in vcf_dict:
                        seqname = chr_dict[r].replace('chr','')
                        bcf = subprocess.run(command.format(fasta=snakein.fasta,
                                                            chr=seqname,
                                                            start=start,
                                                            end=end,
                                                            vcf=vcf_dict[species]),
                                             shell=True, stdout=subprocess.PIPE)
                        ref_seq = bcf.stdout.decode().split('\n')[1]
                        seq_dist[species] = Levenshtein.distance(seq_dict[r], ref_seq)
                    species_list = list(vcf_dict.keys())
                    if seq_dist[species_list[0]] > seq_dist[species_list[1]]:
                        calls[r] = species_list[1]
                    elif seq_dist[species_list[0]] < seq_dist[species_list[1]]:
                        calls[r] = species_list[0]
                    else:
                        calls[r] = 'unknown'
                print('{barcode}\t{rev}\t{fwd}'.format(barcode=barcode,
                                                       rev=calls['rev'],
                                                       fwd=calls['fwd']), file=fout)
