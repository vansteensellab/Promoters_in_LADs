import pysam
import gzip
import os
import re
from Bio import SeqIO
from collections import Counter


def parse_sam(sam_file, starcode_set, mut_dict, max_soft_clip=5,
              min_first_match=10, remap_soft_clip=17):
    remap_list = []
    map_dict = {}
    map_stat_dict = {}
    length_dict = {}
    mut_sum = 0
    for line in pysam.AlignmentFile(sam_file):
        if re.match(r'.*[:_]([A-Z]+)$', line.query_name) is None:
            print('test:')
            print(sam_file)
            print(line.query_name)
        bc_this = re.match(r'.*[:_]([A-Z]+)$', line.query_name).groups(1)[0]
        if bc_this in starcode_set or bc_this in mut_dict:
            if bc_this in mut_dict:
                bc_this = mut_dict[bc_this]
                is_mut = 1
            else:
                is_mut = 0
            if bc_this not in map_stat_dict:
                map_stat_dict[bc_this] = [0, 0, 0, 0]
                length_dict[bc_this] = [0, 0, 0, 0]
            if line.is_reverse:
                start_pos = line.reference_end
                ori = '-'
            else:
                start_pos = line.reference_start
                ori = '+'
            add_mapping = False
            if line.cigarstring is None:
                mapping = ('*', ori, 0)
                map_stat_dict[bc_this][1] += 1
                length_dict[bc_this][1] += len(line.query_sequence)
            else:
                if line.is_reverse:
                    first_cigar = line.cigartuples[-1]
                else:
                    first_cigar = line.cigartuples[0]
                mapping = (line.reference_name, ori, start_pos)
                if first_cigar[0] == 4:
                    if first_cigar[1] <= max_soft_clip:
                        add_mapping = True
                    elif first_cigar[1] >= remap_soft_clip:
                        quality_dict = {'phred_quality':
                                        line.query_qualities[:first_cigar[1]]}
                        query_seq = line.query_sequence[:first_cigar[1]]
                        seq = SeqIO.SeqRecord(query_seq, line.query_name,
                                              description=line.query_name,
                                              letter_annotations=quality_dict)
                        remap_list.append(seq)
                        map_stat_dict[bc_this][3] += 1
                        length_dict[bc_this][3] += len(line.query_sequence)
                    else:
                        map_stat_dict[bc_this][2] += 1
                        length_dict[bc_this][2] += len(line.query_sequence)
                if first_cigar[0] == 0 and first_cigar[1] >= min_first_match:
                    add_mapping = True
            mapping_quality = line.mapping_quality
            if add_mapping:
                mut_sum += is_mut
                map_stat_dict[bc_this][0] += 1
                length_dict[bc_this][0] += len(line.query_sequence)
                seq = line.query_alignment_sequence
                if bc_this in map_dict:
                    if mapping in map_dict[bc_this]:
                        map_dict[bc_this][mapping][0] += 1
                        map_dict[bc_this][mapping][1] += mapping_quality
                        map_dict[bc_this][mapping][2].append(seq)
                    else:
                        map_dict[bc_this][mapping] = [1, mapping_quality, [seq]]
                else:
                    map_dict[bc_this] = {mapping: [1, mapping_quality, [seq]]}
    print(mut_sum)
    return(map_dict, remap_list, map_stat_dict, length_dict)


def write_bed(map_dict_in, out_file, max_dist):
    def sort_keys(key, map_dict_in):
        return (map_dict_in[key][0])
    with open(out_file, 'w') as f_out:
        for bc in map_dict_in:
            this_map_dict = map_dict_in[bc]
            sorted_key_list = sorted(this_map_dict.keys(),
                                     key=lambda elem: sort_keys(elem,
                                                                this_map_dict),
                                     reverse=True)
            total_reads = sum(value[0] for value in this_map_dict.values())
            while len(sorted_key_list) > 0:
                top_key, \
                    top_reads, \
                    top_mapq, \
                    sorted_key_list = refine_map(this_map_dict,
                                                 sorted_key_list,
                                                 max_dist)
                reference_name, ori, start_pos = top_key
                if ori == '+':
                    start = start_pos - 4
                    end = start_pos
                elif ori == '-':
                    start = start_pos
                    end = start_pos + 4
                f_out.write('%s\t%i\t%i\t%s\t%i\t%i\n' % (reference_name,
                                                          start,
                                                          end,
                                                          bc,
                                                          top_reads,
                                                          total_reads))


def top_map(map_dict_in, max_dist):

    def sort_keys(key, map_dict_in):
        return (map_dict_in[key][0])

    map_dict_out = {}
    for bc in map_dict_in:
        this_map_dict = map_dict_in[bc]
        if len(this_map_dict) == 1:
            key = next(iter(this_map_dict))
            reference_name, ori, start_pos = key
            total_reads, mapq_sum, seq_list = this_map_dict[key]
            top_mapq = mapq_sum
            top_reads = total_reads
            second_mapq = 0
            second_reads = 0
            seq = Counter(seq_list).most_common(1)[0][0]
        else:
            total_reads = sum(value[0] for value in this_map_dict.values())
            sorted_key_list = sorted(this_map_dict.keys(),
                                     key=lambda elem: sort_keys(elem,
                                                                this_map_dict),
                                     reverse=True)
            seq_list = this_map_dict[sorted_key_list[0]][2]
            seq = Counter(seq_list).most_common(1)[0][0]
            top_key, \
                top_reads, \
                top_mapq, \
                sorted_key_list = refine_map(this_map_dict, sorted_key_list,
                                             max_dist)
            reference_name, ori, start_pos = top_key
            if len(sorted_key_list) > 0:
                second_key, \
                    second_reads, \
                    second_mapq, \
                    sorted_key_list = refine_map(this_map_dict,
                                                 sorted_key_list,
                                                 max_dist)
            else:
                second_reads = 0
                second_mapq = 0
            sorted_key_list = sorted(this_map_dict.keys(),
                                     key=lambda elem: sort_keys(elem,
                                                                this_map_dict),
                                     reverse=True)
        map_dict_out[bc] = [reference_name, ori, start_pos, total_reads,
                            top_mapq, top_reads, second_mapq, second_reads, seq]
    return(map_dict_out)


def refine_map(this_map_dict, sorted_key_list, max_dist):
    top_key = sorted_key_list.pop(0)
    this_reads = this_map_dict[top_key][0]
    this_mapq = this_map_dict[top_key][1]
    i = 0
    while i < len(sorted_key_list):
        if (sorted_key_list[i][0] == top_key[0] and
                sorted_key_list[i][1] == top_key[1] and
                abs(sorted_key_list[i][2] - top_key[2]) < max_dist):
            close_key = sorted_key_list.pop(i)
            this_reads += this_map_dict[close_key][0]
            this_mapq += this_map_dict[close_key][1]
        # elif (abs(sorted_key_list[i][2] - top_key[2]) < max_dist):
        #     print(sorted_key_list[i][0])
        #     print(top_key[0])
        #     print(sorted_key_list[i][1])
        #     print(top_key[1])
        i += 1
    return(top_key, this_reads, this_mapq, sorted_key_list)


if __name__ == '__main__':
    starcode_set = set()
    snakein = snakemake.input
    snakeout = snakemake.output
    snakeparam = snakemake.params
    threads = snakemake.threads

    with open(snakein.count) as cf:
        for line in cf.readlines():
            barcode = line.split('\t')[0]
            if barcode not in starcode_set:
                starcode_set.add(barcode)
    mut_dict = {}
    if 'mutated' in snakein.keys():
        with open(snakein.mutated) as mf:
            for line in mf.readlines():
                line_split = line.strip().split('\t')
                mut_dict[line_split[0]] = line_split[2]

    (map_dict, remap_list,
        map_stat_dict, length_dict) = parse_sam(snakein.bam, starcode_set,
                                                mut_dict)
    if len(remap_list) > 0:
        with gzip.open(snakeout.remap_fq, "wt") as fqfile:
            SeqIO.write(remap_list, fqfile, 'fastq')
        options = snakeparam.options
        align_command = ("bowtie2 -p %i -t %s"
                         " -x %s -U %s | "
                         "samtools view -Sb - > %s" % (threads,
                                                       ' '.join(options),
                                                       snakeparam.bowtie_index,
                                                       snakeout.remap_fq,
                                                       snakeout.remap))
        os.system(align_command)
        remap_sam_out = parse_sam(snakeout.remap, starcode_set,
                                  mut_dict)
        for bc in remap_sam_out[0]:
            if bc in map_dict:
                this_remap_dict = remap_sam_out[0][bc]
                this_map_dict = map_dict[bc]
                for loc in this_remap_dict:
                    if loc in this_map_dict:
                        this_map_dict[loc] = [a + b for a, b in
                                              zip(this_remap_dict[loc],
                                                  this_map_dict[loc])]
                    else:
                        this_map_dict[loc] = this_remap_dict[loc]
            else:
                map_dict[bc] = remap_sam_out[0][bc]

    else:
        os.system('touch %s; touch %s;' % (snakeout.remap_fq,
                                           snakeout.remap))



    write_bed(map_dict, snakeout.bed, snakeparam.max_dist)

    top_map_dict = top_map(map_dict, snakeparam.max_dist)
    with open(snakeout.table, 'w') as fout:
        fout.write('\t'.join(('barcode', 'seqname', 'ori', 'start_pos',
                              'total_mapped', 'mapq_sum1', 'reads1', 'mapq_sum2',
                              'reads2', 'seq')))
        fout.write('\n')
        for bc in top_map_dict:
            fout.write('\t'.join([bc] + [str(item) for item in
                                         top_map_dict[bc]]))
            fout.write('\n')

    with open(snakeout.stats, 'w') as f_stats:
        with open(snakeout.length, 'w') as f_length:
            f_stats.write('\t'.join(('barcode', 'total_reads', 'aligned',
                                     'aligned_correct', 'realigned',
                                     'realigned_correct')))
            f_stats.write('\n')
            f_length.write('\t'.join(('barcode', 'aligned_correct',
                                      'not_aligned', 'aligned_incorect',
                                      'realigned', 'realigned_correct')))
            f_length.write('\n')
            line = '{}\t{}\t{}\t{}\t{}\t{}\n'
            for bc in map_stat_dict:
                stat_list = map_stat_dict[bc]
                total = sum(stat_list)
                aligned = total - stat_list[1]
                if stat_list[3] > 0:
                    stat_list.append(remap_sam_out[2][bc][0])
                else:
                    stat_list.append(0)
                f_stats.write(line.format(bc, total, aligned, stat_list[0],
                                          stat_list[3], stat_list[4]))

                length_list = length_dict[bc]
                if length_list[3] > 0:
                    length_list.append(remap_sam_out[3][bc][0])
                else:
                    length_list.append(0)
                avg_length = [length_list[i] / stat_list[i]
                              if stat_list[i] > 0 else 0
                              for i in range(0, len(length_list))]
                f_length.write('\t'.join([bc] + [str(round(l, 2))
                                                 for l in avg_length]))
                f_length.write('\n')
