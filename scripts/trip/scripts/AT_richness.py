import sys
import subprocess

if __name__ == '__main__':
    bed = sys.argv[1]
    genome = sys.argv[2]
    command = 'bedtools getfasta -bed %s -fi %s -fo /dev/stdout' % (bed, genome)
    bedtools = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    outs, errs = bedtools.communicate()
    fa_dict = {}
    for line in outs.decode('UTF-8').split('\n'):
        if line.startswith('>'):
            region = line[1:]
            fa_dict[region] = ''
        else:
            fa_dict[region] = ''.join((fa_dict[region], line))
    with open(bed) as bed_file:
        sys.stdout.write('barcode\tAT_ratio\tcount\ttotal\n')
        for line in bed_file.readlines():
            line_split = line.strip().split()
            region = '%s:%s-%s' % tuple(line_split[:3])
            sequence = fa_dict[region]
            AT = sum(sequence.count(n) for n in ('T', 'A', 't', 'a'))
            CG = sum(sequence.count(n) for n in ('C', 'G', 'c', 'g'))
            if (AT+CG) > 0:
                ratio = AT/(AT + CG)
            else:
                ratio = 'NaN'

            sys.stdout.write(line_split[3])
            if ratio != 'NaN':
                sys.stdout.write('\t%f\t' % ratio)
            else:
                sys.stdout.write('\tNaN\t')
            sys.stdout.write('\t'.join(line_split[4:]))
            sys.stdout.write('\n')
