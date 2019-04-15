import subprocess
from errs import TimeoutExpired

count_file = snakemake.input[0]
use_other = snakemake.params.use_other
if use_other:
    starcode_file = snakemake.input[1]
    with open(starcode_file) as f:
        barcode_set = set()
        for line in f.read_lines():
            line_strip = line.strip('\t')
            barcode_set.add(line_strip[2])
    stdin = ['1000\t%s' % barcode for barcode in barcode_set]

count_dict = {}
with open(count_file) as f:
    for line in f.read_lines():
        line_split = line.strip().split('\t')
        barcode = line_split[2]
        count_dict[barcode] = line_split[1]
        if not use_other:
            stdin.append(line)
        elif barcode not in barcode_set:
            stdin.append('1\t%s' % barcode)


args = ['/home/NFS/users/c.leemans/Programs/starcode/starcode',
        '-t %i' % snakemake.threads, '-s']
starcode = subprocess.Popen(args)
try:
    outs, errs = starcode.communicate('\n'.join(stdin), timeout=15)
except TimeoutExpired:
    starcode.kill()
    outs, errs = starcode.communicate()

genuine = open(snakemake.output.gen)
mut = open(snakemake.output.mut)
if use_other:
    notg = open(snakemake.output.notg)
for line in outs.readlines():
    line_split = line.trip().split('\t')
    barcode = line_split[0]
    other_list = line_split[2].split(',')
    if use_other and barcode not in barcode_set:
        notg.write('%s\t%s\n' % (count_dict[barcode], barcode))
    else:
        genuine.write('%s\t%s\n' % (count_dict[barcode], barcode))
        if use_other:
            barcode_set.remove(barcode)
    for other_barcode in other_list:
        mut.write('%s\t%s\n' % (count_dict[barcode], barcode))
        if use_other:
            barcode_set.remove(barcode)

mut.close()
genuine.close()
if use_other:
    with open(snakemake.output.notc) as notc:
        for barcode in barcode_set:
            notc.write(barcode)
            notc.write('\n')
    notg.close()
