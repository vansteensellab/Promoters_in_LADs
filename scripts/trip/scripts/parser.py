read_type = snakemake.params.read_type
structure = snakemake.params.structure
name = snakemake.params.name
fastq_name = snakemake.input
log = snakemake.log
outdir = snakemake.params.outdir

structure_file = '%s/%s.structure.txt' % (outdir, name)

with open(structure_file, 'w') as f:
    f.write(structure)


if type(fastq_name) != str and len(fastq_name)==2:
    paired = '-p %s' % fastq_name[1]
    first = fastq_name[0]
else:
    paired = ''

shell('~t.v.schaik/modules/read-parsing/read_parser.py -r -l %s %s '
      '-b %s %s %s %s' % (log, paired, name, fastq_name, structure_file,
                          outdir))
