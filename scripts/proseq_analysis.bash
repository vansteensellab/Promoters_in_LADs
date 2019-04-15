fastq-dump SRR1554311 SRR1554312

nice -19 fastp -a TGGAATTCTCGGGTGCCAAGG \
               -o SRR1554311_fastp.fastq -i SRR1554311.fastq

zcat SRR1554311.fastq.gz | \
    fastx_trimmer -m 15 -i /dev/stdin | \
    cutadapt -m 15 -a TGGAATTCTCGGGTGCCAAGG -g GTTCAGAGTTCTACAGTCCGACGATC \
             -o SRR1554311_cutadapt.fastq - &

nice -19 bowtie2 -p 10 --very-sensitive -U SRR1554311_cutadapt.fastq \
                 -x /DATA/usr/c.leemans/data/bowtie2/human_rdna_U13369.1 | \
    samtools view -Sbf 4 | \
    bamToFastq -i /dev/stdin -fq /dev/stdout | \
    bowtie2 -p 10 --very-sensitive -U - -x /DATA/usr/c.leemans/data/hg38/hg38_1-22_XYM | \
    sambamba view -S -f bam -F "mapping_quality >= 42 and [XM] < 3" /dev/stdin | \
    sambamba sort -o proseq_rep1_K562_hg38_SRR1554311.bam /dev/stdin &

sambamba index -t 20 proseq_rep1_K562_hg38_SRR1554311.bam &

nice -19 fastp -a TGGAATTCTCGGGTGCCAAGG \
               -o SRR1554312_fastp.fastq -i SRR1554312.fastq

zcat SRR1554312.fastq.gz | \
    fastx_trimmer -m 15 -i /dev/stdin | \
    cutadapt -m 15 -a TGGAATTCTCGGGTGCCAAGG -g GTTCAGAGTTCTACAGTCCGACGATC \
             -o SRR1554312_cutadapt.fastq -

nice -19 bowtie2 --very-sensitive -p 10 -U SRR1554312_cutadapt.fastq \
                 -x /DATA/usr/c.leemans/data/bowtie2/human_rdna_U13369.1 | \
    samtools view -Sbf 4 | \
    bamToFastq -i /dev/stdin -fq /dev/stdout | \
    bowtie2 --very-sensitive -p 10 -U - -x /DATA/data/human/hg38/index/bowtie2/hg38_1-22_XYM | \
    sambamba view -S -f bam -F "mapping_quality >= 42 and [XM] < 3" /dev/stdin | \
    sambamba sort -o proseq_rep2_K562_hg38_SRR1554312.bam /dev/stdin

sambamba index -t 20 proseq_rep2_K562_hg38_SRR1554312.bam &


bamCoverage -p 20 --Offset 1 -bs 1 --filterRNAstrand reverse -b proseq_rep1_K562_hg38_SRR1554311.bam -o proseq_rep1_plus_K562_hg38_SRR1554311.bw
bamCoverage -p 20 --Offset 1 -bs 1 --filterRNAstrand forward -b proseq_rep1_K562_hg38_SRR1554311.bam -o proseq_rep1_min_K562_hg38_SRR1554311.bw

bamCoverage -p 20 --Offset 1 -bs 1 --filterRNAstrand reverse -b proseq_rep2_K562_hg38_SRR1554312.bam -o proseq_rep2_plus_K562_hg38_SRR1554312.bw
bamCoverage -p 20 --Offset 1 -bs 1 --filterRNAstrand forward -b proseq_rep2_K562_hg38_SRR1554312.bam -o proseq_rep2_min_K562_hg38_SRR1554312.bw
