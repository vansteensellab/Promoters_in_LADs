fastq-dump SRR1552480

nice -19 fastx_trimmer -t 20 -i SRR1552480.fastq | \
    bowtie2 -p 15 -U - -x /DATA/usr/c.leemans/data/bowtie2/human_rdna_U13369.1 | \
    samtools view -Sbf 4 | \
    bamToFastq -i /dev/stdin -fq /dev/stdout | \
    bowtie2 -p 15 -U - -x /DATA/data/human/hg38/index/bowtie2/hg38_1-22_XYM | \
    samtools view -Sb | \
    samtools sort -o grocap_wTAP_K562_hg38_SRR1552480.bam

samtools index grocap_wTAP_K562_hg38_SRR1552480.bam



bamCoverage --minMappingQuality 30 --Offset 1 -bs 1 --filterRNAstrand reverse -b grocap_wTAP_K562_hg38_SRR1552480.bam -o grocap_wTAP_plus_K562_hg38_SRR1552480.bw
bamCoverage --minMappingQuality 30 --Offset 1 -bs 1 --filterRNAstrand forward -b grocap_wTAP_K562_hg38_SRR1552480.bam -o grocap_wTAP_min_K562_hg38_SRR1552480.bw
