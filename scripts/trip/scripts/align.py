options = ' '.join(params[0]['options'][params[1]])
index = params[0]['index']
gunzip = "gunzip -c {input}"
## filter for read length
awk = ("awk '{{"
       "       step=NR%4;"
       "       if (step==0 && length(a[2])>6){{"
       "           for (i in a){{"
       "               print a[i]"
       "           }}"
       "           print $0"
       "       }} else if (step!=0){{"
       "           a[step]=$0;"
       "       }}"
       "}}'")
bowtie = 'bowtie2 -p {threads} %s -x %s -U - > {output}' % (options,
                                                            index)
flagstat = 'samtools flagstat {output} > {log}'
print(bowtie)
shell('%s | %s | %s; %s' % (gunzip, awk, bowtie, flagstat))