
######## Snakemake header ########
import sys; sys.path.insert(0, "/home/c.leemans/mydata/miniconda3/envs/trip/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X.\x00\x00\x00cl20170929_TTRIP_K562/mapping/pool2_rep2.2.bamq\x06X7\x00\x00\x00cl20170929_TTRIP_K562/mapping/pool2_rep2.starcode.countq\x07e}q\x08(X\x06\x00\x00\x00_namesq\t}q\n(X\x03\x00\x00\x00bamq\x0bK\x00K\x01\x86q\x0cX\x05\x00\x00\x00countq\rK\x01K\x02\x86q\x0euh\x0bcsnakemake.io\nNamedlist\nq\x0f)\x81q\x10h\x06a}q\x11h\t}q\x12sbh\rh\x0f)\x81q\x13h\x07a}q\x14h\t}q\x15sbubX\x06\x00\x00\x00outputq\x16csnakemake.io\nOutputFiles\nq\x17)\x81q\x18(X.\x00\x00\x00cl20170929_TTRIP_K562/mapping/pool2_rep2.2.bedq\x19X0\x00\x00\x00cl20170929_TTRIP_K562/mapping/pool2_rep2.2.tableq\x1aX;\x00\x00\x00cl20170929_TTRIP_K562/mapping/pool2_rep2.2.parse_stat.tableq\x1bX7\x00\x00\x00cl20170929_TTRIP_K562/mapping/pool2_rep2.2.length.tableq\x1cX9\x00\x00\x00cl20170929_TTRIP_K562/mapping/pool2_rep2.2.remap.fastq.gzq\x1dX4\x00\x00\x00cl20170929_TTRIP_K562/mapping/pool2_rep2.2.remap.bamq\x1ee}q\x1f(h\t}q (X\x03\x00\x00\x00bedq!K\x00N\x86q"X\x05\x00\x00\x00tableq#K\x01N\x86q$X\x05\x00\x00\x00statsq%K\x02N\x86q&X\x06\x00\x00\x00lengthq\'K\x03N\x86q(X\x08\x00\x00\x00remap_fqq)K\x04N\x86q*X\x05\x00\x00\x00remapq+K\x05N\x86q,uh!h\x19h#h\x1ah%h\x1bh\'h\x1ch)h\x1dh+h\x1eubX\x06\x00\x00\x00paramsq-csnakemake.io\nParams\nq.)\x81q/(X"\x00\x00\x00/DATA/data/bowtie2/hg19_ch1-22_XYMq0ccollections\nOrderedDict\nq1)Rq2(X\x01\x00\x00\x001q3]q4X\x10\x00\x00\x00--very-sensitiveq5aX\x01\x00\x00\x002q6]q7X\x16\x00\x00\x00--very-sensitive-localq8auh1)Rq9(h3M\xf4\x01h6K\x14uh6e}q:(h\t}q;(X\x0c\x00\x00\x00bowtie_indexq<K\x00N\x86q=X\x07\x00\x00\x00optionsq>K\x01N\x86q?X\x08\x00\x00\x00max_distq@K\x02N\x86qAX\x03\x00\x00\x00numqBK\x03N\x86qCuh<h0h>h2h@h9hBh6ubX\t\x00\x00\x00wildcardsqDcsnakemake.io\nWildcards\nqE)\x81qF(X\x15\x00\x00\x00cl20170929_TTRIP_K562qGX\n\x00\x00\x00pool2_rep2qHh6e}qI(h\t}qJ(X\x06\x00\x00\x00outdirqKK\x00N\x86qLX\x04\x00\x00\x00nameqMK\x01N\x86qNX\x03\x00\x00\x00numqOK\x02N\x86qPuX\x06\x00\x00\x00outdirqQhGX\x04\x00\x00\x00nameqRhHhBh6ubX\x07\x00\x00\x00threadsqSK\nX\t\x00\x00\x00resourcesqTcsnakemake.io\nResources\nqU)\x81qV(K\nK\x01e}qW(h\t}qX(X\x06\x00\x00\x00_coresqYK\x00N\x86qZX\x06\x00\x00\x00_nodesq[K\x01N\x86q\\uhYK\nh[K\x01ubX\x03\x00\x00\x00logq]csnakemake.io\nLog\nq^)\x81q_}q`h\t}qasbX\x06\x00\x00\x00configqb}qc(X\x06\x00\x00\x00outdirqdX\x15\x00\x00\x00cl20170929_TTRIP_K562qeX\x06\x00\x00\x00groupsqf]qg(]qh(X\x03\x00\x00\x00POIqiX\x03\x00\x00\x00G9aqjX\x04\x00\x00\x00KRABqkX\x03\x00\x00\x00CBXqle]qm(X\t\x00\x00\x00conditionqnX\x08\x00\x00\x00GAL4-POIqoX\x03\x00\x00\x00POIqpX\x04\x00\x00\x00GAL4qqe]qr(X\x03\x00\x00\x00dayqsK\x02K\tK\x0bK\x0cK\x0ee]qt(X\t\x00\x00\x00replicatequK\x01K\x02eeX\n\x00\x00\x00input_fileqvh1)Rqw(X\x04\x00\x00\x00gDNAqxh1)Rqy(X\x0c\x00\x00\x00G9a_GAL4_2_1qz]q{(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_13.fqq|K\neX\x10\x00\x00\x00G9a_GAL4-POI_2_1q}]q~(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_14.fqq\x7fK\neX\x0b\x00\x00\x00G9a_POI_2_1q\x80]q\x81(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_15.fqq\x82K\neX\x0c\x00\x00\x00G9a_GAL4_2_2q\x83]q\x84(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_16.fqq\x85K\neX\x10\x00\x00\x00G9a_GAL4-POI_2_2q\x86]q\x87(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_17.fqq\x88K\neX\x0b\x00\x00\x00G9a_POI_2_2q\x89]q\x8a(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_18.fqq\x8bK\neX\r\x00\x00\x00G9a_GAL4_12_1q\x8c]q\x8d(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_19.fqq\x8eK\neX\x11\x00\x00\x00G9a_GAL4-POI_12_1q\x8f]q\x90(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_20.fqq\x91K\neX\x0c\x00\x00\x00G9a_POI_12_1q\x92]q\x93(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_21.fqq\x94K\neX\r\x00\x00\x00G9a_GAL4_12_2q\x95]q\x96(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_22.fqq\x97K\neX\x11\x00\x00\x00G9a_GAL4-POI_12_2q\x98]q\x99(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_23.fqq\x9aK\neX\x0c\x00\x00\x00G9a_POI_12_2q\x9b]q\x9c(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_24.fqq\x9dK\neX\x10\x00\x00\x00CBX_GAL4-POI_2_1q\x9e]q\x9f(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_13.fqq\xa0K\x00eX\x0c\x00\x00\x00CBX_GAL4_2_1q\xa1]q\xa2(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_14.fqq\xa3K\x00eX\x0b\x00\x00\x00CBX_POI_2_1q\xa4]q\xa5(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_15.fqq\xa6K\x00eX\x10\x00\x00\x00CBX_GAL4-POI_2_2q\xa7]q\xa8(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_16.fqq\xa9K\x00eX\x0c\x00\x00\x00CBX_GAL4_2_2q\xaa]q\xab(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_17.fqq\xacK\x00eX\x0b\x00\x00\x00CBX_POI_2_2q\xad]q\xae(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_18.fqq\xafK\x00eX\x10\x00\x00\x00CBX_GAL4-POI_9_1q\xb0]q\xb1(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_19.fqq\xb2K\x00eX\x0c\x00\x00\x00CBX_GAL4_9_1q\xb3]q\xb4(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_20.fqq\xb5K\x00eX\x0b\x00\x00\x00CBX_POI_9_1q\xb6]q\xb7(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_21.fqq\xb8K\x00eX\x10\x00\x00\x00CBX_GAL4-POI_9_2q\xb9]q\xba(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_22.fqq\xbbK\x00eX\x0c\x00\x00\x00CBX_GAL4_9_2q\xbc]q\xbd(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_23.fqq\xbeK\x00eX\x0b\x00\x00\x00CBX_POI_9_2q\xbf]q\xc0(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_24.fqq\xc1K\x00eX\x11\x00\x00\x00CBX_GAL4-POI_12_1q\xc2]q\xc3(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_48_r2_D12_gDNA.fqq\xc4K\neX\x11\x00\x00\x00CBX_GAL4-POI_12_2q\xc5]q\xc6(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_48_r3_D12_gDNA.fqq\xc7K\neX\r\x00\x00\x00CBX_GAL4_12_1q\xc8]q\xc9(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_50_r2_D12_gDNA.fqq\xcaK\neX\r\x00\x00\x00CBX_GAL4_12_2q\xcb]q\xcc(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_50_r3_D12_gDNA.fqq\xcdK\neX\x0c\x00\x00\x00CBX_POI_12_1q\xce]q\xcf(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_51_r2_D12_gDNA.fqq\xd0K\neX\x0c\x00\x00\x00CBX_POI_12_2q\xd1]q\xd2(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_51_r3_D12_gDNA.fqq\xd3K\neX\r\x00\x00\x00KRAB_GAL4_2_1q\xd4]q\xd5(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_13.fqq\xd6K\neX\x11\x00\x00\x00KRAB_GAL4-POI_2_1q\xd7]q\xd8(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_14.fqq\xd9K\neX\x0c\x00\x00\x00KRAB_POI_2_1q\xda]q\xdb(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_15.fqq\xdcK\neX\r\x00\x00\x00KRAB_GAL4_2_2q\xdd]q\xde(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_16.fqq\xdfK\neX\x11\x00\x00\x00KRAB_GAL4-POI_2_2q\xe0]q\xe1(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_17.fqq\xe2K\neX\x0e\x00\x00\x00KRAB_GAL4_11_1q\xe3]q\xe4(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_19.fqq\xe5K\neX\x0c\x00\x00\x00KRAB_POI_2_2q\xe6]q\xe7(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_18.fqq\xe8K\neX\x12\x00\x00\x00KRAB_GAL4-POI_11_1q\xe9]q\xea(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_20.fqq\xebK\neX\r\x00\x00\x00KRAB_POI_11_1q\xec]q\xed(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_21.fqq\xeeK\neX\x0e\x00\x00\x00KRAB_GAL4_11_2q\xef]q\xf0(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_22.fqq\xf1K\neX\x12\x00\x00\x00KRAB_GAL4-POI_11_2q\xf2]q\xf3(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_23.fqq\xf4K\neX\r\x00\x00\x00KRAB_POI_11_2q\xf5]q\xf6(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_24.fqq\xf7K\neX\x0e\x00\x00\x00KRAB_GAL4_14_1q\xf8]q\xf9(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_50_r1_D14_gDNA.fqq\xfaK\neX\x0e\x00\x00\x00KRAB_GAL4_14_2q\xfb]q\xfc(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_50_r2_D14_gDNA.fqq\xfdK\neX\x12\x00\x00\x00KRAB_GAL4-POI_14_1q\xfe]q\xff(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_74_r1_D14_gDNA.fqr\x00\x01\x00\x00K\neX\x12\x00\x00\x00KRAB_GAL4-POI_14_2r\x01\x01\x00\x00]r\x02\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_74_r2_D14_gDNA.fqr\x03\x01\x00\x00K\neX\r\x00\x00\x00KRAB_POI_14_1r\x04\x01\x00\x00]r\x05\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_76_r1_D14_gDNA.fqr\x06\x01\x00\x00K\neX\r\x00\x00\x00KRAB_POI_14_2r\x07\x01\x00\x00]r\x08\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_76_r2_D14_gDNA.fqr\t\x01\x00\x00K\neuX\x04\x00\x00\x00cDNAr\n\x01\x00\x00h1)Rr\x0b\x01\x00\x00(X\x0c\x00\x00\x00G9a_GAL4_2_1r\x0c\x01\x00\x00]r\r\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_1.fqr\x0e\x01\x00\x00K\neX\x10\x00\x00\x00G9a_GAL4-POI_2_1r\x0f\x01\x00\x00]r\x10\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_2.fqr\x11\x01\x00\x00K\neX\x0b\x00\x00\x00G9a_POI_2_1r\x12\x01\x00\x00]r\x13\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_3.fqr\x14\x01\x00\x00K\neX\x0c\x00\x00\x00G9a_GAL4_2_2r\x15\x01\x00\x00]r\x16\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_4.fqr\x17\x01\x00\x00K\neX\x10\x00\x00\x00G9a_GAL4-POI_2_2r\x18\x01\x00\x00]r\x19\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_5.fqr\x1a\x01\x00\x00K\neX\x0b\x00\x00\x00G9a_POI_2_2r\x1b\x01\x00\x00]r\x1c\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_6.fqr\x1d\x01\x00\x00K\neX\r\x00\x00\x00G9a_GAL4_12_1r\x1e\x01\x00\x00]r\x1f\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_7.fqr \x01\x00\x00K\neX\x11\x00\x00\x00G9a_GAL4-POI_12_1r!\x01\x00\x00]r"\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_8.fqr#\x01\x00\x00K\neX\x0c\x00\x00\x00G9a_POI_12_1r$\x01\x00\x00]r%\x01\x00\x00(Xe\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_9.fqr&\x01\x00\x00K\neX\r\x00\x00\x00G9a_GAL4_12_2r\'\x01\x00\x00]r(\x01\x00\x00(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_10.fqr)\x01\x00\x00K\neX\x11\x00\x00\x00G9a_GAL4-POI_12_2r*\x01\x00\x00]r+\x01\x00\x00(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_11.fqr,\x01\x00\x00K\neX\x0c\x00\x00\x00G9a_POI_12_2r-\x01\x00\x00]r.\x01\x00\x00(Xf\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160406_fastqs_G9a/3893_1_BarcodedPool_NoIndex_TRIP_K562_G9a_12.fqr/\x01\x00\x00K\neX\x10\x00\x00\x00CBX_GAL4-POI_2_1r0\x01\x00\x00]r1\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_1.fqr2\x01\x00\x00K\x00eX\x0c\x00\x00\x00CBX_GAL4_2_1r3\x01\x00\x00]r4\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_2.fqr5\x01\x00\x00K\x00eX\x0b\x00\x00\x00CBX_POI_2_1r6\x01\x00\x00]r7\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_3.fqr8\x01\x00\x00K\x00eX\x10\x00\x00\x00CBX_GAL4-POI_2_2r9\x01\x00\x00]r:\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_4.fqr;\x01\x00\x00K\x00eX\x0c\x00\x00\x00CBX_GAL4_2_2r<\x01\x00\x00]r=\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_5.fqr>\x01\x00\x00K\x00eX\x0b\x00\x00\x00CBX_POI_2_2r?\x01\x00\x00]r@\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_6.fqrA\x01\x00\x00K\x00eX\x10\x00\x00\x00CBX_GAL4-POI_9_1rB\x01\x00\x00]rC\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_7.fqrD\x01\x00\x00K\x00eX\x0c\x00\x00\x00CBX_GAL4_9_1rE\x01\x00\x00]rF\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_8.fqrG\x01\x00\x00K\x00eX\x0b\x00\x00\x00CBX_POI_9_1rH\x01\x00\x00]rI\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_9.fqrJ\x01\x00\x00K\x00eX\x10\x00\x00\x00CBX_GAL4-POI_9_2rK\x01\x00\x00]rL\x01\x00\x00(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_10.fqrM\x01\x00\x00K\x00eX\x0c\x00\x00\x00CBX_GAL4_9_2rN\x01\x00\x00]rO\x01\x00\x00(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_11.fqrP\x01\x00\x00K\x00eX\x0b\x00\x00\x00CBX_POI_9_2rQ\x01\x00\x00]rR\x01\x00\x00(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_12.fqrS\x01\x00\x00K\x00eX\x11\x00\x00\x00CBX_GAL4-POI_12_1rT\x01\x00\x00]rU\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_48_r2_D12_cDNA.fqrV\x01\x00\x00K\neX\x11\x00\x00\x00CBX_GAL4-POI_12_2rW\x01\x00\x00]rX\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_48_r3_D12_cDNA.fqrY\x01\x00\x00K\neX\r\x00\x00\x00CBX_GAL4_12_1rZ\x01\x00\x00]r[\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_50_r2_D12_cDNA.fqr\\\x01\x00\x00K\neX\r\x00\x00\x00CBX_GAL4_12_2r]\x01\x00\x00]r^\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_50_r3_D12_cDNA.fqr_\x01\x00\x00K\neX\x0c\x00\x00\x00CBX_POI_12_1r`\x01\x00\x00]ra\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_51_r2_D12_cDNA.fqrb\x01\x00\x00K\neX\x0c\x00\x00\x00CBX_POI_12_2rc\x01\x00\x00]rd\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_51_r3_D12_cDNA.fqre\x01\x00\x00K\neX\r\x00\x00\x00KRAB_GAL4_2_1rf\x01\x00\x00]rg\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_1.fqrh\x01\x00\x00K\neX\x11\x00\x00\x00KRAB_GAL4-POI_2_1ri\x01\x00\x00]rj\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_2.fqrk\x01\x00\x00K\neX\x0c\x00\x00\x00KRAB_POI_2_1rl\x01\x00\x00]rm\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_3.fqrn\x01\x00\x00K\neX\r\x00\x00\x00KRAB_GAL4_2_2ro\x01\x00\x00]rp\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_4.fqrq\x01\x00\x00K\neX\x11\x00\x00\x00KRAB_GAL4-POI_2_2rr\x01\x00\x00]rs\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_5.fqrt\x01\x00\x00K\neX\x0c\x00\x00\x00KRAB_POI_2_2ru\x01\x00\x00]rv\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_6.fqrw\x01\x00\x00K\neX\x0e\x00\x00\x00KRAB_GAL4_11_1rx\x01\x00\x00]ry\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_7.fqrz\x01\x00\x00K\neX\x12\x00\x00\x00KRAB_GAL4-POI_11_1r{\x01\x00\x00]r|\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_8.fqr}\x01\x00\x00K\neX\r\x00\x00\x00KRAB_POI_11_1r~\x01\x00\x00]r\x7f\x01\x00\x00(Xl\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_9.fqr\x80\x01\x00\x00K\neX\x0e\x00\x00\x00KRAB_GAL4_11_2r\x81\x01\x00\x00]r\x82\x01\x00\x00(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_10.fqr\x83\x01\x00\x00K\neX\x12\x00\x00\x00KRAB_GAL4-POI_11_2r\x84\x01\x00\x00]r\x85\x01\x00\x00(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_11.fqr\x86\x01\x00\x00K\neX\r\x00\x00\x00KRAB_POI_11_2r\x87\x01\x00\x00]r\x88\x01\x00\x00(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160331_fastqs_TRIP_KRAB/3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_12.fqr\x89\x01\x00\x00K\neX\x0e\x00\x00\x00KRAB_GAL4_14_1r\x8a\x01\x00\x00]r\x8b\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_50_r1_D14_cDNA.fqr\x8c\x01\x00\x00K\neX\x0e\x00\x00\x00KRAB_GAL4_14_2r\x8d\x01\x00\x00]r\x8e\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_50_r2_D14_cDNA.fqr\x8f\x01\x00\x00K\neX\x12\x00\x00\x00KRAB_GAL4-POI_14_1r\x90\x01\x00\x00]r\x91\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_74_r1_D14_cDNA.fqr\x92\x01\x00\x00K\neX\x12\x00\x00\x00KRAB_GAL4-POI_14_2r\x93\x01\x00\x00]r\x94\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_74_r2_D14_cDNA.fqr\x95\x01\x00\x00K\neX\r\x00\x00\x00KRAB_POI_14_1r\x96\x01\x00\x00]r\x97\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_76_r1_D14_cDNA.fqr\x98\x01\x00\x00K\neX\r\x00\x00\x00KRAB_POI_14_2r\x99\x01\x00\x00]r\x9a\x01\x00\x00(Xr\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160815_fastqs_CBX5_D12_KRAB_D14/4088_1_BarcodedPool_NoIndex_76_r2_D14_cDNA.fqr\x9b\x01\x00\x00K\neuX\x07\x00\x00\x00mappingr\x9c\x01\x00\x00h1)Rr\x9d\x01\x00\x00(X\n\x00\x00\x00pool1_rep1r\x9e\x01\x00\x00]r\x9f\x01\x00\x00(]r\xa0\x01\x00\x00(X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R1_001_smplIdx_09.fastq.gzr\xa1\x01\x00\x00X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R2_001_smplIdx_09.fastq.gzr\xa2\x01\x00\x00eK\neX\n\x00\x00\x00pool1_rep2r\xa3\x01\x00\x00]r\xa4\x01\x00\x00(]r\xa5\x01\x00\x00(X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R1_001_smplIdx_10.fastq.gzr\xa6\x01\x00\x00X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R2_001_smplIdx_10.fastq.gzr\xa7\x01\x00\x00eK\neX\n\x00\x00\x00pool1_rep3r\xa8\x01\x00\x00]r\xa9\x01\x00\x00(]r\xaa\x01\x00\x00(X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R1_001_smplIdx_11.fastq.gzr\xab\x01\x00\x00X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R2_001_smplIdx_11.fastq.gzr\xac\x01\x00\x00eK\neX\n\x00\x00\x00pool2_rep1r\xad\x01\x00\x00]r\xae\x01\x00\x00(]r\xaf\x01\x00\x00(X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R1_001_smplIdx_12.fastq.gzr\xb0\x01\x00\x00X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R2_001_smplIdx_12.fastq.gzr\xb1\x01\x00\x00eK\neX\n\x00\x00\x00pool2_rep2r\xb2\x01\x00\x00]r\xb3\x01\x00\x00(]r\xb4\x01\x00\x00(X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R1_001_smplIdx_13.fastq.gzr\xb5\x01\x00\x00X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R2_001_smplIdx_13.fastq.gzr\xb6\x01\x00\x00eK\neX\n\x00\x00\x00pool2_rep3r\xb7\x01\x00\x00]r\xb8\x01\x00\x00(]r\xb9\x01\x00\x00(X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R1_001_smplIdx_14.fastq.gzr\xba\x01\x00\x00X@\x00\x00\x00raw_data/3354_1_iPCR_laura_eva_altndx_R2_001_smplIdx_14.fastq.gzr\xbb\x01\x00\x00eK\neuX\x05\x00\x00\x00spiker\xbc\x01\x00\x00]r\xbd\x01\x00\x00(Xm\x00\x00\x00/DATA/usr/l.brueckner/TTRIP_K562/lb20160318_fastqs_TRIP_CBX5/3870_1_BarcodedPool_NoIndex_TRIP_K562_CBX5_25.fqr\xbe\x01\x00\x00K\x00euX\x06\x00\x00\x00configr\xbf\x01\x00\x00X \x00\x00\x00cl20160816_config_K562_TTRIP.txtr\xc0\x01\x00\x00X\t\x00\x00\x00intersectr\xc1\x01\x00\x00h1)Rr\xc2\x01\x00\x00(X\x06\x00\x00\x00repeatr\xc3\x01\x00\x00X\x14\x00\x00\x00{outdir}/repeats.bedr\xc4\x01\x00\x00X\x03\x00\x00\x00ladr\xc5\x01\x00\x00XF\x00\x00\x00/DATA/usr/c.leemans/data/carolineLADs/LAD_K562_continuous_cl160714.bedr\xc6\x01\x00\x00X\x05\x00\x00\x00chromr\xc7\x01\x00\x00X@\x00\x00\x00/DATA/usr/c.leemans/data/tracks/hg19/wgEncodeBroadHmmK562HMM.bedr\xc8\x01\x00\x00uX\x07\x00\x00\x00extractr\xc9\x01\x00\x00h1)Rr\xca\x01\x00\x00X\x06\x00\x00\x00timingr\xcb\x01\x00\x00Xo\x00\x00\x00/DATA/usr/c.leemans/data/tracks/hg19/GSM923448/GSM923448_hg19_wgEncodeUwRepliSeqK562{state}PctSignalRep1.bigWigr\xcc\x01\x00\x00sX\x07\x00\x00\x00nearestr\xcd\x01\x00\x00h1)Rr\xce\x01\x00\x00X\x03\x00\x00\x00cpgr\xcf\x01\x00\x00XD\x00\x00\x00/DATA/usr/c.leemans/data/tracks/hg19/cpgIslandExtUnmasked_140601.bedr\xd0\x01\x00\x00sX\x0c\x00\x00\x00repeatMaskerr\xd1\x01\x00\x00XJ\x00\x00\x00/DATA/usr/c.leemans/data/tracks/hg19/repeatMasker_hg19_fa_out_20140131.bedr\xd2\x01\x00\x00X\t\x00\x00\x00structurer\xd3\x01\x00\x00h1)Rr\xd4\x01\x00\x00(X\x07\x00\x00\x00mappingr\xd5\x01\x00\x00X[\x02\x00\x00ID      5\'      3\'      type    req     second-read     pos   keep-bases\nindex   %i      -       const   present False   fixed   -\nmap_pat1a       GTCACAAGGGCCGGCCACAAC   -       const   present False   fixed   -\nconst_bar       TCGAG\\{16\\}TGATC  -       const_bar       present False   fixed   -\nrev_map_complement      -       TTAACCCTAGAAAGATAATCATATTGTGACGTAC      const   -       False   var   -\nrev_map GTACGTCACAATATGATTATCTTTCTAGGGTTAA      -       const   present True    fixed   -\nfwd_map_complement      -       GATCA[BC]CTCGAGTTGTGGCCGGCCCTTGTGAC     const_bar_comp  -       True    var   -\nr\xd6\x01\x00\x00X\x04\x00\x00\x00gDNAr\xd7\x01\x00\x00Xx\x01\x00\x00ID      5\'      3\'      type    req     second-read     pos                   keep-bases\nindex   %i      -       const   present False   fixed                            -\npat1    GTCACAAGGGCCGGCCACAACTCGAG      -       const   present False   fixed   -\nbarcode 16      -       barcode present False   fixed   -\npat2    TGATCCTGCAGTG   -       const   present False   var   -\nr\xd8\x01\x00\x00X\x04\x00\x00\x00cDNAr\xd9\x01\x00\x00XQ\x01\x00\x00ID      5\'      3\'      type    req     second-read     pos     keep-bases\nindex   %i      -       const   present False   fixed   -\npat1    GTCACAAGGGCCGGCCACAACTCGAG      -       const   present False   fixed   -\nbarcode 16      -       barcode present False   fixed   -\npat2    TGATCCTGCAGTG   -       const   present False   var   -\nr\xda\x01\x00\x00X\x05\x00\x00\x00spiker\xdb\x01\x00\x00XI\x01\x00\x00ID      5\'      3\'      type    req     second-read     pos     keep-bases\nindex   %i      -       const   present False   fixed   -\npat1    GTCACAAGGGCCGGCCACAA    -       const   present False   fixed   -\nbarcode 16      -       barcode present False   fixed   -\npat2    GATCGGTACCCA    -       const   present False   var   -\nr\xdc\x01\x00\x00uX\x08\x00\x00\x00lev_distr\xdd\x01\x00\x00K\x02X\t\x00\x00\x00min_countr\xde\x01\x00\x00h1)Rr\xdf\x01\x00\x00(X\x05\x00\x00\x00spiker\xe0\x01\x00\x00M\xe8\x03X\x04\x00\x00\x00gDNAr\xe1\x01\x00\x00K\x05X\x04\x00\x00\x00cDNAr\xe2\x01\x00\x00K\x00X\x03\x00\x00\x00mapr\xe3\x01\x00\x00K\x03uX\x06\x00\x00\x00bowtier\xe4\x01\x00\x00h1)Rr\xe5\x01\x00\x00(X\x05\x00\x00\x00indexr\xe6\x01\x00\x00h0X\x07\x00\x00\x00optionsr\xe7\x01\x00\x00h2uX\x08\x00\x00\x00max_distr\xe8\x01\x00\x00h9uX\x04\x00\x00\x00ruler\xe9\x01\x00\x00X\t\x00\x00\x00parse_samr\xea\x01\x00\x00ub.'); from snakemake.logging import logger; logger.printshellcmds = False
######## Original script #########
import pysam
import gzip
import os
import re
from Bio import SeqIO
from collections import Counter


def parse_sam(sam_file, starcode_set, max_soft_clip=5, min_first_match=10,
              remap_soft_clip=17):
    remap_list = []
    map_dict = {}
    map_stat_dict = {}
    length_dict = {}
    for line in pysam.AlignmentFile(sam_file):
        bc_this = re.match(r'.*[:_]([A-Z]+)$', line.query_name).groups(1)[0]
        if bc_this in starcode_set:
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
            av_mapq = mapq_sum / total_reads
            freq1 = 1.0
            freq2 = 0.0
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
            av_mapq = top_mapq / top_reads
            freq1 = top_reads / total_reads
            if len(sorted_key_list) > 0:
                second_key, \
                    second_reads, \
                    second_mapq, \
                    sorted_key_list = refine_map(this_map_dict,
                                                 sorted_key_list,
                                                 max_dist)
                freq2 = second_reads / total_reads
            else:
                freq2 = 0
            sorted_key_list = sorted(this_map_dict.keys(),
                                     key=lambda elem: sort_keys(elem,
                                                                this_map_dict),
                                     reverse=True)
        map_dict_out[bc] = [reference_name, ori, start_pos, total_reads,
                            av_mapq, freq1, freq2, seq]
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

    for count_file in snakein.count:
        with open(count_file) as cf:
            for line in cf.readlines():
                barcode = line.split('\t')[0]
                if barcode not in starcode_set:
                    starcode_set.add(barcode)
    (map_dict, remap_list,
        map_stat_dict, length_dict) = parse_sam(snakein.bam[0], starcode_set)
    if len(remap_list) > 0:
        print(remap_list[0])
        with gzip.open(snakeout.remap_fq[0], "wt") as fqfile:
            SeqIO.write(remap_list, fqfile, 'fastq')
        options = snakeparam.options[snakeparam.num]
        align_command = ("bowtie2 -p %i -t %s"
                         " -x %s -U %s | "
                         "samtools view -Sb - > %s" % (threads,
                                                       ' '.join(options),
                                                       snakeparam.bowtie_index,
                                                       snakeout.remap_fq[0],
                                                       snakeout.remap[0]))
        os.system(align_command)
        remap_sam_out = parse_sam(snakeout.remap[0], starcode_set)
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
        os.system('touch %s; touch %s;' % (snakeout.remap_fq[0],
                                           snakeout.remap[0]))



    write_bed(map_dict, snakeout.bed[0], snakeparam.max_dist[snakeparam.num])

    top_map_dict = top_map(map_dict, snakeparam.max_dist[snakeparam.num])
    with open(snakeout.table[0], 'w') as fout:
        fout.write('\t'.join(('barcode', 'seqname', 'ori', 'start_pos',
                              'total_mapped', 'av_mapq', 'freq1', 'freq2',
                              'seq')))
        fout.write('\n')
        for bc in top_map_dict:
            fout.write('\t'.join([bc] + [str(item) for item in
                                         top_map_dict[bc]]))
            fout.write('\n')

    with open(snakeout.stats[0], 'w') as f_stats:
        with open(snakeout.length[0], 'w') as f_length:
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
