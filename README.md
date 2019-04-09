# Promoters_in_LADs
Code used in Leemans et al 2019: "Promoter-intrinsic and local chromatin features determine gene repression in LADs"


## Dependencies

In order to run the code successfully, a couple of Dependencies are required.
Easiest way of installing them would be to use conda:

*conda env create -f config/conda_environment.yaml*

Dependencies:
cutadapt v1.11
chipseq-greylist v1.0.1
deeptools v3.2.0
snakemake v4.0.0
sambamba v0.6.8
samtools v1.9

R Dependencies (separate from conda):


## setting paths
Before generating the data, there are a couple of paths that will need to
be set in each of the config files.

The ChIP and TRIP data is downloaded automatically from SRA.
PROseq and GROcap data can be obtained by running the proseq_analysis.bash and
grocap_analysis.bash scripts respectively.

For SuRE data re-processed to hg38, please contact:
b.v.steensel@nki.nl


## preparing data

In order to generate the data following code can be used.

for the TRIP experiment:

*snakemake -s scripts/trip/trip.snake \
          --configfile config/TRIP_config.yaml \
          --use-conda*


Epigenetic effect on promoter expression (SuRE vs GROcap):

*snakemake -s scripts/promoters_sure/hg38_lad_repression.snakemake \
          --configfile config/hg38_lad_repression.yaml*

When multiple cores are available *"-j [# of cores]"* can be used to increase
the number of cores used.
For TF de-novo motif analysis and affinity calculations, code snippets are in
the respecitive R markdown file.



## creating figures

Reports generating figures can be found in *Promoters_in_LADs/scripts/reports*
Which file creates which figures can be found below:

cl20181031_plasmid_mixes.Rmd
-   Figure 4

cl20181220_TRIP_feature_plots.Rmd
-   Figure S4

cl20190109_lad_detachement.R
-   Figure 2, S2

cl20190109_prom_SuRE_vs_GROcap.R
-   Figure 1, 6I, S1

cl20190128_TRIP_figures.Rmd
-   Figure 3, 6, S3

fc190129_trip_modeling.Rmd
-   Figure 5, S5

cl20190304_enh_affinity.Rmd
-   Figure 7
