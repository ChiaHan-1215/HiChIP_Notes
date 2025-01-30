### Notes for arima MAPS pipeline installation steps
#### install in biowulf conda env
#### the newer version of MAPS is located at https://github.com/HuMingLab/MAPS since Mar, 2023

login Biowulf and setup env,then activate conda

Clone the MAPS from github

`git clone https://github.com/HuMingLab/MAPS.git`

Create conda environment with python 

`conda create -n MAPS_env python=3.10`

Activate the conda environment

`conda activate MAPS_env`

Once activate the conda environment, start install tools

`mamba install -c bioconda deeptools`

`mamba install pandas`

`mamba install  pysam` 

`mamba install -c conda-forge -c bioconda pybedtools` 

`mamba  install -c bioconda macs2`

`mamba install R=4.2`

`$ R` then in R, install package

`install.packages("argparse")`

`install.packages("data.table")`

The error when running huge data, need VGAM older version ...

https://github.com/ijuric/MAPS/issues/29

```r
# in R 
packageurl <- "https://cran.r-project.org/src/contrib/Archive/VGAM/VGAM_1.1-3.tar.gz"

install.packages(packageurl, repos=NULL, type="source")
```


- As for other tools, it's already in Biowulf:

bedtools
	
 samtools
	
 HTSLIB
	
 bcftools
	
 bwa
	
`ml samtools bedtools bwa`

- Once installed, move the `Arima-MAPS_v2.0.sh`  from `MAPS/Arima_Genomics/` directory to `MAPS/bin/`
- It need **absolute pathway** for output for the command will have error. 
- for file name **abcd_**_R1.fastq.gz, the `-I`  inupt in command should  **abcd** 


```
# code example
# the newer version of folder in biowulf is:
# /data/leec20/HiChIP_arima_miseq_novaseq_project_091223

sh /data/leec20/HiChIP_Miseq_0320/MAPS/bin/Arima-MAPS_v2.0.sh\
 -C 1 -p broad -I /data/leec20/HiChIP_Miseq_0320/test_file/Arima-MAPS-test\
 -O /data/leec20/HiChIP_Miseq_0320/test_file/output\
 -o hg38 -b /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa -t 6 -f 1

```

**Note:**
-f 1/0 is design for `-f 0` is for shallow seq and `-f 1` is for Deep sequence

```
#!/bin/bash

source myconda
conda activate MAPS_new
ml samtools bedtools bwa

sh /data/leec20/HiChIP_arima_miseq_novaseq_project_091223/MAPS/bin/Arima-MAPS_v2.0.sh -C 1 -p broad -I /data/leec20/HiChIP_arima_miseq_novaseq_project_091223/ELENTA_fastqs/Sample-2 -O /data/leec20/HiChIP_arima_miseq_novaseq_project_091223/Elenta_Output -o hg38 -b /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa -t 6 -f 0

echo all Done

```

---------


#### Using HiC-pro to generate the files for other downstream tools

- HiC-pro is most used tools https://github.com/nservant/HiC-Pro
- Biowulf has HiCpro (ver: 3.1.0_v2) https://hpc.nih.gov/apps/hicpro.html
- script location in Biowulf for HiC-pro utils: 
`/usr/local/apps/hicpro/3.1.0_v2/HiC-Pro_3.1.0/bin/utils/`
- to run the tools, first need to have `config files`, `genome.size files`, and `enzyme digested file`

**Making enzyme digest file for Armia HiChIP kit**
Since Arima use multiple enzyme digest, need to adjust the script. 

According to github: https://github.com/nservant/HiC-Pro/issues/202 and https://github.com/nservant/HiC-Pro/blob/master/doc/FAQ.md
The script look like this:
`digest_genome.py -r ^GATC G^ANTC -o arima.digest.bed /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa`

**Setting config files for running HiC-pro**
- Config files are located in `/[usr/local/apps/hicpro/config_example.txt](https://github.com/nservant/HiC-Pro/blob/master/config-hicpro.txt)`
- Copy and edit the configuration file 'config-hicpro.txt' in your local folder.
- Here's the example for Arima config files, **NOTE** need to change the cpu mem use, Bowtie path depend on hg19 or hg38, and the path of arima enzyme digest file, the genome.size file

```
#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM - PBS - Start Editing Here !!
#######################################################################
N_CPU = 24
LOGFILE = hicpro_arima.hg19.log

JOB_NAME = arima_mai_test
JOB_MEM = 40G
JOB_WALLTIME = 24:00:00
JOB_QUEUE = norm
JOB_MAIL = 

#########################################################################
## Data
#########################################################################

PAIR1_EXT = _R1
PAIR2_EXT = _R2

#######################################################################
## Alignment options
#######################################################################

FORMAT = phred33
MIN_MAPQ = 30

BOWTIE2_IDX_PATH = /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = genome
GENOME_SIZE = /data/leec20/hichip_CHiC_project/HiC_pro/hg19.chrom.sizes

#######################################################################
## Allele specific
#######################################################################

ALLELE_SPECIFIC_SNP = 

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = /data/leec20/hichip_CHiC_project/HiC_pro/arima_hg19_hicpro_digest.bed
LIGATION_SITE = GATCGATC,GANTGATC,GANTANTC,GATCANTC
MIN_FRAG_SIZE = 100
MAX_FRAG_SIZE = 100000
MIN_INSERT_SIZE = 100
MAX_INSERT_SIZE = 600

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 1
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 2500 5000 10000 25000 500000 1000000
MATRIX_FORMAT = upper

#######################################################################
## ICE Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1

```










