# Notes for arima MAPS pipeline installation steps
## install in biowulf conda env
## the newer version of MAPS is located at https://github.com/HuMingLab/MAPS since Mar, 2023
## Updated installation protocol

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

