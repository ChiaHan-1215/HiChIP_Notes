# HiChIP_HiC_capture_stuff

### this is notes/tools used for HiC_cap/HiChIP project
use **Arima** kit as major tool

## HiC_pro , FitHiChip for Aqua_hichip

### Biowulf has HiCpro (ver: 3.1.0)

### script location: 
/usr/local/apps/hicpro/3.1.0_conda/bin/utils

### config text
/usr/local/apps/hicpro/config_example.txt

### first digest genome using .py 
### the cutting site for Arima : 
https://github.com/nservant/HiC-Pro/blob/master/doc/FAQ.md

### digest script example: 
`digest_genome.py -r ^GATC G^ANTC -o arima.test.bed /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa`

### according to: 
https://github.com/nservant/HiC-Pro/issues/202

### modify config.txt 

### Run HiC-pro
`HiC-Pro -i data -o output_sbatch -c config.txt`

### paper:
https://www.nature.com/articles/s41467-021-21867-0?elqTrackId=07572bf0c6e54475b75496bb5c1d030f

### running in separate mode: -s option

https://nservant.github.io/HiC-Pro/MANUAL.html

### bwahaha in HiC-pro problem
https://github.com/nservant/HiC-Pro/issues/124


###
# Auqa_HiC_pro 

http://py4e-data.dr-chuck.net/known_by_Fikret.html

### convert HiC-pro output to .hic file using juicer:
### /usr/local/apps/juicer/juicer-1.6
### /usr/local/apps/juicer/juicer-1.6/scripts/juicer_tools_1.22.01.jar


`$ module load juicer hicpro`
`$ tar=/usr/local/apps/hicpro/3.1.0_conda/bin/utils`
`$ $tar/hicpro2juicebox.sh -i ../../hic_results/data/rh4_hg19_d/rh4_hg19_d.allValidPairs -g hg19 -j /usr/local/apps/juicer/juicer-1.6/scripts/juicer_tools.jar`


### dump 

`juicer_tools dump observed NONE hic_file.hic 11:17600000:17800000 11:17600000: 17800000 BP 5000 output.txt`

### whats dump value? 
is interaction freq


### The mainly normalization in Aqua
Convert to contacts per million (CPM), then AQuA normalize
https://github.com/GryderArt/AQuA-HiChIP/blob/master/extractAPA_plotAQuA-APA.R

`CPM.result = Matrix * 1000000 / sum of valid_interaction_rmdup of hg19 and mm10 
Final normalization =  CPM.result * Aqua_factor`



### Conver HiC-pro to BEDOE format: the package is inside the HiCcompare/ the guy recent interview 
https://github.com/dozmorovlab/HiCcompare/blob/master/R/hicpro2bedpe.R

### HiCcompare tutorial
https://www.bioconductor.org/packages/devel/bioc/vignettes/HiCcompare/inst/doc/HiCcompare-vignette.html


### file in the /output/hic_results/matrix/sample/raw/resulation/


### need .bed and .matrix

### arima_MAPS pipeline resolution is 5000
`Max_range = 2000000 (in .sh script)`

### .allvalidpairs in HicPro, probably can be covert as bedpe file??


The format of allValidpair

read name / chr_reads1 / pos_reads1 / strand_reads1 / chr_reads2 / pos_reads2 / strand_reads2 / fragment_size [/ allele_specific_tag]

In HiC-Pro, this is the size of the DNA fragments before sequencing. Basically, the distance between the first bases of R1/R2.

To convert the allValidpair to gi_list


### 4C -plot to loop plot 
### Allvalidparis to beep

https://groups.google.com/g/hic-pro/c/qhHNbQVJfYY

https://github.com/nservant/HiC-Pro/issues/362

### cloop to convert dump to bedpe 
https://github.com/YaqiangCao/cLoops/issues/14


### FitHiCHip

GOAL: 
Plan A : scale aquq-facor in the matrix IF, int , done!
Plan B: scale ICEd ???

# Bias is for Bias regression model

Target:
`FitHiChIP_HiCpro.sh`

### after the operation, a file with 8 columns is produced 
### chr start end coverage ispeak mappability GCcontent NoCutSites

xxx.norm.Contact.Matrix'

FitHiChIP.AllBin_CompleteFeat_ICE.bed
Chr	Start End Coverage IsPeak “Bias”


``$RScriptExec ./src/Significance_Features.r -I $Interaction_File -E $AllFeatFile -O $IntFileALLtoALL -C $ChrSizeFile``

### 

Src/Interaction_Sort_Genomic_Distance.r

Cerated a file called “Interactions.sortedGenDist.bed”

## cc= contact count between corresponding interacting regions == HiC-pro raw matrix 3rd cloumn

## coverage also scaled according to aqua-factor

## The bias value won’t change regardless scale aqua-factor or not


#####################
# FitHichip Diff modify ####
#####################

#  the code it has problem dealing with bedgraph

https://github.com/ChenfuShi/HiChIP_peaks/issues/15

# Tool he provide:
https://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html#bedmap

# Bedops is in Biowulf:
https://hpc.nih.gov/apps/bedops.html


# background:

Since the narrow peak from MACS2 can be used as bedgraph,
The “SingleValue” is tag density, we could use that as the “Map element set in bedops tool”  the value will be 7th colum



The “reference set” use the 5000 bin, in the code 789:
bedtools makewindows -g", opt$ChrSizeFile, "-w", BinSize, ">", TargetBinnedChrFile

# format example:

“map element set, the value is density from narrowpeak, at 7th column”
#####
The format need 5 column to run, so need to modify the sample to this:


“Reference set, in 5000 bin size ” (OPTIONAL)
 
# code example:

bedmap --echo --mean  5000_bin_with_bin.bed test_with_name.bedgraph

Or 

bedmap --mean  --skip-unmapped  5000_bin_with_bin.bed test_with_name.bedgraph

Or 

#output


Probably need NAN to 0 since it is according to order of bin 
Can it accept decimal point??

# line 641 end function

# line 802-819 
If i==1, else function, both will be excuated, the code run if==1 first, second round just run else code since the second round i == 2 
 
# separate column by “|” to ‘\t’ 

awk '{sub(/\|/,"\t"); print $0}' result.txt

bedmap --echo --mean standard_bin test_with_name.bedgraph | awk '{sub(/\|/,”\t”); print $0}' > Mergedfile


# replace NAN to 0

sed 's/NAN/0/g' file > finalfile

￼

In R: wondering 

How to escape some regex in R, the “\” is escape sign :
How to escape “\” ==> \\

# no need escape in paste():
/ 
|
# need to escape 
\
“
‘

system(paste("bedmap --echo --mean", f1 , f2, "| awk","'{sub(/\\|/,\"\t\"); print $0}'",">", "MERGED.txt" ))




#######################
### Chip-sed bam file:######
#######################

For searching FitHiCHip diff input chip-seq bam file, found MAPS pipeline has some insight.

# they use this bam file as MACS2 input for calling peak

feather_output/SC829846_CGATGT_L001_20211013_103654//tempfiles/SC829846_CGATGT_L001.shrt.bam

# find this place how to generate this bam file. 
SC829846_CGATGT_L001.shrt.bam


# Command line:
 
MACS2 callpeak -t SC829846_CGATGT_L001.shrt.bam -n SC829846_CGATGT_L001 -g hs --broad --nolambda --broad-cutoff 0.3 --outdir ./OUTPUTFOLDER

-t = input bam file
-n =file name
https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md

Printing with sed or awk a line following a matching pattern

https://stackoverflow.com/questions/17908555/printing-with-sed-or-awk-a-line-following-a-matching-pattern/17914105#17914105

python $cwd/feather/feather_pipe preprocess -o $feather_output -p $dataset_name -f1 $fastq1 -f2 $fastq2 -b $bwa_index -q $mapq -l $length_cutoff -t $threads -c $per_chr -j $generate_hic -d $optical_duplicate_distance	
	

# Excuse the script in the folder with fastq files

python3 py_sc/feather_pipe preprocess -o out/ -p test -f1 t1_R1.fastq.gz -f2 t1_R2.fastq.gz -b $genome -q 30 -l 1000 -c 0 -j 0 -d 2500 -t 12

# HiC-pro output explain:
https://github.com/nf-core/hic/blob/master/docs/output.md
HiC_pro , FitHiChip for Aqua_hichip

# Biowulf has HiCpro (ver: 3.1.0)

# script location: 
/usr/local/apps/hicpro/3.1.0_conda/bin/utils

#config text
/usr/local/apps/hicpro/config_example.txt

# first digest genome using .py 
# the cutting site for Arima : 
https://github.com/nservant/HiC-Pro/blob/master/doc/FAQ.md

# digest script example: 
digest_genome.py -r ^GATC G^ANTC -o arima.test.bed /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa

# according to: 
https://github.com/nservant/HiC-Pro/issues/202
￼
# modify config.txt 

# Run HiC-pro
HiC-Pro -i data -o output_sbatch -c config.txt

# paper:
https://www.nature.com/articles/s41467-021-21867-0?elqTrackId=07572bf0c6e54475b75496bb5c1d030f

# running in separate mode: -s option

https://nservant.github.io/HiC-Pro/MANUAL.html

# bwahaha in HiC-pro problem
https://github.com/nservant/HiC-Pro/issues/124





###
###
Auqa_HiC_pro 

http://py4e-data.dr-chuck.net/known_by_Fikret.html

# convert HiC-pro output to .hic file using juicer:
# /usr/local/apps/juicer/juicer-1.6
#/usr/local/apps/juicer/juicer-1.6/scripts/juicer_tools_1.22.01.jar


$ module load juicer hicpro
$ tar=/usr/local/apps/hicpro/3.1.0_conda/bin/utils
$ $tar/hicpro2juicebox.sh -i ../../hic_results/data/rh4_hg19_d/rh4_hg19_d.allValidPairs -g hg19 -j /usr/local/apps/juicer/juicer-1.6/scripts/juicer_tools.jar


## dump 

juicer_tools dump observed NONE hic_file.hic 11:17600000:17800000 11:17600000: 17800000 BP 5000 output.txt

# dump value? 
is interaction freq


# The mainly normalization in Aqua
Convert to contacts per million (CPM), then AQuA normalize
https://github.com/GryderArt/AQuA-HiChIP/blob/master/extractAPA_plotAQuA-APA.R

CPM.result = Matrix * 1000000 / sum of valid_interaction_rmdup of hg19 and mm10 
Final normalization =  CPM.result * Aqua_factor


###
###

Conver HiC-pro to BEDOE format: the package is inside the HiCcompare/ the guy recent interview 
https://github.com/dozmorovlab/HiCcompare/blob/master/R/hicpro2bedpe.R

# HiCcompare tutorial
https://www.bioconductor.org/packages/devel/bioc/vignettes/HiCcompare/inst/doc/HiCcompare-vignette.html


# file in the /output/hic_results/matrix/sample/raw/resulation/


# need .bed and .matrix

# arima_MAPS pipeline resolution is 5000
Max_range = 2000000 (in .sh script)

# .allvalidpairs in HicPro, probably can be covert as bedpe file??


################

The format of allValidpair

read name / chr_reads1 / pos_reads1 / strand_reads1 / chr_reads2 / pos_reads2 / strand_reads2 / fragment_size [/ allele_specific_tag]

In HiC-Pro, this is the size of the DNA fragments before sequencing. Basically, the distance between the first bases of R1/R2.

To convert the allValidpair to gi_list



#################

# 4C -plot to loop plot 
￼

##
Allvalidparis to beep

https://groups.google.com/g/hic-pro/c/qhHNbQVJfYY

https://github.com/nservant/HiC-Pro/issues/362

# cloop to convert dump to bedpe 
https://github.com/YaqiangCao/cLoops/issues/14



## FitHiCHip###

GOAL: 
Plan A : scale aquq-facor in the matrix IF, int , done!
Plan B: scale ICEd ???

# Bias is for Bias regression model

Target:
FitHiChIP_HiCpro.sh

 # after the operation, a file with 8 columns is produced 
# chr start end coverage ispeak mappability GCcontent NoCutSites

xxx.norm.Contact.Matrix'

FitHiChIP.AllBin_CompleteFeat_ICE.bed
Chr	Start End Coverage IsPeak “Bias”
￼


$RScriptExec ./src/Significance_Features.r -I $Interaction_File -E $AllFeatFile -O $IntFileALLtoALL -C $ChrSizeFile

#############################

Src/Interaction_Sort_Genomic_Distance.r

Cerated a file called “Interactions.sortedGenDist.bed”

# cc= contact count between corresponding interacting regions == HiC-pro raw matrix 3rd cloumn

# coverage also scaled according to aqua-factor

# The bias value won’t change regardless scale aqua-factor or not


#####################
# FitHichip Diff modify ####
#####################

#  the code it has problem dealing with bedgraph

https://github.com/ChenfuShi/HiChIP_peaks/issues/15

# Tool he provide:
https://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html#bedmap

# Bedops is in Biowulf:
https://hpc.nih.gov/apps/bedops.html


# background:

Since the narrow peak from MACS2 can be used as bedgraph,
The “SingleValue” is tag density, we could use that as the “Map element set in bedops tool”  the value will be 7th colum

￼

The “reference set” use the 5000 bin, in the code 789:
bedtools makewindows -g", opt$ChrSizeFile, "-w", BinSize, ">", TargetBinnedChrFile

# format example:

“map element set, the value is density from narrowpeak, at 7th column”
#####
The format need 5 column to run, so need to modify the sample to this:

￼

“Reference set, in 5000 bin size ” (OPTIONAL)
￼
 
# code example:

bedmap --echo --mean  5000_bin_with_bin.bed test_with_name.bedgraph

Or 

bedmap --mean  --skip-unmapped  5000_bin_with_bin.bed test_with_name.bedgraph

Or 

#output, 
￼


Probably need NAN to 0 since it is according to order of bin 
Can it accept decimal point??

# line 641 end function

# line 802-819 
If i==1, else function, both will be excuated, the code run if==1 first, second round just run else code since the second round i == 2 
 
# separate column by “|” to ‘\t’ 

awk '{sub(/\|/,"\t"); print $0}' result.txt

bedmap --echo --mean standard_bin test_with_name.bedgraph | awk '{sub(/\|/,”\t”); print $0}' > Mergedfile


# replace NAN to 0

sed 's/NAN/0/g' file > finalfile

￼

In R: wondering 

How to escape some regex in R, the “\” is escape sign :
How to escape “\” ==> \\

# no need escape in paste():
/ 
|
# need to escape 
\
“
‘

system(paste("bedmap --echo --mean", f1 , f2, "| awk","'{sub(/\\|/,\"\t\"); print $0}'",">", "MERGED.txt" ))




#######################
### Chip-sed bam file:######
#######################

For searching FitHiCHip diff input chip-seq bam file, found MAPS pipeline has some insight.

# they use this bam file as MACS2 input for calling peak

feather_output/SC829846_CGATGT_L001_20211013_103654//tempfiles/SC829846_CGATGT_L001.shrt.bam

# find this place how to generate this bam file. 
SC829846_CGATGT_L001.shrt.bam


# Command line:
 
MACS2 callpeak -t SC829846_CGATGT_L001.shrt.bam -n SC829846_CGATGT_L001 -g hs --broad --nolambda --broad-cutoff 0.3 --outdir ./OUTPUTFOLDER

-t = input bam file
-n =file name
https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md

Printing with sed or awk a line following a matching pattern

https://stackoverflow.com/questions/17908555/printing-with-sed-or-awk-a-line-following-a-matching-pattern/17914105#17914105

python $cwd/feather/feather_pipe preprocess -o $feather_output -p $dataset_name -f1 $fastq1 -f2 $fastq2 -b $bwa_index -q $mapq -l $length_cutoff -t $threads -c $per_chr -j $generate_hic -d $optical_duplicate_distance	
	

# Excuse the script in the folder with fastq files

python3 py_sc/feather_pipe preprocess -o out/ -p test -f1 t1_R1.fastq.gz -f2 t1_R2.fastq.gz -b $genome -q 30 -l 1000 -c 0 -j 0 -d 2500 -t 12

# HiC-pro output explain:
https://github.com/nf-core/hic/blob/master/docs/output.md





