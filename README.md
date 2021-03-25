# 722

## MD5 check ##
````
md5sum - c md5.txt 
````

## Trimmomatic ##

````
#!/bin/bash
# Project directory variable:
# project_dir=/home/katie/dachsProject/WingShapeBSA/20141205_A_DNASeq_PE

# make variable for trimmomatic program location
trim=/usr/local/trimmomatic/trimmomatic-0.36.jar

# make input directory for raw reads
raw_dir=/2/scratch/TylerA/SSD/genomes

# make output directory from trimmomatic outputs
trim_dir=/2/scratch/TylerA/SSD/trimmed

# make path to adapter sequences (to be used with ILLUMINACLIP)
adapter=/home/tylera/scripts/AllAdapters.fa

#list all files to be read (all raw data)
files=(${raw_dir}/*_R1.fastq.gz)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1.fastq.gz`

java -jar ${trim} PE -threads 16 -phred33 \
  -trimlog ${trim_dir}/trimlog.txt \
  ${raw_dir}/${base}_R1.fastq.gz \
  ${raw_dir}/${base}_R2.fastq.gz \
  ${trim_dir}/${base}_R1_PE.fastq.gz \
  ${trim_dir}/${base}_R1_SE.fastq.gz \
  ${trim_dir}/${base}_R2_PE.fastq.gz \
  ${trim_dir}/${base}_R2_SE.fastq.gz \
  ILLUMINACLIP:${adapter}:2:30:10:2 \
  LEADING:5 \
  TRAILING:5 \
  MAXINFO:35:0.5 \
  MINLEN:36 \
  AVGQUAL:20
  
done
````

## bbduk ##

````
#!/bin/bash
# Project directory variable:
# project_dir=/home/katie/dachsProject/WingShapeBSA/20141205_A_DNASeq_PE

# make variable for trimmomatic program location
#trim=/home/tylera/bin/TrimGalore-0.6.6/trim_galore

# make input directory for raw reads
raw_dir=/2/scratch/TylerA/SSD/genomes

# make output directory from trimmomatic outputs
trim_dir=/2/scratch/TylerA/SSD/bbduk

# make path to adapter sequences (to be used with ILLUMINACLIP)
#adapter=/home/tylera/scripts/AllAdapters.fa

#list all files to be read (all raw data)
files=(${raw_dir}/*_R1.fastq.gz)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1.fastq.gz`

/usr/local/BBmap/bbduk.sh \
in1=${raw_dir}/${base}_R1.fastq.gz \
in2=${raw_dir}/${base}_R2.fastq.gz \
out1=${trim_dir}/${base}_R1_trimmed.fastq.gz \
out2=${trim_dir}/${base}_R2_trimmed.fastq.gz \
ref=/home/tylera/scripts/AllAdapters.fa \
threads=8 ftr=149 ktrim=r k=23 mink=6 hdist=1 tpe tbo \
qtrim=rl trimq=20 minlength=36 2> /2/scratch/TylerA/SSD/bbduk/log/${base}.log


  
done
````

## Trim_galore ##

````
#!/bin/bash
# Project directory variable:
# project_dir=/home/katie/dachsProject/WingShapeBSA/20141205_A_DNASeq_PE

# make variable for trimmomatic program location
trim=/home/tylera/bin/TrimGalore-0.6.6/trim_galore

# make input directory for raw reads
raw_dir=/2/scratch/TylerA/SSD/genomes

# make output directory from trimmomatic outputs
trim_dir=/2/scratch/TylerA/SSD/trim_galore

# make path to adapter sequences (to be used with ILLUMINACLIP)
#adapter=/home/tylera/scripts/AllAdapters.fa

#list all files to be read (all raw data)
files=(${raw_dir}/*_R1.fastq.gz)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1.fastq.gz`

${trim} --hardtrim5 150 --paired --cores 4 --retain_unpaired \
${raw_dir}/${base}_R1.fastq.gz ${raw_dir}/${base}_R2.fastq.gz \
--o /2/scratch/TylerA/SSD/trim_galore
  
done


# make variable for trimmomatic program location
trim=/home/tylera/bin/TrimGalore-0.6.6/trim_galore

# make output directory from trimmomatic outputs
trim_dir=/2/scratch/TylerA/SSD/trim_galore

files=(${trim_dir}/150/*150bp_5prime.fq)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1.150bp_5prime.fq`

${trim} --paired --fastqc --cores 4 --retain_unpaired --gzip \
${trim_dir}/150/${base}_R1.150bp_5prime.fq ${trim_dir}/150/${base}_R2.150bp_5prime.fq \
--o /2/scratch/TylerA/SSD/trim_galore
  
done
````

## Index for bbmap ##

````
ref=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta.gz
````

## BBMap ##

````
#!/bin/bash
#
#project directory (to keep each run separate)
#run=trimmed

# directory of processed sequences with trimmomatic
trim_dir=/2/scratch/TylerA/SSD/bbduk

# variable for the reference genome
#refGenome=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta.gz

# make output directory from mapping outputs
output=/2/scratch/TylerA/SSD/bbmap

# make BWA directory path
#dir=/usr/local/BBmap

#list all files to be read (this selects the left end from each PE pair)
#list all files to be read (this selects the left end from each PE pair)
files=(${trim_dir}/*_R1_trimmed_good.fastq.gz)

#echo ${files[@]}

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_trimmed_good.fastq.gz`

/usr/local/BBmap/bbmap.sh t=16 in1=${trim_dir}/${base}_R1_trimmed_good.fastq.gz \
in2=${trim_dir}/${base}_R2_trimmed_good.fastq.gz \
out=${output}/${base}.sam


#ref=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta.gz
#Index commmand ran first


#echo base is
#echo ${base}
#echo name is
#echo ${name}


done
````

## BWA mem ##

````
#!/bin/bash
#
#project directory (to keep each run separate)
#run=trimmed

# directory of processed sequences with trimmomatic
trim_dir=/2/scratch/TylerA/SSD/bbduk

# variable for the reference genome
refGenome=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta.gz

# make output directory from mapping outputs
output=/2/scratch/TylerA/SSD/sam_files

# make BWA directory path
bwa_dir=/usr/local/bwa/0.7.8

cd ${bwa_dir}

#list all files to be read (this selects the left end from each PE pair)
#list all files to be read (this selects the left end from each PE pair)
files=(${trim_dir}/*_R1_trimmed_good.fastq.gz)

#echo ${files[@]}

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_trimmed_good.fastq.gz`

#echo base is
#echo ${base}
#echo name is
#echo ${name}

bwa mem -t 16 -M ${refGenome} \
${trim_dir}/${base}_R1_trimmed_good.fastq.gz \
${trim_dir}/${base}_R2_trimmed_good.fastq.gz \
> ${output}/${base}_bwa_PE.SAM
done

# -M Mark shorter split hits as secondary (for Picard compatibility).
````

## Sam to Bam ##

````
#!/bin/bash

#Specify input directory
sam_dir=.

#Specify output directory
bam_dir=.

files=(${sam_dir}/*.SAM)

for file in ${files[@]} 
do 
name=${file} base=`basename ${name} .SAM`

samtools view -b -@8 ${sam_dir}/${base}.SAM | samtools sort -o ${bam_dir}/${base}.bam 
done 

# -b output to bam
# -q skip alignments with MAP quality less than [int]
# @8 use 8 threads
````

## Merge files ##

````
#!/bin/bash
#Only works if all lanes are L001/L002

#Each sample was run on 2 seperate runs, there is a mathing sequencing
#file in each directory (identical names) Within each directory, each
#library was run on 2 lanes marked as 001 and 002
run1=run1
run2=run2

bam_dir=/2/scratch/TylerA/SSD/bwamap
merged_dir=/2/scratch/TylerA/SSD/bwamap/merged

files=(${bam_dir}/${run1}/*_bwa_PE.bam)

for file in ${files[@]} 
do 
name=${file}
base=`basename ${name} _bwa_PE.bam`
samtools merge ${merged_dir}/${base}_merged_aligned_PE.bam \
${bam_dir}/${run1}/${base}_bwa_PE.bam \
${bam_dir}/${run2}/${base}_bwa_PE.bam

done 
````

## Quality filter ##

````
#for all files

#!/bin/bash
#
#Specify output directory
out_dir=/2/scratch/TylerA/SSD/bwamap/filter30

#temp storage 
in_dir=/2/scratch/TylerA/SSD/bwamap/merged

files=(${in_dir}/*_merged_aligned_PE.bam)

for file in ${files[@]} 
do 
name=${file} 
base=`basename ${name} _merged_aligned_PE.bam`
samtools view -b -q 30 -@ 16 \
${in_dir}/${base}_merged_aligned_PE.bam \
> ${out_dir}/${base}.bam
done 
#-b indicates bam files being used
#-q 30 indicates quality threshold
#-@ 16 indicates threads being used
````

## Calculate coverage ##

````
#!/bin/bash

in_dir=/2/scratch/TylerA/SSD/bwamap/filter30

files=(${in_dir}/*.bam)

for file in ${files[@]} 
do 
name=${file} 
base=`basename ${name} .bam`
samtools depth ${in_dir}/${base}.bam > ${in_dir}/${base}.coverage
done
````

## Call histogram ##

````
#!/bin/bash

in_dir=/2/scratch/TylerA/SSD/bwamap/filter30

files=(${in_dir}/*.bam)

for file in ${files[@]} 
do 
name=${file} 
base=`basename ${name} _merged_aligned_PE.bam`
Rscript /2/scratch/TylerA/SSD/scripts/depth_histogram.R ${in_dir}/${base}.coverage ${base}
done 
````

### Rscript for the histogram ###

````
##coverage_histogram.R
## need next line to call arguments:
args <- commandArgs(trailingOnly = TRUE)

#read in the first argument which should be the file
dat <- read.table(args[1])
#the title should be the second argument (the base name)
title <- args[2]
colnames(dat) <- c("chr","pos","depth")

#make a histogram of the coverage for each file
pdf(paste(title, ".pdf", sep=""))
hist(dat$depth, xlim=c(0,500), breaks=500)
dev.off()
````

## Adding in read groups ##

````
#! /bin/bash

#Variable for project:
project_dir=/2/scratch/TylerA/SSD/bwamap/dedup

#Path to Picard
picdir=/usr/local/picard-tools/picard.jar


files=(${project_dir}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -jar ${picdir} AddOrReplaceReadGroups I=${project_dir}/${base}.bam \
  O=${project_dir}/${base}_RG.bam \
  RGID=1_2 \
  RGLB=library1 \
  RGPL=illumina \
  RGPU=None \
  RGSM=${base}

done
````

## Index for GATK ##

````
#! /bin/bash 

# Index files

#Variable to put indexed files in 
index_dir=/2/scratch/TylerA/SSD/bwamap/gatk/index

#Path to input directory
input=/2/scratch/TylerA/SSD/bwamap/dedup

files=(${input}/*_rmd_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _rmd_RG.bam`

samtools index ${input}/${base}_rmd_RG.bam 

done
````

## Mark Indels ##

````
#! /bin/bash 
#Path to input directory
final_bam=/2/scratch/TylerA/SSD/bwamap/dedup

#Path to output directory
gatk_dir=/2/scratch/TylerA/SSD/bwamap/gatk

#Variable for reference genome (non-zipped)
#index_dir=/home/sarahm/cvl/index_dir
ref_genome=/2/scratch/TylerA/Dmelgenome/gatk/dmel-all-chromosome-r6.23.fa

#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar


files=(${final_bam}/*_rmd_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _rmd_RG.bam`

java -Xmx32g -jar ${gatk} -I ${final_bam}/${base}_rmd_RG.bam \
-R ${ref_genome} \
  -T RealignerTargetCreator \
  -o ${gatk_dir}/${base}.intervals

done
````


