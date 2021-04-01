# 722

### MD5 check ###

My samples were downloaded from the Genome Quebec website using curl on to the McMaster cluster. To ensure they were uploaded without issue I run an MD5 check using the command below. All files print out that they are downloaded 'OK', so I can continue.

````
md5sum -c readSets.md5
````

### Quality check ###

To check te quality of my samples I run `fastqc *.fastq` on them. This will report the number of reads, the proportion of duplication, adapter contamination, and repetetive sequences. I can then run `multiqc .` to compile all the outputs of fastqc to a single summary sheet. The output looks like this:

Screen Shot 2021-04-01 at 8.56.08 AM![image](https://user-images.githubusercontent.com/77504755/113297298-6fb24280-92c8-11eb-94b2-df4574fb5dfe.png)

As expected there is substantial adapter contamination:

Screen Shot 2021-04-01 at 8.59.53 AM![image](https://user-images.githubusercontent.com/77504755/113297492-af792a00-92c8-11eb-93e6-3b6d01d3777a.png)

### Trimmomatic ###

To remove this adapter contamination I use Trimmomatic with the following code:

````
#!/bin/bash

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
  ILLUMINACLIP:${adapter}:2:30:10:2 \ #ILLUMINACLIP:${adapter}:2:30:10:2 Gives the window to search for adapters and the location of my adapter file
  LEADING:5 \ #This tells Trimmomatic to cut off the first and last 5 basepares if they fall below my quality threshold
  TRAILING:5 \ #Same as above
  MAXINFO:35:0.5 \ #MAXINFO:35:0.5 This is the balance between the minimum length I want a read to be and the error rate
  MINLEN:36 \ #This is the minimum length of a read that will be kept
  AVGQUAL:20 #This is the minimum quality to be kept
  
done
````

After this I run the outputs through `fastqc` and `multiqc` again. The adapter contamination has been removed, however there are still a number of 6mer repetitive sequences. I also notice that the aximum sequence length is 151. Although my sequences were supposed to be 150bp in length, Genome Quebec sendan extra basepair on each read. This extra basepar could be a part of the problem.

I decided to try two other trimming programs to see if any will result in no 6mer repetitive sequences. I choose bbduk and trim galore because they are the two I have seen in the literature other than Trimmomatic.

### Trim_galore ###

Trim Galore is a wrapper for cutadapt. Cutadapt trims adapters but dies not quality trim, trim_galore integrates fastqc to trim for quality as well. After a first pass, there is still 151 basepairs for some reads, and 6mer contamination in all reads. So the script below was used to first run it through to hard trim to 150bp, and then screen for adapters and quality. The output still had 6mer contamination and lower quality than Trimmomatic, which makes it an unideal fit for this pipeline.

````
#!/bin/bash

# make variable for trim_galore program location
trim=/home/tylera/bin/TrimGalore-0.6.6/trim_galore

# make input directory for raw reads
raw_dir=/2/scratch/TylerA/SSD/genomes

# make output directory from trimmomatic outputs
trim_dir=/2/scratch/TylerA/SSD/trim_galore

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

#fastqc flag signals that we want to quality trim the reads as well

${trim} --paired --fastqc --cores 4 --retain_unpaired --gzip \
${trim_dir}/150/${base}_R1.150bp_5prime.fq ${trim_dir}/150/${base}_R2.150bp_5prime.fq \
--o /2/scratch/TylerA/SSD/trim_galore
  
done
````

### bbduk ###

Next I tried bbduk. This program is meant for quality as well as adapter trimming. It also has alot of options to customize exactly how you want to trim sequences, so it can be very adaptable. 

````
#!/bin/bash

# make variable for trimmomatic program location
#trim=/home/tylera/bin/TrimGalore-0.6.6/trim_galore

# make input directory for raw reads
raw_dir=/2/scratch/TylerA/SSD/genomes

# make output directory from trimmomatic outputs
trim_dir=/2/scratch/TylerA/SSD/bbduk

#list all files to be read (all raw data)
files=(${raw_dir}/*_R1.fastq.gz)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1.fastq.gz`

/usr/local/BBmap/bbduk.sh \ #Path to bbduk
in1=${raw_dir}/${base}_R1.fastq.gz \
in2=${raw_dir}/${base}_R2.fastq.gz \
out1=${trim_dir}/${base}_R1_trimmed.fastq.gz \
out2=${trim_dir}/${base}_R2_trimmed.fastq.gz \
ref=/home/tylera/scripts/AllAdapters.fa \ #Path to reference adapters
threads=8 ftr=149 ktrim=r k=23 mink=6 hdist=1 tpe tbo \
qtrim=rl trimq=20 minlength=36 2> /2/scratch/TylerA/SSD/bbduk/log/${base}.log

#ftr=149: the basepair to trim after (anchored at 0) so 149 will trim off the 151st bp
#ktrim=r: This flag tells bbduk that once it finds and adapter to trim everything to the right (r) of it
#k=23: This is the expected length of my adapter
#mink=6: This is the smallest possible length of an adapter to trim
#hdist=1: This allows one mismatch in a trimmed region
#tpe: This tells bbduk to trim both reads to the same length if an adapter is found on one pair
#tbo: This trims based on paired overlap detection
#qtrim=rl: Adapters are expected on the right side of reads
#trimq=20: Trim anything qith quality less than 20
#2> : save the log to this folder
  
done
````

BBduk removed far more 6mers and also removed the extra asepair present with the other two trimmers. Quality appears to be high, and although there is still a small amount of 6mer contamination, they do not match any adapter sequences so they may be contamination from yeast or bacteria which are common in fly samples. BBduk also did not removed much more than the other two, which hopefully means I did not accidentally over-trim my reads.



### BWA mem ###

In order to map my reads I download the Drosophila reference genome from flybase using 'curl'
````
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.23_FB2018_04/fasta/dmel-all-chromosome-r6.23.fasta.gz
````
I then need to index it so that BWA knows how to read it.

````
bwa index dmel-all-chromosome-r6.23.fasta.gz
````

Once I have an indexed reference genome I can map my reads to it using 'BWA mem' Using the following script.


````
#!/bin/bash

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
files=(${trim_dir}/*_R1_trimmed_good.fastq.gz)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_trimmed_good.fastq.gz`

bwa mem -t 16 -M ${refGenome} \
${trim_dir}/${base}_R1_trimmed_good.fastq.gz \
${trim_dir}/${base}_R2_trimmed_good.fastq.gz \
> ${output}/${base}_bwa_PE.SAM
done

# -M Mark shorter split hits as secondary (for Picard compatibility).
# -t used to tell BWA to use 16 threads
````


### Sam to Bam ###

I next converted the SAM files outputted from the mapping to BAM files. SAM files are very large, and take up too much space on the cluster, so keeping them in their binary format is less memory intensive. The following script is written to just take any SAM files in the current directory and convert them to BAM. I can then use `rm *SAM` to get rid of the huge SAM files.

During this step I also sort the reads in the samples. This will be necessary in downstream filtering and other steps.

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
# @8 use 8 threads
````

### Merge files ###

My initial samples did not meet the read depth that we paid for, so Genome Quebec re-ran any samples that were less than 105M reads. At this stage I need to merge the samples from the two different runs in to matching files. I use the following script:

````
#!/bin/bash

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

For this to run I create a run1/ and run2/ directory and sort my reads in to them using `mv`. I then need to move the samples that did not need an extra run in to my merged/ folder, because they can cause issues with this script if they do not have a matching pair in the run2/ file. After the files have been merged I can delete my run1/ and run2/ directories.

### Quality filter ###

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
## Deduplicating the bam files ##

````
#! /bin/bash

mapped_dir=/2/scratch/TylerA/SSD/bwamap/filter30
outdir=/2/scratch/TylerA/SSD/bwamap/dedup
picdir=/usr/local/picard-tools/picard.jar

files=(${mapped_dir}/*.bam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
#echo ${name}
java -Xmx2g -jar \
${picdir} MarkDuplicates \
I=${mapped_dir}/${base}.bam \
O=${outdir}/${base}_rmd.bam \
M=${outdir}/dupstat.txt \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES=true
done
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


