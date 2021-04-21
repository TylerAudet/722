# A pipeline for identifying variants of interest in a sexual size dimorphism reversed long term evolve and resequence study

## INTRODUCTION

Body size is a complex polygenic trait with important biological consequences. Physiologically, understanding the underlying genetics of body size variation is important because it influences a number of life-history traits, and provides insight into the evolutionary history of a species (Peters, 1989). Sexual size dimorphism (SSD) is especially interesting because despite sharing much of the autosomal genome, males and females often vary in body size. This suggests that males and females have different optima for body size, and yet must use the same genetic architecture to achieve their own optimal size (Stewart & Rice, 2018). In the model organism D. melanogaster, females are 30% larger on average than males (Rideout et al., 2015). This is likely due to fecundity selection on females creating an optimally larger body size to maximize egg laying (Rideout et al., 2015). Sexually dimorphic body size has been examined in the lab by Rideout et al. (2015) and Millington et al. (2020;2021) and a number of genes of large effect have been identified. Many of these genes are in the insulin signaling pathway such as Drosophila insulin like peptide II (DILPII) and the transcription factor forkhead box O (FOXO) (Millington et al., 2020). The transformer (tra) gene in the sex differentiation pathway has also been shown to be an upstream influencer of sexual size dimorphism (Rideout et al., 2015). These genes have large effects, and were examined in laboratory bred strains, and have yet to be demonstrated in natural populations. Due to the large effect and biological necessity of these genes, they are presumably under stringent purifying selection, and variation in these genes does not seem like a likely explanation for the variability in body size within species.

Experimentally evolved populations combined with whole genome sequencing has created an opportunity to tease apart complex traits such as SSD. Artificially selecting on traits over multiple generations allows the examination of standing genetic variation and how it can respond to natural selection. Combined with this, methods derived for analyzing genomes sequenced from pools of individuals, termed pool-seq, has drastically reduced the cost of sequencing experimentally evolved populations. One such study looking at genomes from experimentally evolved D. melanogaster was conducted by Turner et al. (2011). Turner et al. (2011) selected on body size in D. melanogaster for 100 generations to create a lineage that is larger and a lineage that is smaller in average body size. They then looked for alleles under selection using a program called diffstat to identify peaks of differentiation across the genome. Body size in these lineages were incredibly polygenic, with peaks of selection spread out across the entire genome. The gene ontology for many of these peaks were in post-embryonic development and metamorphosis, and not in large impact genes in the insulin pathway (Turner et al., 2011). They inferred a high heritability for body size because of the rapid response to selection in the large and small selection lines, with a significant difference after 100 generations (Turner et al., 2011).

Parallel to this size concordant selection experiment, another lineage was selected for sex-discordant body size. This sex discordant lineage was selected for small females and large males to reverse the SSD present naturally in D. melanogaster (Stewart & Rice, 2018). Sexually discordant selective pressures create conflict in the genome as each sex diverges toward their optimum and away from the optimum of the other sex, termed intra-locus sexual conflict (Rice, 1992). SSD in particular has been shown to have a strong resistance to change when selected on in outbred populations (Reeve & Fairbairn, 1996; Tigreros & Lewis, 2011). Despite this, after 250 generations Stewart & Rice (2018) successfully reversed the SSD in D. melanogaster. This reversal was due predominantly to a decrease in female body size, and male body size was left unchanged. This is interesting because male body size changed in the lineages where body size was increased or decreased, which suggests a genetic correlation between the sexes may be impacting response (Stewart & Rice 2018). This is an exciting premise, because most traits responsible for sexual dimorphism are autosomal, and therefore this sex specific response may be indicative of an interesting genetic interaction at sexually antagonistic loci (Fairbarn & Roff, 2006). 

Stewart & Rice (2018) conclude that SSD is a highly persistent trait. They then hypothesise that genes responsible for this phenotypic conservation would be extremely polygenic and located far upstream of body size determining loci (Stewart & Rice, 2018). The genomic consequences of this SSD reversal have yet to be examined, and a comparison to the large and small selection lineages have also not been conducted. With this pipeline I plan on identifying SNPs that are present in this sexual size reversal lineage of flies. This will lay the foundation to follow up on Turner et al., (2011) after 378 generations of experimental evolution. This work also allows an examination of the genes responsible for the conservation of SSD in D. melanogaster. I hypothesize that genes under selection in the SSD reversal lineage will be upstream of the genes found by Turner et al. (2011), and that possibly the same alleles may be present in the sex reversal population as the large and small populations. In this examination I plan on establishing a pipeline that can be used to clean and map the genomes and identify SNPs that are under selection. I also plan on looking at Fst as a measure of heterozygosity in the reversal selection lineage to identify locations across the genomes that are under selection to further verify the loci of interest for further analysis.

## Methods

Selection experiment and data acquisition

Selection was performed by A. Stewart and early results published in Turner et al., (2011) and Stewart & Rice (2018). Flies were sorted by size with a novel sieve system and allowed to mate and lay eggs for 24 hours. This was repeated every 14 days for 378 generations. After 378 generations flies were killed in ethanol and stored at -18oC until extraction. DNA was extracted with a Qiagen DNeasy kit. The DNA was sequenced on an Illumina NovaSeq 6000 from a PCR-free shotgun library. Sequencing and library preparation was performed by Genome Quebec.

Pipeline components

Data was first checked for issues from exporting using an MD5 checksum. Next FastQC was used to check quality and size, and BBduk was used to trim adapters and low-quality reads (Bushnell et al., 2017). FastQC was used to verify that adapters were removed (Andrews, 2010). Next reads were mapped to a reference genome using BWA-mem (Li, 2013). SAM files were converted to the smaller BAM files and a second read from the sequencer was merged with the first and all files were filtered to remove anything with read quality under 30 and read-groups were added for GATK compatibility using SAMtools view and merge (Li et al., 2009). Optical duplicates were removed using Picard (Picard Team, 2021). The genome analysis toolkit (GATK) was used to realign around indels (Van der Auwera & O'Connor, 2020). An mpileup file was created using SAMtools for the SSD reversal and control populations to do pairwise analysis. Repetitive regions were removed using Picard and SNPs were called using Varscan (Koboldt et al., 2012). A sync file was also created for population statistics using popoolation2 (Kofler et al., 2011). Population statistics were carried out using R version 4.0.1 (R core team, 2020). To create an R object from a sync file and to calculate Ne Pool-seq was used (Taus et al., 2017). Acer was used to perform a CMH test for selection on allele frequency (Spitzer et al., 2020). Finally, genes of interest found from the CMH were filtered from the SNP calls using vcftools (Danecek et al., 2011). These gees of interest were then annotated using SNPEff (Cingolani et al., 2012).

### MD5 check ###

My samples were downloaded from the Genome Quebec website using `curl` on to the McMaster cluster. To ensure they were uploaded without issue I run an MD5 check using the command below. All files print out that they are downloaded 'OK', so I can continue.

````
md5sum -c readSets.md5

# -c flag denotes that we are checking they are the same as expected
````

### Quality check ###

To check the quality of my samples I run `fastqc *.fastq` on them. This will report the number of reads, the proportion of duplication, adapter contamination, and repetetive sequences. I can then run `multiqc .` to compile all the outputs of fastqc into a single summary sheet. The output looks like this:

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
  ILLUMINACLIP:${adapter}:2:30:10:2 \
  LEADING:5 \
  TRAILING:5 \
  MAXINFO:35:0.5 \
  MINLEN:36 \
  AVGQUAL:20
  
done

#ILLUMINACLIP:${adapter}:2:30:10:2 Gives the window to search for adapters and the location of my adapter file
#LEADING:5 This tells Trimmomatic to cut off the first 5 basepares if they fall below my quality threshold
#TRAILING:5 This tells Trimmomatic to cut off the last 5 basepares if they fall below my quality threshold
#MAXINFO:35:0.5 This is the balance between the minimum length I want a read to be and the error rate
#MINLEN:36 This is the minimum length of a read that will be kept
#AVGQUAL:20 This is the minimum quality to be kept, any reads with quality below 20 are trimmed out

````

After this I run the outputs through `fastqc` and `multiqc` again. Multiqc says the adapter contamination has been removed, however there are still a number of 6mer repetitive sequences that may indicate some adapters were not fully removed. I also notice that the maximum sequence length is 151. Although my sequences were supposed to be 150bp in length, Genome Quebec sends an extra basepair on each read. This extra basepair could be influencing the ability for trimmomatic to identify adapter fragments. For this reason I decided to try two other trimming programs to see if any will result in no 6mer repetitive sequences. I choose bbduk and trim galore because they are the two I have seen in the literature other than Trimmomatic.

### Trim_galore ###

Trim Galore is a wrapper for cutadapt. Cutadapt trims adapters but does not quality trim, trim_galore integrates fastqc to trim for quality as well. After a first pass, there is still 151 basepairs for some reads, and 6mer contamination in all reads. So the script below was used to first run it through to hard trim to 150bp, and then screen for adapters and quality. The output still had 6mer contamination and lower quality than Trimmomatic, which makes it an unideal fit for this pipeline.

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

#--hardtrim5 150 Trimms all reads to be 150bp
#--retain_unpaired Keeps the reads that don't pair after trimming

````

````

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

#--paired denotes we are using paired reads
#--fastqc tells trim_galore that we want quality trimming
#--gzip outputs a gzipped file
#--cores 4 multicore to speed up the process

````

### bbduk ###

Next I tried bbduk. This program is meant for quality as well as adapter trimming. It also has a large number of options to customize exactly how you want to trim sequences, so it can be very adaptable. 

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

BBduk removed far more 6mers and also removed the extra basepair present with the other two trimmers. Quality appears to be high, and although there is still a small amount of 6mer contamination, they do not match any adapter sequences so they may be contamination from yeast or bacteria which are common in fly samples. BBduk also did not removed many more reads than the other two, which hopefully means I did not accidentally over-trim my reads.

### BWA mem ###

In order to map my reads I download the Drosophila reference genome from flybase using 'curl'
````
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.23_FB2018_04/fasta/dmel-all-chromosome-r6.23.fasta.gz
````
I then need to index the reference genome so that BWA knows how to read it.

````
bwa index dmel-all-chromosome-r6.23.fasta.gz
````

Once I have an indexed reference genome I can map my reads to it using 'BWA mem' with the following script.

````
#!/bin/bash

# directory of processed sequences with BBduk
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
# -t used to tell BWA to use 16 threads to speed things up
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

Next I need to filter for quality because accuracy is very important for SNP calling, as low quaslity can result in false positives. To do this I can use `samtools view` to filter out reads below a certain quality score. I ran this code twice, once with a conservative quality of 20, and again with a more stringent quality of 30 as the threshold. I compared the coverage between the two and found that a more stringent quality score removed more reads without impacting the coverage, so I decided to go with a quality filter of 30.

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

### Calculate coverage ###

Next I calculate the coverage of the samples after filtering for each quality threshold using `samtools depth`.

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

The output .coverage file can be visualized with a histogram using the scripts below. First a shell script is created that converts my .coverage files to be compatible with R. Next a script is run which creates a histogram in R. This histogram can be copied to my local computer using `scp` and viewed.

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

### Quality threshold of 20 ###
C1F_20.pdf![image](https://user-images.githubusercontent.com/77504755/113307003-abeaa080-92d2-11eb-9f77-48da9dea49b7.png)

### Quality thresholf of 30 ###

C1F_30.pdf![image](https://user-images.githubusercontent.com/77504755/113307082-c1f86100-92d2-11eb-9aca-cd725e8ecf2f.png)

Coverage was not noticeably reduced by a more stringent filter, so in the interest of calling the fewest false SNPs as possible I filter all samples for quality 30. We expected around 200x coverage using the formula coverage = # of reads * length of reads / length of the genome. My coverage histogram shows coverage over 200x though, which means there could be duplicated reads from PCR optical duplication. It also means we could have copy nuymber variation because some reads have coverage that is double what we expected. Copy number variation can be dealt with downstream in programs by imitting the max coverage of reads for programs that call SNPs or calculate Fst. Optical duplicates however should be removed before moving on.

### Deduplicating the bam files ###

Next I need to remove these optical duplicates. This is done most commonly with Picard. I used the following script:

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

#M=  Where to put summary statistics
#VALIDATION_STRINGENCY=SILENT Improves performance because it doesn't read through unecessary information in the bam file
#REMOVE_DUPLICATES=true Tells picard to not keep duplicates after it finds them
````

Afterwards my coverage looks like this:

***S2F.pdf![image](https://user-images.githubusercontent.com/77504755/113308302-07695e00-92d4-11eb-838f-c6b4c95d68d4.png)

There still appears to be reads over 200x, but fewer. This step also reduced coverage by quite a bit, with far fewer reads close to the 200x mark.

### Combining samples ###

At this stage I combined my replicates and sexes together for each treatment. I also begin only working with my control (C) and SSD reversed (E) treatments. Comparisons between sexes will be done at a later time to look for sexual conflict, and comparison of the large (L) and small (S) flies will also be done at another time. For the purposes of this pipeline I am looking at how SSD-reversal flies compare to the previous finding of Turner et al., (2011). So here I merge my replicates and sexes. I use the same code that I used to merge my two sequences from Genome Quebec above. This time hoever I first put all replicate 1 samples in /run1 and replicate 2 in /run2 and merge. I then put all males in /run1 and all females in /run2 and merge. I now have just two samples: C and E.

### Adding in read groups ###

My next step is to realign around indels. To do this I need read-groups added to my bam files. These can be added using Picard.

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

# The specifications here are just to add bland read-groups because GATK needs read-groups before it can run. I am basically just setting them to defaults just so they are there.
````

## Index for GATK ##

Indel realignment is done using GATK. This is a three step process. With the first script I index all files for GATK so they can be read properly. In the next script, indels are identified and marked. And finally in the third script, the indels are filtered out. Indels can create misalignments in my mapping file. To deal with this mapping error I first need to identify indels with GATK, and then realign around them to correct for mapping alignment errors due to indels.

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
The reference also needs to be indexed for the next steps.

````
samtools faidx dmel-all-chromosome-r6.23.fasta.gz
````

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

````
#! /bin/bash 
#Path to input directory
final_bam=/2/scratch/TylerA/SSD/bwamap/dedup

#Path to output directory
gatk_dir=/2/scratch/TylerA/SSD/bwamap/gatk

#Variable for reference genome (non-zipped)
ref_genome=/2/scratch/TylerA/Dmelgenome/gatk/dmel-all-chromosome-r6.23.fasta


#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar


files=(${final_bam}/*_rmd_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _rmd_RG.bam`

java -Xmx32g -jar ${gatk} -I ${final_bam}/${base}_rmd_RG.bam -R ${ref_genome} \
  -T IndelRealigner -targetIntervals ${gatk_dir}/${base}.intervals \
  -o ${gatk_dir}/${base}_realigned.bam

done
````

### mpileup and sync file ###

A single mpileup will be needed with my sample information for SNP calling and population statistics. I can create an mpileup using `samtools mpileup`. I will also need a sync file to work with for popoolation statistics such as Fst and for statstics testing in R.

I next create my mpileup for my samples:

````
samtools mpileup -B -Q 20 -f \
/2/scratch/TylerA/Dmelgenome/gatk/dmel-all-chromosome-r6.23.fa \
/2/scratch/TylerA/SSD/bwamap/gatk/*_realigned.bam \
> Treatments_combined.mpileup
````
This mpilup need to have the repetetive regions removed because these are regions where SNP false discovery is common. This can be done with popoolation commands.

````
# Create a GTF
 /usr/local/RepeatMasker/RepeatMasker -pa 10 -species drosophila -gff dmel-all-chromosome-r6.23.fasta

# Remove repetetive repetetive

perl /home/tylera/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf /2/scratch/TylerA/Dmelgenome/dmel-all-r6.23.gtf --input /2/scratch/TylerA/SSD/bwamap/Treatment_combined/Treatment_combined.mpileup --output /2/scratch/TylerA/SSD/bwamap/Treatment_combined/Treatment_combined_norepeat.mpileup

````

And I also create sync files for population statistic calculations:

````
java -ea -jar /usr/local/popoolation/mpileup2sync.jar --threads 16 --input Treatment_combined_norepeat.mpileup --output Treatment_combined_norepeat.sync
````

## Popoolation Statistics ##

### Fst for our Experimental population ###

Next, I want to calculate the Fst statistic for my experimental population compared to my control populations to look at structure.

````
perl /home/tylera/bin/popoolation2_1201/fst-sliding.pl \
				--input /2/scratch/TylerA/SSD/bwamap/Treatment_combined/Treatment_combined_norepeat.sync \
				--output /2/scratch/TylerA/SSD/bwamap/Treatment_combined/Treatment_combined.fst \
				--suppress-noninformative \
				--min-count 2 \
				--min-coverage 10 \
				--max-coverage 300 \
				--min-covered-fraction 1 \
				--window-size 1 \
				--step-size 1 \
				--pool-size 100
````

The resulting Fst file can be visualized in R: 

````
library(tidyverse)
library(data.table)
library(ggplot2)


ddat2 <- fread("~/Desktop/SSD/Data/sub.fst")

head(ddat2)

ccol <- ncol(ddat2)

head(ccol)

for (i in 33:ccol){
  ddat2[[i]] <- gsub(".*=","", ddat2[[i]])
}

for (i in 33:ccol){
  ddat2[[i]] <- as.numeric(ddat2[[i]])
}
 
ddat2$meanFst <- rowMeans(subset(ddat2, select = c(33:ccol)), na.rm = TRUE)

ddat2 <- ddat2[ddat2$meanFst!='NaN',]

head(ddat2)

ddat2 <- ddat2[,c(1,2,3,4,5,34)]

colnames(ddat2) <- c('chr', 'window', "num", 'frac', 'meanCov','meanFst')

ddat22L <- ddat2[which(ddat2$chr=='2L'),]
ddat22R <- ddat2[which(ddat2$chr=='2R'),]
ddat23L <- ddat2[which(ddat2$chr=='3L'),]
ddat23R <- ddat2[which(ddat2$chr=='3R'),]
ddat24 <- ddat2[which(ddat2$chr=='4'),]
ddat2X <- ddat2[which(ddat2$chr=='X'),]

ddat2 <- rbind(ddat2X, ddat22L, ddat22R, ddat23L, ddat23R, ddat24)

g <- nrow(ddat2[which(ddat2$chr=='2L'),])
h <- nrow(ddat2[which(ddat2$chr=='2R'),])
i <- nrow(ddat2[which(ddat2$chr=='3L'),])
j <- nrow(ddat2[which(ddat2$chr=='3R'),])
k <- nrow(ddat2[which(ddat2$chr=='4'),])
l <- nrow(ddat2[which(ddat2$chr=='X'),])

#X-2L-2R-3L-3R-4
ddat2$number <-  c((1:l),
                   (l+1):(l+g), 
                   (l+g+1):(l+g+h), 
                   (l+g+h+1):(l+g+h+i),
                   (l+g+h+i+1):(l+g+h+i+j),
                   (l+g+h+i+j+1):(l+g+h+i+j+k))

### PLOTS:

ggplot(ddat2, aes(x=number, y=meanFst, color=chr)) +
  geom_point(size=0.5, show.legend = F) +
  theme(panel.background = element_blank()) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
  xlab("Chromosome") +
  ylab(expression(F[ST])) +
  scale_x_discrete(limits=c(l/2, l+(g/2), (l+g+(h/2)), (l+g+h+(i/2)), (l+g+h+i+(j/2)), (l+g+h+i+j+(k/2))), labels = c("X","2L", "2R", '3L', '3R', "4")) +
  scale_colour_manual(values=c("seagreen", "darkslateblue", 'darkred', 'darkorchid4', 'darkolivegreen', 'darkblue')) +
  theme(text = element_text(size=20),
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15))
````

Single_bp_Fst.png![image](https://user-images.githubusercontent.com/77504755/114544388-8c376e80-9c28-11eb-884c-e6f3455c09d2.png)


## SNP calling ##

SNPs can be indentified using varscan and the following code:

````
java -Xmx32g -jar \
~/bin/VarScan.v2.3.9.jar \
mpileup2cns \
/2/scratch/TylerA/SSD/bwamap/Treatment_combined/Treatment_combined_norepeat.mpileup \
--min-coverage 50 \
--min-reads2 3 \
--p-value 0.1 \
--min-var-freq 0.01 \
--min-freq-for-hom 1 \
--min-avg-qual 20 \
--variants \
--output-vcf 1 \
| bgzip > Treatment_combined_norepeat.vcf

#mpileup2cns: this looks for bothh SNPs and indels
#--min-coverage 50 minimum coverage a read must have to be considered a true variant
#--min-reads 3 minimum reads an allele must have to be considered a true variant
#--p-value 0.1 A less stringet p-value is recomended for pooled data to make sure rare alleles are noyt ruled out
#--min-var-freq 0.01 Essentilly setting it so every single rare allele is examined
#--min-freq-for-hom 1 Saying that to be considered homozygoes a loci must be 100% one bp, to avoid loss of rare alleles
#--min-avg-qual 20 Filtering for quality again
#--variants Denotes that we are looking for variants
#--output-vcf 1 Outputs a vcf
#bgzip


````

Using this code varscan identified ######### SNPs.

## CMH Test ##

````
# R work
rm(list=ls())
#install.packages("/home/tylera/bin/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")
#install.packages("/home/tylera/bin/ACER-1.0.2.tar.gz")

#Loading in the required packages
library(poolSeq)
library(ACER)

#Load in my sync file created by popoolation2
mySync <- read.sync(file="/2/scratch/TylerA/SSD/bwamap/Experimental/sample_indels_repetetive.sync", gen=c(375, 375), repl=c(1, 1),polarization = "minor", keepOnlyBiallelic = TRUE)

#Create a matrix for allele frequences in replicate 1
af<-af(mySync, repl = 1, gen = 0)
af2<-af(mySync, repl = 1, gen = 400)
afMat<-as.matrix(af)
afMat2<-as.matrix(af2)

#Create a matrix of allele frequences for replicate 2
baf<-af(mySync, repl = 2, gen = 0)
baf2<-af(mySync, repl = 2, gen = 400)
bafMat<-as.matrix(baf)
bafMat2<-as.matrix(baf2)

#Verify the genomic positions line up in the replicate 1 and replicate 2 matrices
all(rownames(bafMat) == rownames(bafMat2))
all(rownames(afMat) == rownames(afMat2))

#Give more meaningful names to columns
colnames(afrep1) <- c("R1.C", "R1.E")
colnames(afrep2) <- c("R2.C", "R2.E")

#combine replicate 1 and replicate 2
afrep2<-cbind(bafMat,bafMat2)
afrep1<-cbind(afMat,afMat2)
afMat<-cbind(afrep1,afrep2)
 
#Repeat to create a matrix of coverage
cov<-coverage(mySync, repl = 1, gen = 0)
cov2<-coverage(mySync, repl = 1, gen = 400)
covMat<-as.matrix(cov)
covMat2<-as.matrix(cov2)
covrep1<-cbind(covMat,covMat2)

bcov<-coverage(mySync, repl = 2, gen = 0)
bcov2<-coverage(mySync, repl = 2, gen = 400)
bcovMat<-as.matrix(bcov)
bcovMat2<-as.matrix(bcov2)
covrep2<-cbind(bcovMat,bcovMat2)

colnames(covrep1) <- c("R1.C", "R1.E")
colnames(covrep2) <- c("R2.C", "R2.E")

covMat<-cbind(covrep1,covrep2)


# extracting only the chromosomes I want to look at and removing the random scaffolds
X<-afMat[grep("^X",rownames(afMat)),]
L2<-afMat[grep("^2L",rownames(afMat)),]
R2<-afMat[grep("^2R",rownames(afMat)),]
L3<-afMat[grep("^3L",rownames(afMat)),]
R3<-afMat[grep("^3R",rownames(afMat)),]
LR4<-afMat[grep("^4",rownames(afMat)),]
afdata<-rbind(X,L2,R2,L3,R3,LR4)

#Doing the same for coverage
cX<-covMat[grep("^X",rownames(covMat)),]
cL2<-covMat[grep("^2L",rownames(covMat)),]
cR2<-covMat[grep("^2R",rownames(covMat)),]
cL3<-covMat[grep("^3L",rownames(covMat)),]
cR3<-covMat[grep("^3R",rownames(covMat)),]
cLR4<-covMat[grep("^4",rownames(covMat)),]
covdata<-rbind(cX,cL2,cR2,cL3,cR3,cLR4)

## I got an error because there can not be 0s in coverage, so I use this code to see how many 0s I have in my coverage matrix
#mydatanew=covdata[,-1]                  # first gene name column deleted
#nonzero_row <- covdata[rowSums(covdata) == 0, ]  # filtered row read count above 0 
#dim(nonzero_row) #270

#Omitting NAs
#Also must omit 0s from coverage data, so I convert them to NAs and omit NAs from the coverage data as well.
afdata <- na.omit(afdata)
covdata[covdata==0] <- NA
covdata <- na.omit(covdata)

#dim(afdata)
#dim(covdata)


#create a new dataframe which contains only the names that match
#match <- covdata[rownames(covdata) %in% rownames(afdata),]

#Ensure my allele frequency and coverage matrices still line up
all(rownames(afdata) == rownames(match))

#change the name back to covdata when I've varified I didn't ess the command up
#match <- covdata


#Create the variables I will need for the CMH test
rep<-c(1,2) #Number of replicates
Ne<-500 #Estimated effective population size
tp<-c(0,400) #Generations of evolution for each sample
ps<-ps <- rep(100, 2*length(rep)) #Pool size

#Run the CMH command
pval <- adapted.cmh.test(freq=afdata, coverage=covdata, Ne=rep(Ne,
                         length(rep)), gen=tp, repl=rep, poolSize=ps)

#I get this warning, but the command appeared to run. I will need to verify that it was run correctly.
#Warning messages:
#1: In adapted.cmh.test(freq = afdata, coverage = covdata, Ne = rep(Ne,  :
#  Ne value(s) which are not integer are converted to integer
#2: In adapted.cmh.test(freq = afdata, coverage = covdata, Ne = rep(Ne,  :
#  The counts assuming values 0 or equal the coverage of the considered 
#                      locus are changed to 1 and to coverage-1 respectively.
#3: In adapted.cmh.test(freq = afdata, coverage = covdata, Ne = rep(Ne,  :
#  The counts assuming values 0 or equal the coverage of the considered 
#                      locus are changed to 1 and to coverage-1 respectively.

#Convert the output to a matrix to work with
p<-as.matrix(pval)

#Run a Benjamini Hochberg estimate to control for FDR
padj<-p.adjust(p, method = "hochberg", n = length(p))

#Look at the number of loci with significant p-values
length(which(padj<0.05)) # 30349
length(which(padj<0.01)) # 28694

#Re-attach the SNP locations
afp<-cbind(afdata,padj)
colnames(afp) <- c("R1.C", "R1.E", "R2.C", "R2.E", "padj")

library(dplyr)
df<-as.data.frame(afp)
loci<-filter(df,padj<0.01)

#How many interesting loci do I have?
dim(loci) #28694


````













