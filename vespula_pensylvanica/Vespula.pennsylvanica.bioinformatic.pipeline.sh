# Vespula pennsylvanica - WASP LINKAGE MAP BIOINFORMATIC PIPELINE

# Bioinformatics of linkage mapping of various wasp species 

# In association with the publication "Comparative linkage mapping of wasp species"

# Author: Daniela Zarate, PhD.
________________________________________________________________________________________________________
________________________________________________________________________________________________________
# Important short cuts I should keep handy:
# branch to a new partition and work on an interactive command line 
srun -p short --pty bash -l 
squeue -u danielaz # check on status of jobs 
scancel -u danielaz
________________________________________________________________________________________________________
# full path to reference

REFERENCE=/rhome/danielaz/bigdata/wasps/Vpen_REF/VpenV1.fa

# Find the original raw data here: 
/bigdata/brelsfordlab/abrelsford/purcell/EcoRIrad/wasp_link/Vpen/plate1
/bigdata/brelsfordlab/abrelsford/purcell/EcoRIrad/wasp_link/Vpen/plate2

python sizeMatters.py

The number of 5M < files: 0
The number of 5M <= files < 20M : 0
The number of 20M <= files < 50M : 56
The number of 50M <= files < 100M : 20
The number of files >= 100M : 0

# Current personal working directory:
/rhome/danielaz/bigdata/wasps/vespula_pennsylvanica
________________________________________________________________________________________________________

# compile the plate 1 and plate 2 sample names (n - 161)

cat samples.plate1.txt samples.plate2.txt > all.samples.txt

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.
# Use the pearBwaTool.sh script to process files in tandem with loopsub.sh 
# This will use PEAR to merge paired end reads and then BWA to map with a reference genome

vi pearBwaTool.sh
chmod +x pearBwaTool.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-05:00:00  
#SBATCH --output=PLACEHOLDER.pearBwaTool.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --job-name="PLACEHOLDER-pearBwaTool-job"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short

# Increase limit for jobs 
ulimit -n 9999

# Print current date
date

# Load software
module load pear
module load bwa-mem2
module load samtools/1.16


# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

# Absolute paths to reference genomes
REFERENCE=/rhome/danielaz/bigdata/wasps/Vpen_REF/VpenV1.fa

# keep files compressed, but if they are uncompressed, remember to change the file extensions (remove .gz below)

# use PEAR to merge paired ends and remove the adaptor sequence
# PEAR is an ultrafast, memory-efficient and highly accurate pair-end read merger.

# RE RUN !!!
pear -f PLACEHOLDER.1.fq.gz -r PLACEHOLDER.2.fq.gz -o PLACEHOLDER --threads 8 # -k -q 20 -t 100 

# -f: forward read
# -r: reverse read
# -o: output name 
# --keep-original: Do not reverse and complement the reverse reads when writing \
# the unassembled and discarded reads output.
# Do not need for rads, need for whole genomes (check?)
# outputs assembled, discarded, and unassembled forward and reverse fastq files
# -t: Specify the minimum length of reads after trimming the low quality part
# -q: Specify the quality score threshold for trimming the low quality part of a read. If the quality scores of two consecutive bases are strictly less than the specified threshold, the rest of the read will be trimmed

# combine read2 remnant singles with pear assembled singles
zcat PLACEHOLDER.rem.2.fq.gz >> PLACEHOLDER.assembled.fastq

# For BWA compatibility, edit the names of the Read 2 reads to match read 1
# Unassembled forward reads have a /1 at the end of the read name whereas unassembled reverse have /2 

cat PLACEHOLDER.unassembled.reverse.fastq | sed 's|_2$|_1|' > PLACEHOLDER.converted.fastq

# assenbled have _1 
# unassembled have _2 
# converted should have _1 
# unassembled_forward is _1 
# unassembled_reverse is _2 

# Map unassembled forward and reverse reads to reference genome use BWA and then sort with Samtools
bwa-mem2 mem -t 8 ${REFERENCE} PLACEHOLDER.unassembled.forward.fastq PLACEHOLDER.converted.fastq | \
samtools sort - > PLACEHOLDER.pair.bam

# Map assembled reads to reference genome use BWA and then sort with Samtools
bwa-mem2 mem -t 8 ${REFERENCE} PLACEHOLDER.assembled.fastq | \
samtools sort - > PLACEHOLDER.single.bam

# -S : ignored for compatability with previous samtools versions, Previously requiered if input was in SAM format, but now \
# format is automatically detected by examining first characters of input
# --with-header: -h : include the header in the output
# --uncompressed : -u : uncompressed

# merge the two bams
samtools merge PLACEHOLDER.all.bam PLACEHOLDER.pair.bam PLACEHOLDER.single.bam -@ 8
# -@ : number of threads 

#index the final bam
samtools index PLACEHOLDER.all.bam

# Print name of node
hostname
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

# Then use loopsub on command line to submit all the individual scripts to hpcc 
# alternatively, use bash loopsub.sh on the command line to submit as script

while read i ; do sed "s/PLACEHOLDER/$i/g" pearBwaTool.sh  > pearBwaTool.$i.sh; sbatch pearBwaTool.$i.sh ; done < V.pennsylvanica.samples.txt

The average assembled reads rate is: 34.690074534161475%
The total number of assembled reads is: 43311375
The average discarded reads rate is: 0.0%
The total number of discarded reads is: 0
The average assembled reads across individuals is 269014.75155279506
The average discarded reads across individuals is 0.0

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

cd bams
# Compile all bam file names in one txt file:
ls *bam > V.pennsylvanica.bam.names.txt

vi mpileup.sh 
chmod +x mpileup.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-30:00:00     # 6 hours; adjust as needed
#SBATCH --output=V.pennsylvanica.mpile.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=END
#SBATCH --job-name="V.pennsylvanica.mpileup.log"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# This script used in the linkage map pipeline analysis after creating BAM files \
# using pearBwaTool.sh scripts. This script takes samtools and bcftools to call SNPS \
# and generate VCF files.

# Print current date
date

# Load software
module load samtools/1.16
module load bcftools

# Use the command if running on more than 500 samples:
# If not, leave it hashtagged out 
# ulimit -n 9999

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR
REFERENCE=/rhome/danielaz/bigdata/wasps/Vpen_REF/VpenV1.fa
BAMLIST=/rhome/danielaz/bigdata/wasps//vespula_pennsylvanica/bams/V.pennsylvanica.bam.names.txt

# run bcftools mpileup on whole genome. Requires a file bamlist.txt with paths to all bam files.

bcftools mpileup \
 --output-type u \
  --min-MQ 20 \
  --annotate DP  \
  --fasta-ref ${REFERENCE} \
  --bam-list ${BAMLIST} | \
  bcftools call \
  --output-type z \
  --multiallelic-caller \
  --variants-only \
  --format-fields GQ \
  --output V.pennsylvanica.vcf.gz

# bcftools mpileup commands: 
# Generate VCF or BCF containing genotype likelihoods for one or multiple alignment files 
# --output-type :-O : option (u) : -Ou: -> uncompressed BCF ; USE this option when piping between BCF subcommands \
# in order to speed up performance 
# --min-MQ INT : -q :  Skip alignments with mapQ smaller than INT [0]
# --annotate : -a : Comma-seperated list of FORMAT and INFO tags to output \
# DP = FORMAT/DP = Number of high-quality bases (Number=1,Type=Integer)
# --fasta-ref : -f : reference genome
# --regions: -r : provide region to specify genomic region to focus on, leave out if running on whole genome
# --bam-list : -b : list of bam files

# bcftools call commands: 
# --output-type: -O : b|u|z|v : -Oz: compressed VCF (z)
# --multiallelic-caller : -m : Alternative model for multiallelic and rare-variant calling (conflicts with -c
# --variants-only: -v:  Output variant sites only
# --format-fields: -f : comma-seperated list of format fields to output for each sample, currently GQ and GP \
# fields are supported
# --output : -o: output file name  
# --skip-variants indels 

# print node name
hostname
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.     

# Filter for minDP and maxMissing (n=161)

vcftools --gzvcf V.pennsylvanica.vcf.gz  --minDP 1 --max-missing 0.8 --out V.pennsylvanica --recode 

# --minDP = Includes only genotypes greater than or equal to the "--minDP" value; requires that the "DP" FORMAT tag is specified for all sites.
# --max-missing = Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).

# After filtering,  kept 161 out of 161 Individuals
# After filtering, kept 22,609 out of a possible 93,124 Sites
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Remove individuals with >25% missing data 
# check out the statistics for missing data for each individual 

vcftools --vcf V.pennsylvanica.recode.vcf --missing-indv --out V.pennsylvanica

# Generates a file reporting the missingness on a per-individual basis. The file has the suffix ".imiss".
# Cut-off for exclusind individuals with large chunks of missing data: 25%

# Use the missingIndv.filter.py script to output data that is above and below the threshold:

python missingIndv.filter.py --maxMISS 0.25 --path . 
# outputs : goodData.txt (data of samples that pass threshold)
# outputs : badData.txt (data of samples that do NOT pass threshold)
# outputs : badIndv.txt (list of samples that do NOT pass threshold) -> to be piped into next command

less V.pennsylvanicaq.imiss.badData.txt
# No individuals removed 
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

module load plink
plink --vcf V.pennsylvanica.recode.vcf --pca --aec --out V.pennsylvanica
# --pca = principal component analysis
# --out = specifies output name 
# Eigenvectors written to polyergus.A.eigenvec & eigenvalues written to polyergus.A.eigenval 
# Eigenval file: 10 lines 
# Eigenvec file: 96 lines of PC 1-10 
# run time: < 1 min. 

awk -v OFS='\t' '{ print $2, $3, $4, $5 }' V.pennsylvanica.v2.eigenvec > V.pennsylvanica.v2.eigenvec.cut 
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/wasps/vespula_pennsylvanica/vcfs/V.pennsylvanica.v2.eigen*' .

M148.all.bam
vcftools --vcf V.pennsylvanica.recode.vcf --remove-indv M148.all.bam --out V.pennsylvanica.plink.indv --recode 

plink --vcf V.pennsylvanica.plink.indv.recode.vcf  --pca --aec --out V.pennsylvanica.v2

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# filter for minor allele frequency 

# Allele frequency is defined as the number of times an allele appears over all individuals at that site, \
# divided by the total number of non-missing alleles at that site.
# --maf <float> : Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value 
# --mac : minor allele count only use for small families 
# --maf of 0.15 or .1 is fine. 

vcftools --vcf V.pennsylvanica.plink.indv.recode.vcf --out V.pennsylvanica.maf --maf 0.15 --recode 
# After filtering, kept 8,944 out of a possible 22,609 Sites
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# remove Indels and keep only SNPs as variants 

vcftools --vcf V.pennsylvanica.maf.recode.vcf --out V.pennsylvanica.snps --remove-indels --recode 
# After filtering, kept 7804 out of a possible 8944 Sites
# Don't do this because we're trying to keep as many as possible 
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# In order to identify heterozygous sites, we used VCFtools --hardy flag :
vcftools --vcf V.pennsylvanica.maf.recode.vcf --hardy --out V.pennsylvanica

# --hardy: 
# Reports a p-value for each site from a Hardy-Weinberg Equilibrium test (as defined by Wigginton, \
# Cutler and Abecasis (2005)). The resulting file (with suffix ".hwe") also contains the Observed \
# numbers of Homozygotes and Heterozygotes and the corresponding Expected numbers under HWE.
# GT Calls: (0/0) = homozygous for reference , (0/1) = heterozygote, (1/1) = homozygous for alternate, (./.) = missing data

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Start off with males (n = 160)
# Then filter out sites where heterozygote allele count is greater than 5% of total allele count (excessive heterozygosity)
cut -f 1-3 V.pennsylvanica.hwe | tr '/' '\t' | awk '$4>0.08*($3+$4+$5)' | cut -f 1-5 > V.pennsylvanica.0.8.badSNPs.txt
wc -l V.pennsylvanica.*.badSNPs.txt # Bad SNPs = 6105 if its 5% 

# JP's had 3.6 K SNPs at the end. 

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.
# Retain only high quality sites by excluding the shitty sites.
vcftools --vcf V.pennsylvanica.maf.recode.vcf --exclude-positions V.pennsylvanica.1.0.badSNPs.txt --out V.pennsylvanica.filtered --recode 
# After filtering, kept 3,263 out of a possible 8944 Sites

# For the final filtered file, convert all heterozygous sites (0/1) to (./.) for missing data
sed 's+0/1+./.+g' V.pennsylvanica.filtered.SNPs.recode.vcf > V.pennsylvanica.converted.SNPs.recode.vcf
# To check to make sure this was successful, run --hardy again:
vcftools --vcf V.pennsylvanica.converted.SNPs.recode.vcf --hardy --out V.pennsylvanica.converted
# Review the .hwe file, the sites should all be balanced for OBS HOM1 and OBS HOM2 and a 0 for hetero calls: e.g. 36/0/36
less V.pennsylvanica.converted.hwe

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# recode to a 012 format 
vcftools --vcf V.pennsylvanica.converted.SNPs.recode.vcf --012 --out V.pennsylvanica

# This option outputs the genotypes as a large matrix. Three files are produced. \
# The first, with suffix ".012", contains the genotypes of each individual on a separate line.\
# Genotypes are represented as 0, 1 and 2, where the number represent that number of \
# non-reference alleles. Missing genotypes are represented by -1. The second file,\
# with suffix ".012.indv" details the individuals included in the main file. \
# The third file, with suffix ".012.pos" details the site locations included in the main file.
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

vi V.pennsylvanica.header.txt

# INPUT FILES
V.pennsylvanica.012 
V.pennsylvanica.012.indv 
V.pennsylvanica.012.pos

# METHOD
Males

# MCT HEADER
population_type DH
population_name pennsylvanica
distance_function       kosambi
cut_off_p_value 0.000001
no_map_dist     30
no_map_size     1
missing_threshold       0.2
estimation_before_clustering    no
detect_bad_data yes
objective_function      ML

# These parameters no longer needed in header file, as they are automatically stripped from data file:
# number_of_loci   3263
# number_of_individual   160

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Use vcf2mst2.0.py to generate a MST INPUT file from the final VCF:
# Then run MSTmap with the output file, using the mstmap singularity image:
# Last, run LinkageGroupStats.py for scaffold stats: 

module load singularity

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

## RUN 1: p = 0.000001 
singularity run pandas.simg vcf2mst3.0.py --header V.pennsylvanica.header.txt --output V.pennsylvanica.mst.R1.txt 
singularity run mstmap.simg V.pennsylvanica.mst.R1.txt V.pennsylvanica.mst.out.R1.txt  | tee V.pennsylvanica.mst.R1.log
singularity run pandas.simg LinkageGroupsStats.py --input V.pennsylvanica.mst.out.R1.txt --output V.pennsylvanica.mst.R1.stats --threshold 5 --moderate 30 --severe 40 

# Out of 50 linkage groups, 48 had more than 5 scaffolds, out of those 48 groups, 2 had moderate to severe gaps. 
# There are a total of 2 severe gaps, and 0 moderate gaps across all linkage groups
# 25 LGs after splitting at gaps, resolves into expected karyotype


__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

## RUN 2: p = 0.000005 
singularity run pandas.simg vcf2mst3.0.py --header V.pennsylvanica.header.txt --output V.pennsylvanica.mst.R2.txt 
singularity run mstmap.simg V.pennsylvanica.mst.R2.txt V.pennsylvanica.mst.out.R2.txt  | tee V.pennsylvanica.mst.R2.log
singularity run pandas.simg LinkageGroupsStats.py --input V.pennsylvanica.mst.out.R2.txt --output V.pennsylvanica.mst.R2.stats --threshold 5 --moderate 30 --severe 40 

# Out of 48 linkage groups, 46 had more than 5 scaffolds, out of those 46 groups, 4 had moderate to severe gaps. 
# There are a total of 4 severe gaps, and 0 moderate gaps across all linkage groups
# 25 LGs after splitting at gaps, resolves into expected karyotype

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

## RUN 3: p = 0.0000001
singularity run pandas.simg vcf2mst3.0.py --header V.pennsylvanica.header.txt --output V.pennsylvanica.mst.R3.txt 
singularity run mstmap.simg V.pennsylvanica.mst.R3.txt V.pennsylvanica.mst.out.R3.txt  | tee V.pennsylvanica.mst.R3.log
singularity run pandas.simg LinkageGroupsStats.py --input V.pennsylvanica.mst.out.R3.txt --output V.pennsylvanica.mst.R3.stats --threshold 5 --moderate 30 --severe 40 

# Out of 52 linkage groups, 50 had more than 5 scaffolds, out of those 50 groups, 0 had moderate to severe gaps. 
# There are a total of 0 severe gaps, and 0 moderate gaps across all linkage groups

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Run 3 emerges as the most appropriate. Going with this one 

# Review the stats log to determine where the gaps are.
# Download the mst output file to build the linkage map. 
# Choose the best LGs of the duplicates and build in excel. 
# Create a "stripped" version that is ready for use in the LinkageMappingWizard.py script. 

_.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Run the Linkage Map Wizard python script. 

# load singularity so you can run the singularity image with the bioinformatic tools necessary to run the LinkageMapWizard: 
module load singularity  

# run the image, calling on python3, providing the stripped xlsx file and an argument that tells MatPlotLib to only graph scaffolds with a minimum of X markers:
# To run for both scaffolds and contigs use the LinakgeMapWizard.py
# To run for just scaffolds, use the LinkageMapWizard.1.0.py

singularity run BioInforMaticGizmo.simg python3 ./LinkageMapWizard.py --input  v.pennsylvanica.original.linkage.map.xlsx --scaff_thresh 3

# It will output two files, one file the ".condensed" that is the .txt file of the centimorgan x position file used to plot the scatter plot and 
# the .png file which is the scatterplot of the linkage groups
V.pennsylvanica.inverted.linkageMap.xlsx.png
V.pennsylvanica.inverted.linkageMap.xlsx.condensed

scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/wasps/vespula_pennsylvanica/linkageMap/v.pennsylvanica.ori*' .

## Figure out which scaffolds need to be inverted and then flip them in excel and then run again: 
# Use the Inverted Linkage Map Wizard

singularity run BioInforMaticGizmo.simg python3 ./InvLinkageMapWizard.py  --input  V.pennsylvanica.inverted.linkageMap.xlsx  --scaff_thresh 3

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Split the output of continuous centiMorgan and position into seperate files, one per scaffold
# I just do this manually for now. 
# I do it by downloading the .condensed file and then copying and pasting it in editor. 

cd scaffs 

# I also need to add the header to all of the files 
# Again, doing this manually for now. 

Scaffold        centiMorgan     Position


# Now, we need to prune the scaffolds so that only one rad tag per 150 bp is retained. 
# Do this manually for now, but eventually should make this automated via python. 
# Use vi to delete every 2nd + lines of the same rad tag. 
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__. 

# OK, doing it manually is super annoying, here's a nice python script:

while read i ; do sed "s/PLACEHOLDER/$i/g" prunePlusCalc.sh > prunePlusCalc.$i.sh; sbatch prunePlusCalc.$i.sh ; done<Scaffolds.txt

vi prunePlusCalc.sh
chmod +x prunePlusCalc.sh

ls Scaffold* > list.scaffolds
sed 's/\.txt$//' list.scaffolds > Scaffolds.txt

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-01:00:00     # 6 hours; adjust as needed
#SBATCH --output=PLACEHOLDER.std
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --job-name="PLACEHOLDER"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# This script is used to call the prune_rads.py script to prune markers to one per radtag. 
# For use before calculating recombination. 

# Print current date
date

# Load software
module load singularity

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

singularity run BioInforMaticGizmo.simg python3 prune_rads.py --input PLACEHOLDER.txt --output PLACEHOLDER --scaffold PLACEHOLDER
singularity run BioInforMaticGizmo.simg python3 slope_calc.py --input PLACEHOLDER.pruned.txt --output PLACEHOLDER

# print node name
hostname

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__. 

# Download all the pruned scaffolds and recombination rates to local to graph in R. 

scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/wasps/vespula_pensylvanica/linkageMap/scaffs/pruned_scaffolds/Scaffold*.txt' .
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/wasps/vespula_pensylvanica/linkageMap/scaffs/recombination/Scaffold*.recombination.txt' .

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__. 

singularity run BioInforMaticGizmo.simg python3 prune_rads.py --input V.pensylvanica.inv.linkageMap --output V.pensylvanica.inv.linkageMap --scaffold 

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__. 

# Calculate recombination again for reviewer's suggestion:

cd /rhome/danielaz/bigdata/wasps/vespula_pensylvanica/linkageMap/scaffs

# METHOD ONE 
# Calculate all pairwise comparisons without averaging.

singularity run BioInforMaticGizmo.simg python3 melee_slope_calc.py --input test.txt --output test

# Submit in batch, all of the scaffolds: 

while read i ; do sed "s/PLACEHOLDER/$i/g" prunePlusCalc.sh > prunePlusCalc.$i.sh; sbatch prunePlusCalc.$i.sh ; done<Scaffolds.txt

# Copy them to local:
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/wasps/vespula_pensylvanica/linkageMap/scaffs/pruned_scaffolds/Scaffold*.recombination.txt' .


# METHOD TWO
# Calculate a weighted average per window. 

singularity run BioInforMaticGizmo.simg python3 weighted_calc.py --input test.txt --output test
while read i ; do sed "s/PLACEHOLDER/$i/g" prunePlusCalc.sh > prunePlusCalc.$i.sh; sbatch prunePlusCalc.$i.sh ; done<Scaffolds.txt
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/wasps/vespula_pensylvanica/linkageMap/scaffs/pruned_scaffolds/Scaffold*.weighted.recombination.txt' .
