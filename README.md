# Detecting disease associated pathways by whole-genome sequencing data. 
The analysis process mainly includes 6 steps and depending on python or R. As follows, we provide a complete and concise procedure code along with explanations, facilitating researchers to conduct personalized analyses. Copy the code to your script as according to your needs and make the necessary modifications as prompted (e.g., if your genome sequencing data does not require processing, you may ignore the first step. Alternatively, if you adopt an HHG P value < 0.001 as the threshold for significant association between pathways and disease risk, you may skip the permutation calculation in step 5.6).

## Step 1. WGS —> Genotype ####
In order to perform polygenic risk score (PRS) analysis and gene expression prediction, whole genome sequencing data in the format of .bed, .bim, .fam is required.
Here's an outline of the processing whole-genome sequencing (WGS) raw sequencing data into .bed, .bim, and .fam formats

### 1.1	Software and reference genome sequence preparation.
Download and install the following software: BWA, Samtools, GATK4.0, fastp, plink.
Download the reference genome sequence (for example: GRCh38.fa).

### 1.2	Quality control (QC) of raw sequencing data
QC and remove sequencing adapters and low-quality sequences
#Example code
fastp -i Sample1_1.fastq.gz -o Sample1_1.clean.fastq.gz -I Sample1_2.fastq.gz -O Sample1_2.clean.fastq.gz -f 5 -t 5 -F 5 -T 5 -5 -W 5 -M 20 -Q -l 50 -c -w 4
Description of some parameters: -i: the input file name of read1. -o: the output file name of read1. -I: the input file name of read2. -O: the output file name of read2. -f 5 -t 5 -F 5 -T 5: global trimming options. Trim low-quality base number parameters at the beginning and end of the sequence. -5 -W 5 -M 20: per read cutting by quality options. -Q: --disable_quality_filtering. -l, --length_required. -c, --correction. -w, --thread.

