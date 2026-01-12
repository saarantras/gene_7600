**Gene 760 – Problem Set 3**

The purpose of this problem set is to familiarize students with the analysis, manipulation and visualization of ChIP-seq and ATAC-seq data. The datasets we are using were obtained from Inoue et al., Cell Stem Cell 25: 713 (2019). The goal of this study was to characterize gene regulatory dynamics during neural induction. A copy of the manuscript is available in the PS3 directory on Canvas.

By the end of this problem set, you will have learned how to use:
-	The SRA Toolkit to download and extract FASTQ files
-	FastQC and Cutadapt to process reads
-	Bowtie2 to map short reads to a genome
-	SAMtools to manipulate BAM/SAM files
-	MACS2 to call peaks in ChIP-Seq and ATAC-seq data
-	BEDTools to manipulate interval data (e.g., merging and intersecting peak calls from different datasets)

Students are to submit a gzipped tarball called *<Your_NetID>_PS3.tar.gz* containing:
- *<Your_NetID>_PS3_answers.md: Answers to the questions below*
- *<Your_NetID>_H3K27ac_t0_rep1_peaks.gz*
- *<Your_NetID>_H3K27me3_t0_rep1_peaks.gz*
- *<Your_NetID>_ATACseq_t0_rep1_peaks.gz*
- The intersected and merged H3K27me3 bed files in Question 7
- *<Your_NetID>_H3K27ac_shared_t0_t72*
- *<Your_NetID>_ATACseq_shared_t0_t72*

to `/gpfs/gibbs/project/gene760/shared/Dropbox/PS3` on **Sunday, March 9 at 11:59 pm**.

You will be running commands and processing all data in `palmer_scratch` space. This provides more usable – and temporary - disk space for large analyses than the ‘project’ directory. However, keep in mind that _all files in scratch60 are deleted after 60 days, so don’t use it for permanent storage! Your home directory does not have sufficient storage space for all your intermediate files, which may cause your processes to stall.

**This problem set is computationally intensive! Do not run on the login node!**
We highly recommend using Slurm/SimpleQueue to run jobs in parallel and take advantage of more than one core (multithreading). Note, you will need to use specific command line options if you want to run Bowtie2 and SAMtools in “multithread” mode.

**Question 2 will take a long time to run. You may want to consider working on Questions 7-10 in the meantime as the necessary files are provided.**

For this problem set, you will be working with time-series ChIP-seq and ATAC-seq data obtained at specific time points during induction of human neural stem cells from embryonic stem cells. The targets in the ChIP-seq experiments are H3K27ac, which marks active promoters and enhancers, and H3K27me3, which marks repressed regions.  These data are available at GEO under the accession number GSE115046. Please look at this record to familiarize yourself with the GEO format.

To save time and effort, we are only using a subset of the data in the paper. We are also only asking you to call peaks for each dataset at a single time point (t0). We are providing peak calls for all time points in BED format and in hg19 coordinates.

The datasets are located in `/gpfs/gibbs/project/gene760/shared/folders/Data/PS3`-k
In this directory, you will find:

1. fastq files containing raw reads for peak calling:
•	K27ac_rep1_0hr.fastq.gz
•	K27me3_rep1_0hr.fastq.gz
•	Input_rep1_0hr.fastq.gz
•	ATACseq_rep1_0hr_1.fastq.gz
•	ATACseq_rep1_0hr_2.fastq.gz

2. A subdirectory called 'Peaks' with BED files containing peak calls:
•	There are 12 files in this subdirectory, containing processed peak calls for ATAC-seq, H3K27ac and H3K27me3 data at the start of induction and at the end (0 hours and 72 hours). Two replicates (designated rep1 and rep2 in the file name) were carried out for each experimental condition and target.

For alignments, pre-built Bowtie2 index files for hg19 are hosted by YCRC at: 
‘/gpfs/gibbs/data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index’

#Quality control and reads alignment

###Question 1
Check the quality of the reads (all fastq.gz files) via FastQC and apply Cutadapt to process the ATAC-seq reads. Use default parameters for FastQC (-o flags) and Cutadapt (-a, -A, -o, -p flags).

**ANS**: Write down the commands you used and a description of what each command does.

Commands: 

module load FastQC

fastqc -o . *.fastq.gz

Description: 
First command loads the FastQC module
Second uses a wildcard on files ending in .fastq.gz and uses them as input for the fastqc command to check the quality of the reads. The output directory is specified by the argument -o as the current directory (./)

Command: 
module load cutadapt

[can also do with unzipped files]
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o ATACseq_rep1_0hr_1_trimmed.fastq -p ATACseq_rep1_0hr_2_trimmed.fastq ATACseq_rep1_0hr_1.fastq.gz ATACseq_rep1_0hr_2.fastq.gz

Description: 
First command loads the cutadapt module

Second command trims the Nextera Tansposase adaptor sequence (CTGTCTCTTATACACATCT) from the paired-end (-p) read input files: ATACseq_rep1_0hr_1.fastq and ATACseq_rep1_0hr_2.fastq and saves the trimmed output fastq file in the same directory. The remaining input files did not require trimming based on the Fastqc report.



[not required for answer but TA reference of summary output]:
Total read pairs processed:         42,326,401
  Read 1 with adapter:              29,973,583 (70.8%)
  Read 2 with adapter:              29,837,875 (70.5%)
Pairs written (passing filters):    42,326,401 (100.0%)
Total basepairs processed: 8,549,933,002 bp
  Read 1: 4,274,966,501 bp
  Read 2: 4,274,966,501 bp
Total written (filtered):  6,055,987,640 bp (70.8%)
  Read 1: 3,025,244,351 bp
  Read 2: 3,030,743,289 bp



**ANS**: Why are we using Cutadapt and what does it do?
Cutadapt trims input fastqc files with overrepresented adaptor sequences. These sequences need to be removed for quality control purposes as they are not representative of biologically relevant genomic sequences and therefore may interfere with downstream parts of the analysis (e.g. alignment).



###Question 2

Generate alignments for each dataset (H3K27ac, H3K27me3, ATAC-seq and Input) using Bowtie2. Please note that the parameters for ATAC-seq reads alignment should be different from that for H3K27ac AND H3K27me3 (paired end vs single end). Many genomics tools report useful information to “standard error” (stderr). For each alignment, save the stderr to a log file. Use SAMtools to convert the output to a BAM file. Bowtie2 flags to be used (-X 500, --very-sensitive, -N 1, -x, -U or -1/-2, -S, and 2>). SAMtools flags to be used (-b -q 30).

**ANS**: Write the commands you used and a description of what each command does.

Commands (descriptions in comments):

# load alignment module
module load Bowtie2

# alignment command (bowtie2) using 8 cores for parallel processing (-p), increasing the speed of the alignment, -X 500 is the maximum fragment length in bp for paired-end alignemnts, --very-sensitive specifies the alignment type (more accurate but slower alignment option). -N is the number of mismatches allowed during seeding in a multiseed alignment (either 0 or 1); setting this to 1 instead of 0 makes the alignment slower but increases sensitivity. -x specifies the directory to the Bowtie reference genome index; this is a binary representation of the genome from which sequence info can quickly and effecienty be extracted from during the alignment. -U specifies the input for single-end reads and -1 and -2 for paired end input

# for single-end files: H3K27ac, H3K27me3, Input
bowtie2 -p 16 --very-sensitive -N 1 -x /gpfs/gibbs/data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome -U $InputDir/K27ac_rep1_0hr.fastq.gz -S $OutputDir/K27ac_rep1_0hr.sam 2> K27ac_rep1_0hr.log

bowtie2 -p 16 --very-sensitive -N 1 -x /gpfs/gibbs/data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome -U $InputDir/K27me3_rep1_0hr.fastq.gz -S $OutputDir/K27me3_rep1_0hr.sam 2> K27me3_rep1_0hr.log

bowtie2 -p 16 --very-sensitive -N 1 -x /gpfs/gibbs/data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome -U $InputDir/Input_rep1_0hr.fastq.gz -S $OutputDir/Input_rep1_0hr.sam 2> Input_rep1_0hr.log

# for paired-end files: ATAC
bowtie2 -p 8 -X 500 --very-sensitive -N 1 -x /gpfs/gibbs/data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome -1 ATACseq_rep1_0hr_1_trimmed.fastq -2 ATACseq_rep1_0hr_2_trimmed.fastq -S ATACseq_rep1_0hr.sam 2> ATACseq_rep1_0hr.log

# load module to manipulate sam files, including sam to bam conversions
module load SAMtools

# Convert SAM to BAM file 
# -b defined output to be BAM file, -q sets a threshold of 30 to exclude alignment reads with qualities below 30 and -o is the output bam file name corresponding to the .sam input file
samtools view -b -q 30 -o K27ac_rep1_0hr.bam K27ac_rep1_0hr.sam
samtools view -b -q 30 -o K27me3_rep1_0hr.bam K27me3_rep1_0hr.sam
samtools view -b -q 30 -o Input_rep1_0hr.bam Input_rep1_0hr.sam
samtools view -b -q 30 -o ATACseq_rep1_0hr.bam ATACseq_rep1_0hr.sam


**ANS**: Report the number and percent of reads that aligned concordantly exactly 1 time for each alignment.
Input_rep1_0hr: 25804858 (70.52%) aligned exactly 1 time
K27ac_rep1_0hr: 23375853 (75.85%) aligned exactly 1 time
K27me3_rep1_0hr: 22923988 (66.21%) aligned exactly 1 time
ATACseq_rep1_0hr: 15768069 (37.27%) aligned concordantly exactly 1 time


**ANS**: What is the difference between a SAM file and a BAM file? What are the advantages or disadvantages of each?

SAM files are human-readable and store sequence alignment outputs, with each line representing a single alignment. It includes information about the reads mapped such as the read ID, genomic location, and quality. An advantage of them being human readable is that they can be used for debugging and direct inspection. A disadvantage is that they are large files and therefore less space efficient.

BAM files are compressed binary versions of SAM, and are therefore smaller (take up less space) and can also be indexed to access specific regions/positions of the alignment output. A disadvantage of BAM files is that they are not human readable they require specific software such as SAMTools to visualize and manipulate. An advantage is that since they are compressed they can be used for analysis of large amounts of input reads with better space/storage efficiency. 

###Question 3

Use SAMtools to remove the reads mapped to mitochondria from only the ATAC-seq bam file. Hint: open each bam file with SAMtools and then pipe to “grep -v MT”.

**ANS**: Write the commands you used and a description of what each command does.
[write answer here]

module load SAMtools

Command:
samtools view -h ATACseq_rep1_0hr.bam | grep -v "MT" | samtools view -b -o ATAC_seq_rep1_0hr_noMT.bam -

Description: 
The view command from the samtools module allows access to the binary .bam file by presenting it in a human readable format, this output is piped into the grep command that will exclude (-v) any line that includes "MT", while keepin the header (-h), thereby excluding mitochondrial reads, and these nonmitochondrial reads will be saved in a binary (-b) output (-o) file: ATAC_seq_rep1_0hr_noMT.bam

**ANS**: Explain why mitochondria reads needs to be removed in ATAC-seq analysis

ATAC-seq analyzes areas of open chromatin in nuclear/genomic DNA. Mitochondrial DNA (not bound by histones) are considered noise in this case, are not relevant to this analysis and are therefore excluded.

###Question 4

Using SAMtools, sort and create an index for each BAM file. SAMtools flags to be used for sorting (-o), indexing does not require any flags.

**ANS**: Write the SAMtools commands you used to sort and index the BAM files.

samtools sort ATAC_seq_rep1_0hr_noMT.bam -o	ATAC_seq_rep1_0hr_noMT_sorted.bam
samtools sort Input_rep1_0hr.bam -o	Input_rep1_0hr_sorted.bam
samtools sort K27ac_rep1_0hr.bam -o	K27ac_rep1_0hr_sorted.bam
samtools sort K27me3_rep1_0hr.bam -o K27me3_rep1_0hr_sorted.bam

samtools index ATAC_seq_rep1_0hr_noMT_sorted.bam
samtools index Input_rep1_0hr_sorted.bam
samtools index K27ac_rep1_0hr_sorted.bam
samtools index K27me3_rep1_0hr_sorted.bam

**ANS**: What is the advantage of indexing a BAM file?
Indexing a BAM file allows for quick searching and retrieval of sequence information from the BAM file. This is important for time efficiency, especially in cases of random sequence retrieval from a large file (i.e. having an index avoids having to scan through the entire file).


**Note: once you’ve finished these questions, only keep the sorted .bam and .bai files to save disk space!**

#Calling peaks

###Question 5

Use MACS2 to call H3K27ac and H3K27me3 peaks separately against the input control. Use MACS2 to call ATAC peaks. Output your peak calls in a BED file named *<Your_NetID>_<condition>_t0_rep1_peaks.bed* (conditions should be one of “H3K27ac”, “H3K27me3”, “ATAC”). Use the “--broad” flag for H3K27me3 peak calling. You can read about the broadPeaks format at the UCSC Genome Browser FAQ: https://genome.ucsc.edu/FAQ/FAQformat.html#format13

**ANS**: Write down the commands you used and a description of what each command does.

Commands: 

module load MACS2

macs2 callpeak -t K27ac_rep1_0hr_sorted.bam -c Input_rep1_0hr_sorted.bam -f BAM -g hs -n rua2_H3K27ac_t0_rep1_peaks.bed
macs2 callpeak -t K27me3_rep1_0hr_sorted.bam -c Input_rep1_0hr_sorted.bam --broad -f BAM -g hs -n rua2_H3K27me3_t0_rep1_peaks.bed
macs2 callpeak -t ATAC_seq_rep1_0hr_noMT_sorted.bam -f BAMPE -g hs -n rua2_ATAC_t0_rep1_peaks.bed

Description: 
The callpeak command from the macs2 module will take in the treatment file (e.g. K27ac_rep1_0hr_sorted.bam) along with the control input (for ChIP-seq; Input_rep1_0hr_sorted.bam) and call peaks (--broad type for H3K27me3), with the human reference genome size (-g hs). The output is specified as -f BAM for single-end reads and -f BAMPE for paired-end.

**ANS**: Report the total number of peaks called for each file.
wc -l rua2_ATAC_t0_rep1_peaks.bed_peaks.narrowPeak = 89272 

wc -l rua2_H3K27me3_t0_rep1_peaks.bed_peaks.broadPeak =  8004 

wc -l rua2_H3K27ac_t0_rep1_peaks.bed_peaks.narrowPeak = 60919 


**ANS**: What does the --broad flag do? Why would you want to use it for H3K27me3 data?

The --broad flag will generate composite broad regions by combining nearby highly enriched regions into a broad region with a loose cutoff
The broad flag is more appropriate when calling peaks on histone modifications that span entire gene bodies or larger loci. NarrowPeaks is more appropriate for proteins that bind specific motifs, or narrow parts of the genome such as transcription factors 

References: 
https://github.com/crazyhottommy/ChIP-seq-analysis/blob/master/part1_peak_calling.md
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html


**ANS**: Why can we use the same input control for H3K27ac and H3K27me3?

The input sample represents a background control, which goes through the same ChIP processing as the experimental samples (i.e. up to the point of library prep), but without protein enrichment by adding an antibody (e.g. antibody against H3K27ac or H3K27me3). Therefore this same input control can be used for both H3K27ac and H3K27me3 as it represents the nonenriched noise/background, that can be applied to normalize both the H3K27ac and H3K27me3 samples.

**ANS**: Describe the contents of a bed file. 
Bedfiles need three essential columns describing sequence data: chromosome number, start position and stop position. Other columns may be present with additional data (e.g. name/identifier, score, strand, etc.)

**ANS:** What are the outputs from MACS2? How does a narrowPeak file differ from a BED 
file?

The output of MACS2 are a few files: .bed, .xls, .narrowPeak, model.r

The R script is a file that generates a PDF image about the model, summarizing the data.
The .xlsx, .bed, .narrowPeaks file contain similar info on the location of peaks in the genome. summits.bed is specific to the summit/max point of the peak. .narrowPeaks as a few extra columns corresponding to:
- column 5 (score) that specifies how the peak will appear in the browser 
- column 7 (signalValue) that measures the overall, average fold-change enrichment of the 
- column 8 (pValue) a measure of statistical significance 
- column 10 (peak) which is the position of the summit, relative to peak start

Reference: 
https://miawang113.wordpress.com/2019/01/30/output-files-from-macs-2/

#Annotating peak calls

###Question 6

Use BEDTools and the pre-processed peak call files in the ‘Peaks’ directory to intersect overlapping peak calls from each replicate for each condition, using a 1 bp minimum overlap. Name each output file as *<Your_NetID>_<Condition>_intersected_peaks* (e.g., ‘netID_K27me3_0hr_intersected_peaks, etc.).

**ANS**: Write down the commands you used and a description of what each command does.

commands:

bedtools intersect -f 1E-9 -a ${wd}/ATACseq_3TP_0hr_rep1_peaks.narrowPeak -b ${wd}/ATACseq_3TP_0hr_rep2_peaks.narrowPeak > sl2867_ATACseq_3TP_0hr_intersected_peaks
bedtools intersect -f 1E-9 -a ${wd}/ATACseq_3TP_72hr_rep1_peaks.narrowPeak -b ${wd}/ATACseq_3TP_72hr_rep2_peaks.narrowPeak > sl2867_ATACseq_3TP_72hr_intersected_peaks
bedtools intersect -f 1E-9 -a ${wd}/rep1_0hr_K27ac_peaks.narrowPeak -b ${wd}/rep2_0hr_K27ac_peaks.narrowPeak > sl2867_K27ac_0hr_intersected_peaks
bedtools intersect -f 1E-9 -a ${wd}/rep1_72hr_K27ac_peaks.narrowPeak -b ${wd}/rep2_72hr_K27ac_peaks.narrowPeak > sl2867_K27ac_72hr_intersected_peaks
bedtools intersect -f 1E-9 -a ${wd}/rep1_0hr_K27me3_peaks.broadPeak -b ${wd}/rep2_0hr_K27me3_peaks.broadPeak > sl2867_K27me3_0hr_intersected_peaks
bedtools intersect -f 1E-9 -a ${wd}/rep1_72hr_K27me3_peaks.broadPeak -b ${wd}/rep2_72hr_K27me3_peaks.broadPeak > sl2867_K27me3_72hr_intersected_peaks

description: 
for each of these commands (1 command per condition): the bedtool intersect is used to find overlapping reads between the two input replicate files (-a and -b), with a 1 bp minimum overlap (-f 1E-9) and the output is redirected (>) to the given file names 

**ANS**: How many peaks are preserved in the intersected file compared to each of the	individual replicate files?

Command for (1-6): wc -l *intersected_peaks

1)  21812 rua2_ATACseq_3TP_0hr_intersected_peaks
wc -l ATACseq_3TP_0hr*
  30968 ATACseq_3TP_0hr_rep1_peaks.narrowPeak
  39882 ATACseq_3TP_0hr_rep2_peaks.narrowPeak
  
2)  34060 rua2_ATACseq_3TP_72hr_intersected_peaks
wc -l ATACseq_3TP_72hr_rep*
   42881 ATACseq_3TP_72hr_rep1_peaks.narrowPeak
   53941 ATACseq_3TP_72hr_rep2_peaks.narrowPeak

3)  50791 rua2_K27ac_0hr_intersected_peaks
wc -l *0hr_K27ac*
   64737 rep1_0hr_K27ac_peaks.narrowPeak
   71407 rep2_0hr_K27ac_peaks.narrowPeak

4)  41050 rua2_K27ac_72hr_intersected_peaks
wc -l *72hr_K27ac*
   48764 rep1_72hr_K27ac_peaks.narrowPeak
   53218 rep2_72hr_K27ac_peaks.narrowPeak

5)  6698 rua2_K27me3_0hr_intersected_peaks
wc -l *0hr_K27me3*
   9382 rep1_0hr_K27me3_peaks.broadPeak
   8838 rep2_0hr_K27me3_peaks.broadPeak
   
6)  6230 rua2_K27me3_72hr_intersected_peaks
wc -l *72hr_K27me3*
   6801 rep1_72hr_K27me3_peaks.broadPeak
  10509 rep2_72hr_K27me3_peaks.broadPeak


Now, use BEDTools to merge H3K27ac or H3K27me3 peaks in each intersected file generated above that are 100 bp apart or less. Name each output file as *<Your_NetID>_<Condition>_merged_peaks* (e.g., ‘netID_K27me3_0hr_merged_peaks, etc.).

Commands to merge:

#-k1,1 means sort by first field and restrict sorting to that field
#-k2,2n means sort by second field and restrict sorting to that field, numerical sort

sort input bed first:
sort -k1,1 -k2,2n sl2867_ATACseq_3TP_0hr_intersected_peaks > sl2867_ATACseq_3TP_0hr_intersected_peaks.sorted.bed
sort -k1,1 -k2,2n sl2867_ATACseq_3TP_72hr_intersected_peaks > sl2867_ATACseq_3TP_72hr_intersected_peaks.sorted.bed
sort -k1,1 -k2,2n sl2867_K27ac_0hr_intersected_peaks > sl2867_K27ac_0hr_intersected_peaks.sorted.bed
sort -k1,1 -k2,2n sl2867_K27ac_72hr_intersected_peaks > sl2867_K27ac_72hr_intersected_peaks.sorted.bed
sort -k1,1 -k2,2n sl2867_K27me3_0hr_intersected_peaks > sl2867_K27me3_0hr_intersected_peaks.sorted.bed
sort -k1,1 -k2,2n sl2867_K27me3_72hr_intersected_peaks > sl2867_K27me3_72hr_intersected_peaks.sorted.bed

Merge:
bedtools merge -d 100 -i sl2867_ATACseq_3TP_0hr_intersected_peaks.sorted.bed > sl2867_ATACseq_3TP_0hr_merged_peaks.bed
bedtools merge -d 100 -i sl2867_ATACseq_3TP_72hr_intersected_peaks.sorted.bed > sl2867_ATACseq_3TP_72hr_merged_peaks.bed
bedtools merge -d 100 -i sl2867_K27ac_0hr_intersected_peaks.sorted.bed > sl2867_K27ac_0hr_merged_peaks.bed
bedtools merge -d 100 -i sl2867_K27ac_72hr_intersected_peaks.sorted.bed > sl2867_K27ac_72hr_merged_peaks.bed
bedtools merge -d 100 -i sl2867_K27me3_0hr_intersected_peaks.sorted.bed > sl2867_K27me3_0hr_merged_peaks.bed
bedtools merge -d 100 -i sl2867_K27me3_72hr_intersected_peaks.sorted.bed > sl2867_K27me3_72hr_merged_peaks.bed

**ANS**: How many peaks are in each merged file? How does this compare to the number of peaks in each intersected file?

wc -l *merged_peaks.bed
  21808 sl2867_ATACseq_3TP_0hr_merged_peaks.bed
  34053 sl2867_ATACseq_3TP_72hr_merged_peaks.bed
  44702 sl2867_K27ac_0hr_merged_peaks.bed
  36079 sl2867_K27ac_72hr_merged_peaks.bed
   6698 sl2867_K27me3_0hr_merged_peaks.bed
   6230 sl2867_K27me3_72hr_merged_peaks.bed
   
This is very similar to most of the intersected files, with K27ac (0hr and 72hr) showing the most variability: 

wc -l *intersected_peaks
   21812 rua2_ATACseq_3TP_0hr_intersected_peaks
   34060 rua2_ATACseq_3TP_72hr_intersected_peaks
   50791 rua2_K27ac_0hr_intersected_peaks
   41050 rua2_K27ac_72hr_intersected_peaks
    6698 rua2_K27me3_0hr_intersected_peaks
    6230 rua2_K27me3_72hr_intersected_peaks
    
**ANS**: Why might you want to merge histone modification peaks that are close to each other?

Histone modifications may stretch along large regions of the genome (e.g. across entire gene bodies), representing an area of gene activation/repression. Peaks close to each other, although appearing separate and distinct in reads/sequencing output, may actually represent a biologically continuous stretch of regulatory marks.

**Note:** we are only asking you to submit the H3K27me3 intersected and merged peak files in your PS3 tarball.

###Question 7	 

HOMER is a suite of tools for analyzing multiple types of functional genomics data.
If you get stuck, see: http://homer.ucsd.edu/homer/ngs/index.html (tutorials 1-4 and 6). The HOMER utilities are available in the *Tools/HOMER/bin* folder in our course directory. Adding Adding:
**PATH=$PATH:/gpfs/gibbs/project/gene760/shared/folders/Tools/HOMER/bin** to your **.bashrc** will allow you to call these tools directly from the command line rather than executing them via file path. 

Use HOMER to annotate each intersected file obtained in Question 6. HOMER uses RefSeq annotations to classify each peak into one of the following categories:
•	‘exon’ (3’ UTR, 5' UTR, coding, non-coding)
•	‘intron’
•	‘intergenic’
•	‘promoter-TSS’ (defined from -1kb to +100bp)
•	‘TTS’ (defined from -100 bp to +1kb)

For each set of annotated peaks, extract the lines corresponding to promoters and intronic/intergenic regions and save them to separate files. Convert each of these files to BED format. Hint: Use the **awk** command line tool, column 8 from the Homer output contains the annotation classification. 

**ANS**: Write the commands you used and a description of what each command does.

[not required for full marks] To add HOMER to .bashrc:

echo 'export PATH=$PATH:/gpfs/gibbs/project/gene760/shared/folders/Tools/HOMER/bin' >> ~/.bashrc
source ~/.bashrc

[required]

1. Annotate:

# Use for loop to loop through intersected peak files for each timepoint (6 files to annotate)
#use hg19, this was specified in the manucript that the sequences were aligned to hg19
Command:
for file in *intersected_peaks.sorted.bed; do
    annotatePeaks.pl $file hg19 > ${file%intersected_peaks.sorted.bed}annotated.bed 2> ${file%intersected_peaks.sorted.bed}annotation_stats.log
done

Description: 
The command will loop over each of the intersected, sorted bed files in the directory and apply the annotatePeaks.pl script from Homer, with the output bed file containing the annotation of each read in column 8. The 2> exports the standard error output to the .log file


2. Extract promoter and intronic/intergenic regions and save them to output bed file: 

for file in *annotated.bed; do
    awk '{if($8=="promoter-TSS") print $0}' $file > ${file%annotated.bed}promoters
    awk '{if($8=="intron" || $8=="intergenic") print $0}' $file > ${file%annotated.bed}intronic_intergenic
done

Description: 
The command will loop over each of the annotated bed files in the directory and apply the awk command twice, once to filter for lines with matching phrases to "promoter-TSS" in column 8 and once for "intron" or "intergnic". If it finds a match, it will print the entire line ($0) too an output file with an appropriate name.

3. Convert to bed file: 

for file in *promoters; do
    cut -f2-4 $file > ${file%}.bed
done

for file in *intronic_intergenic; do
    cut -f2-4 $file > ${file%}.bed
done

Description:
These for loops will go through the promoter files and the intronic/intergenic converting them to bed files by extracting the three required columns of a bed file (chromosome number, start position and end positions), corresponding to columns 2-4 from the annotated files.

**ANS**: Report the number of peaks in each BED file.

Promoters:
wc -l *promoters.bed
    9079 rua2_ATACseq_3TP_0hr_promoters.bed
   11206 rua2_ATACseq_3TP_72hr_promoters.bed
   11890 rua2_K27ac_0hr_promoters.bed
   11266 rua2_K27ac_72hr_promoters.bed
     760 rua2_K27me3_0hr_promoters.bed
     730 rua2_K27me3_72hr_promoters.bed

Intronic/Intergenic:
wc -l *intronic_intergenic.bed
    5156 rua2_ATACseq_3TP_0hr_intronic_intergenic.bed
   10144 rua2_ATACseq_3TP_72hr_intronic_intergenic.bed
   20588 rua2_K27ac_0hr_intronic_intergenic.bed
   16089 rua2_K27ac_72hr_intronic_intergenic.bed
    2255 rua2_K27me3_0hr_intronic_intergenic.bed
    2080 rua2_K27me3_72hr_intronic_intergenic.bed


**ANS**: Consider the biological functions of H3K27ac and H3K27me3. How would you interpret the presence of H3K27ac enriched regions at promoters? H3K27me3 peaks at promoters?

H3K27ac is an activating histone mark, associated with active enhancers and therefore its presence at promoters might indicate high transcription and expression of that gene and that it may be marking a cis regulatory element. H3K27me3 is a repressive histone mark and its peaks at promoters might indicate low transcription/expression of that gene.

**ANS**: What about H3K27ac at intergenic/intronic regions? H3K27me3 at intergenic/intronic regions?

At intronic/intergenic regions, these marks may be associated with enhancing and repressing trans regulatory elements (for H3K27ac and H3K27me3, respectively). This is because they are further from the transcription start sites and may contact distant promoters of one or more genes to serve their regulatory roles.

**ANS**: Using BEDTools, intersect your H3K27ac and H3K27me3 peak calls at both time points. How many sites overlap?

Used the intersect.bed files (representing both replicates):

At 0hr:
bedtools intersect -a sl2867_K27ac_0hr_intersected_peaks -b sl2867_K27me3_0hr_intersected_peaks > 0hr_K27ac_K27me3_overlap.bed

wc -l 0hr_K27ac_K27me3_overlap.bed
372 0hr_K27ac_K27me3_overlap.bed

At 72hr:
bedtools intersect -a sl2867_K27ac_72hr_intersected_peaks -b sl2867_K27me3_72hr_intersected_peaks > 72hr_K27ac_K27me3_overlap.bed

wc -l 72hr_K27ac_K27me3_overlap.bed
323 72hr_K27ac_K27me3_overlap.bed

372 sites overlap at 0hrs and 323 overlap at 72 hrs

**ANS**: Using BEDTools, intersect your H3K27ac and ATAC-seq calls. What fraction of ATAC-seq calls land in H3K27ac regions at each time point?

At 0hr:
bedtools intersect -a sl2867_ATACseq_3TP_0hr_intersected_peaks -b sl2867_K27ac_0hr_intersected_peaks > 0hr_ATAC_K27ac_overlap.bed

At 72hr:
bedtools intersect -a sl2867_ATACseq_3TP_72hr_intersected_peaks -b sl2867_K27ac_72hr_intersected_peaks > 72hr_ATAC_K27ac_overlap.bed

wc -l 0hr_ATAC_K27ac_overlap.bed
10084 0hr_ATAC_K27ac_overlap.bed

wc -l 72hr_ATAC_K27ac_overlap.bed 
13729 72hr_ATAC_K27ac_overlap.bed

At 0hrs around 46.2% (10084/21812) of ATAC sites land in H3K27ac regions 
At 72 hrs around 40.3% (13729/34060) of ATAC sites land in H3K27ac regions 


**ANS**: How would you interpret the presence of both H3K27ac and H3K27me3 at promoters or intergenic/intronic regions, particularly at time point 0?

The presence of both H3K27ac and H3K27me3 at promoters or intergenic/intronic regions likely represents a state of bivalency, in which both activating and repressive histone marks are located together. These are "poised" or "primed" for either activation or repression. At a later time point, the repressive H3K27me3 may be lost, leading to gene expression


###Question 8

Using BEDTools and the BED files you generated in Question 7 above, identify promoter and intergenic and intronic H3K27ac and ATAC-seq peaks shared between the two time points, and peaks that are unique to each. Save your intergenic/intronic H3K27ac and ATAC-seq peaks in BED format as follows:
		*<Your_NetID>_H3K27ac_shared_t0_t72*
		*<Your_NetID>_ATACseq_shared_t0_t72*

**ANS**:  Write the commands you used and a description of what each command does.

COMMANDS: 

SHARED peaks:
0hr and 72hr shared H3K27ac at promoters:
bedtools intersect -a sl2867_K27ac_0hr_promoters.bed -b sl2867_K27ac_72hr_promoters.bed > sl2867_H3K27ac_shared_promoter_t0_t72.bed

0hr and 72hr shared H3K27ac at intergenic/intronic
bedtools intersect -a sl2867_K27ac_0hr_intronic_intergenic.bed -b sl2867_K27ac_72hr_intronic_intergenic.bed > sl2867_H3K27ac_shared_t0_t72.bed

0hr and 72hr shared ATAC at promoters 
bedtools intersect -a sl2867_ATACseq_3TP_0hr_promoters.bed -b sl2867_ATACseq_3TP_72hr_promoters.bed > sl2867_ATAC_shared_promoter_t0_t72.bed

0hr and 72hr shared ATAC at intergenic/intronic
bedtools intersect -a sl2867_ATACseq_3TP_0hr_intronic_intergenic.bed -b sl2867_ATACseq_3TP_72hr_intronic_intergenic.bed > sl2867_ATACseq_shared_t0_t72.bed


Unique Peaks:
0hr and 72hr unique H3K27ac at promoters:
bedtools intersect -v -a sl2867_K27ac_0hr_promoters.bed -b sl2867_K27ac_72hr_promoters.bed > sl2867_H3K27ac_unique_promoter_t0_t72.bed
bedtools intersect -v -a sl2867_K27ac_72hr_promoters.bed -b sl2867_K27ac_0hr_promoters.bed > sl2867_H3K27ac_unique_promoter_t72_t0.bed

0hr and 72hr unique H3K27ac at intergenic/intronic
bedtools intersect -v -a sl2867_K27ac_0hr_intronic_intergenic.bed -b sl2867_K27ac_72hr_intronic_intergenic.bed > sl2867_H3K27ac_unique_t0_t72.bed
bedtools intersect -v -a sl2867_K27ac_72hr_intronic_intergenic.bed -b  sl2867_K27ac_0hr_intronic_intergenic.bed > sl2867_H3K27ac_unique_t72_t0.bed

0hr and 72hr unique ATAC at promoters 
bedtools intersect -v -a sl2867_ATACseq_3TP_0hr_promoters.bed -b sl2867_ATACseq_3TP_72hr_promoters.bed > sl2867_ATAC_unique_promoter_t0_t72.bed
bedtools intersect -v -a sl2867_ATACseq_3TP_72hr_promoters.bed -b sl2867_ATACseq_3TP_0hr_promoters.bed > sl2867_ATAC_unique_promoter_t72_t0.bed


0hr and 72hr unique ATAC at intergenic/intronic
bedtools intersect -v -a sl2867_ATACseq_3TP_0hr_intronic_intergenic.bed -b sl2867_ATACseq_3TP_72hr_intronic_intergenic.bed > sl2867_ATACseq_unique_t0_t72.bed
bedtools intersect -v -a sl2867_ATACseq_3TP_72hr_intronic_intergenic.bed -b sl2867_ATACseq_3TP_0hr_intronic_intergenic.bed > sl2867_ATACseq_unique_t72_t0.bed

DESCRIPTION:
Bedtools intersect is used to find shared peaks between each time point, as well as unique (-v) for each time point 0 and 72.

**ANS**:  How many peaks are shared for each feature? How many are unique?

Number of shared peaks: 
wc -l *shared*
  8316 rua2_ATAC_shared_promoter_t0_t72.bed
  2823 rua2_ATACseq_shared_t0_t72
  9466 rua2_H3K27ac_shared_promoter_t0_t72.bed
  9212 rua2_H3K27ac_shared_t0_t72
   
Number of unique peaks: 
wc -l *unique*
   829 rua2_ATAC_unique_promoter_t0_t72.bed
  2986 rua2_ATAC_unique_promoter_t72_t0.bed
  2337 rua2_ATACseq_unique_t0_t72.bed
  7334 rua2_ATACseq_unique_t72_t0.bed
  3152 rua2_H3K27ac_unique_promoter_t0_t72.bed
  2592 rua2_H3K27ac_unique_promoter_t72_t0.bed
 12248 rua2_H3K27ac_unique_t0_t72.bed
  7853 rua2_H3K27ac_unique_t72_t0.bed

###Question 9	 

The sequence read archive (SRA) is an NCBI repository for high throughput sequencing data and is a place that many researchers publicly provide access to their raw data. 
The sequencing files that you used in this problem set were retrieved for you from this repository. In order to retrieve this data, you can use the SRA-Toolkit on McCleary. 
You will be retrieving only one of these files to familiarize yourself with this process and to save disk space. 
Here is a link to the full dataset from Inoue et al: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115046 

**ANS**: Write down the commands you used to load and configure the SRA-Toolkit.

module load SRA-Toolkit  # load

vdb-config #configure; do when running for the first time 

Use the ‘prefetch’ utility to retrieve the file with SRR identifier: SRR7230423. 
**ANS**: Write down the commands you used and a description of what each command does.

Command:
prefetch SRR7230423 -O /vast/palmer/scratch/gene760/gene760_rua2/pset3/SRA

Description: 
prefetch is used to retrieve the datasets with the given SRR ID to the specified output (-O) directory.

Next, use the ‘fasterq-dump’ utility with the ‘--split-files’ flag to output your fastq files.
**ANS**: Write down the commands you used and a description of what each command does.

Command:

fasterq-dump --split-files --gzip SRR7230423

Description:
fasterq-dump command specifically obtains the fastq files (split into each of the two reads and zipped for space effeciency) from the given SRR dataset.

**ANS**: Why is it necessary to use the --split-files flag when working with paired-end sequencing data?
To obtain two distinct files, R1 and R2, for each of the paired-end reads that will be used in the downstream analysis steps

Run FastQC on the output files and compare the resulting QC report to ATACseq_rep1_0hr_1.fastq.gz and ATACseq_rep2_0hr_2.fastq.gz.
**ANS**: Are these files the same as those retrieved from SRR7230423? How many reads does each file have?

module load FastQC

fastqc -o . *.fastq.gz

[any reasonable interpretation]
The files retrieved via the SRA toolkit seem to have less reads that are longer reads:
171994 total sequences with a sequence length of 150 for SRR7230423 
42326401 total sequences with a length of 101 for ATACseq_rep1_0hr
However the %GC is similar 46%, and 47% for SRR7230423 and ATAC_seq_rep1, respectively. 