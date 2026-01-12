**GENE 760 -- Problem Set 1**

The purpose of this problem set is to familiarize you with the analysis of mRNA-seq data. **Note: DataCamp also offers a tutorial on RNA-seq analysis using DESeq2 in R which will be very helpful for this Problem Set:** <https://app.datacamp.com/learn/courses/rna-seq-with-bioconductor-in-r>

By the end of this problem set, you will have learned how to use:

-   STAR to map spliced RNA-seq reads to a genome

-   featureCounts to quantify gene expression as read counts

-   DESeq2 to perform differential expression analysis

-   Enrichr to perform gene ontology analysis on a list of genes

Please submit a **gzipped tarball** called `<Your_NetID>_PS1.tar.gz` containing:

-   `<Your_NetID>_PS1_answers.md`: Answers to the questions below

-   `<Your_NetID>_DEseq2.R`: Commented R script used to perform differential expression analysis

-   `<Your_NetID>_2wk_heatmap.png`

-   `<Your_NetID>_6wk_heatmap.png`

to `gpfs/gibbs/project/gene760/shared/Dropbox/PS1` by **Friday, February 7 at 11:59 pm.**

**Remember to run commands and process all data in scratch60.**

You will be working with data from Gjoneska et al., Nature **518:**365-369 (2015). This study generated chromatin and transcriptome maps from hippocampus tissue of a mouse model for Alzheimer's disease. In these mice, induced brain-specific overexpression of the protein p25 results in aberrant activation of Cyclin-dependent kinase 5 (Cdk5), leading to accumulation of amyloid-β protein and subsequent neurodegeneration. Accumulation of amyloid-β occurs soon after induction of p25, before overt signs of neurodegeneration are apparent. Chromatin and transcriptome maps were generated from hippocampal tissues in both control (CK) and inducible p25 transgenic (CK-p25) animals 2 weeks and 6 weeks post-induction, and differential analysis was performed. We will focus on the transcriptome maps for this problem set.

Datasets are located in `/gpfs/gibbs/project/gene760/shared/folders/Data/PS1`

In this directory, you will find 12 gzipped FASTQ sequence files:

-   Control samples (2 week time point):

    -   `CK_2wk_1_R1.fastq.gz, CK_2wk_1_R2.fastq.gz`

    -   `CK_2wk_2_R1.fastq.gz, CK_2wk_2_R2.fastq.gz`

    -   `CK_2wk_3_R1.fastq.gz, CK_2wk_3_R2.fastq.gz`

-   Treatment samples (2 week time point):

    -   `CKp25_2wk_1_R1.fastq.gz, CKp25_2wk_1_R2.fastq.gz`

    -   `CKp25_2wk_2_R1.fastq.gz, CKp25_2wk_2_R2.fastq.gz`

    -   `CKp25_2wk_3_R1.fastq.gz, CKp25_2wk_3_R2.fastq.gz`

Pre-built genome index files for mm10 are located in:
`/gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/STAR_index/`

A transcript annotation file in GTF format is located at:

`/gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/gencode.vM23.annotation.gtf`

**STAR and featureCounts**

Note: Gjoneska et al. generate strand-specific RNA-seq reads using Illumina's TruSeq Stranded Total RNA prep kit (see Methods)--i.e. only the strand from first strand synthesis is sequenced.

1)  Use STAR to align RNA-seq reads for the treatment and control samples from the 2-week time point. Report alignments with no more than 2 mismatches. Generate the output as a sorted BAM file.

**ANS:** Write the commands you used and a description of what each command does, including the meaning of any arguments.

commands:
module load STAR
STAR --genomeDir /gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/STAR_index/ --readFilesCommand zcat --readFilesPrefix /gpfs/gibbs/project/gene760/shared/folders/Data/PS1/ --readFilesIn CK_2wk_1_R1.fastq.gz CK_2wk_1_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --runThreadN 20 --outFileNamePrefix CK_2wk_1
STAR --genomeDir /gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/STAR_index/ --readFilesCommand zcat --readFilesPrefix /gpfs/gibbs/project/gene760/shared/folders/Data/PS1/ --readFilesIn CK_2wk_2_R1.fastq.gz CK_2wk_2_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --runThreadN 20 --outFileNamePrefix CK_2wk_2
STAR --genomeDir /gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/STAR_index/ --readFilesCommand zcat --readFilesPrefix /gpfs/gibbs/project/gene760/shared/folders/Data/PS1/ --readFilesIn CK_2wk_3_R1.fastq.gz CK_2wk_3_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --runThreadN 20 --outFileNamePrefix CK_2wk_3
STAR --genomeDir /gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/STAR_index/ --readFilesCommand zcat --readFilesPrefix /gpfs/gibbs/project/gene760/shared/folders/Data/PS1/ --readFilesIn CKp25_2wk_1_R1.fastq.gz CKp25_2wk_1_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --runThreadN 20 --outFileNamePrefix CKp25_2wk_1
STAR --genomeDir /gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/STAR_index/ --readFilesCommand zcat --readFilesPrefix /gpfs/gibbs/project/gene760/shared/folders/Data/PS1/ --readFilesIn CKp25_2wk_2_R1.fastq.gz CKp25_2wk_2_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --runThreadN 20 --outFileNamePrefix CKp25_2wk_2
STAR --genomeDir /gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/STAR_index/ --readFilesCommand zcat --readFilesPrefix /gpfs/gibbs/project/gene760/shared/folders/Data/PS1/ --readFilesIn CKp25_2wk_3_R1.fastq.gz CKp25_2wk_3_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --runThreadN 20 --outFileNamePrefix CKp25_2wk_3

description: The STAR command will align the paired input read (--readFilesIn) fastq files to the reference genome directory (STAR index files, with the path provided following the --genomeDir option), generating a bam file with the aligned reads from the input file with their genomic locations. This genome alignment command is designed specifically for raw data files from RNA seq experiments. The --outSAMtype BAM SortedByCoordinate portion specifies that the output should be a sorted BAM file. The --outFilterMismatchNmax 2 option specifies that alignments reported should have no more than 2 mismatches. Finally, --runThreadN 20 means that 20 CPUs/cores should be used to run the command, and -outFileNamePrefix is given to specify the output file names.


**ANS:** What are the contents of the genome index directory, and what does STAR use them for?
Contents:
binary genome sequence, suffix arrays, text chromosome names/lengths, splice junctions coordinates, transcripts/genes information

STAR uses them to map the sequences in the fastq files to the reference genome and provide information about the location of the sequence on the genome.

2)  Use the featureCounts tool in the Subread module to quantify gene expression for each sample. Subread's featureCounts tool uses alignments generated by STAR to calculate transcript abundances.

Use the transcript annotation file provided.

Info on featureCounts is available here: <http://subread.sourceforge.net/featureCounts.html>

We have also provided the Subread documentation in the PS1 directory on Canvas.

Note: For paired-end data, the counting will be faster if the alignment files provided are sorted by name so that aligned pairs are next to each other in sequence.

**ANS:** Write the commands you used and a description of what each command does.

commands:
module load Subread
featureCounts -p --countReadPairs -t exon -g gene_name -a /gpfs/gibbs/project/gene760/shared/folders/Annotation/mm10/gencode.vM23.annotation.gtf -o counts.txt CK_2wk_1Aligned.sortedByCoord.out.bam CK_2wk_2Aligned.sortedByCoord.out.bam CK_2wk_3Aligned.sortedByCoord.out.bam CKp25_2wk_1Aligned.sortedByCoord.out.bam CKp25_2wk_2Aligned.sortedByCoord.out.bam CKp25_2wk_3Aligned.sortedByCoord.out.bam

description: This featureCounts command quantifies gene expression by processing multiple paired-end dataset bam input files attained from the STAR alignment (i.e. sample1.bam...sample6.bam). The -p --countReadPairs options specify that the input files are from paired reads, -t specifies that counting should be done at the exon level, -g indicates that features should be grouped by gene names, -a if the GTF annotation file path, and -o is the output file name where the counts will be stored.

**ANS:** In your own words, describe the logic featureCounts uses to assign reads to genes.

featureCounts compares the coordinates of the reads from the .bam file to the coordinates of the features (exons) in the GTF annotation file.
It will check to see if reads (with their given locations in the .bam files) overlap these features. If at least one read base from the read overlaps the feature/exon, then a count is generated and recorded in the .txt output file, with its corresponding unique gene ID (information in the .gtf file).

**ANS:** What are the advantages of stranded RNA-seq reads for quantifying gene expression?

Stranded RNA-seq helps to accurately quantify gene expression levels in cases when different genes may have overlapping genomic loci but are transcribed from different strands (i.e. forward vs. reverse). Stranded RNA-seq reads preserve information of the strand from which the read was obtained. Quantifying sense and antisense gene expression can give information on potential regulatory interactions (Zhao et al., 2015). Other advantages include greater reproducibility of downstream steps of the analysis (e.g .differential expression analysiss), more accurate mapping, higher splice isoform resolution and allowing for de-novo transcriptome assemblies (Signal et al., 2015)

References: 
Zhao, S., Zhang, Y., Gordon, W. et al. Comparison of stranded and non-stranded RNA-seq transcriptome profiling and investigation of gene overlap. BMC Genomics 16, 675 (2015). https://doi.org/10.1186/s12864-015-1876-7
Signal B, Kahlke T. how_are_we_stranded_here: quick determination of RNA-Seq strandedness. BMC Bioinformatics. 2022 Jan 22;23(1):49. doi: 10.1186/s12859-022-04572-7. PMID: 35065593; PMCID: PMC8783475.


For the following questions, use reads counts from the 2-week time point (generated above) as well as counts from the 6-week time point. Read counts for control and treatment samples from the 6-week time point are located in

`/gpfs/gibbs/project/gene760/shared/folders/Data/PS1/6wk_gencodeM23_featureCounts.txt`

**DEseq2 (R package)**

DEseq2 utilizes read counts from featureCounts to build statistical models of gene expression and identify differentially expressed genes. An extensive description of how DEseq2 works can be found here:

<http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>

3)  Perform differential expression analysis for the 2-week and 6-week time points independently using the filtered read counts generated in question 3.

**ANS**: Save the R script containing all the commands you used to perform differential expression analysis as `<Your_NetID>_DEseq2.R`

**ANS:** How many differentially expressed genes with an adjusted p-value less than 0.01 are identified at 2 weeks? And at 6 weeks?

# If low counts filtered:
2wks: 1678
6wks: 3381

# If not:
2wks: 1648
6wks: 3220

**ANS:** Use R to generate heatmaps showing the top differentially expressed genes for each time point. Save the plots as `<Your_NetID>_2wk_heatmap.png` and `<Your_NetID>_6wk_heatmap.png`

4)  Using the DEseq2 output you generated above, determine the number and identity of genes specifically upregulated in each time point (upregulated at one time point but not the other).

**ANS:** How many genes are significantly upregulated at a log~2~FC <span class="underline">\></span> 1? Save these to a text file for each time point (for use in question 5).
# If low counts filtered:
2wks: 831
6wks: 1137

# If not:
2wks: 826
6wks: 1106

**ANS:** Write the commands you used to determine this and a description of what the command does.
Commands:
sig_up_2wk <- res_2wk[res_2wk$log2FoldChange >= 1 & res_2wk$padj < 0.01, ]
nrow(sig_up_2wk)
write.table(sig_up_2wk, file = "./deseq2_log2FC1_padj0.01_upreg_2wk.txt", sep = "\t", quote = F, row.names = T)

sig_up_6wk <- res_6wk[res_6wk$log2FoldChange >= 1 & res_6wk$padj < 0.01, ]
nrow(sig_up_6wk)
write.table(sig_up_6wk, file = "./deseq2_log2FC1_padj0.01_upreg_6wk.txt", sep = "\t", quote = F, row.names = T)

Description:
The first line subsets the res_2wk results object to those that have log2FC greater than or equal to 1 and adjusted p-value less than 0.01.
The second line prints the number of rows that meet that criteria.
The last line saves the subsetted results as a tab-delimited file.

**ANS:** How many genes are specifically upregulated at 2 weeks? And at 6 weeks?
# If low counts filtered:
2wks: 198
6wks: 504

# If not:
2wks: 207
6wks: 487

**ANS:** DESeq reports the probability that a gene is differentially expressed based on fold-change and level of expression. Why is incorporating the level of expression into the model important?

The level of expression may provide more insight into relevant biological information (e.g. pathways active, etc) that may not be captured by fold change alone.
For example, GeneA and GeneB may both have a fold change of 2, seemingly having similar differences between treatment and control. However, the level of expression of GeneA may have been 2 in the control and 4 in the treatment, while GeneB may have been 200 in the control and 400 in the treatment (arbitrary/hypothetical values). Based on expression levels, GeneB may have a broader or more relevant biological impact in the treatment.
Additionally, expression levels can allow for comparisons of different transcripts in the same group (i.e. treatment or control) relative to each other.

**Functional annotation and gene ontology analysis with DAVID**

5) Enrichr (https://maayanlab.cloud/Enrichr/) is an intuitive, easy to use tool for functionally annotating lists of genes. These gene lists can be derived from differential expression analyses such as the one carried out by DESeq2 here, or from other functional genomic studies. Use Enrichr to identify GO Biological Process enrichments for genes upregulated at each time point (Question 4). 

**ANS:** What are the top 15 GO Biological Process terms from Enrichr?

2wks:
Defense Response To Virus (GO:0051607)
Mitotic Sister Chromatid Segregation (GO:0000070)
Defense Response To Symbiont (GO:0140546)
DNA-templated DNA Replication (GO:0006261)
DNA Metabolic Process (GO:0006259)
Positive Regulation Of Cytokine Production (GO:0001819)
Positive Regulation Of Cell Cycle Process (GO:0090068)
Sister Chromatid Segregation (GO:0000819)
DNA Replication (GO:0006260)
Positive Regulation Of Tumor Necrosis Factor Superfamily Cytokine Production (GO:1903557)
Response To Cytokine (GO:0034097)
Negative Regulation Of Mitotic Metaphase/Anaphase Transition (GO:0045841)
Positive Regulation Of Tumor Necrosis Factor Production (GO:0032760)
Mitotic Nuclear Division (GO:0140014)
Microtubule Cytoskeleton Organization Involved In Mitosis (GO:1902850)

6wks:
Positive Regulation Of Cytokine Production (GO:0001819)
Positive Regulation Of Tumor Necrosis Factor Superfamily Cytokine Production (GO:1903557)
Positive Regulation Of Tumor Necrosis Factor Production (GO:0032760)
Regulation Of Tumor Necrosis Factor Production (GO:0032680)
Defense Response To Symbiont (GO:0140546)
Defense Response To Virus (GO:0051607)
Positive Regulation Of Intracellular Signal Transduction (GO:1902533)
Inflammatory Response (GO:0006954)
Positive Regulation Of Response To External Stimulus (GO:0032103)
Cellular Response To Lipopolysaccharide (GO:0071222)
Response To Cytokine (GO:0034097)
Positive Regulation Of Interleukin-6 Production (GO:0032755)
Regulation Of Interleukin-6 Production (GO:0032675)
Pattern Recognition Receptor Signaling Pathway (GO:0002221)
Cellular Response To Molecule Of Bacterial Origin (GO:0071219)

**ANS:** What biological processes do your results implicate in Alzheimer's disease? (5-6 sentences) 

Based on my results, it seems like in the initial stages of Alzheimer's disease (AD) induction (week 2), inflammatory responses are occurring.
These include positive regulation of cytokine production and Tumor Necrosis Factor Superfamily cytokine production. There are also processes involved in cell division, such as mitotic sister chromatid segregation (top GO term), DNA Replication, and mitosis-related GO terms.
In the later stages, cell division and cell cycle progression terms become less prominent, and immune processes (both production of immune components and response to them) become more prominent. For example, positive regulation of cytokines and tumor necrosis factors (top 3 GO terms). Continuous immune activation is a known component of AD pathology (Kinnery et al., 2018).

References: 
Kinney JW, Bemiller SM, Murtishaw AS, Leisgang AM, Salazar AM, Lamb BT. Inflammation as a central mechanism in Alzheimer's disease. Alzheimers Dement (N Y). 2018 Sep 6;4:575-590. doi: 10.1016/j.trci.2018.06.014. PMID: 30406177; PMCID: PMC6214864.

**ANS:** How do your results compare to those found in Gjoneska et al.? (1-2 sentences)

Similar to Gjoneska et al. (Fig 1), my results show a transient increase in cell cycle progression (GO seen in week 2 but not in week 6). I also see an upregulation of immune response genes (e.g. cytokine, tumor necrosis and Interleukin-6 production, cytokine and inflammatory response in week 6).
