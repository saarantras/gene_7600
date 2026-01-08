The goal of this problem set is to familiarize students with the manipulation and analysis of human genetics data. In this assignment, you will profile genome variation information and attempt to answer biologically relevant questions. 

You’ve been provided with a VCF file (*Trio-NA12891-NA12892-NA12878.vcf.gz*). This file contains SNPs and indels from three individuals. Individuals NA12891, NA12892, and NA12878 form a trio (dad, mom, and child) that you will be analyzing in depth later. NA12878 has been sequenced and analyzed countless times and has been awarded the notorious title of “The Most well studied Genome in the World”.

The variant call format (VCF) is a generic text file format for storing genome variation data such as SNPs, indels, and structural variants, together with rich annotations. It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome. The format can also contain genotype information for multiple samples for each position, which are stored as additional columns. A VCF file is usually stored in a compressed manner and can be indexed for fast data retrieval of variants from a range of positions on the reference genome.

***

The files for PS4 are available in `/gpfs/gibbs/project/gene760/shared/folders/Data/PS4/`, including:

* Trio VCF: (*Trio-NA12891-NA12892-NA12878.vcf.gz*)

* Pre-built gnomAD gnotate file: (*gnomad.genomes.v3.1.sites.chr1-22_XYM.zip*)

* Genome in a Bottle difficult regions: (*GRCh38_alldifficultregions.bed*)

***

You will use the following software:
* BCFtools (https://samtools.github.io/bcftools/bcftools.html)

* Slivar (https://github.com/brentp/slivar)

* BEDtools (https://bedtools.readthedocs.io/en/latest/)

* VEP (https://useast.ensembl.org/Homo_sapiens/Tools/VEP)

 _Note that Slivar is not available as a module. The software binary is available under */gpfs/gibbs/project/gene760/shared/folders/Data/PS4/slivar*_
    
***

Students should submit:
> *<Your_NetID>_PS4_answers.md: Answers to the questions below including code; your commands should be commented with their (intended) purpose!*

> *<Your_NetID>_mendel.R: R script for Question 6*

> *<Your_NetID>_denovo_variants.R: R script for Question 7*

to `/gpfs/gibbs/project/gene760/shared/Dropbox/PS4` on **Friday, April 5 at 11:59 pm**.

*** 

### **1.** Let’s get some basic statistics of variants in the VCF file. Use BCFtools to answer following questions.
  
#### a. What fraction of variants are SNPs or indels?

```
Examining the output we see:

SN      0       number of samples:      3
SN      0       number of records:      5837075
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 4736243
SN      0       number of MNPs: 0
SN      0       number of indels:       1100832
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0

So:
fraction SNPs = 4736243/5837075=0.8114...
fraction indels = 1100832/5837075=0.1885...
```

```bash
srun --pty -p day -t 4:00:00 bash #not required
module load BCFtools
bcftools stats  /gpfs/gibbs/project/gene760/shared/folders/Data/PS4/Trio-NA12891-NA12892-NA12878.vcf.gz > out # or similar
```

#### b. Count the number of SNPs for each substitution type. Report it in the form of a table with the six possibilities and their corresponding counts. 

(Hint: Considering that DNA is double stranded, A -> C is same as T->G and can be written as A/T -> C/G).

| **Substitution type** | **Count** | **BP change** | **Sum** |
| --------------------- | --------- | ------------- | ------- |
| A>C                   | 193353    | A/T -> C/G    | 386267  |
| T>G                   | 192914    |               |         |
| A>G                   | 758952    | A/T -> G/C    | 1518878 |
| T>C                   | 759926    |               |         |
| A>T                   | 170007    | A/T -> T/A    | 341427  |
| T>A                   | 171420    |               |         |
| C>A                   | 205816    | C/G -> A/T    | 411424  |
| G>T                   | 205608    |               |         |
| C>G                   | 204678    | C/G -> G/C    | 410383  |
| G>C                   | 205705    |               |         |
| C>T                   | 833683    | C/G -> T/A    | 1667864 |
| G>A                   | 834181    |               |         |


```
Same command as above...
```

#### c. Among SNPs, what substitution is the most prevalent? Why is this substitution type the most prevalent?

```
Any reasonable answer. Key points are:
1. C/G -> T/A is the most common
2. This is because transitions are more common than transversions (optionally : chemically similar -> more likely to be tolerated / not recognized by repair enzymes)
3. Mutations away from C/G pairs (instead of towards, as 'A/T -> G/C') are more frequent due to CpG methylation (optionally: consequent deamination)
```

```
No command
```

***

### **2.** Since the human genome is diploid, a variant can be found on either (heterozygous) or both (homozygous alternate) chromosomes. Use BCFtools to extract data for the individual NA12878 and perform the following tasks.

#### a. Count the total number of genotype calls that are homozygous reference, heterozygous, and homozygous alternate for that individual.

```

Output of below commands is:

1539360 0|0
1330707 0|1
1345315 1|0
1621693 1|1

Then, recalling that 1 means alt and 0 means ref:

Homozygous ref = 0|0 = 1539360
Homozygous alt = 1|1 = 1621693
Heterozygous = 0|1 + 1|0 = 1330707 + 1345315 = 2676022

```

```
cd palmer_scratch #optional

bcftools view -O z -s NA12878 /gpfs/gibbs/project/gene760/shared/folders/Data/PS4/Trio-NA12891-NA12892-NA12878.vcf.gz -o NA12878.vcf.gz 

# Explanations 
# view : subcommand to extract information
# -O z : compressed VCF output (Optional)
# -s NA12878 : extract that sample
# Then input & subsq. output filepaths. 

bcftools query -f '[%GT]' NA12878.vcf.gz | sort | uniq -c

# I anticipate students may try to just use grep or similar text-processing tools instead of bcftools. I think this could be fine. For example they could match three characters after `GT\t` and get the same information, you don't necessarially need to use bcftools here. 

```


***

### **3.** The allele frequency indicates how often a variant can be observed if we sequence many individuals, and it can be obtained from resources such as gnomAD or the 1000 Genomes Project (1KGP). Use Slivar and the pre-built annotation file (*gnomAD gnotate file*: gnomad.genomes.v3.1.sites.chr1-22_XYM.zip) to annotate the trio VCF file with allele frequencies from gnomAD; store it as Trio-NA12891-NA12892-NA12878.gnomAD_annotated.vcf.gz. 

### The VCF file is already pre-annotated with allele frequencies (INFO field, ID=AF) from 1KGP. ***The allele frequency from 1KGP is to be used for Q3 in case there are troubles with gnomAD annotation***. Write a script (use BCFtools if needed) to perform the following tasks.

#### a. Write the code to annotate the trio VCF file with gnomAD allele frequency using Slivar:

```bash
srun --pty -p day --mem-per-cpu=64G -t 4:00:00 bash # not required
dat_root="/gpfs/gibbs/project/gene760/shared/folders/Data/PS4/" # not required

${dat_root}slivar expr --vcf ${dat_root}Trio-NA12891-NA12892-NA12878.vcf.gz --gnotate ${dat_root}gnomad.genomes.v3.1.sites.chr1-22_XYM.zip --out-vcf Trio-NA12891-NA12892-NA12878.gnomAD_annotated.vcf.gz  > out 2>&1

# ${dat_root}slivar : calls slivar
# expr : subcommand to filter / annotate
# --vcf Trio... : pick the VCF file
# --gnotate ...gnomad.genomes.v3.1.si... annotate w/ this file

```

#### b. Count the fraction of variants that are common (>5%), low allele frequency (1-5%) and rare (<1%) per 1KGP allele frequency:

```
common=5328038/5837075=0.91279245
low=316153/5837075=0.05416292
rare=192884/5837075=0.03304463
```

```bash
bcftools view -i 'AF>0.05' Trio-NA12891-NA12892-NA12878.gnomAD_annotated.vcf.gz | bcftools view -H | wc -l > common &

# low allele frequency 
bcftools view -i 'AF>=0.01 && AF<=0.05' Trio-NA12891-NA12892-NA12878.gnomAD_annotated.vcf.gz | bcftools view -H | wc -l > low &

#rare 
bcftools view -i 'AF<0.01' Trio-NA12891-NA12892-NA12878.gnomAD_annotated.vcf.gz | bcftools view -H | wc -l > rare &

# -H to avoid miscount due to header !
```

#### c. Count the number of variants that are rare in *at least one* of the five super-populations in 1KGP (AFR, AMR, EAS, EUR, SAS):

```
815366
```

```bash
bcftools view -H -i 'AF_AFR<0.01 | AF_AMR<0.01 | AF_EAS<0.01 | AF_EUR<0.01 | AF_SAS<0.01'  Trio-NA12891-NA12892-NA12878.gnomAD_annotated.vcf.gz | wc -l
```
***

### **4.** The VCF file is annotated with predicted functional consequences (INFO field, ID=CSQ) of alternate alleles (non-reference) by VEP (https://useast.ensembl.org/Homo_sapiens/Tools/VEP). Write a script to answer the following questions (use BCFtools if needed).

#### a. Briefly explain why many variants are annotated with multiple consequences for a single gene:

```
Any reasonable answer. Key points:
1. Multiple isoforms for one gene means multiple potential conseqeunces
2. Multiple alt alleles for tri+ allelic sites
3. Non mutually exclusive consequence codes (e.g. _intron_variant_ & _splice_region_variant_.)
```

#### b. Count the number of variants in the following categories: synonymous, missense, and probable loss of function (pLOF: stop gained, frameshift, splice_donor/acceptor). Does this count make sense biologically? Why?

```
synonymous_variant:16001
missense_variant:15692
pLOF:1409

Calcs (not required):
pLOF=stop_gained+frameshift_variant+splice_acceptor_variant+splice_donor_variant
pLOF=185+383+348+493

For the text answer: key point is negative selection reduces frequency of deleterious mutations.

```

```bash
bcftools query -f '%CSQ\n' Trio-NA12891-NA12892-NA12878.gnomAD_annotated.vcf.gz > csq

declare -a arr=("synonymous_variant" "missense_variant" "stop_gained" "frameshift_variant" "splice_acceptor_variant" "splice_donor_variant")

for i in "${arr[@]}"
do
   echo "$i"
   cat csq | grep $i | wc -l
done
```

***


### **5.** Mapping short reads is difficult in certain genomic regions, particularly due to repeats. The Genome in a Bottle consortium (GIAB) identified such regions of the human genome (GRCh38_alldifficultregions.bed; BED file format). This file can be used with other benchmarking tools to identify true/false positive variant calls. Use BEDtools to calculate the fraction of variants that overlap with the GIAB difficult regions. Briefly, how do you interpret these results?

```
1943107/5837075=0.33289053

That fraction are in difficult-to-map regions, and so may be incorrect. 

```

```bash
dat_root="/gpfs/gibbs/project/gene760/shared/folders/Data/PS4/" # not required

bedtools intersect -a ${dat_root}Trio-NA12891-NA12892-NA12878.vcf.gz -b ${dat_root}GRCh38_alldifficultregions.bed | wc -l
```

***

### **6.** Write a R script to count the number of variants that clearly violate the rules of Mendelian inheritance, given the trio’s relationships to one another. (See diagram in the HTML version on canvas.)

[![`Mendelian Inheritance`](https://i.imgur.com/lMdkj5b.jpeg)](https://i.imgur.com)

(Including R script here for convenience):

```R
# load library
library(vcfR)

# read vcf
vcf <- read.vcfR("/gpfs/gibbs/project/gene760/shared/folders/Data/PS4/Trio-NA12891-NA12892-NA12878.vcf.gz", verbose = FALSE)

# set the counter
vio_count <- 0

# FOR LOOP START
for (i in 1:nrow(vcf@gt)){
  # what is genotype of the child
  child_genotype = vcf@gt[i,"NA12878"]

  # what is genotype of the mom
  mom_genotype = vcf@gt[i,"NA12892"]

  # what is genotype of the dad
  dad_genotype = vcf@gt[i,"NA12891"]

  # Run if-else conditions depending on child_genotype
  # in what conditions, child genotype and parent genotypes violate the rules of Mendelian inheritance
  ## feel free to rewrite these if else commands in your own way
  if (child_genotype == "1|1"){
    if ((mom_genotype == "0|0") | (dad_genotype == "0|0")) {
      #no source of alt (1) allele
      vio_count <- vio_count + 1
    }
  } else if (child_genotype == "0|0"){
    if ((mom_genotype == "1|1") | (dad_genotype == "1|1")) {
▶     # No source of ref (0) allele
      vio_count <- vio_count + 1
    }
  } else {
  ▶ #child is het
    if ((mom_genotype == "1|1") & (dad_genotype == "1|1") |
        (mom_genotype == "0|0") & (dad_genotype == "0|0")) {
▶     #parents homo : no source of one of the two alleles in child
      vio_count <- vio_count + 1
    }
  } # IF-ELSE STATMENT END

} # FOR LOOP END

print(vio_count)
```
#### a. How many variants violate the rules of Mendelian inheritance?

```
2125
```

#### b. Describe some potential reasons that could explain these violations:

```
Students should touch on 2+ ('some') of the following:
1. De-novo mutation
2. Loss of heterozygosity due to DNA repair in meiosis
3. Errors in sequencing
4. Errors in mapping
5. Other reasonable answer

Students might see "non-mendelian" and think "LD" (at least that's what my brain auto-completed) but that's not relevant here where we consider only one variant at at time.
```

***

### **7.** De novo mutations are those that are present in a child and absent in both parents, potentially occurring as a mutation during gametogenesis during the parent. Write a R script to find _de novo_ variants in NA12878. To do this, you should identify variants that are heterozygous in the child but homozygous reference in both parents. The _de novo_ variants should be written in a new VCF file named Trio-NA12891-NA12892-NA12878.denovo.vcf. Use your preferred method to answer the following questions.

#### a. How many _de novo_ variants are observed?

```
669
```


```R
# load library
library(vcfR)

# read vcf
vcf <- read.vcfR("/gpfs/gibbs/project/gene760/shared/folders/Data/PS4/Trio-NA12891-NA12892-NA12878.vcf.gz", verbose = FALSE)

# store the row index which is a de novo variant
de_novo_index <- c()

for (i in 1:nrow(vcf@gt)){
  # what is genotype of the child
  child_genotype = vcf@gt[i,"NA12878"]

  # what is genotype of the mom
  mom_genotype = vcf@gt[i,"NA12892"]

  # what is genotype of the dad
  dad_genotype = vcf@gt[i,"NA12891"]

  # in which conditions, the variant is a de novo variant
  if (((child_genotype == "0|1") | (child_genotype == "1|0")) &
      (mom_genotype == "0|0") &
      (dad_genotype == "0|0")){
    de_novo_index <- c(de_novo_index, i)
  }
}

# write the subset of vcf object using de_novo_index
write.vcf(vcf[de_novo_index,], file = "Trio-NA12891-NA12892-NA12878.denovo.vcf.gz")
```

Then

```bash
zcat Trio-NA12891-NA12892-NA12878.denovo.vcf | grep -vE '^#' | wc -l
```
#### b. How many _de novo_ variants are present in gnomAD? How many are in the GIAB difficult regions of the genome? Should we be more confident or less confident that these variants (present in gnomAD or overlap with GAB difficult regions) are true _de novo_ variants? Why?

(You can assume that those variants observed in gnomad have passed quality control filters and are compatible with life + not sequencing errors. )

```
De-novo in gnomad = 606
De-novo in GIAB difficult regions = 604

Those in difficult regions are less trustworthy. Those in gnomad are more trustworthy.
```

```bash
# Students might run gnomad annotation again on resulting de-novo variants, or may have called de-novo variants on gnomad annotated VCF to begin with. 

dat_root="/gpfs/gibbs/project/gene760/shared/folders/Data/PS4/" # optional

#re-annotation, if required
${dat_root}slivar expr --vcf Trio-NA12891-NA12892-NA12878.denovo.vcf.gz --gnotate  ${dat_root}gnomad.genomes.v3.1.sites.chr1-22_XYM.zip --out-vcf Trio-NA12891-NA12892-NA12878.denovo.gnomAD_annotated.vcf.gz  > out_second 2>&1
## end of re-annotation.

#assuming re-annotation was performed:
bcftools view -i 'gnomad_af >0' Trio-NA12891-NA12892-NA12878.denovo.gnomAD_annotated.vcf.gz | grep -vE '^#' | wc -l
#If re-annotation wasn't performed, replace filename with one from earlier step : just denovo
#produces 606

module load BEDTools

bedtools intersect -a Trio-NA12891-NA12892-NA12878.denovo.gnomAD_annotated.vcf.gz -b ${dat_root}GRCh38_alldifficultregions.bed > difficult.bed
#quickly checking that there are no header lines... (shouldn't be but...)

wc -l difficult.bed
#604
```

