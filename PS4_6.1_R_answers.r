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