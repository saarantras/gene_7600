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