#rm(list=ls())

#Formatting allele groups
hwe <- data_final

hwe$allele1 <-
  gsub("B\\*", "", as.character(hwe$allele1))
hwe$allele2 <-
  gsub("B\\*", "", as.character(hwe$allele2))

hwe_split <- split(hwe, hwe$SamplesGroupYear)

alleles20 <- as.data.frame(paste("2020_sample", ","," ", hwe_split$`2020`$allele1, hwe_split$`2020`$allele2, sep = ""))
alleles21 <- as.data.frame(paste("2021_sample", ","," ", hwe_split$`2021`$allele1, hwe_split$`2021`$allele2, sep = ""))


#Build the Genepop format
sink("genepop_format.txt")
cat("Title: HMPGL UTI 2021 2020 REDOME HLA INFO \n")
cat("HLAB \n")
cat("Pop \n")
invisible(apply(alleles20, 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(alleles21, 1,function(x) cat(x,"\n")))
sink()


outfile <- test_HW("results/genepop_format.txt",which='Proba', outputFile = "results/output_genepop.txt",batches = 10000,iterations = 10000, dememorization = 10000)



