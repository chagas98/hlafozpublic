rm(list=ls())

#Formatting allele groups
rawtable20$B_AlleleGroup1 <-
  gsub("B\\*", "", as.character(rawtable20$B_AlleleGroup1))
rawtable20$B_AlleleGroup2 <-
  gsub("B\\*", "", as.character(rawtable20$B_AlleleGroup2))
rawtable21$B_AlleleGroup1 <-
  gsub("B\\*", "", as.character(rawtable21$B_AlleleGroup1))
rawtable21$B_AlleleGroup2 <-
  gsub("B\\*", "", as.character(rawtable21$B_AlleleGroup2))

alleles20 <- as.data.frame(paste("2020_sample", ","," ", rawtable20$B_AlleleGroup1, rawtable20$B_AlleleGroup2, sep = ""))
alleles21 <- as.data.frame(paste("2021_sample", ","," ", rawtable21$B_AlleleGroup1, rawtable21$B_AlleleGroup2, sep = ""))

#Build the Genepop format
sink("genepop_format.txt")
cat("Title: HMPGL UTI 2021 2020 REDOME HLA INFO \n")
cat("HLAB \n")
cat("Pop \n")
invisible(apply(alleles20, 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(alleles21, 1,function(x) cat(x,"\n")))
sink()


outfile <- test_HW("genepop_format.txt",which='Proba', outputFile = "output/output_genepop.txt",batches = 100000,iterations = 100000, dememorization = 100000)



