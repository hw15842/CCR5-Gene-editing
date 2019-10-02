
##########################
## MR-phewas_six_traits ##
##########################

# pull in args
args  <-  commandArgs(trailingOnly = TRUE)
Trait   <-  toString(args[1]) #Phenotypic exposure or outcome of interest
WD <- toString(args[2])

setwd(WD)

library (devtools)
library(plyr)
library(data.table)
library(TwoSampleMR)


### Read in the clumped SNPs for the traits we are going to use as the exposures ##
		## Specified in the .sh file which trait we are using for the exposure


read.file <- paste(Trait, "clumped.txt", sep="_")

exposure_trait <- read.table(read.file)

######################################
### perform MR on the MRbase traits ##
######################################

## MRbase outcome traits ##

outcome_traits_mrbase <- read.csv("list_of_395_traits.csv", header=T)

#remove NAs 
outcome_traits_mrbase <- outcome_traits_mrbase[!is.na(outcome_traits_mrbase$MR.Base.id),]

## Run the MR ##


run_MR <- function(outcome_ID){
  
  outcome_1 <- extract_outcome_data(snps = exposure_trait$SNP, outcome = outcome_ID, proxies=FALSE)
  dat1 <- harmonise_data(exposure_trait, outcome_1)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

results_MR_MRbase <- lapply(outcome_traits_mrbase$MR.Base.id, run_MR)
results_MR_MRbase <- ldply(results_MR_MRbase, data.frame)#
savefile1 <- paste (Trait, "results_MRbase.txt", sep="_")#
write.table(results_MR_MRbase, file = savefile1, quote=F, sep="\t")


######################################
### perform MR on the OTHER traits ##
######################################


#### Have already extracted out the clumped SNPs from the other traits to use as the outcome for each of the 6 exposures ###

## read in the SNPs for the exposure extracted from the outcomes



## First read in the SNPs from broad institute as have same column names so can format 

read_in_outcome_traits_1 <- function(outcome_trait) {

	pastefile <- paste("/newhome/hw15842/PhD/CCR5/MR-PheWAS_six_traits", Trait, outcome_trait, sep="/")
	X <- fread(pastefile)
	Y <- format_data(X, type="outcome", snp_col="SNP", beta_col="Beta", se_col="se", eaf_col="EAF", effect_allele_col="A1", other_allele_col="A2", pval_col="P", pos_col="pos", chr_col="CHR")
	return(Y)
}

broad_institute  = list.files(path=paste("/newhome/hw15842/PhD/CCR5/MR-PheWAS_six_traits", Trait, sep="/"), pattern="*sumstats.gz_.txt")

broad_institute 

Z1 <- lapply(broad_institute, read_in_outcome_traits_1)
names(Z1) <- broad_institute

str(Z1)
summary(Z1)




## Then read in the non broad institute ones so can format those 

read_in_outcome_traits_2 <- function(outcome_trait) {

	pastefile <- paste("/newhome/hw15842/PhD/CCR5/MR-PheWAS_six_traits", Trait, outcome_trait, sep="/")
	X <- fread(pastefile)
	Y <- format_data(X, type="outcome", snp_col="SNP", beta_col="b", se_col="se", eaf_col="freq", effect_allele_col="A1", other_allele_col="A2", pval_col="p")
	return(Y)
}

not_broad_institute = list.files(path=paste("/newhome/hw15842/PhD/CCR5/MR-PheWAS_six_traits", Trait, sep="/"), pattern="*.assoc.tsv.gz_.txt")
not_broad_institute

Z2 <- lapply(not_broad_institute, read_in_outcome_traits_2)
names(Z2) <- not_broad_institute

str(Z2)
summary(Z2)


### read in the AD outcome SNPs as different format again ### https://ctg.cncr.nl/documents/p1651/AD_sumstats.readme shows what each column is

read_in_outcome_traits_3 <- function(outcome_trait) {

	pastefile <- paste("/newhome/hw15842/PhD/CCR5/MR-PheWAS_six_traits", Trait, outcome_trait, sep="/")
	X <- fread(pastefile)
	Y <- format_data(X, type="outcome", snp_col="SNP", beta_col="BETA", se_col="SE", eaf_col="MAF", effect_allele_col="BP       A1", other_allele_col="A2", pval_col="P")
	return(Y)
}

AD = list.files(path=paste("/newhome/hw15842/PhD/CCR5/MR-PheWAS_six_traits", Trait, sep="/"), pattern="*AD_sumstats_Jansenetal.txt.gz_.txt")
AD


Z3 <- lapply(AD, read_in_outcome_traits_3)
names(Z3) <- AD

str(Z3)
summary(Z3)


outcome_list <- c(Z1, Z2, Z3)

str(outcome_list)
summary(outcome_list)


## Run the MR

run_MR <- function(outcome_ID){
  
  #outcome_1 <- extract_outcome_data(snps = exposure_trait$SNP, outcome = outcome_ID, proxies=FALSE)
  dat1 <- harmonise_data(exposure_trait, outcome_ID)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

results_MR_other <- lapply(outcome_list, run_MR)
results_MR_other <- ldply(results_MR_other, data.frame)

savefile2 <- paste (Trait, "results_other.txt", sep="_")

write.table(results_MR_other, file = savefile2, quote=F, sep="\t")



### merge together 


head(results_MR_MRbase)
head(results_MR_other)

results_MR_six_traits_all <- rbind(results_MR_MRbase, results_MR_other)

write.table(results_MR_six_traits_all, file="results_MR_six_traits_all.txt", quote=F, sep="\t", row.names=F)

























