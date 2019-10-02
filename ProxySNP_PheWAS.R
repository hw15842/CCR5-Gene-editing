library (devtools)
library(plyr)
library(data.table)

install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)

ao <- available_outcomes()

## Outcome traits that are in MRbase
outcome_traits_mrbase <- read.csv("list_of_395_traits.csv", header=T)

## Remove NAs from mr.base.ids column in outcome_traits_mrbase ##

outcome_traits_mrbase_na.rm <- outcome_traits_mrbase[!is.na(outcome_traits_mrbase$MR.Base.id),]



## Check which outcomes are giving nulls and remove them from the dataset
outcomes_check_nulls_ <- function(outcome_ID){
  extract_outcome_data(snps = exposure_CCR5$SNP, outcome = outcome_ID )
}

outcomes_with_nulls <- lapply(outcome_traits_mrbase_na.rm$MR.Base.id, outcomes_check_nulls_)
outcomes_with_nulls_2 <- ldply(outcomes_with_nulls, data.frame)

df_new <- sapply(outcomes_with_nulls, is.null)
table(df_new) 
## rows 6 and 9 are null so remove those from the dataset

outcome_traits_mrbase_na.rm <- outcome_traits_mrbase_na.rm[-c(6,9),]


### Running PheWAS for the proxy  ###

## MR base traits ##

proxy_phewas <- function(outcome_ID){
  
  outcome_1 <- extract_outcome_data(snps = "rs113010081", outcome = outcome_ID)
  return(outcome_1)
}

results_proxy_phewas <- lapply(outcome_traits_mrbase$MR.Base.id, proxy_phewas)

results_proxy_phewas <- ldply(results_proxy_phewas, data.frame)

write.table(results_proxy_phewas_MRbase, file="results_proxy_phewas_MRbase_rs113010081.txt", sep="\t")


