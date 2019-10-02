###################################
##### Clumping_all_six_traits #####
###################################



library (devtools)
library(plyr)
library(data.table)
library(pkgcond)
library(TwoSampleMR)

setwd("/newhome/hw15842/PhD/CCR5/ALL_OTHER_TRAITS")

### Read in the summary stats (sig pvalues only) for the 6 sig traits 

temp = list.files(path="/newhome/hw15842/PhD/CCR5/ALL_OTHER_TRAITS", pattern="*sig_pvals.txt")
temp
for (i in 1:length(temp)) assign(temp[i], read.table(temp[i], header=T))

names(blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT_sig_pvals.txt) <- c("SNP", "CHR", "position", "effect_allele", "other_allele", "reference_allele", "eaf", "beta", "se", "pval", "N", "info")
names(blood_LYMPHOCYTE_COUNT_sig_pvals.txt) <- c("SNP", "CHR", "position", "effect_allele", "other_allele", "reference_allele", "eaf", "beta", "se", "pval", "N", "info")
names(blood_MEAN_CORPUSCULAR_HEMOGLOBIN_sig_pvals.txt) <- c("SNP", "CHR", "position", "effect_allele", "other_allele", "reference_allele", "eaf", "beta", "se", "pval", "N", "info")
names(blood_MEAN_SPHERED_CELL_VOL_sig_pvals.txt) <- c("SNP", "CHR", "position", "effect_allele", "other_allele", "reference_allele", "eaf", "beta", "se", "pval", "N", "info")
names(body_HEIGHTz_sig_pvals.txt) <- c("SNP", "CHR", "position", "effect_allele", "other_allele", "reference_allele", "eaf", "beta", "se", "pval", "N", "info")
names(ulcerative_colitis_sig_pvals.txt) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "N")

HLSRC <- blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT_sig_pvals.txt
LC <- blood_LYMPHOCYTE_COUNT_sig_pvals.txt
MCH <- blood_MEAN_CORPUSCULAR_HEMOGLOBIN_sig_pvals.txt
MSCV <- blood_MEAN_SPHERED_CELL_VOL_sig_pvals.txt
height <- body_HEIGHTz_sig_pvals.txt
UC <- ulcerative_colitis_sig_pvals.txt

traits_combined <- list(HLSRC, LC, MCH, MSCV, height, UC)
names(traits_combined) <- c("HLSRC", "LC", "MCH", "MSCV", "height", "UC")
lapply(traits_combined, head)



format_data_function <- function (trait){

	X <- format_data(trait, type="exposure")
	return(X)
}

traits_formatted <- lapply(traits_combined, format_data_function)


lapply(traits_formatted, head)

## CLUMP



clumping_by_chunk <- function(trait){

	data_split <- split(trait, (as.numeric(rownames(trait))-1) %/% 1000)
	section_numbers <- c(1:(length(data_split)))
		
		data_by_chunks <- function(section){
  		clumped_data_out <- try(suppress_messages((clump_data(data_split[[section]])), pattern="rs"))
	  	return(clumped_data_out)
		}

	X <- lapply(section_numbers, data_by_chunks)
	X <- ldply(X, data.frame)
	X <- clump_data(X) # Clump the data again so fully clumped 
	return(X)

}


try(HLSRC_clumped <- clumping_by_chunk(traits_formatted$HLSRC))
try(LC_clumped <- clumping_by_chunk(traits_formatted$LC))
try(MCH_clumped <- clumping_by_chunk(traits_formatted$MCH))
try(MSCV_clumped <- clumping_by_chunk(traits_formatted$MSCV))
try(height_clumped <- clumping_by_chunk(traits_formatted$height))
try(UC_clumped <- clumping_by_chunk(traits_formatted$UC))

try(write.table(HLSRC_clumped, file="HLSRC_clumped.txt", quote=F, sep="\t"))
try(write.table(LC_clumped, file="LC_clumped.txt", quote=F, sep="\t"))
try(write.table(MCH_clumped, file="MCH_clumped.txt", quote=F, sep="\t"))
try(write.table(MSCV_clumped, file="MSCV_clumped.txt", quote=F, sep="\t"))
try(write.table(height_clumped, file="height_clumped.txt", quote=F, sep="\t"))











