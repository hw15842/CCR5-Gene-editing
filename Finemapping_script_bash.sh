#!/bin/bash


###################################################################
#### This is all the code for making the finemapping data work ####
###################################################################

#just keeps the SNPs 500kb either side of proxy
cat <(sed -n 1p blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT.sumstats.gz_reordered_for_finemapping.txt) <(awk -F "\t" '{ if($2 == 3 && $3 >= 45915921 && $3 <= 46915921) { print }}' blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT.sumstats.gz_reordered_for_finemapping.txt) > HLSRC_500kb_either_side_proxy.txt

# removes duplicate IDS (duplicate due to being triallelic)
awk '!seen[$1]++' HLSRC_500kb_either_side_proxy.txt > HLSRC_500kb_either_side_proxy_triallelic_IDs_removed.txt

# keeps just the SNPs, no header
awk '{print $1'} HLSRC_500kb_either_side_proxy_triallelic_IDs_removed.txt | sed '1d' > HLSRC_SNPs.txt

# makes the LD matrix
plink --r2 square --bfile eur --extract HLSRC_SNPs.txt --out HLSRC_matrix

# to find out what SNPs were not removed, some removed sue to not being in the eur.bim file
plink --bfile eur --extract HLSRC_SNPs.txt --write-snplist --out HLSRC_SNPs_kept

# make sure SNPs in matrix same as ones in snp stats file
grep -Fwf HLSRC_SNPs_kept.snplist HLSRC_500kb_either_side_proxy_triallelic_IDs_removed.txt > HLSRC_SNP_data_in_matrix.txt


#### To check ### all should be the same length

wc -l HLSRC_SNPs_kept.snplist 
wc -l HLSRC_SNP_data_in_matrix.txt
wc -l HLSRC_matrix.ld

## Change to space deleminated file and add .z extension

sed 's/\t/ /g' HLSRC_SNP_data_in_matrix.txt  > HLSRC.z

sed 's/\t/ /g' HLSRC_matrix.ld > HLSRC.ld

## add headers back in 

sed  -i '1i rsid chromosome position noneff_allele eff_allele maf beta se' HLSRC.z

## LC ##

cat <(sed -n 1p blood_LYMPHOCYTE_COUNT.sumstats.gz_reordered_for_finemapping.txt) <(awk -F "\t" '{ if($2 == 3 && $3 >= 45915921 && $3 <= 46915921) { print }}' blood_LYMPHOCYTE_COUNT.sumstats.gz_reordered_for_finemapping.txt) > LC_500kb_either_side_proxy.txt          

awk '!seen[$1]++' LC_500kb_either_side_proxy.txt > LC_500kb_either_side_proxy_triallelic_IDs_removed.txt     

awk '{print $1'} LC_500kb_either_side_proxy_triallelic_IDs_removed.txt | sed '1d' > LC_SNPs.txt                                

plink --r2 square --bfile eur --extract LC_SNPs.txt --out LC_matrix                

plink --bfile eur --extract LC_SNPs.txt --write-snplist --out LC_SNPs_kept       

grep -Fwf LC_SNPs_kept.snplist LC_500kb_either_side_proxy_triallelic_IDs_removed.txt > LC_SNP_data_in_matrix.txt  


wc -l LC_SNPs_kept.snplist 
wc -l LC_SNP_data_in_matrix.txt
wc -l LC_matrix.ld


sed 's/\t/ /g' LC_SNP_data_in_matrix.txt  > LC.z
sed  -i '1i rsid chromosome position noneff_allele eff_allele maf beta se' LC.z

sed 's/\t/ /g' LC_matrix.ld > LC.ld

### MCH ###

cat <(sed -n 1p blood_MEAN_CORPUSCULAR_HEMOGLOBIN.sumstats.gz_reordered_for_finemapping.txt) <(awk -F "\t" '{ if($2 == 3 && $3 >= 45915921 && $3 <= 46915921) { print }}' blood_MEAN_CORPUSCULAR_HEMOGLOBIN.sumstats.gz_reordered_for_finemapping.txt) > MCH_500kb_either_side_proxy.txt          

awk '!seen[$1]++' MCH_500kb_either_side_proxy.txt > MCH_500kb_either_side_proxy_triallelic_IDs_removed.txt     

awk '{print $1'} MCH_500kb_either_side_proxy_triallelic_IDs_removed.txt | sed '1d' > MCH_SNPs.txt                                

plink --r2 square --bfile eur --extract MCH_SNPs.txt --out MCH_matrix                

plink --bfile eur --extract MCH_SNPs.txt --write-snplist --out MCH_SNPs_kept       

grep -Fwf MCH_SNPs_kept.snplist MCH_500kb_either_side_proxy_triallelic_IDs_removed.txt > MCH_SNP_data_in_matrix.txt  


wc -l MCH_SNPs_kept.snplist 
wc -l MCH_SNP_data_in_matrix.txt
wc -l MCH_matrix.ld

sed 's/\t/ /g' MCH_SNP_data_in_matrix.txt  > MCH.z
sed  -i '1i rsid chromosome position noneff_allele eff_allele maf beta se' MCH.z

sed 's/\t/ /g' MCH_matrix.ld > MCH.ld



### MSCV ###

cat <(sed -n 1p blood_MEAN_SPHERED_CELL_VOL.sumstats.gz_reordered_for_finemapping.txt) <(awk -F "\t" '{ if($2 == 3 && $3 >= 45915921 && $3 <= 46915921) { print }}' blood_MEAN_SPHERED_CELL_VOL.sumstats.gz_reordered_for_finemapping.txt) > MSCV_500kb_either_side_proxy.txt          

awk '!seen[$1]++' MSCV_500kb_either_side_proxy.txt > MSCV_500kb_either_side_proxy_triallelic_IDs_removed.txt     

awk '{print $1'} MSCV_500kb_either_side_proxy_triallelic_IDs_removed.txt | sed '1d' > MSCV_SNPs.txt                                

plink --r2 square --bfile eur --extract MSCV_SNPs.txt --out MSCV_matrix                

plink --bfile eur --extract MSCV_SNPs.txt --write-snplist --out MSCV_SNPs_kept       

grep -Fwf MSCV_SNPs_kept.snplist MSCV_500kb_either_side_proxy_triallelic_IDs_removed.txt > MSCV_SNP_data_in_matrix.txt  


wc -l MSCV_SNPs_kept.snplist 
wc -l MSCV_SNP_data_in_matrix.txt
wc -l MSCV_matrix.ld

sed 's/\t/ /g' MSCV_SNP_data_in_matrix.txt  > MSCV.z
sed  -i '1i rsid chromosome position noneff_allele eff_allele maf beta se' MSCV.z


sed 's/\t/ /g' MSCV_matrix.ld > MSCV.ld


### Height ###

cat <(sed -n 1p body_HEIGHTz.sumstats.gz_reordered_for_finemapping.txt) <(awk -F "\t" '{ if($2 == 3 && $3 >= 45915921 && $3 <= 46915921) { print }}' body_HEIGHTz.sumstats.gz_reordered_for_finemapping.txt) > height_500kb_either_side_proxy.txt          

awk '!seen[$1]++' height_500kb_either_side_proxy.txt > height_500kb_either_side_proxy_triallelic_IDs_removed.txt     

awk '{print $1'} height_500kb_either_side_proxy_triallelic_IDs_removed.txt | sed '1d' > height_SNPs.txt                                

plink --r2 square --bfile eur --extract height_SNPs.txt --out height_matrix                

plink --bfile eur --extract height_SNPs.txt --write-snplist --out height_SNPs_kept       

grep -Fwf height_SNPs_kept.snplist height_500kb_either_side_proxy_triallelic_IDs_removed.txt > height_SNP_data_in_matrix.txt  


wc -l height_SNPs_kept.snplist 
wc -l height_SNP_data_in_matrix.txt
wc -l height_matrix.ld

sed 's/\t/ /g' height_SNP_data_in_matrix.txt  > height.z
sed  -i '1i rsid chromosome position noneff_allele eff_allele maf beta se' height.z

sed 's/\t/ /g' height_matrix.ld > height.ld


### UC ##

cat <(sed -n 1p UC_chr3_with_positions_FMorder.txt) <(awk -F "\t" '{ if($2 == 3 && $3 >= 45915921 && $3 <= 46915921) { print }}' UC_chr3_with_positions_FMorder.txt) > UC_500kb_either_side_proxy.txt          

awk '!seen[$1]++' UC_500kb_either_side_proxy.txt > UC_500kb_either_side_proxy_triallelic_IDs_removed.txt     

awk '{print $1'} UC_500kb_either_side_proxy_triallelic_IDs_removed.txt | sed '1d' > UC_SNPs.txt                                

plink --r2 square --bfile eur --extract UC_SNPs.txt --out UC_matrix                

plink --bfile eur --extract UC_SNPs.txt --write-snplist --out UC_SNPs_kept       

grep -Fwf UC_SNPs_kept.snplist UC_500kb_either_side_proxy_triallelic_IDs_removed.txt > UC_SNP_data_in_matrix.txt  


wc -l UC_SNPs_kept.snplist 
wc -l UC_SNP_data_in_matrix.txt
wc -l UC_matrix.ld

sed 's/\t/ /g' UC_SNP_data_in_matrix.txt  > UC.z
sed  -i '1i rsid chromosome position noneff_allele eff_allele maf beta se' UC.z


sed 's/\t/ /g' UC_matrix.ld > UC.ld

