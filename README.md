# EGYREF_Diabetes_Obesity
A workflow descriptive steps for the comparison of AFs and Genotypes among the Egyptian population genome database and the 1000 genome different subpopulations.

# Genes Data Pre-processing
**1- 1000 genome 16 genes files concatenation, Deduplication, and normalization to Biallelic**

```
$ cd EgyRef/2022.metabolic_elhadidi
$ bcftools concat -D -o analysis/genes_1000_DEDUP -O z *_1000g.vcf.gz && bcftools norm -m-any -o analysis/genes_1000_DEDUP_biallelic.vcf.gz -O z analysis/genes_1000_DEDUP
```
**2- EGYREF 16 genes files concatenation, Deduplication, and normalization to Biallelic**
```
$ bcftools concat -D -o analysis/genes_EGYREF_DEDUP -O z *_egyptians.vcf.gz && bcftools norm -m-any -o analysis/genes_EGYREF_DEDUP_biallelic.vcf.gz -O z analysis/genes_EGYREF_DEDUP
```
**3- Creation of the intermediate files needed for the extraction of the Allele frequenceies for each subpopulation variant in the 1000 genome**

_In this step we aim to create a tab separated subfile from the orginial VCF file to be used as a grep target for the next step_
```
 $ zcat analysis/genes_1000_DEDUP_biallelic.vcf.gz| awk ' BEGIN{FS=OFS="\t"} gsub (";","\t",$8)’ | cut -f 1-17 > analysis/1000_all_tab
```
_Column 8 contains the info for the AFs of the individual subpopulations so we aim to change its separator to make each subpopulation in a separate columns and then cut all the needed information only_

**4- Same Creation of the intermediate file needed for the extraction of the Allele frequenceies for the variants in the EGYREF genome**
```
zcat analysis/genes_EGYREF_DEDUP_biallelic.vcf.gz | sed '/^#/d' |sed 's/;/\t/g'| cut -f 1-9 | sed 's/chr//g' > analysis/EGYREF_all_tab
```
_EGYREF VCF file chromosome field was normalized like the previous 1000g file by removing the "chr" from the begining of each entry at the first column_

**6- creating a list that contains all the chr_position from the 2 concatenated VCF files to be used as the base of the joining process in the next step**
```
$ (cat analysis/1000_all_tab | cut -f 1,2 | sed 's/\t/_/g' && cat analysis/EGYREF_all_tab | cut -f 1,2 | sed 's/\t/_/g')|sort | uniq > analysis/EGYREF_1000g_all_chr_positions  
```
**7- creation of the subpopulations' AFs files based on the existence of an AF value for such a chromosome_position in each subpopulation**

_To do so we have to join the 2 files [EGYREF_1000g_all_chr_positions] and [1000_all_tab] based on the chromosome and position and then check for the state of AF if it euqals zero or not. By the end we'll have 5 files for each subpopulation that contains the variants that exist in each of them resembled by chr_position_Alternative allele to be used for the intersection and comparison in the next step_

Here we Will be using a bash script like follows: 
```
$ cd analysis
$ nano script
#! /bin/bash
readarray -t pop < pop;
readarray -t colnum < colnum; 
for i in {0..4};
do
	awk -v j="${colnum[i]}" 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $1,$2,$5,$j,$3, a[$1,$2]}' EGYREF_1000g_all_chr_positions  1000_all_tab|
		 sed "s|${pop[i]}_AF=||g" |awk '$4 !=0 {print}' > "${pop[i]}"_sub_SharedPositions_AF &&
			sed 's/ /\t/g' "${pop[i]}"_sub_SharedPositions_AF|cut -f 1,2,3| sed 's/\t/_/g' > "${pop[i]}"_sub_SharedPositions;
done

$ bash script
```
_"pop" is a list that contains names of the subpopulations [EAS,AMR,AFR,EUR,SAS] each in a line. While "colnum" is another list with the same length that contains the column numbers where each subpopulation's AF exists [13,14,15,16,17] repectively._

**8- creation of the the EGYREF AF list**
```
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $1,$2,$5,$9,$3, a[$1,$2]}' EGYREF_1000g_all_chr_positions EGYREF_all_tab > EGYREF_SharedPositions_AF && cut -f 1,2,3 EGYREF_SharedPositions_AF| sed 's/\t/_/g' > EGYREF_SharedPositions
```
# Populations genotype principal component analysis (PCA) extraction

```
$ mkdir plink && cd plink
```
**1- Merging all the samples in one BCF file with unified annotation**
```
$ bcftools merge -0 -m all -o genes_all.vcf.gz -O z ../genes_1000_DEDUP_biallelic.vcf.gz genes_EGYREF_DEDUP_biallelic.vcf.gz
$ bcftools norm -m-any genes_all.vcf.gz | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' > genes_all.bcf; bcftools index genes_all.bcf
```
**2- making the bed file required for the plink**

_Here we're using plink v2_
```
$ plink2 --bcf genes_all.bcf --vcf-idspace-to _ --const-fid --allow-extra-chr 0  --make-bed --out genes_all --vcf-half-call r
```
**3- Variants Pruning with MAF=0.05 and indep-pairwise 50 5 0.5 and exclusion of the pruned varianted**
```
$ mkdir pruned && cd pruned
$ plink2 --bfile ../genes_all --maf 0.05 --indep-pairwise 50 5 0.5 --out genes_all; plink2 --bfile ../genes_all --extract pruned/genes_all.prune.in --make-bed --out pruned/genes_all
```
**4- extraction of the eigen values and eigen vectors of the PCA**
```
plink2 --bfile genes_all --pca
```
# R analysis and visualization
