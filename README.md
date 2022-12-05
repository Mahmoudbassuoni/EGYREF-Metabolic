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
$ nano script

```
