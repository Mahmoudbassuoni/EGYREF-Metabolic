# EGYREF_Diabetes_Obesity
A workflow descriptive steps for the comparison of AFs and Genotypes among the Egyptian population genome database and the 1000 genome different subpopulations.

# Genes Data Pre-processing

`$ cd EgyRef/2022.metabolic_elhadidi`
 '$ bcftools concat -D -o new/genes_1000_DEDUP -O z *_1000g.vcf.gz'
