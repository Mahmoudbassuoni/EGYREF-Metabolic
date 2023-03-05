# EGYREF-Diabetes 
A workflow descriptive pipline for the comparison of AFs and Genotypes among the Egyptian population genome database and the 1000 genome different subpopulations for the genes responsible for the 4 metabolic diseases [Diabetes,Obesity,Hypertension,hyperlipidemia].Total of 19 genes were found to be sharing variants of the 4 diseases [ALDH2, APOE, CDKAL1, FTO, GIPR, HECTD4, LPL, MIR642A, MIR642B, MIR6761, OAS1, PDILT, POC5, PPARG, RNU6-680P, RPL7AP60, SLC39A8, TOMM40, UMOD].The aim for this project is to step on the genetic distance among the 1000 genome populations and the Egyptian one in those speicif diseases. Data for this project was obtained from [Wohlers, I., Künstner, A., Munz, M. et al. An integrated personal and population-based Egyptian genome reference. Nat Commun 11, 4719 (2020). https://doi.org/10.1038/s41467-020-17964-1].

# Genes Data Pre-processing
**1- 1000 genome 16 genes files indexing, concatenation, Deduplication, and normalization to Biallelic**

```
mkdir analysis
```
```
cd ~/EgyRef/2022.metabolic_elhadidi \ 
;readarray -t genes < genes \
; for i in ${genes[@]}; do bcftools index ${i}/${i}_1000g.vcf.gz; done \ 
;bcftools concat -D -a -o analysis/genes_1000_DEDUP -O z */*_1000g.vcf.gz && bcftools norm -m-any -o analysis/genes_1000_DEDUP_biallelic.vcf.gz -O z analysis/genes_1000_DEDUP
```
**2- EGYREF 16 genes files concatenation, Deduplication, and normalization to Biallelic**
```
for i in ${genes[@]}; do bcftools index ${i}/${i}_egyptians.vcf.gz; done \
; bcftools concat -D -a -o analysis/genes_EGYREF_DEDUP -O z */*_egyptians.vcf.gz && bcftools norm -m-any -o analysis/genes_EGYREF_DEDUP_biallelic.vcf.gz -O z analysis/genes_EGYREF_DEDUP
```
**3- Creation of the intermediate files needed for the extraction of the Allele frequenceies for each subpopulation variant in the 1000 genome**

_In this step we aim to create a tab separated subfile from the orginial VCF file to be used as a grep target for the next step_
```
zcat analysis/genes_1000_DEDUP_biallelic.vcf.gz| awk 'BEGIN{FS=OFS="\t"} gsub (";","\t",$8)' | cut -f 1-17 > analysis/1000_all_tab
```
_Column 8 contains the info for the AFs of the individual subpopulations so we aim to change its separator to make each subpopulation in a separate columns and then cut all the needed information only_

**4- Same Creation of the intermediate file needed for the extraction of the Allele frequenceies for the variants in the EGYREF genome**
```
zcat analysis/genes_EGYREF_DEDUP_biallelic.vcf.gz | sed '/^#/d' |sed 's/;/\t/g'| cut -f 1-9 | sed 's/chr//g' > analysis/EGYREF_all_tab
```
_EGYREF VCF file chromosome field was normalized like the previous 1000g file by removing the "chr" from the begining of each entry at the first column_

**6- creating a list that contains all the chr_position from the 2 concatenated VCF files to be used as the base of the joining process in the next step**
```
(cat analysis/1000_all_tab | cut -f 1,2 | sed 's/\t/_/g' && cat analysis/EGYREF_all_tab | cut -f 1,2 | sed 's/\t/_/g')|sort | uniq > analysis/EGYREF_1000g_all_chr_positions  
```
**7- creation of the subpopulations' AFs files based on the existence of an AF value for such a chromosome_position in each subpopulation**

_To do so we have to join the 2 files [EGYREF_1000g_all_chr_positions] and [1000_all_tab] based on the chromosome and position and then check for the state of AF if it euqals zero or not. By the end we'll have 5 files for each subpopulation that contains the variants that exist in each of them resembled by chr_position_Alternative allele to be used for the intersection and comparison in the next step_

Here we Will be using a bash script like follows: 
```
cd analysis ; \
nano script
```
```
#! /bin/bash
readarray -t pop < pop;
readarray -t colnum < colnum; 
sed -i 's/_/\t/g' EGYREF_1000g_all_chr_positions;
for i in {0..4};
do
	awk -v j="${colnum[i]}" 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $1,$2,$5,$j,$3, a[$1,$2]}' EGYREF_1000g_all_chr_positions  1000_all_tab|
		 sed "s|${pop[i]}_AF=||g" |awk '$4 !=0 {print}' > "${pop[i]}"_sub_SharedPositions_AF &&
			sed 's/ /\t/g' "${pop[i]}"_sub_SharedPositions_AF|cut -f 1,2,3| sed 's/\t/_/g' > "${pop[i]}"_sub_SharedPositions;
done
```
```
bash script
```
_"pop" is a list that contains names of the subpopulations [EAS,AMR,AFR,EUR,SAS] each in a line. While "colnum" is another list with the same length that contains the column numbers where each subpopulation's AF exists [13,14,15,16,17] repectively._

**8- creation of the the EGYREF AF list**
```
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $1,$2,$5,$9,$3, a[$1,$2]}' EGYREF_1000g_all_chr_positions EGYREF_all_tab > EGYREF_SharedPositions_AF && cut -f 1,2,3 EGYREF_SharedPositions_AF| sed 's/\t/_/g' > EGYREF_SharedPositions
```
# Populations genotype principal component analysis (PCA) extraction

```
mkdir plink ; cd plink
```
**1- Samples sorting ,indexing and finally Annotating the IDs both EGYREF and 1000g**
```
bcftools sort ../genes_EGYREF_DEDUP_biallelic.vcf.gz -o genes_EGYREF_DEDUP_biallelic_sorted.vcf.gz -O z ;bcftools index genes_EGYREF_DEDUP_biallelic_sorted.vcf.gz ;bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' genes_EGYREF_DEDUP_biallelic_sorted.vcf.gz  -o EGYREF.vcf.gz -O z
```
```
bcftools sort ../genes_1000_DEDUP_biallelic.vcf.gz -o genes_1000_DEDUP_biallelic_sorted.vcf.gz -O z ;bcftools index genes_1000_DEDUP_biallelic_sorted.vcf.gz ;bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' genes_1000_DEDUP_biallelic_sorted.vcf.gz -o 1000g.vcf.gz -O z 
```
**2- Files preparation for the merging process using R**
```
zcat EGYREF.vcf.gz | sed '/^##/d' > EGY.tsv; sed -i 's/chr//g' EGY.tsv
```
```
zcat 1000g.vcf.gz | sed '/^##/d' > g1000.tsv
```
**3- Joining the 2 files to get the common IDs with their samples genotypes in one file**
```
library(plyr)
merged <-join(x= EGY,y = g1000, by="ID",type ="inner",match="all")
write.table(merged,quote = FALSE,row.names = FALSE,sep = "\t","~/EgyRef/2022.metabolic_elhadidi/analysis/plink/merged.vcf")
```
**4- Clean the merge output from the extra columns and adjust the header for the bcftools to recognize it**
```
cat merged.vcf | cut -f 120-127 --complement| sed '1s/^/#/' > merged_complete.vcf
```
```
zcat ../genes_EGYREF_DEDUP_biallelic.vcf.gz | grep "##" > header; cat header merged_complete.vcf > temp; mv temp merged_complete.vcf
```
**5- bcftools final sorting and indexing before going to the plink**
```
bgzip merged_complete.vcf; tabix merged_complete.vcf.gz
```
```
bcftools sort merged_complete.vcf.gz -o merged_complete_sorted.vcf.gz -O z ; bcftools index merged_complete_sorted.vcf.gz
```
**6- making the bed file required for the plink**

_Here we're using plink v2_
```
plink2 --vcf merged_complete_sorted.vcf.gz --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --make-pgen --sort-vars --out genes_all  --vcf-half-call r
```
```
plink2 --pfile genes_all --make-bed --out genes_all 
```
**7- Variants Pruning with MAF=0.05 and indep-pairwise 50 5 0.5 and exclusion of the pruned varianted**
```
mkdir pruned && cd pruned
```
```
plink2 --bfile ../genes_all --maf 0.05 --indep-pairwise 50 5 0.5 --out genes_all; plink2 --bfile ../genes_all --extract genes_all.prune.in --make-bed --out genes_all
```
**8- extraction of the eigen values and eigen vectors of the PCA**
```
mkdir PCA && cd PCA
```
```
plink2 --bfile ../genes_all --pca
```
# R analysis and visualization
**1- Subpopulation lists intersection with the Egyptian list**

_UpSetR package was used to draw the intersection of the six populations’ lists_
```
library(UpSetR)
x <- list('AFR'=AFR_sub_SharedPositions$V1, 'EGYREF'=EGYREF_SharedPositions$V1, 
      'AMR'= AMR_sub_SharedPositions$V1,'EUR'= EUR_sub_SharedPositions$V1,
      'EAS'=EAS_sub_SharedPositions$V1,
      'SAS' = SAS_sub_SharedPositions$V1 )

upset(fromList(x),sets = c("EGYREF","AFR","AMR","EAS","EUR","SAS"),sets.x.label="Sets Size",
      show.numbers = "yes",line.size = 0.5,set_size.numbers_size=10,set_size.scale_max= 70000
      ,order.by = "freq",matrix.color = "gray32",main.bar.color = "brown",
      set_size.show=TRUE,sets.bar.color="blue4",point.size = 2,text.scale=c(1.3, 1.4, 1.3, 1, 1.3, 1.1))
```
**2- Allele Frequencies heatmap and PCA Visualization**

_(1) Join function in R “plyr” package was used on the 1000 and the Egyptian genome files
based on the chromosome, position, and alternative allele columns_
```
library(plyr)
joined <- join(x= all_tab_1000,y = all_tab_Egyref, by= c("V1"="V1","V2"="V2","V5"="V5"),type ="inner",match="all")
write.table(joined,sep = "\t","~/EgyRef/2022.metabolic_elhadidi/analysis/joined")
```
_(2) AFs [joined] table is Prepared for a heatmap to be drawn_ 

```
cd ~/EgyRef/2022.metabolic_elhadidi/analysis
```
```
sed -i 's/"//g' joined; \
sed 's/NA/0/g' joined| sed 's/EAS_AF=//g'| sed 's/AMR_AF=//g'|sed 's/EUR_AF=//g' | sed 's/AFR_AF=//g'| sed 's/SAS_AF=//g'|cut -f 2,3,6,14-18,24 |sed 's/AF=//g'|sed '1d' |sed '1i\Chr\tPosition\tALT_Allele\tEAS_AF\tAMR_AF\tAFR_AF\tEUR_AF\tSAS_AF\tEGYREF_AF' > joined_full
```
_(3) heatmap Visualization_
```
All_var <- joined_full[,-3][,-2][,-1]
matrix <- as.matrix(All_var)
pheatmap::pheatmap(mat = matrix,scale = "column", cellheight = 0.02, cellwidth = 60 ,main = "Common positions Allele Frequencies",
                   labels_col =c("EAS","AMR","AFR","EUR","SAS","EGYREF"),angle_col = 315,
                   fontsize_col= 13,cluster_cols= TRUE, cluster_rows= TRUE, color=colorRampPalette(c("navy", "white", "red"))(50) )

```
_(4) calculation and visualization the AFs' PCA_

```
#Load needed libraries 
library(plotly)
library("RColorBrewer")
#Load the main file AF
joined_full <- data.frame(read.table("~/EgyRef/2022.metabolic_elhadidi/analysis/joined_full"
                                     , header=TRUE, skip=0,))
#Convert the data frame to a matrix and clean the first 3 columns to retain only numbers  
pop_mat <- as.matrix(joined_full[,-3][,-2][,-1])
#Transpose the matrix
pop_mat_t <- t(pop_mat)
#extract populations' names and set the populations' names as to be used as metadata 
x <- colnames(pop_mat)
pop<- data.frame(x)
#Calculate the PCA
PCA <- prcomp(pop_mat_t,rank. = 3)
#check the Cumlative proportions for each PC
PCA_Summary <- summary(PCA)
PCA_importance <- PCA_Summary$importance
#extract the PCA components
components <- PCA[['x']]
components <- data.frame(components)
# PCs proportion of variances
PC1_proportion <- PCA_importance[2,1]
PC2_proportion <- PCA_importance[2,2]
PC3_proportion <- PCA_importance[2,3]
#Plotting the PCA in 2D brewer.pal(n = 8, name = "Set1")
tit = 'PCA Populations total variations'
fig <- plot_ly(components, x=components$PC1,y=components$PC2, color = pop$x, colors = brewer.pal(n = 8, name = "Set1")  , name = pop$x) %>%
  add_markers(size = 20)
fig <- fig %>%
  layout(
    title = tit,
    xaxis=list(title=paste0("PC1 (", scales::percent(PC1_proportion),")")),
    yaxis=list(title=paste0("PC2 (", scales::percent(PC2_proportion),")")),
    scene = list(bgcolor = "#e5ecf6")
  )

fig
```
**3- Genotype PCA Visualization**

_** samples_names were extracted from the merged file and then joined with the meta data file downloaded from the 1000g website (Auton et al., 2015)to annotate the 1000g samples, while the rest 110 samples were the Egyptian samples obtained from  (Wohlers et al., 2020)._

```
#Loading needed libraries 
library(plotly)
library(plyr)
library(readxl)
library(RColorBrewer)
#Loading the EigenValues and EigenVectors
eigenvectors <- data.frame(read.table("~/EgyRef/2022.metabolic_elhadidi/analysis/plink/pruned/plink2.eigenvec"
                               , header=FALSE, skip=0,))
eigenvalues <- data.frame(read.table("~/EgyRef/2022.metabolic_elhadidi/analysis/plink/pruned/plink2.eigenval"
                                     , header=FALSE, skip=0,))
#Data_Cleaning 
rownames(eigenvectors) <- eigenvectors[,2]
eigenvectors <- eigenvectors[,3:ncol(eigenvectors)]
#Check Summary 
summary(eigenvectors)
#Getting the Proportion of Variances in different PCs
proportionvariances <- data.frame(PC = 1:10, proportionvariances = eigenvalues/sum(eigenvalues)*100)
#Visualization of the proportions of variances
a <- ggplot(proportionvariances,aes(x = PC,y=V1))+ geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
# calculate the cumulative sum of the percentage variance explained 
cumsum(proportionvariances$V1)
#2D Basic plot 
plot(eigenvectors[,1], eigenvectors[,2])
#Loading the samples metadata (populations)
samples <- data.frame(read.table("~/EgyRef/2022.metabolic_elhadidi/analysis/plink/samples_pop_meta.csv"
                                 , header=FALSE, skip=0,sep = "\t"))
samples_names <- data.frame(read.table("~/EgyRef/2022.metabolic_elhadidi/analysis/plink/samples_names"
                                 , header=FALSE, skip=0,sep = "\t"))
#sorting the samples based on the first columns(names) by joining the 2 files from the 1000g website and the samples we had in our file
samples <- join(samples_names,samples,by="V1")

#Better looking Visualization with colors annotation 
names(eigenvectors)[1] <- "PC1"
names(eigenvectors)[2] <- "PC2"
names(eigenvectors)[3] <- "PC3"
names(eigenvectors)[4] <- "PC4"
names(eigenvectors)[5] <- "PC5"
PC1_proportion <- round(proportionvariances[1,2])/100 
PC2_proportion <- round(proportionvariances[2,2])/100 
PC3_proportion <- round(proportionvariances[3,2])/100 
tit = 'PCA genotypes'
#axx <- list(title=paste0("PC1 (", scales::percent(PC1_proportion),")"))
#axy <- list(title=paste0("PC2 (", scales::percent(PC2_proportion),")"))
#axz <- list(title=paste0("PC3 (", scales::percent(PC3_proportion),")"))

fig <- plot_ly(eigenvectors,x= ~PC1,y= ~PC2,z=~PC3, color = samples_annotation$X.2, colors = brewer.pal(n = 8, name = "Set1")) %>%
  add_markers(size = 20)
fig <- fig %>%
  layout(
    title = tit,
    xaxis=list(title=paste0("PC1 (", scales::percent(PC1_proportion),")")),
    yaxis=list(title=paste0("PC2 (", scales::percent(PC2_proportion),")")),
    zaxis=list(title=paste0("PC3 (", scales::percent(PC3_proportion),")")),
    scene = list(bgcolor = "white",xaxis=axz,yaxis=axy,zaxis=axx)
  )

fig
```
