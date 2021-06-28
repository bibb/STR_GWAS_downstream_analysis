# STR imputation GWAS meta-analysis in Parkinson's Disease
#### Author: Bernabe Bustos. Northwestern University. Chicago. Il.
## Downstream analysis
### Contents:
1. GCTA conditional-joint analysis on meta-analisys results (STRs   and SNPs+STRs)
2. Hudson plots for STRs and SNPs
3. LD analysis with surrounding 1+-MB STRs and meta-5 risk PD loci
4. eSTR analysis with NABEC and GTEx brain tissues
5. GCTA-LDMS heritability calculation for STRs, SNPs+STRs and SNPs
6. MAGMA gene-wise, gene-set and pathway enrichment analyses

## 1. GCTA conditional-joint analysis on meta-analisys results (STRs and SNPs+STRs)

```
#Input file used from meta-analysis including variant ID, alleles, frequency in meta-analysis, beta, se, p value and sample size: 

head SNPs_STR_I80.wRealSample_size.wHeader.v3.txt
SNP	A1	A2	Freq	b	se	p	N
1_74790_C_G	C	G	0.9854	0.0439	0.2309	0.8494	4705
1_74792_G_A	A	G	0.01456	-0.0439	0.2309	0.8494	4705
1_649192_A_T	A	T	0.9027	-0.0597	0.1362	0.6613	3721
1_657788_G_C	C	G	0.8904	-0.0681	0.1318	0.6052	3721
1_662029_G_A	A	G	0.886	-0.0910	0.1191	0.4446	4705
1_662201_G_A	A	G	0.8807	-0.0821	0.1174	0.4844	4705
1_662414_C_T	T	C	0.107	0.0924	0.1331	0.4872	3721
1_662622_G_A	A	G	0.09904	0.0686	0.1239	0.5798	4705
1_662857_G_A	A	G	0.8802	-0.0134	0.1027	0.8959	7349

#data was also filtered for I2 < 80%

#genotypes used:  chr1-22.split.SNPs_STRs.PDGWAS.uniqueID.DR2_03.merged. It is a GWAS subset from the meta-analysis (~5K PD cases and ~10K controls).

#Run GCTA COJO for STRs only

gcta64 --bfile chr1-22.split.SNPs_STRs.PDGWAS.uniqueID.DR2_03.merged \
--maf 0.01 \
--extract STRs_only_RawID.txt \
--cojo-file SNPs_STR_I80.wRealSample_size.wHeader.v3.txt \
--cojo-slct \
--cojo-p 5.34e-6 \
--out STR_only.gcta.cojo.v3_wLocalGWAS

#Run GCTA COJO for STRs+SNPs

gcta64 --bfile chr1-22.split.SNPs_STRs.PDGWAS.uniqueID.DR2_03.merged \
--maf 0.01 \
--cojo-file SNPs_STR_I80.wRealSample_size.wHeader.v3.txt \
--cojo-slct \
--cojo-p 5.34e-6 \
--out SNPs_STRs.gcta.cojo.v3_wLocalGWAS
```

## 2. Hudson plot for STRs and SNPs

```
#STR meta-analysis data

head Full.STR.meta.summary.stats.I80.Chr_BP_P.wID.txt
SNP	CHR	POS	pvalue
2	19	17030013	0.6479
3	3	187113276	0.921
4	20	13321137	0.283
5	4	169891985	0.1468
6	9	18812966	0.07399
7	14	99804514	0.4961
8	17	70546828	0.0879
9	15	47417999	0.6157
10	4	187874150	0.7323

#SNP meta-analysis data

head SNPs.Rdata.wID.txt
SNP	CHR	POS	pvalue
2	1	74792	0.8494
3	1	649192	0.6613
4	1	657788	0.6052
5	1	662029	0.4446
6	1	662201	0.4844
7	1	662414	0.4872
8	1	662622	0.5798
9	1	662857	0.8959
10	1	663097	0.9599

# hudson plot in R/3.6.3

library(hudson)

gwas.t <-read.table("Full.STR.meta.summary.stats.I80.Chr_BP_P.wID.txt",header=T)
gwas.b <-read.table("SNPs.Rdata.wID.txt",header=T)
gmirror(top=gwas.t, bottom=gwas.b, tline=5.34E-6, bline=5E-8,
  chrcolor1 = "dodgerblue4", chrcolor2 = "Deepskyblue4",
  background = "grey84", file = "STR_SNPs_HudsonPlot.v6", hgt = 14, hgtratio = 0.5, wi = 24, res = 300)

```
## 3. LD analysis with surrounding 1+-MB STRs and meta-5 risk PD loci

```
#filtered the 34 top STRs obtained after GCTA-COJO

wc -l candidate_TopSTRs.txt
34

#Obtain variants with high LD (r2>0.5) within 1 mb up- and downstream

plink --bfile chr1-22.split.SNPs_STRs.PDGWAS.uniqueID.DR2_03.merged \
--extract STRs_only_RawID.txt \
--r2 \
--ld-snp-list candidate_TopSTRs.txt \
--ld-window-kb 1000 \
--ld-window-r2 0.5 \
--out TOPSTRs_vs_allSTRs_LD05

#Obtain variants with high LD (r2>0.5) with variants from meta-5 PD GWAS

#file GW_STR_wMETA5_top.RawID.txt > all genome-wide STRs from meta + Meta5 PD risk SNPs

wc -l GW_STR_wMETA5_top.RawID.txt 
301

#file meta_5_90_SNPs.RawID_inIPDGCGWAS.txt > SNPs from meta-5 90 loci present in genotype data

wc -l meta_5_90_SNPs.RawID_inIPDGCGWAS.txt
87

#run LD

plink --bfile chr1-22.split.SNPs_STRs.PDGWAS.uniqueID.DR2_03.merged \
--extract GW_STR_wMETA5_top.RawID.txt \
--r2 \
--ld-snp-list meta_5_90_SNPs.RawID_inIPDGCGWAS.txt \
--ld-window-kb 1000 \
--ld-window-r2 0.5 \
--out TOPSTRs_wTOPMETA5_LD05

```
#### Run local Locuszoom plots for the 4 independent STRs
```
# This software requires a file in "metal" format, with 2 columns: SNPID and p value. SNPIDs must be in rsID format or chr:pos
# For STRs there is no currently available LD reference panel to download, therefore a custom file was made with PLINK, using the PD GWAS cohort of ~5K cases and 10K controls.
# A personalized dataset for gene annotations and variat localization was made based on the STR positions. Gene coordinates were obtained from gencode19
# The following code is to generate a plot for the candidate STR in chromosome 17.

locuszoom \
--plotonly \
--verbose \
--metal STR_I80.ID_P_N.ForLocusZoom.colonID.wCHR.sorted.v2.SNP_P.txt \
--refsnp chr17:15941750:TTTTTTTTTT \
--flank 400kb \
--build hg19 \
--ld 17_15941750_TTTTTTTTTT.localPDGWAS_LD_forLocuszoom.wDprime.ld.formatted.final.wCHR.txt \
--delim tab \
--pvalcol P \
--markercol SNP \
--snpset NULL \
geneFontSize=.8 \
smallDot=.3 \
largeDot=.9 \
format=pdf \
ymax=50 \
legend=auto \
metalRug='Plotted SNPs' \
rfrows=10 \
weightCol=Weight \
warnMissingGenes=T \
showAnnot=FALSE \
--db STR_database.db \
--prefix 17_15941750_TTTTTTTTTT.localPDGWAS
```

## 4. eSTR analysis with NABEC and GTEx brain tissue RNASEq data

### NABEC data
#### NABEC imputed VCF filtering and formatting steps

```
# Index imputed VCFs

tabix -p vcf ${i}.nabec.snp.str.imputed.vcf.gz

# Print variants with "STR" in ID

bcftools view -H ${i}.nabec.snp.str.imputed.vcf.gz |\
cut -f3 |\
awk '\$1 ~ /STR/ {print}' > chr${i}.STRs_IDs.txt

# Extract STR variants

bcftools view -i 'ID=@chr${i}.STRs_IDs.txt' -Oz -o ${i}.nabec.imputed.STROnly.bcftools.vcf.gz ${i}.nabec.snp.str.imputed.vcf.gz

# Modify header for multiallelic spliting with VT

zcat ${i}.nabec.imputed.STROnly.bcftools.vcf.gz  |\
sed 's/##INFO=<ID=DR2,Number=1/##INFO=<ID=DR2,Number=A/g' |\
bgzip -c > ${i}.nabec.imputed.STROnly.bcftools.DR2A.vcf.gz

vt decompose -s ${i}.nabec.imputed.STROnly.bcftools.DR2A.vcf.gz -o VTSplit.${i}.nabec.imputed.STROnly.bcftools.DR2A.vcf.gz

# Filter sites with DR < 0.3, recode variant ID to chr_pos_alt and write IDs to a file

bcftools view -e 'INFO/DR2 < 0.3' VTSplit.${i}.nabec.imputed.STROnly.bcftools.DR2A.vcf.gz -Ou |\
bcftools annotate --set-id '%CHROM\_%POS\_%ALT' |\
bcftools view -H |\
cut -f3 > chr${i}.STRs_DR2_03.RawID.txt

# Second multiallelic spliting now with bcftools on non-split vcf file. Write "raw" variant ID, remove variants with DR2 < 0.3 according to VT multiallelic splitting and filtering. Filter variants with missingness rate > 5%, HWE P < 0.0001 and MAF < 1%

bcftools norm -m-both ${i}.nabec.imputed.STROnly.bcftools.DR2A.vcf.gz -Ou |\
bcftools annotate --set-id '%CHROM\_%POS\_%ALT' -Ou |\
bcftools view -i 'ID=@chr${i}.STRs_DR2_03.RawID.txt' |\
bcftools +fill-tags -Ou |\
bcftools filter -e 'F_MISSING > 0.05' -e 'INFO/HWE < 0.0001' -Ou |\
bcftools view -q 0.01 -Oz -o BCFTSplit.${i}.nabec.imputed.STROnly.bcftools.DR2A_03.MAF01_GENO05_HWE0001.vcf.gz

#merge into a single file

ls -ltrh BCFTSplit.*.nabec.imputed.STROnly.bcftools.DR2A_03.MAF01_GENO05_HWE0001.vcf.gz | awk '{print $9}' > VCF_files.txt

bcftools merge -a `cat VCF_files.txt` -Oz -o chr1-22.split.onlySTRs.NABEC.uniqueID.DR2_03.Maf01_hwe0001_CR90.FixedIDs.vcf.gz
```
#### NABEC RNASeq data 
```
# Starting file -> nabec.wd_autosomal_qtnorm.csv

# Transpose file so gene names are in the first column and samples on each column

awk -f transposer.awk nabec.wd_autosomal_qtnorm.csv > nabec.wd_autosomal_qtnorm.TR.csv

# With downloaded gencode.v37.annotation.gtf.gz from gencode website and extract chromosome, start & end position and gene ID

zcat gencode.v37.annotation.gtf.gz |\
tail -n+6 |\
awk '$3 == "gene" {print $10 "\t" $1 "\t" $4 "\t" $5}' |\
sed 's/"//g' |\
sed 's/;//g' > gene_start_end_ENSID.v37.txt

# Select gene ID column from genecode file, then annotate it with all gene expression values

cut -f1 gene_start_end_ENSID.v37.txt |\
tail -n2+ > gene_start_end_ENSID.v37.GeneID.txt

awk 'FNR==NR {a[$1]=$0; next} {if ($1 in a) {print a[$1]} else {print $1 "\t" "NA"}}' nabec.wd_autosomal_qtnorm.TR.csv gene_start_end_ENSID.v37.GeneID.txt > nabec.gene_start_end_ENSID.v37.IDwGeneExpression.txt

# complete formatted gene expression file with chr start end geneID samples.

paste nabec.gene_start_end_ENSID.v37.txt nabec.gene_start_end_ENSID.v37.IDwGeneExpression.txt |\
cut -f2 | bgzip -c > nabec.gene_start_end_ENSID.v37.hg38.txt.gz

## Liftover hg38 coordinates to hg19

# Join columns with colons in order to create a large 4th column for the bedfile containing all samples' gene-expression values, so they don't get lost during the liftover process

zcat nabec.gene_start_end_ENSID.v37.hg38.txt.gz |\
sed 's/ //g' |\
sed 's/\t/:/g' > nabec.gene_start_end_ENSID.v37.hg38.colons.bed

zcat nabec.gene_start_end_ENSID.v37.hg38.txt.gz |\
sed 's/ //g' |\
cut -f1-3 |\
paste - nabec.gene_start_end_ENSID.v37.hg38.colons.bed |\
bgzip -c > nabec.gene_start_end_ENSID.v37.hg38.colons.bed.gz

zcat nabec.gene_start_end_ENSID.v37.hg38.colons.bed.gz |\
head -1 > nabec.gene_start_end_ENSID.v37.hg38.colons.bed.header

# Run Liftover
liftOver -bedPlus=4 nabec.gene_start_end_ENSID.v37.hg38.colons.bed.gz hg38ToHg19.over.chain nabec.gene_start_end_ENSID.v37.colons.hg19.bed nabec.gene_start_end_ENSID.v37.hg38.colons.UnMapped.bed


# Format output

sed 's/chr//g' nabec.gene_start_end_ENSID.v37.colons.hg19.bed |\
cat nabec.gene_start_end_ENSID.v37.hg38.colons.bed.header - |\
(sed -u 1q; sort -V -k1,1 -k2,2) | sed 's/:/\t/g' | cut -f1,2,3,7- |\
bgzip -c > nabec.gene_start_end_ENSID.v37.colons.hg19.bed.gz

tabix -p bed nabec.gene_start_end_ENSID.v37.colons.hg19.bed.gz
```
#### PEER factors calculations
```

# file NABEC.genes_individuals.tr.tab -> it's the gene expression file nabec.gene_start_end_ENSID.v37.colons.hg19.bed.gz transposed so individual IDs are in the first column and gene expression values are on each column vertically, like the original .csv file

# file covars.tr.tab > contains individual ID's in the first column and the other columns are SEX, AGE, PC1, PC2, PC3, PC4 and PC5. Data was obtained from "nabec.apr.2018.sample_info.txt" and PCA was calculated from the imputed VCF with "plink --pca" function

peertool -f NABEC.genes_individuals.tr.tab --has_header -c covars.tr.tab -n 45

```

#### Extract the candidate STRs from VCF file and their corresponding chromosomes from the gene expression file. Without the latter, FASTQTL won't work.
```

# Variant extraction from VCF

# file all_nominated_STRs.ld.txt -> all STR candidate variants obtained by LD analysis

bcftools view -i 'ID=@all_nominated_STRs.ld.txt' chr1-22.split.onlySTRs.NABEC.uniqueID.DR2_03.Maf01_hwe0001_CR90.FixedIDs.vcf.gz -Oz -o STR_candidates.split.onlySTRs.NABEC.uniqueID.DR2_03.Maf01_hwe0001_CR90.FixedIDs.vcf.gz

tabix -p vcf STR_candidates.split.onlySTRs.NABEC.uniqueID.DR2_03.Maf01_hwe0001_CR90.FixedIDs.vcf.gz

# Chromosome extraction from expression data

for i in 1 2 3 4 5 7 8 10 12 15 16 17 18 20
do

zcat nabec.gene_start_end_ENSID.v37.colons.hg19.bed.gz |\
awk -v var="$i" '$1 == var {print}' >> nabec.gene_start_end_ENSID.v37.colons.hg19.candidateVars.bed

done

bgzip -c nabec.gene_start_end_ENSID.v37.colons.hg19.candidateVars.bed > nabec.gene_start_end_ENSID.v37.colons.hg19.candidateVars.bed.gz

tabix -p bed nabec.gene_start_end_ENSID.v37.colons.hg19.candidateVars.bed.gz

```
#### Run FastQTL
```
# Enviroment requirements: anaconda3, samtools, R/3.6.3
# Indidivual IDs in VCF file should match with gene expression dataset. VCF file can have more individuals than bed file, but not otherwise.
# Number of chunks should be >= number of chromosomes in data
# Adaptive permutation between 1000 and 10000 is recomended by developer and it also generates a file with FDR corrected p values (q values).

unset PYTHONPATH

python run_FastQTL_threaded.py \
STR_candidates.split.onlySTRs.NABEC.uniqueID.DR2_03.Maf01_hwe0001_CR90.FixedIDs.vcf.gz \
nabec.gene_start_end_ENSID.v37.colons.hg19.candidateVars.bed.gz frontal_cortex.nabec_results \
--covariates Nabec.covars.Sex_Age.PCA_wPEER.txt \
--permute 1000 10000 \
--window 1e6 \
--chunks 14 \
--threads 14 \
frontal_cortex.nabec_results
```
### GTEx v.8 data
#### GTEx v.8 genotypes were downloaded through dbGaP authorized access (dbGaP Accession phs000424.v8.p2)
#### GTEx genotype data liftover (hg38 to hg19), STR imputation and QC.
```
# Raw data QC

vcftools --gzvcf ../dbGaP-25827/77345/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/GenotypeFiles/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz \
--maf 0.01 \
--hwe 0.0001 \
--minDP 20 \
--minGQ 20 \
--max-missing 0.9 \
--recode \
--recode-INFO-all \
--stdout |\
bgzip -c > GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.MAF01_HWE4_DP20_GQ20_CR90.vcf.gz

# Separate file in individual chromosomes

for i in {1..22}
do

vcftools --gzvcf GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.MAF01_HWE4_DP20_GQ20_CR90.vcf.gz \
--chr chr${i} \
--recode \
--recode-INFO-all \
--stdout |\
bgzip -c > chr${i}.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.MAF01_HWE4_DP20_GQ20_CR90.vcf.gz

done

# Liftover VCF from hg38 to hg19 with GATK

for i in {1..22}
do

gatk LiftoverVcf \
--INPUT chr${i}.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.MAF01_HWE4_DP20_GQ20_CR90.vcf.gz \
--OUTPUT chr${i}.hg19.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.MAF01_HWE4_DP20_GQ20_CR90.vcf.gz \
--CHAIN hg38ToHg19.over.chain \
--REJECT chr${i}.rejected_variants.vcf \
--REFERENCE_SEQUENCE hg19.fa

done

# Data harmonization with STR reference panel (using conform-gt.24May16.cee.jar from Beagle)

for i in {1..22}
do

java -jar -Xmx20g conform-gt.24May16.cee.jar \
gt=chr${i}.hg19.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.MAF01_HWE4_DP20_GQ20_CR90.vcf.gz \
ref=1kg.snp.str.chr${i}.vcf.gz \
chrom=${i} \
match=POS \
out=chr${i}.harmonized.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.MAF01_HWE4_DP20_GQ20_CR90.vcf.gz

done 

# STR imputation with Beagle

java -jar -Xmx110g beagle.28Sep18.793.jar \
gt=chr${i}.harmonized.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.MAF01_HWE4_DP20_GQ20_CR90.vcf.gz.vcf.gz \
ref=1kg.snp.str.chr${i}.vcf.gz \
map=plink.chr${i}.GRCh37.map \
nthreads=22 \
window=30 \
overlap=3.0 \
out=GTEx.imputed.snp.str.chr${i}.vcf.gz

# STR extraction and QC (same steps as done with NABEC)

tabix -p vcf GTEx.imputed.snp.str.chr${i}.vcf.gz

# Print variants with "STR" in ID

bcftools view -H GTEx.imputed.snp.str.chr${i}.vcf.gz |\
cut -f3 |\
awk '\$1 ~ /STR/ {print}' > chr${i}.STRs_IDs.txt

# Extract STR variants

bcftools view -i 'ID=@chr${i}.STRs_IDs.txt' -Oz -o GTEx.imputed.STROnly.chr${i}.vcf.gz GTEx.imputed.snp.str.chr${i}.vcf.gz

# Modify header for multiallelic spliting with VT

zcat GTEx.imputed.STROnly.chr${i}.vcf.gz  |\
sed 's/##INFO=<ID=DR2,Number=1/##INFO=<ID=DR2,Number=A/g' |\
bgzip -c > GTEx.imputed.STROnly.chr${i}.DR2A.vcf.gz

vt decompose -s GTEx.imputed.STROnly.chr${i}.DR2A.vcf.gz -o VTSplit.GTEx.imputed.STROnly.chr${i}.DR2A.vcf.gz

# Filter sites with DR < 0.3, recode variant ID to chr_pos_alt and write IDs to a file

bcftools view -e 'INFO/DR2 < 0.3' VTSplit.GTEx.imputed.STROnly.chr${i}.DR2A.vcf.gz -Ou |\
bcftools annotate --set-id '%CHROM\_%POS\_%ALT' |\
bcftools view -H |\
cut -f3 > chr${i}.STRs_DR2_03.RawID.txt

# Second multiallelic splitting now with bcftools on non-split vcf file. Write "raw" variant ID, remove variants with DR2 < 0.3 according to VT multiallelic splitting and filtering. Filter variants with missingness rate > 5%, HWE P < 0.0001 and MAF < 1%

bcftools norm -m-both GTEx.imputed.STROnly.chr${i}.DR2A.vcf.gz -Ou |\
bcftools annotate --set-id '%CHROM\_%POS\_%ALT' -Ou |\
bcftools view -i 'ID=@chr${i}.STRs_DR2_03.RawID.txt' |\
bcftools +fill-tags -Ou |\
bcftools filter -e 'F_MISSING > 0.05' -e 'INFO/HWE < 0.0001' -Ou |\
bcftools view -q 0.01 -Oz -o BCFTSplit.GTEx.imputed.STROnly.chr${i}.DR2A_03.MAF01_GENO05_HWE0001.vcf.gz

#merge into a single file

ls -ltrh BCFTSplit.GTEx.imputed.STROnly.chr*.DR2A_03.MAF01_GENO05_HWE0001.vcf.gz | awk '{print $9}' > VCF_files.txt

bcftools merge -a `cat VCF_files.txt` -Oz -o chr1-22.BCFTSplit.GTEx.imputed.STROnly.DR2A_03.MAF01_GENO05_HWE0001.vcf.gz
```
#### GTEx v.8 RNASeq data processing
##### Fully processed, filtered and normalized gene expression matrices (in BED format) for each tissue, where downloaded from here -> https://gtexportal.org/home/protectedDataAccess (file GTEx_Analysis_v8_eQTL_expression_matrices.tar)
##### Covariates file was obtained from this file -> GTEx_Analysis_v8_eQTL_covariates.tar.gz
#### GTEx v.8 expression matrix liftover from hg38 to hg19 for the 13 brain tissues
```

cat Brain_regions.txt
Amygdala
Anterior_cingulate_cortex_BA24
Caudate_basal_ganglia
Cerebellar_Hemisphere
Cerebellum
Cortex
Frontal_Cortex_BA9
Hippocampus
Hypothalamus
Nucleus_accumbens_basal_ganglia
Putamen_basal_ganglia
Spinal_cord_cervical_c-1
Substantia_nigra

for i in `cat Brain_regions.txt`
do

# format file to keep expression values in 4th column

zcat Brain_${i}.v8.normalized_expression.bed.gz |\
sed 's/\t/:/g' > Brain_${i}.v8.normalized_expression.colons.bed

zcat Brain_${i}.v8.normalized_expression.bed.gz |\
cut -f1-3 |\
paste - Brain_${i}.v8.normalized_expression.colons.bed |\
bgzip -c > Brain_${i}.v8.normalized_expression.wcolons.bed.gz

zcat Brain_${i}.v8.normalized_expression.wcolons.bed.gz |\
head -1 > Brain_${i}.v8.normalized_expression.wcolons.bed.header

# lift over coordinates

liftOver -bedPlus=4 \
Brain_${i}.v8.normalized_expression.wcolons.bed.gz \
hg38ToHg19.over.chain \
Brain_${i}.v8.normalized_expression.wcolons.hg19.bed \
Brain_${i}.v8.normalized_expression.UnMapped.bed

# expanding columns joined by colons to tabs and reformat file for usage with FASTQTL

sed 's/chr//g' Brain_${i}.v8.normalized_expression.wcolons.hg19.bed |\
cat Brain_${i}.v8.normalized_expression.wcolons.bed.header - |\
(sed -u 1q; sort -V -k1,1 -k2,2) |\
sed 's/:/\t/g' |\
cut -f1,2,3,7- |\
bgzip -c > Brain_${i}.v8.normalized_expression.wcolons.hg19.wHeader.final.bed.gz

# Index bed file

tabix -p bed Brain_${i}.v8.normalized_expression.wcolons.hg19.wHeader.final.bed.gz

#remove intermediate files

rm Brain_${i}.v8.normalized_expression.colons.bed
rm Brain_${i}.v8.normalized_expression.wcolons.bed.header
rm Brain_${i}.v8.normalized_expression.wcolons.hg19.bed
rm Brain_${i}.v8.normalized_expression.wcolons.bed.gz

done
```

#### Extract the candidate STRs from VCF file and their corresponding chromosomes from the gene expression file. Without the latter, FASTQTL won't work.

```
bcftools view -i 'ID=@all_nominated_STR.ld' chr1-22.split.onlySTRs.GTEx.uniqueID.DR2_03.Maf01_hwe0001_CR90.vcf.gz.vcf.gz |\
-Oz |\
-o STR_candidates.split.onlySTRs.GTEx.uniqueID.DR2_03.Maf01_hwe0001_CR90.vcf.gz

tabix -p vcf STR_candidates.split.onlySTRs.GTEx.uniqueID.DR2_03.Maf01_hwe0001_CR90.vcf.gz


# Filter expression datasets

for tissue in `cat Brain_regions.txt`
do
for i in 2 3 4 5 7 8 10 12 15 16 17 18 20
do

zcat Brain_${tissue}.v8.normalized_expression.wcolons.bed.gz | awk -v var="$i" '$1 == var {print}' >> Brain_${tissue}.v8.normalized_expression.wcolons.candidateVars.bed

done

bgzip -c Brain_${tissue}.v8.normalized_expression.wcolons.candidateVars.bed.gz > Brain_${tissue}.v8.normalized_expression.wcolons.candidateVars.bed..gz

tabix -p bed Brain_${tissue}.v8.normalized_expression.wcolons.candidateVars.bed..gz


rm Brain_${tissue}.v8.normalized_expression.wcolons.candidateVars.bed

done
```
#### Run FastQTL
```
# Enviroment requirements: anaconda3, samtools, R/3.6.3
# Indidivual IDs in VCF file should match with gene expression dataset. VCF file can have more individuals than bed file, but not otherwise.
# Number of chunks should be >= number of chromosomes in data
# Adaptive permutation between 1000 and 10000 is recomended by developer and it also generates a file with FDR corrected p values (q values).

unset PYTHONPATH

for i in `cat Brain_regions.txt`
do

unset PYTHONPATH

run_FastQTL_threaded.py STR_candidates.split.onlySTRs.GTEx.uniqueID.DR2_03.Maf01_hwe0001_CR90.vcf.gz \
Brain_${i}.v8.normalized_expression.wcolons.candidateVars.bed.gz \
NormalizedGTEx_eSTR_${i}_results_CandidateVars_PERMUTE_PASS \
--covariates Brain_${i}.v8.covariates.txt \
--permute 1000 10000 \
--window 1e6 \
--chunks 13 \
--threads 13 \
--ma_sample_threshold 10 \
--maf_threshold 0.01 \
-o NormalizedGTEx_${i}_results_CandidateVars_PERMUTE_PASS
done
```
#### Significant gene selection, qq-plot and box-plots
```
#For NABEC results

zcat frontal_cortex.nabec_results/frontal_cortex.nabec_results.genes.txt.gz | head -1 | awk '{print $1 "\t" $7 "\t" $8 "\t" $11 "\t" $17 "\t" $18 "\t" $14 "\t " $15}' > frontal_cortex.nabec_results.genes.SignificantResults.txt

zcat frontal_cortex.nabec_results/frontal_cortex.nabec_results.genes.txt.gz | awk '$(NF-1) < 0.05 {print $1 "\t" $7 "\t" $8 "\t" $11 "\t" $17 "\t" $18 "\t" $14 "\t " $15}' >> frontal_cortex.nabec_results.genes.SignificantResults.txt

# For GTEx results
for i in `cat Brain_regions.txt`
do
zcat NormalizedGTEx_eSTR_${i}_results_CandidateVars_PERMUTE_PASS.genes.txt.gz |\
head -1 |\
awk '{print $1 "\t" $7 "\t" $8 "\t" $11 "\t" $17 "\t" $18 "\t" $14 "\t " $15}' > NormalizedGTEx_eSTR_${i}_results_CandidateVars_PERMUTE_PASS.genes.SignificantResults.txt

zcat NormalizedGTEx_eSTR_${i}_results_CandidateVars_PERMUTE_PASS.genes.txt.gz |\
awk '$(NF-1) < 0.05 {print $1 "\t" $7 "\t" $8 "\t" $11 "\t" $17 "\t" $18 "\t" $14 "\t " $15}' >> NormalizedGTEx_eSTR_${i}_results_CandidateVars_PERMUTE_PASS.genes.SignificantResults.txt

done

# QQ-plot in R

# load each data
data1 = read.table("frontal_cortex.nabec_results_CandidateVars_PERMUTATION_PASS_newGeneIDs.genes.txt.gz", hea=T, stringsAsFactors=F)
data2 = read.table("NormalizedGTEx_eSTR_Amygdala_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data3 = read.table("NormalizedGTEx_eSTR_Anterior_cingulate_cortex_BA24_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data4 = read.table("NormalizedGTEx_eSTR_Caudate_basal_ganglia_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data5 = read.table("NormalizedGTEx_eSTR_Cerebellar_Hemisphere_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data6 = read.table("NormalizedGTEx_eSTR_Cerebellum_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data7 = read.table("NormalizedGTEx_eSTR_Cortex_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data8 = read.table("NormalizedGTEx_eSTR_Frontal_Cortex_BA9_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data9 = read.table("NormalizedGTEx_eSTR_Hippocampus_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data10 = read.table("NormalizedGTEx_eSTR_Hypothalamus_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data11 = read.table("NormalizedGTEx_eSTR_Nucleus_accumbens_basal_ganglia_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data12 = read.table("NormalizedGTEx_eSTR_Putamen_basal_ganglia_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data13 = read.table("NormalizedGTEx_eSTR_Spinal_cord_cervical_c-1_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)
data14 = read.table("NormalizedGTEx_eSTR_Substantia_nigra_results_CandidateVars_PERMUTE_PASS.genes.txt.gz", hea=T, stringsAsFactors=F)

# sort p values for each data

p1.sorted <- sort(data1$pval_beta)
p2.sorted <- sort(data2$pval_beta)
p3.sorted <- sort(data3$pval_beta)
p4.sorted <- sort(data4$pval_beta)
p5.sorted <- sort(data5$pval_beta)
p6.sorted <- sort(data6$pval_beta)
p7.sorted <- sort(data7$pval_beta)
p8.sorted <- sort(data8$pval_beta)
p9.sorted <- sort(data9$pval_beta)
p10.sorted <- sort(data10$pval_beta)
p11.sorted <- sort(data11$pval_beta)
p12.sorted <- sort(data12$pval_beta)
p13.sorted <- sort(data13$pval_beta)
p14.sorted <- sort(data14$pval_beta)

# generate color ramp

fc <- colorRampPalette(c("lightblue", "darkblue"))

# qq plot generation

pdf(file="qqplot.pdf")
plot( -log10(ppoints(p11.sorted )), -log10(p11.sorted) , col="#00008B", pch = 19, cex= 0.8, xlim=c(0, 3.5), ylim=c(0,70))
points( -log10(ppoints(p4.sorted )), -log10(p4.sorted) , col="#3C008B", cex= 0.8,pch = 19)
points( -log10(ppoints(p6.sorted )), -log10(p6.sorted) , col="#77008B", cex= 0.8, pch = 19)
points( -log10(ppoints(p7.sorted )), -log10(p7.sorted) , col="#8B0063", cex= 0.8, pch = 19)
points( -log10(ppoints(p1.sorted )), -log10(p1.sorted) , col="#8B0028", cex= 0.8, pch = 19)
points( -log10(ppoints(p8.sorted )), -log10(p8.sorted) , col="#8B1400", cex= 0.8, pch = 19)
points( -log10(ppoints(p10.sorted )), -log10(p10.sorted) , col="#8B4F00", cex= 0.8, pch = 19)
points( -log10(ppoints(p5.sorted )), -log10(p5.sorted) , col="#8B8B00", cex= 0.8, pch = 19)
points( -log10(ppoints(p9.sorted )), -log10(p9.sorted) , col="#4F8B00", cex= 0.8, pch = 19)
points( -log10(ppoints(p12.sorted )), -log10(p12.sorted) , col="#148B00", cex= 0.8, pch = 19)
points( -log10(ppoints(p3.sorted )), -log10(p3.sorted) , col="#008B28", cex= 0.8, pch = 19)
points( -log10(ppoints(p13.sorted )), -log10(p13.sorted) , col="#008B63", cex= 0.8, pch = 19)
points( -log10(ppoints(p2.sorted )), -log10(p2.sorted) , col="#00778B", cex= 0.8, pch = 19)
points( -log10(ppoints(p14.sorted )), -log10(p14.sorted) , col="#003C8B", cex= 0.8, pch = 19)

# add lines for better resolution of points trayectory

lines( -log10(ppoints(p11.sorted )), -log10(p11.sorted) , col=alpha(#00008B", 0.5) col=alpha(rgb(0,0,0), 0.5)"#00008B", cex =0.8)
lines( -log10(ppoints(p4.sorted )), -log10(p4.sorted) , col="#3C008B", cex= 0.8)
lines( -log10(ppoints(p6.sorted )), -log10(p6.sorted) , col="#77008B", cex= 0.8)
lines( -log10(ppoints(p7.sorted )), -log10(p7.sorted) , col="#8B0063", cex= 0.8)
lines( -log10(ppoints(p1.sorted )), -log10(p1.sorted) , col="#8B0028", cex= 0.8)
lines( -log10(ppoints(p8.sorted )), -log10(p8.sorted) , col="#8B1400", cex= 0.8)
lines( -log10(ppoints(p10.sorted )), -log10(p10.sorted) , col="#8B4F00", cex= 0.8)
lines( -log10(ppoints(p5.sorted )), -log10(p5.sorted) , col="#8B8B00", cex= 0.8)
lines( -log10(ppoints(p9.sorted )), -log10(p9.sorted) , col="#4F8B00", cex= 0.8)
lines( -log10(ppoints(p12.sorted )), -log10(p12.sorted) , col="#148B00", cex= 0.8)
lines( -log10(ppoints(p3.sorted )), -log10(p3.sorted) , col="#008B28", cex= 0.8)
lines( -log10(ppoints(p13.sorted )), -log10(p13.sorted) , col="#008B63", cex= 0.8)
lines( -log10(ppoints(p2.sorted )), -log10(p2.sorted) , col="#00778B", cex= 0.8)
lines( -log10(ppoints(p14.sorted )), -log10(p14.sorted) , col="#003C8B", cex= 0.8)

# Box plots for chr17 STR in GTEx data
# First, genotypes are extracted from VCF file, then gene-expression data for the specific gene is extracted from bed file. Format look like this:

head GTEx_IIDs.genotypes.expression_chr17Var.txt
GTEx_genotypes	Nucleus_accumbens_basal_ganglia_ENSG00000187688	Hypothalamus_ENSG00000187688	Frontal_Cortex_BA9_ENSG00000187688	Anterior_cingulate_cortex_BA24_ENSG00000187688	Hippocampus_ENSG00000141027	Anterior_cingulate_cortex_BA24_ENSG00000275413	Amygdala_ENSG00000275413	Substantia_nigra_ENSG00000275413	
BB	-1.089661981946107	-0.18427060915128232	0.7291513430522927	0.9319713123431904	0.10589879442185142	-1.9264031529639818	-0.6151411045959735	-1.5475149163347155
AB	-2.3319288684462447	-0.33603814037182317	-1.904706898156732	-1.6716437380157467	-0.04531601592937276	-1.4941549086209047	0.2533471031357997	-0.43871291859610834
AA	0.26867818583276143	1.739929801420972	0.4727891209922672	0.39995564109547443	-1.0338531513302158	0.1701847241123691	-1.4815444904261041	-1.7116753065097285
BB	-1.0248823169739363	-0.4795056533309502	0.5540487914423362	-0.1359106806464805	-2.2555889554631765	-0.01693748732987151	-0.4377900815764892	-0.6141284723406751
AA	-0.7338269926088062	0.9328891732059703	1.3351777361189363	1.4941549086209052	-1.2031567637736242	-0.5321897309193041	0.019282950895712157	-1.0088564614556432
BB	-0.10515079053017912	0.43072729929545744	0.24453282622639816	1.2738891541514752	-0.6278221076294491	0.9584461410107845	1.23889437958136	0.7517957253750428
AB	-0.8882910058929048	1.1903263992957545	0.6744897501960817	0.11883585040329893	0.5735607771007943	0.6533766109715351	0.13538473551751684	1.0458026143123607
AA	-0.7500875614194734	0.3672262071417851	-1.335177736118937	-1.2366522415199066	-0.5558502170512656	2.2111272410853284	0.7112600919034947	1.7116753065097288
AB	0.9255501634629518	0.8885053006369341	0.44116864743261897	1.8278799714620912	-0.4197079794617776	0.3097425315385302	-0.3135719665290585	0.8109553641377003

# R code for boxplot

library(reshape2)
library(ggplot2)

data <- read.table("GTEx_IIDs.genotypes.expression_chr17Var.txt", header = TRUE)

DF <- melt(data, id.var = "GTEx_genotypes")  

pdf(file="GTEx_Chr17Var.boxplots.pdf", width=9, height=3)
ggplot (data= DF, aes(x=variable, y=value)) + geom_boxplot(aes(fill=GTEx_genotypes)) + theme_classic() + scale_fill_grey()
dev.off()
```
## 5. GCTA-LDMS heritability calculation for STRs, SNPs+STRs and SNPs
### Method taken from here -> https://cnsgenomics.com/software/gcta/#GREMLinWGSorimputeddata
```
# File STRs_fromMeta.txt  -> all STRs from meta-analysis
# File SNPs_fromMeta.txt  -> all SNPs from meta-analysis
# File SNPS_STRs_fromMeta.txt -> all SNPs and STRs from meta-analysis

# The following steps were performed identically for each file above. The code below is for the STRs only analysis.

# Run LDScore regression

gcta64 --bfile chr1-22.split.SNPs_STRs.PDGWAS.uniqueID.DR2_03.merged.MAF01_GENO05_HWE0001 --extract STRs_fromMeta.txt --ld-score-region 200 --out str_only --thread-num 27

# Stratify snps in 4 quantiles in R

lds_seg = read.table("str_only.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, "str_only_group1.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "str_only_group2.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "str_only_group3.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "str_only_group4.txt", row.names=F, quote=F, col.names=F)

# Run GRM for each group

for i in {1..4}
do

gcta64 --bfile chr1-22.split.SNPs_STRs.PDGWAS.uniqueID.DR2_03.merged --extract str_only_group${i}.txt --make-grm --out str_only_group${i} --thread-num 10

done

# Run REML analysis, assuming a global prevalence of PD of 0.5% and using covariates such as sex and PC1-10

gcta64 --reml --mgrm STRs_only_multi_GRMs.txt --pheno UK.MF.US.GER.FRcas_strand.lifted_remake_RawID.QC_10072019.modifiedIDs.pheno --prevalence 0.005 --covar sex.covar --qcovar PC1_5.modifiedIDs.covar --out STRs_only_sig_h2 --thread-num 10
```
## 6. MAGMA gene-wise, gene-set and pathway enrichment analyses
#### Gene-wise STR enrichment analysis
```
# Annotate variants to genes using a window of 35 kb upstream and 10 kb downstream, in order to include regulatory regions (PMID:32152537)
# Since FUMA uses magma along with the 1000 genomes project reference panel, we used here the same STR dataset used for imputation as reference panel.

magma --annotate window=35,10 --snp-loc newID.split.bcftools.STR_only.1kg.snp.str.chr1-22.bim --gene-loc NCBI37.3.gene.loc --out vars_annotated_35_10kb.txt

# Generate MAGMA raw file
# File STR_I80.ID_P_N.txt -> contains all STRs (I2 < 80%) from meta-analysis with ID and p value

magma --bfile newID.splitted.bcftools.STR_only.1kg.snp.str.chr1-22 --pval STR_I80.ID_P.txt N=39087 --gene-annot vars_annotated_35-10kb.txt.genes.annot --out MAGMA_STR_1kg3_35-10kb.TotalN
```
#### Gene-set and pathway STR enrichment analysis
```
# Run gene-set enrichment using gene-sets and pathways from the Molecular Signature database (MSigDB)
# File C2_all_C5_BP_MF_CC.symbols.formatted.gmt -> contains curated pathways from the C2 collection and C5 gene-ontology collection from MSigDB.

magma --gene-results MAGMA_STR_1kg3_35-10kb.TotalN.genes.raw  --set-annot C2_all_C5_BP_MF_CC.symbols.formatted.gmt --out MAGMA_STR_1kg3_35-10kb.TotalN.genes.C2_all_C5_BP_MF_CC.enrichment

# Run gene-property analysis with overexpressed genes from GTEx general and specific tissues

# File gtex_v8_ts_general_avg_log2TPM.wGeneNames.unique.txt was obtained from FUMA developers directly by email. Details about how the data was generated is here -> https://fuma.ctglab.nl/tutorial#magma
# Tutorial explain data generation for GTEx v6 and v7 but they sent me the v8 analysis as I requested it.

# For 30 general tissues

magma --gene-results MAGMA_STR_1kg3_35-10kb.genes.raw --gene-covar gtex_v8_ts_general_avg_log2TPM.wGeneNames.unique.txt --model direction=greater condition=Average --out MAGMA_STR_1kg3_35-10kb.GTEx_General

# For 54 specific tissues

magma --gene-results MAGMA_STR_1kg3_35-10kb.genes.raw --gene-covar gtex_v8_ts_avg_log2TPM.wGeneNames.unique.txt --model direction=greater condition=Average --out MAGMA_STR_1kg3_35-10kb.genes.GTEx_54Tissues
```
#### For the Dropviz 88 single-cell RNASeq expression enrichment analysis, we submitted the MAGMA_STR_1kg3_35-10kb.genes.raw file to the FUMA online server (https://fuma.ctglab.nl/celltype)