
# KTN1 GWAS and whole genome analysis
## 1. GWAS (IPDGC genotype data)
#### 1. Subset PLINK Binaries + Convert to VCFs
#### 2. Output Frequency with PLINK
#### 3. Logistic Regression with PLINK
#### 4. Annotate VCFs with Annovar 

    #!/bin/bash

    #hg19
    
    ##set up environment
    module load plink samtools rvtests R annovar

    mkdir $KTN1_DIR/
    cd $KTN1_DIR/

    mkdir hardcallsNoNeuroX hardcallsNoNeuroX/bin 
    mkdir hardcallsNoNeuroX/freq hardcallsNoNeuroX/score 
    mkdir hardcallsNoNeuroX/burden hardcallsNoNeuroX/vcf 
    mkdir hardcallsNoNeuroX/annotation hardcallsNoNeuroX/burden 
    mkdir hardcallsNoNeuroX/logistic hardcallsNoNeuroX/freq

    # ==================== 1. Subset PLINK Binaries to cis-KTN1 region ====================
    ## Remove NeuroX + keep males + genotype quality of 15% or less missingness
    
    ##KTN1 region 56025790 - 56168244 +/- MB
    ##hg19 coordinates

    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --chr 14 --geno 0.15 --from-bp 55025790 --to-bp 57168244 \
    --make-bed --out hardcallsNoNeuroX/bin/KTN1.GWAS.MB
        # Done 11.01.2020
        # PLINK output:
            # Total genotyping rate in remaining samples is 0.461474.
            # 15399 variants removed due to missing genotype data (--geno).
            # 6188 variants and 32338 people pass filters and QC.
            # Among remaining phenotypes, 14671 are cases and 17667 are controls.


    plink --bfile hardcallsNoNeuroX/bin/KTN1.GWAS.MB --recode vcf --out hardcallsNoNeuroX/vcf/KTN1.GWAS.MB

    cd hardcallsNoNeuroX/vcf
    bgzip KTN1.GWAS.MB.vcf
    tabix -f -p vcf KTN1.GWAS.MB.vcf.gz

    # ==================== 2. Output Frequency via PLINK ==================== 
    cd $KTN1_DIR/
    # Overall
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --chr 14 --from-bp 55025790 --to-bp 57168244 \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --geno 0.15 --freq \
    --out hardcallsNoNeuroX/freq/KTN1.MB
        # Done 11.01.2021


    # Case-control
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --chr 14 --from-bp 55025790 --to-bp 57168244 \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --geno 0.15 --freq case-control  \
    --out hardcallsNoNeuroX/freq/KTN1.MB
        # Done 11.01.2021

    # ==================== 3. Output Logistic Regression via PLINK ====================
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --chr 14 --from-bp 55025790 --to-bp 57168244 \
    --geno \
    --covar $FILE_DIR/IPDGC_all_samples_covariates.tab \
    --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE \
    --assoc \
    --out hardcallsNoNeuroX/logistic/KTN1
        # Done 12.01.2020

    # ==================== 4. Annotate VCFs via Annovar ====================

    # Downloading databases 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
    # annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome humandb/

    table_annovar.pl $KTN1_DIR/hardcallsNoNeuroX/vcf/KTN1.GWAS.MB.vcf.gz $FILE_DIR/annovar/humandb/ -buildver hg19 \
    --thread 16 \
    -out $KTN1_DIR/hardcallsNoNeuroX/annotation/KTN1.MB.annovar \
    -remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
    -operation f,g,g,f \
    -nastring . \
    -vcfinput

    cd $KTN1_DIR/hardcallsNoNeuroX/annotation
    head -1 KTN1.MB.annovar.hg19_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1.MB.annovar.hg19_multianno.txt > KTN1.MB.trimmed.annotation.txt
        # Done 12.01.2020


## 2. Whole Genome Analysis (AMP-PD)
#### 0. Prep file
#### 1. Annotate with Annovar
#### 2. Burden analysis with RVTESTS
 
    #!/bin/bash
    
    ##set up environment
    module load plink samtools rvtests annovar

    ## setting up directories
    cd $KTN1_DIR/

    mkdir AMPPD

    ## ======================= 0. Prep file  ========================= 

    plink --bfile /data/CARD/PD/AMP-PD/AMPv2.5_sampleQC_EURO_no_cousins_SEXPHENO --chr 14 --from-bp 55580207 --to-bp 55684584 --recode vcf-fid --out /data/CARD/PD/AMP-PD/PETER_KTN1/NCBI/pheno_KTN1
	plink --bfile /data/CARD/PD/AMP-PD/AMPv2.5_sampleQC_EURO_no_cousins_SEXPHENO --chr 14 --from-bp 55580207 --to-bp 55684584 --make-bed --out /data/CARD/PD/AMP-PD/PETER_KTN1/NCBI/pheno_KTN1
	plink --bfile /data/CARD/PD/AMP-PD/AMPv2.5_sampleQC_EURO_no_cousins_SEXPHENO --chr 14 --from-bp 55559072 --to-bp 55701526 --recode vcf-fid --out /data/CARD/PD/AMP-PD/PETER_KTN1/ENSEMBL/pheno_KTN1
	plink --bfile /data/CARD/PD/AMP-PD/AMPv2.5_sampleQC_EURO_no_cousins_SEXPHENO --chr 14 --from-bp 55559072 --to-bp 55701526 --make-bed --out /data/CARD/PD/AMP-PD/PETER_KTN1/ENSEMBL/pheno_KTN1


	## Format files
	bgzip pheno_KTN1.vcf
	tabix -f -p vcf pheno_KTN1.vcf.gz

	## Extract variants from Mao et al.
	plink --bfile AMPv2.5_sampleQC_EURO_no_cousins_SEXPHENO --chr 14 --extract range Mao.txt --make-bed --out pheno_KTN1_Mao_vars
	plink --bfile pheno_KTN1_Mao_vars --assoc --out pheno_KTN1_Mao_vars

    ## =======================  1. Annotation via Annovar ========================= 

    
    # vim job.sh and paste the following:
	# vim job.sh and paste the following:
	#!/bin/bash
	table_annovar.pl pheno_KTN1.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
	--thread 16 \
	--out KTN1.annovar \
	-arg '-splicing 15',,, \
	-remove \
	-protocol refGene,avsnp147,ljb26_all,gnomad312_genome \
	-operation g,f,f,f \
	-nastring . \
	--vcfinput
	## sbatch --cpus-per-task=24 --mem=64g --mail-type=BEGIN,TIME_LIMIT_50,END --time=02:00:00 job.sh

	## Remove header of annotated file and remove carriers columns for simplicity
	cut -f1-52 KTN1.annovar.hg38_multianno.txt > KTN1.annovar.trimmed.annotation.txt

	## Extract variants according to variant category
	## Coding
	awk '$6=="exonic" {print}' KTN1.annovar.trimmed.annotation.txt > KTN1.annovar.trimmed.annotation.coding.txt
	awk '{print $1" "$2" "$2" "$7}' KTN1.annovar.trimmed.annotation.coding.txt > KTN1.annovar.trimmed.annotation.coding.SNPs.txt

	## Coding - non-synonymous
	awk '$9=="nonsynonymous" {print}' KTN1.annovar.trimmed.annotation.txt > KTN1.annovar.trimmed.annotation.nonsynonymous.txt

	## Splicing
	awk '$6=="splicing" {print}' KTN1.annovar.trimmed.annotation.txt > KTN1.annovar.trimmed.annotation.splicing.txt
	awk '{print $1" "$2" "$2" "$7}' KTN1.annovar.trimmed.annotation.splicing.txt > KTN1.annovar.trimmed.annotation.splicing.SNPs.txt

	## Extract variant categories
	plink --vcf pheno_KTN1.vcf.gz --extract range KTN1.annovar.trimmed.annotation.coding.SNPs.txt --recode vcf --out pheno_KTN1_coding
	plink --vcf pheno_KTN1.vcf.gz --extract range KTN1.annovar.trimmed.annotation.nonsynonymous.txt --recode vcf --out pheno_KTN1_coding_nonsyn
	plink --vcf pheno_KTN1.vcf.gz  --extract range KTN1.annovar.trimmed.annotation.splicing.SNPs.txt --recode vcf --out pheno_KTN1_splicing

	bgzip pheno_KTN1_coding.vcf
	tabix -f -p vcf pheno_KTN1_coding.vcf.gz

	bgzip pheno_KTN1_coding_nonsyn.vcf
	tabix -f -p vcf pheno_KTN1_coding_nonsyn.vcf.gz

	bgzip pheno_KTN1_splicing.vcf
	tabix -f -p vcf pheno_KTN1_splicing.vcf.gz  

    ## =======================  2. Burden analysis on all variants via RVTESTS  ========================= 
 
   	## Burden all the variants
	rvtest --noweb --inVcf pheno_KTN1.vcf.gz --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --covar /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --out AMP_PD_BURDEN.KTN1

	## Burdens by frequency bins
	rvtest --noweb --inVcf pheno_KTN1.vcf.gz --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --covar /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --out AMP_PD_BURDEN.KTN1_MAX_MAF_0.05_0.03 --freqUpper 0.05 --freqLower 0.03
	rvtest --noweb --inVcf pheno_KTN1.vcf.gz --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --covar /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --out AMP_PD_BURDEN.KTN1_MAX_MAF_0.03_0.01 --freqUpper 0.03 --freqLower 0.01
	rvtest --noweb --inVcf pheno_KTN1.vcf.gz --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --covar /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --out AMP_PD_BURDEN.KTN1_MAX_MAF_0.01 --freqUpper 0.01

	## Burdens by variant category
	rvtest --noweb --inVcf pheno_KTN1_coding.vcf.gz --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --covar /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --out AMP_PD_BURDEN.KTN1_CODING
	rvtest --noweb --inVcf pheno_KTN1_coding_nonsyn.vcf.gz --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --covar /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --out AMP_PD_BURDEN.KTN1_CODING_NONSYN
	rvtest --noweb --inVcf pheno_KTN1_splicing.vcf.gz --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --covar /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --out AMP_PD_BURDEN.KTN1_SPLICING

	######## ASSOCIATION ADJUSTED BY COVARIATES ########
	plink --bfile pheno_KTN1 --logistic --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --covar /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --out pheno_KTN1
	plink --bfile pheno_KTN1 --assoc --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --out pheno_KTN1
	awk '{print $1, $3, $5, $6}' pheno_KTN1.assoc > frequencies.txt

	######## FREQUENCY ########
	plink --bfile pheno_KTN1 --freq --pheno /data/CARD/PD/AMP-PD/PETER_KTN1/EUR_encoded_AMP_covs_SEPT2021_allPCs_PD_PHENO.txt --pheno-name PD_PHENO --out pheno_KTN1
	awk '$5<="0.05" {print}' pheno_KTN1.frq > RARE

	######## PREPARE ST2 ########
	%%bash
	awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $50}' KTN1.annovar.trimmed.annotation.txt > ANNOTATIONS.txt

	R
	library(data.table)
	LOGISTIC <- fread("pheno_KTN1.assoc.logistic", header =T)
	LOGISTIC$Start <- LOGISTIC$BP
	ANNOT <- fread("ANNOTATIONS.txt", header =T)
	temp <- merge(LOGISTIC, ANNOT, by="Start")
	FREQUENCIES <- fread("frequencies.txt", header =T)
	FREQUENCIES$Start <- FREQUENCIES$BP
	merged <- merge(FREQUENCIES, temp, by="Start") 
	subsetted <- subset(merged, select=c("CHR.x", "SNP", "BP.x", "Ref", "Alt", "Func.refGene", "Gene.refGene","ExonicFunc.refGene", "AAChange.refGene", "OR", "P", "gnomad312_AF_nfe", "F_A", "F_U", "TEST"))
	fwrite(subsetted, file = "temp.txt", na = "NA", quote = F, row.names = F, sep = "\t")

	%%bash
	head -1 temp.txt > header.txt
	grep "ADD" temp.txt > final_pheno_KTN1.assoc.logistic
	cat header.txt final_pheno_KTN1.assoc.logistic > ST2.txt

# Summary Statistics 
### 1. GWAS Summary Statistics (Nalls 2019)
    # ==================== 1. Extract gene from GWAS sum stats + Zoom plots ====================
    
    #bash
    # KTN1 positions on hg19 (14: 56025790 - 56168244) 
    # cis-KTN1 (14: 55025790 - 57168244) --> KTN1.MB.txt

    ##making from $KTN1_DIR/annotation/KTN1.trimmed.annotation/txt
    #awk '{ print $1,$2,$3,$8 }' $KTN1_DIR/hardcallsNoNeuroX/annotation/KTN1.MB.trimmed.annotation.txt > $KTN1_DIR/sum_stats/KTN1.MB.txt

    R
    ##run in biowulf
    ## Nalls 2019
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(plyr)
    #data <- fread("Nalls2019_no23andMe.tbl", header = T)
    data <- fread("$FILE_DIR/nallsEtAl2019_excluding23andMe_allVariants.tab", header = T)
    data$CHR <- ldply(strsplit(as.character(data$SNP), split = ":"))[[1]]
    data$CHR <- sub("chr", "", data$CHR)
    data$BP <- ldply(strsplit(as.character(data$SNP), split = ":"))[[2]]
    genesTemp <- read.table("KTN1.MB.txt", sep = " ", header = T)
    colnames(genesTemp) <- c("CHR","START","STOP","GENE")
    genes <- subset(genesTemp, CHR != "X")
    for(i in 1:length(genes$GENE))
    {
      thisGene <- genes$GENE[i]
      thisChr <- genes$CHR[i]
      lower <- genes$START[i]
      upper <- genes$STOP[i]
      output <- subset(data, CHR == thisChr & BP >= lower & BP <= upper)
      file_name <- paste(thisGene,"_risk_variants_mb.tab", sep = "")
      fwrite(output, file = file_name, na = "NA", quote = F, row.names = F, sep = "\t")
    }
    
    ######## Zoom plot ##########
    module load locuszoom
    locuszoom --metal $KTN1_DIR/sum_stats/cisKTN1_variants.Nalls.tab --pvalcol p --markercol SNP --chr 14 --start 55025790 --end 57168244  --pop EUR --build hg19 --source 1000G_March2012 --plotonly &


### 2. Age of Onset GWAS Summary Statistics (Blauwendraat 2019)
    
    #R
    ## Blauwendraat 2019
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(plyr)
    #data <- fread("Blauwendraat2019_no23andMe.tbl", header = T)
    data <- fread("IPDGC_AAO_GWAS_sumstats_april_2018.txt", header = T)
    data$CHR <- ldply(strsplit(as.character(data$MarkerName), split = ":"))[[1]]
    data$BP <- ldply(strsplit(as.character(data$MarkerName), split = ":"))[[2]]
    chr14 <-subset(data, CHR == 'chr14')

    genesTemp <- read.table("KTN1.MB.txt", sep = "", header = F)
    colnames(genesTemp) <- c("CHR","START","STOP","GENE")
    genes <- subset(genesTemp, CHR != "X")

    for(i in 1:length(genes$GENE))
    {
      thisGene <- genes$GENE[i]
      thisChr <- genes$CHR[i]
      lower <- genes$START[i]
      upper <- genes$STOP[i]
      output <- subset(chr14, CHR == thisChr & BP >= lower & BP <= upper)
      file_name <- paste(thisGene,"_aao_variants_mb.tab", sep = "")
      fwrite(output, file = file_name, na = "NA", quote = F, row.names = F, sep = "\t")
    }
    
    ######### Zoom plot ###########
   
    locuszoom --metal $KTN1_DIR/sum_stats/cisKTN1_variants.Blauwendraat.tab --pvalcol P-value --markercol SNP --chr 14 --start 55025790 --end 57168244 --flank 1000kb --pop EUR --build hg19 --source 1000G_March2012 --plotonly &

# KTN1 expression analysis
## 1. eQTL (GTEx v8) & BRAINEAC
## 2. Expression covariate correction + T-test (AMP-PD)
    all_ktn1_expression_notebook_github.ipynb
    all_ktn1_expression_notebook.ipynb
	all_ktn1_wgs_sqc.ipynb
	BL_all_idio_ktn1_expression_notebook.ipynb
	BL_all_ktn1_expression_notebook.ipynb
	BL_all_PCs_ktn1_expression_notebook.ipynb
	BL_pdbp_ktn1_expression_notebook.ipynb
	BL_ppmi_ktn1_expression_notebook.ipynb
	ktn1_expression_notebook.ipynb
	m24_all_PCs_ktn1_expression_notebook.ipynb
    
