#!/bin/bash

ver=$1  # ver4
region=$2 # start31.7end32.2_highLDfiltering_50_5_0.8QUAL20
haplonetver=$3 # ver1_sujung
crossval_N=$4 #1

path="/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/$ver/${region}_${haplonetver}"

cd "$path"

# 1. sampleQC (not 2copy, not A+B=C4, not C4>HERV after phasing)
vcftools --gzvcf whole.eagle.phased.vcf.gz --remove rmsamples.txt --recode --recode-INFO-all --stdout | bgzip > whole.eagle.phased.sampleQC.vcf.gz
tabix -p vcf whole.eagle.phased.sampleQC.vcf.gz

# AC field add for IMPUTE5 running
bcftools +fill-tags whole.eagle.phased.sampleQC.vcf.gz -Oz -o whole.eagle.phased.sampleQC.AC.vcf.gz -- -t AN,AC
tabix -p vcf whole.eagle.phased.sampleQC.AC.vcf.gz


bcftools +fill-tags whole.eagle.phased.vcf.gz -Oz -o whole.eagle.phased.AC.vcf.gz -- -t AN,AC
tabix -p vcf whole.eagle.phased.AC.vcf.gz

# 2. subsampling reference panel for leaveoneout crossvalidation
# + 3. masking(removing) C4 + HLA region of test samples for imputation
# 생각해보니 ref 마다 FRQ.frq 파일을만들어야...? 아니면 그냥 beagle5로 독자적으로?

mkdir "crossval$crossval_N"
LIST="$path/survived_sams.txt"

# i=1
# END=$(wc -l < "$LIST")


LIST="$path/survived_sams.txt"

cat "$LIST" | xargs -I {} -P 100 bash -c '
    FILE=$(echo "{}" | tr -d "\n")
    echo "Now [{}]-th processing"
    vcftools --gzvcf whole.eagle.phased.sampleQC.AC.vcf.gz --remove-indv "$FILE" --recode --recode-INFO-all --stdout | bgzip > "./crossval'"$crossval_N"'/ref_panel{}.vcf.gz"
    vcftools --gzvcf whole.eagle.phased.sampleQC.AC.vcf.gz --indv "$FILE" --recode --recode-INFO-all --stdout | gzip -c > "./crossval'"$crossval_N"'/test_panel{}.phased.vcf.gz"
    vcftools --gzvcf whole.eagle.phased.sampleQC.AC.vcf.gz --indv "$FILE" --exclude-positions masked_alleles.txt --recode --recode-INFO-all --stdout | bgzip > "./crossval'"$crossval_N"'/test_panel{}.vcf.gz"
    tabix -p vcf "./crossval'"$crossval_N"'/ref_panel{}.vcf.gz"
    tabix -p vcf "./crossval'"$crossval_N"'/test_panel{}.vcf.gz"
'


gzip -d ./crossval$crossval_N/*.phased.vcf.gz