#!/bin/bash

ver=$1 # ver3 or ver4

# combine kih and unist vcf (unphased ? ) for phasing 
vcfpath="/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm/0.vcf_to_plink/$ver"
cd /kimlab_wd/yuo1996/C4_HLA_ref/haplonet_proximal_with_Eagle/combine_kih/1.mk_vcf_for_phasing/$ver
gzip -dc $vcfpath/newsnps.vcf.gz  >  newsnps.vcf
gzip -dc $vcfpath/newsnps.kih.vcf.gz  >  newsnps.kih.vcf
bgzip ./newsnps.vcf
bgzip ./newsnps.kih.vcf
bcftools sort newsnps.vcf.gz -Oz -o newsnps.sorted.vcf.gz
bcftools sort newsnps.kih.vcf.gz -Oz -o newsnps.kih.sorted.vcf.gz

tabix -p vcf newsnps.sorted.vcf.gz 
tabix -p vcf newsnps.kih.sorted.vcf.gz 

bcftools merge newsnps.sorted.vcf.gz newsnps.kih.sorted.vcf.gz -Ov -o merged.vcf 
gzip merged.vcf



#start !! /haplonet_proximal_with_Eagle/combine_kih/code/0.vcf_to_plink.sh 에서 위에 코드 돌린 것 그냥 가져옴 
cd /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver}
cp /kimlab_wd/yuo1996/C4_HLA_ref/haplonet_proximal_with_Eagle/combine_kih/1.mk_vcf_for_phasing/${ver}/merged.vcf.gz ./1stQC_with_vcftools.vcf.gz


# convert vcf.gz to plink 
plink --vcf 1stQC_with_vcftools.vcf.gz \
        --double-id \
        --make-bed --out 1st_plink


# plink basic qcs
plink --bfile 1st_plink \
        --list-duplicate-vars suppress-first --out dup  # 없음걍.

# plink --bfile 1st_plink \
#         --geno 0.05 \
#         --maf 0.005 \
#         --hwe 1e-7 \
#         --make-bed \
#         --out 2nd_plinkQC

plink --bfile 1st_plink \
        --geno 0.1 \
        --maf 0.005 \
        --hwe 1e-7 \
        --make-bed \
        --out 2nd_plinkQC


cp /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/code/0.vcf_to_plink.sh ./

