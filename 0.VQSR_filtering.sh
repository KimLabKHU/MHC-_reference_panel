#!/bin/bash

ver=$1 # ver3


# 2. VQSR filtering

cd /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver}

vcftools --gzvcf /kimlab_wd/yuo1996/C4_analysis/KIH/real1K/hla_1k_chr6.jointed.hard-filtered.recal.pass.vcf.gz \
        --remove-filtered-geno-all \
        --min-alleles 2 \
        --max-alleles 2 \
        --min-meanDP 15 \
        --minDP 15 \
        --minGQ 20 \
        --minQ 20 \
        --recode \
        --recode-INFO-all \
        --stdout | gzip -c > 1stQC_with_vcftools.vcf.gz


# it take some times (several miniutes)
# wc -l 171,457 (minDP 5-10, minGQ 15-20, min-meanDP 15-20)
# min-meanDP 30은 너무 쎄고 25까진 다버팀


# 2. hg38 convert

CrossMap.py vcf --compress \
        /kimlab_wd/yuo1996/tools/hg19ToHg38.over.chain \
        1stQC_with_vcftools.vcf.gz \
        /kimlab_wd/yuo1996/C4_analysis/2.Genomestrip_preprocessing_reference/Homo_sapiens_assembly38/mybwa_inx/Homo_sapiens_assembly38.fasta \
        1stQC_with_vcftools.hg38.vcf


vcftools --gzvcf 1stQC_with_vcftools.hg38.vcf.gz \
        --chr chr6 \
        --from-bp 25000000 \
        --to-bp 35000000 \
        --recode \
        --recode-INFO-all \
        --stdout | gzip -c >  chr6.MHC25to35mb.vcf.gz


# convert vcf.gz to plink 
plink --vcf chr6.MHC25to35mb.vcf.gz \
        --double-id \
        --make-bed --out 1st_plink

# plink basic qcs
plink --bfile 1st_plink \
        --list-duplicate-vars suppress-first --out dup  # 없음걍.

plink --bfile 1st_plink \
        --geno 0.1 \
        --maf 0.005 \
        --hwe 1e-7 \
        --make-bed \
        --out 2nd_plinkQC




# 3. extracting MAF 

plink --bfile 2nd_plinkQC --freq --missing --hardy

cd /kimlab_wd/yuo1996/C4_HLA_ref/try6/0.vcf_to_plink/${ver}/  # ref 
plink --bfile 2nd_plinkQC --freq --missing --hardy

# 4. compare with UNIST 1K SNPs (from WGS)
Rscript /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/code/0_1.compare_MAF.R --snpver ${ver}


# 5. extract newly selected (both-satisfied) SNPs in UNIST

cd /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver}/  # ref 

vcftools --gzvcf /kimlab_wd/yuo1996/C4_HLA_ref/try6/0.vcf_to_plink/${ver}/1stQC_with_vcftools.vcf.gz \
        --positions /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver}/newsnps.txt \
        --recode \
        --recode-INFO-all \
        --stdout | bgzip > newsnps.vcf.gz

# convert vcf.gz to plink 
plink --vcf newsnps.vcf.gz \
        --double-id \
        --make-bed --out new_plink


cp /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/code/0.VQSR_filtering.sh /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver}/


# 6. extract newly selected SNPs in KIH (for haplotype clusteirng)

cd /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver} # ref 


vcftools --gzvcf ./chr6.MHC25to35mb.vcf.gz \
        --positions /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver}/newsnps.txt \
        --recode \
        --recode-INFO-all \
        --stdout | bgzip > ./newsnps.kih.vcf.gz


# 7. Merge them 
bcftools sort newsnps.vcf.gz -Oz -o newsnps.sorted.vcf.gz
bcftools sort newsnps.kih.vcf.gz -Oz -o newsnps.kih.sorted.vcf.gz

tabix -p vcf newsnps.sorted.vcf.gz 
tabix -p vcf newsnps.kih.sorted.vcf.gz 

bcftools merge newsnps.sorted.vcf.gz newsnps.kih.sorted.vcf.gz -Ov -o merged.vcf 


# 8. "." snps ID change

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'ID="."' merged.vcf > missing_snps.txt
awk '{print $1"\t"$2"\t"$1"_"$2"_"$3"_"$4}' missing_snps.txt > update_list.txt
sort -k1,1 -k2,2n update_list.txt > sorted_update_list.txt
bgzip sorted_update_list.txt
tabix -s 1 -b 2 -e 2 sorted_update_list.txt.gz

bcftools annotate --annotations sorted_update_list.txt.gz --columns CHROM,POS,ID merged.vcf -o merged.updated.vcf

gzip merged.updated.vcf



# 9. last plink change

mkdir ./formerge 
mv * ./formerge

# convert vcf.gz to plink 
plink --vcf ./formerge/merged.updated.vcf.gz \
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


cp /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/code/0.VQSR_filtering.sh ./0.VQSR_filtering.${ver}.sh

