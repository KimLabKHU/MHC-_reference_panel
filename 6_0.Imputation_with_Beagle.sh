#!/bin/bash

ver=$1  # ver4
region=$2 # start31.7end32.2_highLDfiltering_50_5_0.8QUAL20
haplonetver=$3 # ver1_sujung
crossval_N=$4 #1

path=/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/${ver}/${region}_${haplonetver}/crossval${crossval_N}

mkdir /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/${ver}
mkdir /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/${ver}/${region}_${haplonetver}
mkdir /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/${ver}/${region}_${haplonetver}/crossval${crossval_N}

outpath=/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/${ver}/${region}_${haplonetver}/crossval${crossval_N}


LIST="/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/${ver}/${region}_${haplonetver}/survived_sams.txt"


cat "$LIST" | xargs -I {} -P 200 bash -c '
  FILE={}
  /kimlab_wd/yuo1996/tools/impute5/impute5_v1.2.0/impute5_v1.2.0_static \
  --h '"${path}"'/ref_panel${FILE}.vcf.gz \
  --g '"${path}"'/test_panel${FILE}.vcf.gz \
  --r 6:24000000-36000000 \
  --buffer-region 6:24000000-36000000 \
  --m /kimlab_wd/yuo1996/tools/shapeit4/maps/chr6.b38.gmap.gz \
  --threads 70 \
  --o '"${outpath}"'/test_panel${FILE}.imputed.vcf.gz \
  --l '"${outpath}"'/test_panel${FILE}.imputed.log
'


# cd /kimlab_wd/yuo1996/tools/beagle5.4
# cat "$LIST" | xargs -I {} -P 100 bash -c '
#   FILE={}
#   java -jar beagle.22Jul22.46e.jar gt='"${path}"'/test_panel${FILE}.vcf.gz \
#   ref='"${path}"'/ref_panel${FILE}.vcf.gz \
#   impute=true \
#   ap=true \
#   gp=true \
#   nthreads=1 \
#   chrom=6 \
#   map=/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/dependency/plink.chr6.GRCh38.map \
#   out='"${outpath}"'/test_panel${FILE}.imputed
# '


# toolpath="/kimlab_wd/yuo1996/tools/total_impute/rd-imputation-accuracy-develop/"
# cat "$LIST" | xargs -I {} -P 100 bash -c '
#   FILE={}
#   bash '"${toolpath}"'/Imputation_score_onlyfor_impute5.sh -i '"${path}"'/test_panel${FILE}.vcf.bgzip.gz \
#                                       -r '"${path}"'/ref_panel${FILE}.vcf.bgzip.gz \
#                                       -t 1 \
#                                       -o '"${outpath}"'/test_panel${FILE}.imputed.vcf.gz \
#                                       -c 6 
# '

####
gzip -d ${outpath}/*

#용량넘차지해서 삭제
#rm /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/${ver}/${region}_${haplonetver}/crossval${crossval_N}/ref*

paste <(echo $ver $region $haplonetver $crossval_N try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine) /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/code/6_0.Imputation_with_Beagle.sh >  ${outpath}/../6_0.Imputation_with_Beagle.log
