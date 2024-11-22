#!/bin/bash

ver=$1  # snp ver (try1 folder name)
field=$2 # 2field or ggroup

cd /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/${ver}

if [ "$field" == "2field" ]
then
    #cut -f 1-2 HLA.ped | tr "\t" " "  > ref.list
    cut -f 1-2 ../HLA_2field.ped | cut -f 2 > tmp
    paste tmp tmp > ref.list

    # sample selection !!!
    plink --bfile /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver}/2nd_plinkQC --keep ref.list --make-bed --out 3.whole_samples

    mkdir /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/${ver}/field2


    # Makereference v2
    cd /kimlab_wd/yuo1996/tools/HLA-TAPAS-master

    python -m MakeReference \
        --variants /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/${ver}/3.whole_samples \
        --chped /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/HLA_2field.ped \
        --hg 38 \
        --out /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/${ver}/field2/output \
        --dict-AA /kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/kih_combine/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.not2f \
        --dict-SNPS /kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/kih_combine/HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.not2f \
        --phasing \
        --nthreads 200 \
        --save-intermediates

#        --dict-AA /kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange_pri\
#        --dict-SNPS /kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange_pri \


elif [ "$field" == "ggroup" ]
then

    cut -f 1-2 ../HLA_ggroup.ped | cut -f 2 > tmp
    paste tmp tmp > ref.list

    # sample selection !!!
    plink --bfile /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/${ver}/2nd_plinkQC --keep ref.list --make-bed --out 3.whole_samples

    cd /kimlab_wd/yuo1996/tools/HLA-TAPAS-master

    mkdir /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/${ver}/ggroup

    python -m MakeReference \
        --variants /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/${ver}/3.whole_samples \
        --chped /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/HLA_ggroup.ped \
        --hg 38 \
        --out /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/${ver}/ggroup/output \
        --dict-AA MakeReference/data/hg38/kih_combine/HLA_DICTIONARY_AA.hg38.Ggroup_4field \
        --dict-SNPS MakeReference/data/hg38/kih_combine/HLA_DICTIONARY_SNPS.hg38.Ggroup_4field \
        --phasing \
        --nthreads 200 \
        --save-intermediates

fi

cp /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/code/2.MakeReference_with_HLA.sh /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/${ver}/

#        --chped /kimlab_wd/yuo1996/C4_HLA_ref/try6/2.MakeReference_with_HLA/${ver}/HLA.ped \


#    --dict-AA MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.Ggroup \
#    --dict-SNPS MakeReference/data/hg38/HLA_DICTIONARY_SNPS.hg38.imgt3320.Ggroup \

#    --dict-AA /kimlab_wd/yuo1996/tools/SNP2HLA_package_v1.0.3/MakeReference/HLA_DICTIONARY_AA \
#    --dict-SNPS /kimlab_wd/yuo1996/tools/SNP2HLA_package_v1.0.3/MakeReference/HLA_DICTIONARY_SNPS \