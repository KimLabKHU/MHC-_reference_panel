#!/bin/bash
ver=$1  # snp ver (try1 folder name)
region=$2 # region=$2 #start31.7end32.2QUAL30
haplonetver=$3 # ver2_best

folderpath="/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/"
codepath="/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/code/"

# mkdir ${folderpath}/0.vcf_to_plink
# mkdir ${folderpath}/1.HLA_for_NomenCleaner
# mkdir ${folderpath}/2.MakeReference_with_HLA
# mkdir ${folderpath}/3.phasing_with_Beagle_C4_HLA
# mkdir ${folderpath}/5.Imputation_C4HLA_using_Beagle
# mkdir ${folderpath}/7.Accuracy

mkdir ${folderpath}/0.vcf_to_plink/${ver}
#mkdir ${folderpath}/1.HLA_for_NomenCleaner/${ver} # 필요없음 합쳐서 가는거로 수정.
mkdir ${folderpath}/2.MakeReference_with_HLA/${ver}
mkdir ${folderpath}/3.phasing_with_Beagle_C4_HLA/${ver}
mkdir ${folderpath}/5.Imputation_C4HLA_using_Beagle/${ver}
mkdir ${folderpath}/7.Accuracy/${ver}

mkdir ${folderpath}/3.phasing_with_Beagle_C4_HLA/${ver}/${haplonetver}


# start !
### from try6
#bash ${codepath}/0.vcf_to_plink.sh ${ver} # 각각 큐시하고 합쳐서 플링크 큐시한번더하는 형태로.
bash ${codepath}/0.VSQR_filtering.sh ${ver} 

# < HLA관련해선 snpver상관없이 한번만 하면됨.
Rscript ${codepath}/1.HLA_for_NomenCleaner.R --snpver ${ver} # 여기서 샘플큐씨, HLA (8개모두 cov > 15) + C4 (HERV > 9)typing 관련된. 
cd /kimlab_wd/yuo1996/tools/HLA-TAPAS-master
python -m NomenCleaner --hped /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/nomen_input.ped --out /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/nomen_2field --2field
# ggroup ref은 그냥 HLA-LA 결과그대로 쓰면되서 nomencleaner할필요없음; 하면 오히려 ggroup없는 얘들은 0처리되버려서안됨.
# ggroup도 하긴 해야됨.. 3필드중에 ggroup으로 바뀌는 얘들이있어서..ㅇㄴ..
python -m NomenCleaner --hped /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/nomen_input.ped --out /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/nomen_ggroup --Ggroup
Rscript ${codepath}/1_2.check_Nomencleaner_result.R 
Rscript ${codepath}/1_3.makeref_compare.R  # ref 다시만듬

#cd ${folderpath}/1.HLA_for_NomenCleaner
#cat no2f.txt nogg.txt  | sort | uniq > nodict
#dictionary에 없는 allele있는 샘플 제거.  -> 걍 앞 코드에서 해버림. 1.HLA_for_NomenCleaner.R 에서

# 여까지.한번만하면됨. >

bash ${codepath}/2.MakeReference_with_HLA.sh ${ver} 2field
bash ${codepath}/2.MakeReference_with_HLA.sh ${ver} ggroup

cp ${folderpath}/2.MakeReference_with_HLA/${ver}/ref.list ${folderpath}/2.MakeReference_with_HLA/

# 여까지하고 C4_genomestrip 코드로 넘어감. -> 한번만하면됨. (haplonet 말고 그전단계까지.)
# haplonet 코드 돌림.

#Rscript ${codepath}/3_0.Merge_C4_into_bgl_forphaisng.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver}
Rscript ${codepath}/3_0.Merge_C4_into_bgl_forphaisng_sampleQC.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver} --qual ${qual}
Rscript ${codepath}/3_0.Merge_C4_into_bgl_forphaisng_sampleQC_sujung.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver} --qual ${qual}
#cp ${codepath}/3_0.Merge_C4_into_bgl_forphaisng.R ${folderpath}/3.phasing_with_Beagle_C4_HLA/${ver}/  # for checking haplonet version

###
# 만약에 Sample QC 버전이면 region/haplonet QUAL~ 붙여서 수정해서 돌리기.
bash ${codepath}/3_1.Phasing_wholesample_with_Beagle_for_missingC4.sh ${ver} ${region} ${haplonetver}

Rscript ${codepath}/3_2.Remove_not2copy_samples_after_phasing.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver}
Rscript ${codepath}/3_2.Remove_not2copy_samples_after_phasing_sujung.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver}

bash ${codepath}/3_3.Subsampling_for_crossval.sh ${ver} ${region} ${haplonetver} 1

bash ${codepath}/6_0.Imputation_with_Beagle.sh ${ver} ${region} ${haplonetver} 1

# 이거 돌려야 phased 된 결과 (C4) 랑 imputation결과랑 비교가능? cf. HLA는 0없었나???
Rscript ${codepath}/After3_0.justforcheck_phasedC4_and_typedC4.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver} --crossval crossval1
(걍안돌려도됨)Rscript ${codepath}/After3_0.justforcheck_phasedC4_and_typedC4_sujung.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver} --crossval crossval1

Rscript ${codepath}/7_1.Accuracy_check.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver} --crossval crossval1
Rscript ${codepath}/7_1.Accuracy_check_sujung.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver} --crossval crossval1

Rscript ${codepath}/7_2.Accuracy_check_in_C4allele.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver} --crossval crossval1
Rscript ${codepath}/7_2.Accuracy_check_in_C4allele_sujung.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver} --crossval crossval1


# additional analysis for haploid C4 structure and HLA-DRB1 1501
Rscript ${codepath}/9.LD_bt_haploid_C4_and_HLA_DRB1.R --snpver ${ver} --region ${region} --haplonetver ${haplonetver}



# two residue variants ID correction

Rscript ${codepath}/10.ref_ID_change_file.R
cd /kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/ver4/start31.7end32.2_highLDfiltering_50_5_0.8QUAL20_ver1_sujung

bcftools view -H whole.eagle.phased.sampleQC.AC.vcf.gz | \
awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$2; next} $3 in a {$3=a[$3]}1' ID_change.txt - > modified_variants.txt

bcftools view -h whole.eagle.phased.sampleQC.AC.vcf.gz > header.txt
cat header.txt modified_variants.txt | bgzip -c > whole.eagle.phased.sampleQC.AC.variant_corrected.vcf.gz
tabix -p vcf whole.eagle.phased.sampleQC.AC.variant_corrected.vcf.gz