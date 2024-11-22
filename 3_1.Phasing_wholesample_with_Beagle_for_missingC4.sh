#!/bin/bash

ver=$1  # ver4
region=$2 # start31.7end32.2_highLDfiltering_50_5_0.8QUAL20
haplonetver=$3 # ver1_sujung

haplonet_folder="${region}_${haplonetver}"

# conver bgl to vcf
cd /kimlab_wd/yuo1996/tools/HLA-TAPAS-master/dependency
path="/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/"${ver}/${haplonet_folder}/


# rm intron region
grep -v intron ${path}/whole.bgl > ${path}/intronrm.whole.bgl
grep -v intron ${path}/whole.markers > ${path}/intronrm.whole.markers

cp ${path}/intronrm.whole.bgl ${path}/whole.bgl 
cp ${path}/intronrm.whole.markers ${path}/whole.markers 


java -jar beagle2vcf.jar 6 ${path}/intronrm.whole.markers ${path}/intronrm.whole.bgl 0 > ${path}/whole.bgl.vcf


########################################## eagle
tabix -p vcf ${path}/whole.bgl.vcf
/kimlab_wd/yuo1996/tools/Eagle_v2.4.1/eagle --geneticMapFile /kimlab_wd/yuo1996/tools/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz  \
                                            --numThreads 200 \
                                            --vcf ${path}/whole.bgl.vcf \
                                            --outPrefix ${path}/whole.eagle.phased

# change vcf to beagle (just for quality check)
cd ${path}
zcat whole.eagle.phased.vcf.gz | java -jar /kimlab_wd/yuo1996/tools/HLA-TAPAS-master/dependency/vcf2beagle.jar 0 whole.eagle.phased
gzip -kd whole.eagle.phased.bgl.gz
