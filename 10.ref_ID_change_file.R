rm(list=ls())
library(dplyr);library(tidyr)

setwd("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/ver4/start31.7end32.2_highLDfiltering_50_5_0.8QUAL20_ver1_sujung")
markers_pre=read.table("./whole.markers", header=F)

binaryaa=read.table("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/ver4/field2/output.bim", header=F)
ba=filter(binaryaa, grepl("^AA_", V2)) %>% filter( V5 != "p") %>% rename( ID=V2, REF=V6, ALT=V5)
ba$newID=paste0(ba$ID, "_", ba$REF, "_", ba$ALT)

res = rename(ba, CHROM=V1, POS=V4, OLD_ID=ID, NEW_ID=newID) %>% select( CHROM, POS, OLD_ID, NEW_ID) 
head(res)
res=select(res, OLD_ID, NEW_ID)

write.table( res, "ID_change.txt", col.names = FALSE, row.names = FALSE, quote=F, sep="\t" )
