rm(list=ls())
library(dplyr);library(tidyr)
library(bigreadr)
setwd("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/ver4/start31.7end32.2_highLDfiltering_50_5_0.8QUAL20_ver1_sujung")


bgl=fread2("./whole.phased.sampleQC.bgl", header=T)
c4_bgl=bgl[grep("copy", bgl$id), ]
hlaid=bgl[grep("HLA", bgl$id), "id"]
hla2field=hlaid[grep("^[^:]+:[^:]+$", hlaid)] %>% .[-grep("gg_" ,.)]
HLA_bgl=bgl[which(bgl$id %in% hla2field), ]
HLA_bgl[1:10,1:10]

i=1
bgl[1:10,1:10]
res=matrix(NA, ncol=9, nrow=(dim(bgl)[2]-2))
rownames(res)=colnames(bgl)[-c(1:2)]
colnames(res)=c("c4", "HLA_A", "HLA_B", "HLA_C", "HLA_DRB1", "HLA_DQA1", "HLA_DQB1", "HLA_DPA1", "HLA_DPB1")
head(res)

for(i in 1:(dim(bgl)[2]-2)){

print(i)
    res[i, "c4"]=c4_bgl[c4_bgl[, i+2] == "T", "id"] %>%  gsub("copy", "" ,. ) %>% gsub("_", "" , .) %>% gsub("C4", "", .) %>% gsub("HERV", "L", .) %>% paste0( . , collapse = "_")

    res[i, "HLA_A"]=    ifelse( length(HLA_bgl[grepl("^HLA_A", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"])==1,  HLA_bgl[grepl("^HLA_A", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"], NA)
    res[i, "HLA_B"]=    ifelse( length(HLA_bgl[grepl("^HLA_B", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"])==1,  HLA_bgl[grepl("^HLA_B", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"], NA)
    res[i, "HLA_C"]=    ifelse( length(HLA_bgl[grepl("^HLA_C", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"])==1,  HLA_bgl[grepl("^HLA_C", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"], NA)
    res[i, "HLA_DRB1"]=    ifelse( length(HLA_bgl[grepl("^HLA_DRB1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"])==1, HLA_bgl[grepl("^HLA_DRB1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"], NA)
    res[i, "HLA_DQA1"]=    ifelse( length(HLA_bgl[grepl("^HLA_DQA1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"])==1, HLA_bgl[grepl("^HLA_DQA1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"], NA)
    res[i, "HLA_DQB1"]=    ifelse( length(HLA_bgl[grepl("^HLA_DQB1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"])==1, HLA_bgl[grepl("^HLA_DQB1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"], NA)
    res[i, "HLA_DPA1"]=    ifelse( length(HLA_bgl[grepl("^HLA_DPA1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"])==1, HLA_bgl[grepl("^HLA_DPA1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"], NA)
    res[i, "HLA_DPB1"]=    ifelse( length(HLA_bgl[grepl("^HLA_DPB1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"] )==1, HLA_bgl[grepl("^HLA_DPB1", HLA_bgl$id) & HLA_bgl[, i+2] == "T", "id"], NA)

}

res=as.data.frame(res)
head(res)

sort( table(res$c4) )
filter(res, c4=="A0_B2_L1")
filter(res, c4=="A0_B1_L1")
filter(res, HLA_DRB1=="HLA_DRB1*15:01")$c4 %>% table()
filter(res, HLA_DRB1=="HLA_DRB1*13:02")$c4 %>% table()
