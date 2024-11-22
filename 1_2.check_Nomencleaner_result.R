# #알아냄. map 파일 순서가다름. 안에서도 뒤짚혀있는듯?>? 확인했나..??????? -> 수정? 
# #한번만 했으면됨.
# rm(list=ls())
# library(dplyr); library(tidyr); library(stringr)
# setwd("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference")

# f4_AA <- read.table("./data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.map", header=F)
# f4_SNP <- read.table("./data/hg38/HLA_DICTIONARY_SNPS.hg38.imgt3320.map", header=F)

# f4_AA <- separate( f4_AA , V2, into=c("a", "b", "c", "d", "e"), sep="_" ) 
# f4_SNP <- separate( f4_SNP , V2, into=c("a", "b", "c", "d", "e"), sep="_" ) 

# head(f4_AA)
# f4_AA$b %>% unique()
# f4_SNP$b %>% unique()
# filter(f4_SNP, b=="SNPS") %>% head

# f2_AA <- read.table("./2field_dic_from_SNP2HLA/HLA_DICTIONARY_AA.hg38.imgt3320.map", header=F)
# f2_SNP <- read.table("./2field_dic_from_SNP2HLA/HLA_DICTIONARY_SNPS.hg38.imgt3320.map", header=F)
# head(f2_SNP)
# head(f2_AA)

# f2_AA <- separate(f2_AA, V2, into=c("a", "b", "c", "d"), sep="_", remove=FALSE) %>% select( !c("a","d") ) 
# f2_SNP <- separate(f2_SNP, V2, into=c("a", "b", "c"), sep="_", remove=FALSE) %>% select( !c("a", "c") ) 

# f2_AA$c <- as.numeric(f2_AA$c)
# f2_AA$b %>% unique()
# f2_SNP$b %>% unique()

# head(f2_AA)
# head(f2_SNP)

# f2_AA_update = arrange( f2_AA, b, c)  %>% select(V1, V2, V3, V4)
# head(f2_AA_update)
# nrow(f2_AA_update)

# f2_SNP_update = arrange( f2_SNP, b ) 
# head(f2_SNP)
# f2_SNP_update$b %>% unique()
# # f2_SNP_update = rbind(f2_SNP_update[which(f2_SNP_update$b=="A"), ] ,
# #                 f2_SNP_update[which(f2_SNP_update$b=="B"), ][ order( f2_SNP_update[which(f2_SNP_update$b=="B"), "V4"] , decreasing = T) , ] ,
# #                 f2_SNP_update[which(f2_SNP_update$b=="C"), ][ order( f2_SNP_update[which(f2_SNP_update$b=="C"), "V4"] , decreasing = T) , ] ,
# #                 f2_SNP_update[which(f2_SNP_update$b=="DPA1"), ][ order( f2_SNP_update[which(f2_SNP_update$b=="DPA1"), "V4"] , decreasing = T) , ] ,
# #                 f2_SNP_update[which(f2_SNP_update$b=="DPB1"), ] ,
# #                 f2_SNP_update[which(f2_SNP_update$b=="DQA1"), ] ,
# #                 f2_SNP_update[which(f2_SNP_update$b=="DRB1"), ][ order( f2_SNP_update[which(f2_SNP_update$b=="DRB1"), "V4"] , decreasing = T) , ] 
# #             ) %>% select(V1, V2, V3, V4)

# head(f2_SNP_update)
# head(f2_SNP_update)

# nrow(f2_AA); nrow(f2_SNP)
# nrow(f2_AA_update); nrow(f2_SNP_update)

# #

# write.table(f2_AA_update, "./2field_dic_from_SNP2HLA/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.map", row.names=F, col.names=F, sep="\t" , quote=F)
# write.table(f2_SNP_update, "./2field_dic_from_SNP2HLA/HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.map", row.names=F, col.names=F, sep="\t" , quote=F)

rm(list=ls())
library(dplyr); library(tidyr)

setwd("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/")
res_2f <- read.table( paste0( "./1.HLA_for_NomenCleaner/nomen_2field.chped"), header=F, sep="\t")[,-c(2:6)]
res_gg <- read.table( paste0( "./1.HLA_for_NomenCleaner/nomen_ggroup.chped"), header=F, sep="\t")[,-c(2:6)]
input <- read.table( paste0( "./1.HLA_for_NomenCleaner/nomen_input.ped"), header=F, sep="\t")[,-c(2:6)]
nrow(res_2f);nrow(res_gg);nrow(input)
# colnames(res_2f) <- c("ID", rep("A",2), rep("B",2), rep("C",2), rep("DPA1",2), rep("DPB1",2), rep("DQA1",2), rep("DQB1",2), rep("DRB1",2))
# colnames(res_gg) <- c("ID", rep("A",2), rep("B",2), rep("C",2), rep("DPA1",2), rep("DPB1",2), rep("DQA1",2), rep("DQB1",2), rep("DRB1",2))
# colnames(input) <- c("ID", rep("A",2), rep("B",2), rep("C",2), rep("DPA1",2), rep("DPB1",2), rep("DQA1",2), rep("DQB1",2), rep("DRB1",2))


head(res_2f)
head(res_gg)
head(input)

table(res_2f$V1==res_gg$V1)
table(input$V1==res_gg$V1)
table( res_2f==0 )
table( res_gg==0 )

# 2field 0 없애기
head(res_2f)
head(res_gg)
res_2f[which(res_2f==0, arr.ind=TRUE)[,1], ] 
input[which(res_2f==0, arr.ind=TRUE) ] 

# res_2f[which(res_2f==0, arr.ind=TRUE) ]  <- "DPA1*02:02"
res_2f[which(res_2f=="DPA1*02:07", arr.ind=TRUE) ]  <- "DPA1*02:07:01G"

# ggroup 0 없애기
res_gg[which(res_gg==0, arr.ind=TRUE)[,1], ] %>% head
input[which(res_gg==0, arr.ind=TRUE) ] %>% unique()
res_gg[which(res_gg==0, arr.ind=TRUE) ]  <- input[which(res_gg==0, arr.ind=TRUE) ]

head(res_2f)
head(res_gg)

table( res_2f==0 )
table( res_gg==0 )

#### allele들이 dictionary 에 다 있는지 함 확인.

# 2field
unique( res_2f[,-1]  %>% as.matrix() %>% as.vector() )
dict_2f <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.txt", header=F, sep="\t")
not2f=setdiff( unique( res_2f[,-1]  %>% as.matrix() %>% as.vector() ) , dict_2f$V1 ) # 얘네들은 여기서 안풀리니까 AA level에서 분석불가? 
not2f
write.table(not2f, paste0("./1.HLA_for_NomenCleaner/not2f.txt"), col.names=F , row.names=F , sep="\t" , quote=F)

# ggroup
unique( res_gg[,-1]  %>% as.matrix() %>% as.vector() )
dict_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.Ggroup.txt", header=F, sep="\t")
dict_4f <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.txt", header=F, sep="\t")

notgg=setdiff( unique( res_gg[,-1]  %>% as.matrix() %>% as.vector() ) , dict_gg$V1 )
notgg_from4f = intersect(notgg,  dict_4f$V1)
intersect(notgg,  dict_2f$V1)
setdiff( notgg , dict_4f$V1 ) # 뒤에 두개는 아예 없어서 못함. 
notgg_from2f = intersect( dict_2f$V1,  setdiff( notgg , dict_4f$V1 )  )

# 2f/4f어디에도 없는 gg allele
nogg=setdiff( notgg, union(dict_4f$V1, dict_2f$V1) )
nogg
write.table(nogg, paste0("./1.HLA_for_NomenCleaner/nogg.txt"), col.names=F , row.names=F , sep="\t" , quote=F)

write.table(notgg_from4f, paste0("./1.HLA_for_NomenCleaner/notgg_from4f.txt"), col.names=F , row.names=F , sep="\t" , quote=F)
write.table(notgg_from2f, paste0("./1.HLA_for_NomenCleaner/notgg_from2f.txt"), col.names=F , row.names=F , sep="\t" , quote=F)

# dictionary에 없는 allele들 frequency 확인 -> 뒤에서 어떻게 coding 되는지도 확인할 것. 

hla_result_all <- read.table("/kimlab_wd/yuo1996/HLA/HLALA_working_DIR/kih/hla_result_all.txt", header=T )
head(hla_result_all)
filter(hla_result_all, Allele=="A*11:32")
filter(hla_result_all, Allele=="DPA1*02:02:01")
filter(hla_result_all, Allele=="DQB1*06:146")


# # union
# setdiff( c( unique( res_2f[,-1]  %>% as.matrix() %>% as.vector() ) ,   unique( res_gg[,-1]  %>% as.matrix() %>% as.vector() )       )    , unique( c(dict_gg$V1, dict_2f$V1 ) ) )  # 얘넨 어디에도없음 둘다에서. 분석불가. 

# hla <- read.table(paste0( "./1.HLA_for_NomenCleaner/", opt$snpver, "/hla_sampleQC.txt"), header=T )
# filter(hla, Allele=="DQB1*06:146" | Allele=="DPA1*02:02:01")
# filter(hla, Allele %in% not2f)$Allele %>% table()

# setdiff( setdiff( unique( res_2f[,-1]  %>% as.matrix() %>% as.vector() ) , dict_2f$V1 ) , dict_gg$V1 ) 


### make ped file for Makereference

#res_gg$V1 <- gsub("KOREA1K", "KOREA1K-", res_gg$V1)
#res_2f$V1 <- gsub("KOREA1K", "KOREA1K-", res_2f$V1)

res_gg <- res_gg %>% mutate(V2=V1, V3=0, V4=0, V5=0, V6=0) %>% relocate( paste0("V", 2:6) , .after=V1)
res_2f <- res_2f %>% mutate(V2=V1, V3=0, V4=0, V5=0, V6=0) %>% relocate( paste0("V", 2:6) , .after=V1)


nrow(res_gg);nrow(res_2f)
nrow(res_gg)
nrow(res_2f)

head(res_gg)
head(res_2f)

res_gg$V1=gsub("KOREA1K", "KOREA1K-", res_gg$V1)
res_2f$V1=gsub("KOREA1K", "KOREA1K-", res_2f$V1)
res_gg$V2=gsub("KOREA1K", "KOREA1K-", res_gg$V2)
res_2f$V2=gsub("KOREA1K", "KOREA1K-", res_2f$V2)


write.table( res_2f, paste0( "./2.MakeReference_with_HLA/HLA_2field.ped") , col.names=F, row.names=F, sep="\t", quote=F)
write.table( res_gg, paste0( "./2.MakeReference_with_HLA/HLA_ggroup.ped") , col.names=F, row.names=F, sep="\t", quote=F)
