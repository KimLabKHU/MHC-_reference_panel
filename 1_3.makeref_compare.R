# # # compare position bewteen 2field and 4field ref
# library(dplyr); library(tidyr)
# # AA level
# dict_2f <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.txt", header=F, sep="\t")
# dict_2f <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.not2f.txt", header=F, sep="\t")
# map_2f <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.map", header=F, sep="\t")

# dict_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.Ggroup_4field.txt", header=F, sep="\t")
# map_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.Ggroup_4field.map", header=F, sep="\t")

# dict_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.txt", header=F, sep="\t")
# map_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.map", header=F, sep="\t")

# head(dict_gg)
# head(map_2f)
# head(map_gg)

# map2f = separate( map_2f , V2, sep="_", into=c("AA", "type", "AApos", "bp"))
# mapgg = separate( map_gg , V2, sep="_", into=c("AA", "type", "AApos", "bp", "exon")) 

# dict2f = separate( dict_2f, V1, sep="\\*", into=c("type", "allele") ) 
# dictgg = separate( dict_gg, V1, sep="\\*", into=c("type", "allele") ) 

# head(mapgg)
# head(dictgg)

# dict2f %>% mutate( n = nchar(V2)) %>% select(type, n) %>% table()
# nchar( dict2f$V2 )

# h="B"
# setdiff( filter(map2f, type==h)$AApos , filter(mapgg, type==h)$AApos )
# setdiff( filter(mapgg, type==h)$AApos , filter(map2f, type==h)$AApos )


# hla_class <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
# h="C"
# for( h in hla_class ){
# #print( grep( paste0("AA_", h, "_") , map_2f$V2, value=T) %>% length() == grep(paste0("AA_", h, "_") , map_gg$V2, value=T) %>% length() )
# a= setdiff( filter(map2f, type==h)$AApos , filter(mapgg, type==h)$AApos ) %>% length()
# b =setdiff( filter(mapgg, type==h)$AApos , filter(map2f, type==h)$AApos ) %>% length()
# print(paste( h, "_", a,b) )

# }


# head(mapgg)
# head(dictgg)

# for( h in hla_class ){
#     print(h)
#     filter(map2f, type==h) %>% nrow() %>% print()
#     nchar( filter(dict2f, type==h)$V2 )  %>% table() %>% print()


# }



# # snp level  생각해보니 2field는 AA만가져올거라서 SNP은 필요없음.
# dict_2f <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.txt", header=F, sep="\t")
# map_2f <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.map", header=F, sep="\t")
# dict_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_SNPS.hg38.Ggroup_4field.txt", header=F, sep="\t")
# map_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_SNPS.hg38.Ggroup_4field.map", header=F, sep="\t")

# dict_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_SNPS.hg38.imgt3320.txt", header=F, sep="\t")
# map_gg <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_SNPS.hg38.imgt3320.map", header=F, sep="\t")

# head(dict_gg)
# head(map_2f)
# head(map_gg)

# map2f = separate( map_2f , V2, sep="_", into=c( "SNP", "AA", "bp"))
# mapgg = separate( map_gg , V2, sep="_", into=c("SNP", "AA", "type", "bp", "exon")) 

# #dict_2f = combined
# dict2f = separate( dict_2f, V1, sep="\\*", into=c("type", "allele") ) 
# dictgg = separate( dict_gg, V1, sep="\\*", into=c("type", "allele") ) 

# head(map2f)
# head(mapgg)
# filter(mapgg, AA=="DRB1") %>% nrow()
# filter(map2f, AA=="A")
# hla_class <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")

# for( h in hla_class ){
# #print( grep( paste0("AA_", h, "_") , map_2f$V2, value=T) %>% length() == grep(paste0("AA_", h, "_") , map_gg$V2, value=T) %>% length() )
# a= setdiff( filter(map2f, AA==h)$bp , filter(mapgg, AA==h)$bp ) %>% length()
# b= setdiff( filter(mapgg, AA==h)$bp , filter(map2f, AA==h)$bp ) %>% length()
# print(paste( a,b) )

# }

# h="A"
# setdiff( filter(mapgg, AA==h)$bp , filter(map2f, AA==h)$bp )
# filter(map2f, AA=="A")  %>% nrow()
# filter(mapgg, AA=="A") %>% nrow
# filter(mapgg, bp=="29945456")

# nchar( filter(dict2f, type=="DRB1")$V2 )  %>% table()
# nchar( filter(dictgg, type=="DRB1")$V2 )  %>% table()


# nchar( filter(dict2f, type=="A")$V2 )  %>% table()
# nchar( filter(dictgg, type=="DQB1")$V2 )  %>% table()


# for( h in hla_class ){
#     print(h)
#     filter(mapgg, AA==h) %>% nrow() %>% print()
#     nchar( filter(dictgg, type==h)$V2 )  %>% table() %>% print()


# }




# new ggroup ref combined 2/3 from 4field dict.
rm(list=ls())
library(dplyr); library(tidyr)

##### 1. 2+4
#AA
setwd("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA")
f2_aa.txt <- read.table("./HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.txt", header=F, sep="\t")
f4_aa.txt <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.txt", header=F, sep="\t")
gg_aa.txt <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.Ggroup.txt", header=F)


not2f <- read.table(paste0("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/not2f.txt") , header=F)
#not2f <- separate(not2f, V1, into=c("type", "alle") , sep="\\*", remove=FALSE) %>% .[-which(.$type %in% c("C", "DQB1")) ,] %>% select(V1)
#not2f = not2f[-c(8,11:12),] %>% as.data.frame()

# not2f 가 f4_aa.txt에 있는지 확인.
no2f=setdiff(not2f[,1], c(gg_aa.txt$V1, f4_aa.txt$V1) ) # DQB1*06:146 빼고 다있음. 
no2f
write.table(no2f, "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/no2f.txt", col.names=F, row.names=F, quote=F, sep="\t")

head(f4_aa.txt)
head(f2_aa.txt)
head(gg_aa.txt)

f2_aa.txt$V3 %>% head()
combined = rbind(f2_aa.txt, f4_aa.txt[which(f4_aa.txt$V1 %in% not2f[,1]) , ] %>% mutate(V3="") )
combined = rbind(combined, gg_aa.txt[which(gg_aa.txt$V1 %in% not2f[,1]) , ] %>% mutate(V3="") )

combined = arrange(combined, V1)
which( combined$V1=="DPA1*02:07:01G" )

write.table(combined, "./kih_combine/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.not2f.txt", col.names=F, row.names=F, quote=F, sep="\t")
system("cp HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.map ./kih_combine/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.not2f.map")


# SNP
setwd("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA")
f2_snp.txt <- read.table("./HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.txt", header=F, sep="\t")
head(f2_snp.txt)

key = cbind( unlist( lapply( f2_snp.txt$V1, function(x) strsplit( x , "\\*")[[1]][1] ) ),
unlist( lapply( f2_snp.txt$V2, function(x) nchar(x) ) )   ) %>% unique() %>% as.data.frame() 
repeat_x <- function(n) {paste(rep("x", n), collapse = "")}
key$V3 <- sapply(key$V2, repeat_x)

not2f <- read.table(paste0("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/not2f.txt") , header=F) %>%.[,1]
#not2f <- separate(not2f, V1, into=c("type", "alle") , sep="\\*", remove=FALSE) %>% .[-which(.$type %in% c("C", "DQB1")) ,] %>% select(V1)

not2f=cbind(not2f , unlist( lapply( not2f, function(x) strsplit( x , "\\*")[[1]][1] ) ) ) %>% as.data.frame() %>% rename(V1=V2)
not2f=merge(not2f, key, by="V1") %>% select(not2f, V3) %>% rename(V1=not2f, V2=V3)

combined=rbind(f2_snp.txt, not2f)
combined = arrange(combined, V1)
head(combined)
which( combined$V1=="DPA1*02:07:01G" )

write.table(combined, "./kih_combine/HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.not2f.txt", col.names=F, row.names=F, quote=F, sep="\t")
system("cp HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.map ./kih_combine/HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.not2f.map")


##### 2. g+4
setwd("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38")

gg_aa.txt <- read.table("./HLA_DICTIONARY_AA.hg38.imgt3320.Ggroup.txt", header=F)

f2_aa.txt <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.txt", header=F, sep="\t")
f4_aa.txt <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.txt", header=F, sep="\t")

#f4_aa.txt <- read.table("./HLA_DICTIONARY_AA.hg38.imgt3320.txt", header=F)
head(gg_aa.txt)
head(f4_aa.txt)

# g group이 없는 2or3 추출
notgg_from4f <- read.table(paste0("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/notgg_from4f.txt") , header=F)
notgg_from2f <- read.table(paste0("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/1.HLA_for_NomenCleaner/notgg_from2f.txt") , header=F)

combined = rbind(gg_aa.txt, f2_aa.txt[which(f2_aa.txt$V1 %in% notgg_from2f[,1]) , ][1:2] ) %>% arrange(V1) 
combined = rbind(combined, f4_aa.txt[which(f4_aa.txt$V1 %in% notgg_from4f[,1]) , ][1:2] ) %>% arrange(V1) 


#f4_aa.txt_2or3 = f4_aa.txt[ sapply( f4_aa.txt$V1 , function(x) length( strsplit(x, ":")[[1]] ) ) < 4 , ] 
#combined = rbind(gg_aa.txt, f4_aa.txt_2or3) %>% arrange(V1) 
#combined$type = sapply( combined$V1, function(x) strsplit(x, "*", fixed=TRUE)[[1]][1])
#head(combined)
#combined$type %>% unique()
write.table(combined, "./kih_combine/HLA_DICTIONARY_AA.hg38.Ggroup_4field.txt", col.names=F, row.names=F, quote=F, sep="\t")

gg_snp.txt <- read.table("./HLA_DICTIONARY_SNPS.hg38.imgt3320.Ggroup.txt", header=F)

f2_snp.txt <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/2field_dic_from_SNP2HLA/HLA_DICTIONARY_SNPS.hg38.imgt3320.orderchange.txt", header=F, sep="\t")
f4_snp.txt <- read.table("/kimlab_wd/yuo1996/tools/HLA-TAPAS-master/MakeReference/data/hg38/HLA_DICTIONARY_SNPS.hg38.imgt3320.txt", header=F, sep="\t")

#f4_snp.txt <- read.table("./HLA_DICTIONARY_SNPS.hg38.imgt3320.txt", header=F)

# tmp <- f2_snp.txt[which(f2_snp.txt$V1 %in% notgg_from2f[,1]) , ][1:2] 
# tmp$V2 <- paste0(tmp$V2, ".")
# combined = rbind(gg_snp.txt, tmp) %>% arrange(V1) 
# combined = rbind(combined, f4_snp.txt[which(f4_snp.txt$V1 %in% notgg_from4f[,1]) , ][1:2] ) %>% arrange(V1) 

combined = rbind(gg_snp.txt, f4_snp.txt[which(f4_snp.txt$V1 %in% notgg_from4f[,1]) , ][1:2] ) %>% arrange(V1) 

#f4_snp.txt_2or3 = f4_snp.txt[ sapply( f4_snp.txt$V1, function(x) length( strsplit(x, ":")[[1]] ) ) < 4 , ] 
#combined = rbind(gg_snp.txt, f4_snp.txt_2or3) %>% arrange(V1) 
#combined$type = sapply( combined$V1, function(x) strsplit(x, "*", fixed=TRUE)[[1]][1])
#head(combined)
#combined$type %>% unique()
write.table(combined, "./kih_combine/HLA_DICTIONARY_SNPS.hg38.Ggroup_4field.txt", col.names=F, row.names=F, quote=F, sep="\t")

system("cp HLA_DICTIONARY_AA.hg38.imgt3320.Ggroup.map ./kih_combine/HLA_DICTIONARY_AA.hg38.Ggroup_4field.map")
system("cp HLA_DICTIONARY_SNPS.hg38.imgt3320.Ggroup.map ./kih_combine/HLA_DICTIONARY_SNPS.hg38.Ggroup_4field.map")

# At linux
#cp HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.map HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.not2f.map
#cp HLA_DICTIONARY_AA.hg38.imgt3320.Ggroup.map HLA_DICTIONARY_AA.hg38.imgt3320.Ggroup_4field.map
#cp HLA_DICTIONARY_SNPS.hg38.imgt3320.Ggroup.map HLA_DICTIONARY_SNPS.hg38.imgt3320.Ggroup_4field.map
