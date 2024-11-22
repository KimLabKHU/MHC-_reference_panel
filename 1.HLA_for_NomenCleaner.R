rm(list=ls())
library(dplyr);library(tidyr)

hla_result_all_uni <- read.table("/kimlab_wd/yuo1996/HLA/HLALA_working_DIR/hla_result_all.txt", header=T )
hla_result_all_kih <- read.table("/kimlab_wd/yuo1996/HLA/HLALA_working_DIR/kih/hla_result_all.txt", header=T )

#c4kih= read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/kih/total_table.txt", header=T)
#intersect( unique(hla_result_all_kih$ID) , c4kih$SAMPLE ) %>% length()

hla_result_all_kih$ID %>% unique()
# hla_result_all=hla_result_all_uni ; hla_result_all=hla_result_all_kih 
hla_result_all=rbind( hla_result_all_uni , hla_result_all_kih) %>% arrange(Locus, ID, Chromosome)

hla_result_all$ID %>% unique() %>% length()
head(hla_result_all)
nrow(hla_result_all)

#### select HLA locus
unique(hla_result_all$Locus)
hla_class <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
hla <- filter( hla_result_all, Locus %in% hla_class)
hla$ID %>% unique() %>% length()

nrow(hla)

#######  sample QC for relateness (only for UNIST) 
setwd("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/")

# 1. unrelted
unrelated <- read.table(paste0( "./0.vcf_to_plink/ver3/2nd_plinkQC.fam"), header=F, sep=" ")
unrelated <- gsub("-", "", unrelated$V2, fixed=TRUE)
length(unrelated)

#kihsnpsam = read.table("/kimlab_wd/yuo1996/C4_analysis/KIH/sam", header=F)
#kihsnpsam = as.vector(as.matrix( kihsnpsam[1,] ))
#length(kihsnpsam)
#intersect( kihsnpsam, unique( hla$ID) ) %>% length()
#intersect(unrelated, unique( hla$ID) ) %>% length()

head(hla)
hla$ID %>% unique() %>% length()  
hla <- filter( hla, ID %in% unrelated )
hla$ID %>% unique() %>% length()  

# 2. sample which have coverage >= 15 in all of eight HLA genes
n=length(hla_class)
n=8
rm_ID <- hla %>% select(Locus, ID, AverageCoverage) %>% unique()  %>% 
            mutate( inx = ifelse( AverageCoverage >=15 ,1,0)) %>% group_by(ID)  %>% mutate( inx_sum = sum(inx)) %>% 
            select(ID, inx_sum)%>% unique() %>% filter( inx_sum < n)  %>%  .$ID# %>% nrow

rm_ID
if( length(rm_ID) != 0 ){
hla <- hla[-which((hla$ID) %in% rm_ID ),  ]
}
hla$ID %>% unique() %>% length()  # 1000
hla$ID %>% unique() %>% grep("KOREA", .) %>% length()


# 3. sample QC with C4 typing (herv > 9 filtering)
total_table_uni <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/3rd_ver/total_table.txt", header=T )
total_table_kih <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/kih/total_table.txt", header=T )
total_table = rbind( total_table_uni, total_table_kih)
total_table$SAMPLE= gsub("-", "", total_table$SAMPLE)

nrow(total_table)

hla <- filter( hla, ID %in% total_table$SAMPLE)
hla$ID %>% unique() %>% length()  
hla = filter(hla, !(ID %in% filter(hla, Allele=="DQB1*06:146")$ID  ))


# 4. sampleQC with batch effect alleles  (DPA1*01:11)
nrow(hla)
hla$ID %>% unique() %>% length()  
hla = filter(hla, !(ID %in% filter(hla,  Allele=="DPA1*01:11" )$ID  ))
hla$ID %>% unique() %>% length()  


write.table(hla, paste0( "./1.HLA_for_NomenCleaner/hla_sampleQC.txt"), col.names=T, row.names=F, sep="\t", quote=F)

#### 1. no HLA allle QC ver



#### 2. HLA allele QC ver
nrow(hla) # 7136

hla$QC_result <- "out"
hla[which(hla$Q1 >0.9 & hla$AverageCoverage >20 & hla$proportionkMersCovered ==1), "QC_result"] <- "in"
filter(hla, QC_result =="in") %>% nrow
#hla_QC <- filter( hla, Q1 >0.9 , AverageCoverage >20 , proportionkMersCovered==1 ) #LocusAvgColumnError < 0.3 
#nrow(hla_QC)  # 6809

head(hla)


########## 추가!!! 
head(hla)

# 1. 2개가 동시에 typing된 경우 처리
hla[ which( sapply(strsplit(hla$Allele, ";"), function(x) length(x) ) ==2 ) ,]  %>% .$Allele %>% table()
hla[ which( sapply(strsplit(hla$Allele, ";"), function(x) length(x) ) ==2 ) ,]  %>% filter(Locus=="B")
hla[ which( sapply(strsplit(hla$Allele, ";"), function(x) length(x) ) ==2 ) ,]  %>% filter(Locus=="DPA1", QC_result=="out")  #%>% .$LocusAvgColumnError %>% summary()

hla[which(hla$Allele=="B*07:05:01;B*07:06"), "Allele"] <- "B*07:05:01G"
hla[which(hla$Allele=="DPA1*02:01:01;DPA1*02:01:08"), "Allele"] <- "DPA1*02:01:01G"
hla[which(hla$Allele=="DPA1*02:02:01"), "Allele"] = "DPA1*02:07:01G"

hla$resolution <- "notG"
hla[ which( substr( hla$Allele  ,nchar(hla$Allele),nchar(hla$Allele)) =="G" ), "resolution"] <- "G"

hla[which(hla$resolution=="notG"), "resolution"] <- ifelse( sapply(hla[which(hla$resolution=="notG"), "Allele"]  , function(x) length(gregexpr(":", x)[[1]])) ==2, "3field", "2field" )

head(hla)
select(hla, perfectG, resolution) %>% table() # !!!!
select(hla, Locus, Allele, perfectG, resolution ) %>% unique %>% select( perfectG, resolution) %>% table()

filter(hla, perfectG==1, resolution!="G") %>% .$Allele %>% table() # 이런얘들은 뭐지?  Ggroup이 애초에 없는 경우인가??
filter(hla, perfectG==1, resolution=="2field") %>% .$Allele %>% table() # 이런얘들은 뭐지?  Ggroup이 애초에 없는 경우인가??
filter(hla, perfectG==1, resolution=="3field") %>% .$Allele %>% table() # 이런얘들은 뭐지?  Ggroup이 애초에 없는 경우인가??

filter(hla, perfectG==0, resolution=="G") %>% .$Allele %>% table() # 전부 위에서 바꾼 두가지 경우
filter(hla, perfectG==0, resolution!="G") %>% head # 전부 G를 붇여줄까?
filter(hla, perfectG==0, resolution!="G") %>% .$Allele %>% table() 
filter(hla, perfectG==0, resolution=="3field") %>% head



# # just for checking
# # grep -v '#' hla_nom_g.txt | tr "*;" "\t" > ginof
# ginfo <- read.table("/kimlab_wd/yuo1996/HLA/ginfo", header=F, sep="\t")
# ginfo <- select(ginfo , !V2)
# filter(ginfo, V4=="")$V3  # Ggroup이 없는 2/3field
# setdiff( unique(filter(hla, perfectG==1, resolution!="G")$Allele ),
# paste0( ginfo[which(ginfo$V4==""), "V1"], "*", ginfo[which(ginfo$V4==""), "V3"] )
# )
# setdiff( unique(filter(hla, perfectG==1, resolution=="2field")$Allele ),
# paste0( ginfo[which(ginfo$V4==""), "V1"], "*", ginfo[which(ginfo$V4==""), "V3"] )
# )

# setdiff( unique( filter(hla, perfectG==0, resolution!="G")$Allele ) ,
# paste0( ginfo[which(ginfo$V4==""), "V1"], "*", ginfo[which(ginfo$V4==""), "V3"] )
# ) # -> 전부다임 (25) : 전부 해당 ggroup이 있다는 소리 전부 있는데 왜 ggroupd으로 안나왔지? ㅋ특히 3필드.
# setdiff( unique( filter(hla, resolution=="2field")$Allele ) ,
# paste0( ginfo[which(ginfo$V4==""), "V1"], "*", ginfo[which(ginfo$V4==""), "V3"] )
# ) # ggroup이 없는 얘들을 뺐으니 있는 얘들.? -> 전부다임 (25) : 전부 해당 ggroup이 있다는 소리 전부 있는데 왜 ggroupd으로 안나왔지? ㅋ특히 3필드.
# filter(hla, Allele %in% unique( filter(hla, perfectG==0, resolution!="G")$Allele ))$Allele %>% table()


# # 3field -> g group 로 매뉴얼리 바꿔줌.!! 이거내가할필요없지않나? 내가할피료없음 Nomencleaner가 알아서해준다... (이때 ggroup이 없는 3field는 그냥냅둠?)
# for( a in 1:length( hla[which(hla$resolution=="3field"),"Allele"] )){
#     i=hla[which(hla$resolution=="3field"),"Allele"][a]
#     h = strsplit(i, "*", fixed=TRUE)[[1]][1]
#     c = strsplit(i, "*", fixed=TRUE)[[1]][2]
    
#     nog = filter(ginfo, V1==h, V4=="")[which( filter(ginfo, V1==h, V4=="")$V3 ==c ), ]
#     g = filter(ginfo, V1==h, V4!="")[ which( sapply( strsplit( filter(ginfo, V1==h, V4!="")$V3 , "/" ), function(x) grep( paste0("^", c), x) %>% length() ) >= 1 ) , ]

#     if( nrow(nog)!=1 & nrow(g)!=0){
#         hla[which(hla$resolution=="3field"),"Allele"][a] <- g$V4
#     }
# }

# # take some time..?
# filter(hla, resolution=="3field")$Allele %>% table() %>% sort(decreasing=T)




##########

#### input for NomenCleaner in HLA-TAPAS
# 샘플단위로는 다빼고 , allele단위는 QC_result에 표시헤놓은 상태. 

# # 1. for ggroup 용 : 2field는 "0"으로 처리. ggroup(higher level)로 변환되지 않게.
# hla_forG <- hla
# head(hla_forG)
# filter(hla_forG, resolution=="2field")$Allele %>% unique() %>% sort()
# filter(hla_forG, resolution=="2field") %>% arrange( ID, Locus ) %>% head

# hla_forG[which(hla_forG$resolution=="2field"), "Allele"] %>% table() # ggroup타이핑이 안된 allele들
# hla_forG[which(hla_forG$resolution=="2field"), "Allele"] <- "0"
# write.table(hla_forG, paste0( "./1.HLA_for_NomenCleaner/", opt$snpver, "/hla_forG.txt"), col.names=T, row.names=F, sep="\t", quote=F)


# # 2/G , 2/3 인경우 억지로 둘다 Ggroup or 3field 로 -> Makereference에서 둘다 gor3로 처리 -> 나중에 한쪽만 "A"로 바꾸기. 
# idx = hla_forG[which(hla_forG$Allele==0), c("Locus", "ID")] 
# a=0
# for(i in 1:nrow(idx)){
#     a=a+1
#     print(a)
#     filter(hla_forG, Locus==idx[i,1], ID==idx[i,2])$resolution %>% sort() %>% print()
# }


# for(i in 1:nrow(idx) ){
#    if( filter(hla_forG, Locus==idx[i,1], ID==idx[i,2])$Allele %>% unique() %>% length() ==2 ){
#     hla_forG[ which(hla_forG$Locus==idx[i,1] & hla_forG$ID==idx[i,2] & hla_forG$Allele==0), "Allele"] <- hla_forG[ which(hla_forG$Locus==idx[i,1] & hla_forG$ID==idx[i,2] & hla_forG$Allele!=0), "Allele"]    
#    }
# }

# # 2/2 인 경우? -> 전부 "0"으로해놔야.(해논것) -> Makereference에서 missing 처리 -> 나중에 "A"로 바꾸기


# #hla_QC[grep(";" , hla_QC$Allele ),]  %>% head  # 92 명
# hla_input <- hla_forG %>% separate(Allele , into = c("allele1", "allele2"),sep=";" ) %>% arrange( ID, Locus)

# #system( paste0( "mkdir ./1.HLA_for_NomenCleaner/", opt$snpver) )
# #write.table(hla_input, paste0( "./1.HLA_for_NomenCleaner/", opt$snpver, "/hla_input" ), col.names=T, row.names=F, sep="\t", quote=F)

# nomen_input <- data.frame()
# for(i in unique(hla_input$ID)  ) {

#     allele_df <- data.frame(i, i, 0,0,0,0)
#     for(h in hla_class){
#         allele <- hla_input[which(hla_input$ID ==i & hla_input$Locus==h), "allele1"]
#         allele_df <- cbind(allele_df , data.frame( allele[1] , allele[2]))    
#     }

#     nomen_input <- rbind(nomen_input , allele_df)
# }


# head(nomen_input)
# write.table(nomen_input, paste0( "./1.HLA_for_NomenCleaner/", opt$snpver, "/nomen_input_ggroup.ped"), col.names=F, row.names=F, sep="\t", quote=F)


# 2. 2field용 : 그냥원래코드대로하면됨.

hla_for2 <- hla
hla_input <- hla_for2  %>% arrange( ID, Locus)


#system( paste0( "mkdir ./1.HLA_for_NomenCleaner/", opt$snpver) )
#write.table(hla_input, paste0( "./1.HLA_for_NomenCleaner/", opt$snpver, "/hla_input" ), col.names=T, row.names=F, sep="\t", quote=F)

nomen_input <- data.frame()
for(i in unique(hla_input$ID)  ) {

    allele_df <- data.frame(i, i, 0,0,0,0)
    for(h in hla_class){
        allele <- hla_input[which(hla_input$ID ==i & hla_input$Locus==h), "Allele"]
        allele_df <- cbind(allele_df , data.frame( allele[1] , allele[2]))    
    }

    nomen_input <- rbind(nomen_input , allele_df)
}


head(nomen_input)
dim(nomen_input)
write.table(nomen_input, paste0( "./1.HLA_for_NomenCleaner/nomen_input.ped"), col.names=F, row.names=F, sep="\t", quote=F)



#cd /kimlab_wd/yuo1996/tools/HLA-TAPAS-master
#python -m NomenCleaner --hped /kimlab_wd/yuo1996/C4_HLA_ref/try6/1.HLA_for_NomenCleaner/ver2/nomen_input.ped --out /kimlab_wd/yuo1996/C4_HLA_ref/try6/1.HLA_for_NomenCleaner/ver2/nomen_2field --2field
#python -m NomenCleaner --hped /kimlab_wd/yuo1996/C4_HLA_ref/try6/1.HLA_for_NomenCleaner/ver2/nomen_input.ped --out /kimlab_wd/yuo1996/C4_HLA_ref/try6/1.HLA_for_NomenCleaner/ver2/nomen_4field --4field
#python -m NomenCleaner --hped /kimlab_wd/yuo1996/C4_HLA_ref/try6/1.HLA_for_NomenCleaner/ver2/nomen_input.ped --out /kimlab_wd/yuo1996/C4_HLA_ref/try6/1.HLA_for_NomenCleaner/ver2/nomen_ggroup --Ggroup



