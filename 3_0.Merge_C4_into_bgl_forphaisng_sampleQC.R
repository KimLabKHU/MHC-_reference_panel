rm(list=ls())
library(bigreadr)
library(dplyr);library(tidyr)
library(optparse)

##################################  arguement setting ########################################
args_list <- list(
    make_option( "--snpver" , type="character" , default= "ver3" , help="snpversion", metavar="character"),
    make_option( "--region" , type="character" , default= "start31.7end32.2" , help="haplonet_region", metavar="character"),
    make_option( "--haplonetver" , type="character" , default= "haplonet_RMSE_under0.9" , help="snpversion", metavar="character"),
    make_option( "--qual" , type="integer" , default= 5 , help="QUAL")
)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)
 
# opt=list(); opt$snpver="ver3"; opt$region="start31.7end32.2"; opt$haplonetver="ver2_best"; opt$qual=30
# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2"; opt$haplonetver="ver1_best"; opt$qual=30
###############################################################################################


##########################################
load("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/C4_genomestrip/vcf_to_GT.RData")
setwd("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/C4_genomestrip/")

system( paste0( "grep '#CHROM' ", "bat123456789.C4genotype.output.vcf", " > " , "C4_samples"  ))
header <- readLines( "C4_samples" )  # C4AB likelihood 고려 안한버전.
samples <- strsplit( header , split="\t")[[1]]
samples <- samples[-c(1:9)]
rownames(C4HERV_GT) <- samples

system( paste0( "grep '#CHROM' ", "bat123456789.C4ABgenotype.output.vcf", " > " , "C4AB_samples"  ))
header <- readLines( "C4AB_samples" )  # C4AB likelihood 고려 안한버전.
samples <- strsplit( header , split="\t")[[1]]
samples <- samples[-c(1:9)]
rownames(C4AB_GT) <- samples

C4_copy <- cbind(C4HERV_GT, C4AB_GT) 
C4_GT <- select(C4_copy , GT_C4,  GT_C4A, GT_C4B , GT_HERV)

head(C4_copy)
head(C4_GT)
nrow(C4_GT)

C4_GT$GT_C4 %>% table()
C4_GT$GT_HERV %>% table()
nrow(C4_GT)
head(C4_GT)

C4_GT[ rownames(C4_GT)=="NIH21D7020955", ]

########################################### !!!!! change by haplonet_proximal !!!!!!!!

C4_GT$SAMPLE <- rownames(C4_GT)
head(C4_GT)
# haplonet_path <- paste0("/kimlab_wd/yuo1996/C4_HLA_ref/haplonet_proximal_with_Eagle/combine_kih/4.haplonet/", opt$snpver, "/", opt$region, "/training/result/CON_genotyping_forcrosval/only_all/") 
haplonet_path <- paste0("/kimlab_wd/yuo1996/C4_HLA_ref/haplonet_proximal_with_Eagle/combine_kih/4.haplonet/", opt$snpver, "/", opt$region, "/training/result/QUAL", opt$qual, "/CON_genotyping_forcrosval/only_all/") 

set = read.table( paste0(haplonet_path, opt$haplonetver) , header=F )

C4_path <- paste0( haplonet_path, grep("_C4_", set[,1]  ,value=T))
C4A_path <- paste0( haplonet_path, grep("_C4A_", set[,1]  ,value=T))
C4B_path <- paste0( haplonet_path, grep("_C4B_", set[,1]  ,value=T))
HERV_path <- paste0( haplonet_path, grep("_HERV_", set[,1]  ,value=T))

# 여기서 실제 copy number에 해당하는 genotype으로 갈아끼움. 
#C4
geno_prob_C4 <- read.table( C4_path , header=T)
head(geno_prob_C4)
C4_GT <- merge(C4_GT, geno_prob_C4[c("SAMPLE", "estimated_geno_matched_typedC4")] , by="SAMPLE", all.x=TRUE) %>% select(!GT_C4) %>% rename( GT_C4=estimated_geno_matched_typedC4) %>% relocate( GT_C4, .after=SAMPLE)  
head(C4_GT)

#C4A
geno_prob_C4 <- read.table( C4A_path , header=T)
head(geno_prob_C4)
C4_GT <- merge(C4_GT, geno_prob_C4[c("SAMPLE", "estimated_geno_matched_typedC4A")] , by="SAMPLE", all.x=TRUE) %>% select(!GT_C4A) %>% rename( GT_C4A=estimated_geno_matched_typedC4A) %>% relocate( GT_C4A, .after=GT_C4)  

#C4B
geno_prob_C4 <- read.table( C4B_path , header=T)
C4_GT <- merge(C4_GT, geno_prob_C4[c("SAMPLE", "estimated_geno_matched_typedC4B")] , by="SAMPLE", all.x=TRUE) %>% select(!GT_C4B) %>% rename( GT_C4B=estimated_geno_matched_typedC4B) %>% relocate( GT_C4B, .after=GT_C4A)  

#HERV
geno_prob_C4 <- read.table( HERV_path , header=T)
C4_GT <- merge(C4_GT, geno_prob_C4[c("SAMPLE", "estimated_geno_matched_typedHERV")] , by="SAMPLE", all.x=TRUE) %>% select(!GT_HERV) %>% rename( GT_HERV=estimated_geno_matched_typedHERV) %>% relocate( GT_HERV, .after=GT_C4B)  

rownames(C4_GT) <- C4_GT$SAMPLE
C4_GT <- select( C4_GT, !SAMPLE)
head(C4_GT)
# filter(C4_GT, GT_C4==0 | GT_C4A==0 | GT_C4B==0 | GT_HERV==0  )

# when only_under0.9
#C4_GT[which(C4_GT == 0, arr.ind=TRUE)] <- NA # when only_under0.9


#C4_GT[which(C4_GT == 0, arr.ind=TRUE)] <- "9/9"
#filter(C4_GT, GT_C4=="9/9" | GT_C4A=="9/9" | GT_C4B=="9/9" | GT_HERV=="9/9"  ) 

# filter(C4_GT, GT_C4=="9/9") %>% nrow()
# filter(C4_GT, GT_C4A=="9/9") %>% nrow()
# filter(C4_GT, GT_C4B=="9/9") %>% nrow()
# filter(C4_GT, GT_HERV=="9/9") %>% nrow()

nrow(C4_GT)
head(C4_GT)

############################################################################

if(!dir.exists(  paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/",  opt$snpver,"/",opt$region, "QUAL",opt$qual, "_", opt$haplonetver )   )){ dir.create(path= paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/",  opt$snpver,"/", opt$region,  "QUAL",opt$qual, "_", opt$haplonetver)  ) } 

############ Merge  HLA.bgl with C4 info  for phasing
bgl_hla <- fread2(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/", opt$snpver, "/", "/integrated.bgl"), sep=" ", header=F)
bgl_hla[1:10,1:10]
dim(bgl_hla)
bgl_hla_markers <- fread2( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/", opt$snpver,  "/integrated.markers_redefined") , sep=" ")
head(bgl_hla_markers)
dim(bgl_hla_markers)

C4_markers <- c( paste0("C4_copy", c(2,1,3,4) ), 
                 paste0("C4A_copy", c(2,0,1,3) ),
                 paste0("C4B_copy", c(2,0,1,3) ), 
                 paste0("HERV_copy", c(2,0,1,3,4) ) ) %>% as.data.frame()

colnames(C4_markers) <- "C4"

# make C4.bgl
bgl_hla[1:10,1:10]
bgl_sam_order <- bgl_hla[2,] %>% as.matrix() %>% as.vector()
bgl_sam_order <- bgl_sam_order[-c(1:2)] # HLA sample QCed sampleds

# HLA sampleQCed 된 샘플들로 C4 table 만듦.
C4_insert_for_bgl <- matrix(NA, ncol= (length(bgl_sam_order)+1), nrow= (nrow(C4_markers) +1) )
C4_insert_for_bgl[1,] <-c("id" , bgl_sam_order)
C4_insert_for_bgl[-1,1] <-C4_markers$C4 
C4_insert_for_bgl <- as.data.frame(C4_insert_for_bgl)
C4_insert_for_bgl[1:18,1:10]

C4_GT[which(rownames(C4_GT) %in% c("KOREA1K-11", "KOREA1K-15")),]

# fill in C4_insert_for_bgl with A(absent)/P(present)

i=1
# samples after QC = 446
for( i in 1:length(unique(bgl_sam_order)) ){
    
    a<- 2*i
    nam <- C4_insert_for_bgl[ 1,a]
    tmp <- t( C4_GT[which(rownames(C4_GT) ==nam ),] ) %>% as.data.frame()
    colnames(tmp) <- "abc"
    tmpG <- tmp %>%  separate( abc, into=c("G1", "G2"), sep="/")
        
        for(j in 1:2){
            C4_insert_for_bgl[ grep("C4_copy", C4_insert_for_bgl$V1), a ] <-  case_when(
                                                                                tmpG[1,j]==2 ~  c( "T", "A", "A", "A" ), # copy2
                                                                                tmpG[1,j]==1 ~  c( "A", "T", "A", "A"), # copy1
                                                                                tmpG[1,j]==3 ~  c( "A", "A", "T", "A"), # copy3
                                                                                tmpG[1,j]==4 ~  c( "A", "A", "A", "T"), #copy4
                                                                                )


            C4_insert_for_bgl[ grep("C4A_copy", C4_insert_for_bgl$V1), a ] <-  case_when(
                                                                                tmpG[2,j]==2 ~  c( "T", "A", "A", "A" ), # copy2
                                                                                tmpG[2,j]==0 ~  c( "A", "T", "A", "A"), # copy0
                                                                                tmpG[2,j]==1 ~  c( "A", "A", "T", "A"), # copy1
                                                                                tmpG[2,j]==3 ~  c( "A", "A", "A", "T"), # copy3
                                                                                )

            C4_insert_for_bgl[grep("C4B_copy", C4_insert_for_bgl$V1), a ] <-  case_when(
                                                                                tmpG[3,j]==2 ~  c( "T", "A", "A", "A" ),
                                                                                tmpG[3,j]==0 ~  c( "A", "T", "A", "A"),
                                                                                tmpG[3,j]==1 ~  c( "A", "A", "T", "A"),
                                                                                tmpG[3,j]==3 ~  c( "A", "A", "A", "T"),
                                                                                )
            C4_insert_for_bgl[ grep("HERV_copy", C4_insert_for_bgl$V1), a ] <-  case_when(
                                                                                tmpG[4,j]==2 ~  c( "T", "A", "A", "A", "A" ),  #copy2
                                                                                tmpG[4,j]==0 ~  c( "A", "T", "A", "A", "A"),   #copy0
                                                                                tmpG[4,j]==1 ~  c( "A", "A", "T", "A", "A"),   #copy1
                                                                                tmpG[4,j]==3 ~  c( "A", "A", "A", "T", "A"),   #copy3
                                                                                tmpG[4,j]==4 ~  c( "A", "A", "A", "A", "T"),   #copy4
                                                                                )
             
            # C4_insert_for_bgl[ grep("C4_copy", C4_insert_for_bgl$V1), a ] <-  case_when(
            #                                                                     tmpG[1,j]==0 ~  c( "P", "A", "A", "A" ), # copy2
            #                                                                     tmpG[1,j]==1 ~  c( "A", "P", "A", "A"), # copy1
            #                                                                     tmpG[1,j]==2 ~  c( "A", "A", "P", "A"), # copy3
            #                                                                     tmpG[1,j]==3 ~  c( "A", "A", "A", "P"), #copy4
            #                                                                     )


            # C4_insert_for_bgl[ grep("C4A_copy", C4_insert_for_bgl$V1), a ] <-  case_when(
            #                                                                     tmpG[2,j]==0 ~  c( "P", "A", "A", "A" ), # copy2
            #                                                                     tmpG[2,j]==1 ~  c( "A", "P", "A", "A"), # copy0
            #                                                                     tmpG[2,j]==2 ~  c( "A", "A", "P", "A"), # copy1
            #                                                                     tmpG[2,j]==3 ~  c( "A", "A", "A", "P"), # copy3
            #                                                                     )

            # C4_insert_for_bgl[grep("C4B_copy", C4_insert_for_bgl$V1), a ] <-  case_when(
            #                                                                     tmpG[3,j]==0 ~  c( "P", "A", "A", "A" ),
            #                                                                     tmpG[3,j]==1 ~  c( "A", "P", "A", "A"),
            #                                                                     tmpG[3,j]==2 ~  c( "A", "A", "P", "A"),
            #                                                                     tmpG[3,j]==3 ~  c( "A", "A", "A", "P"),
            #                                                                     )
            # C4_insert_for_bgl[ grep("HERV_copy", C4_insert_for_bgl$V1), a ] <-  case_when(
            #                                                                     tmpG[4,j]==0 ~  c( "P", "A", "A", "A", "A" ),  #copy2
            #                                                                     tmpG[4,j]==1 ~  c( "A", "P", "A", "A", "A"),   #copy0
            #                                                                     tmpG[4,j]==2 ~  c( "A", "A", "P", "A", "A"),   #copy1
            #                                                                     tmpG[4,j]==3 ~  c( "A", "A", "A", "P", "A"),   #copy3
            #                                                                     tmpG[4,j]==4 ~  c( "A", "A", "A", "A", "P"),   #copy4
            #                                                                     )

            a<-a+1
        }
}


C4_insert_for_bgl[1:18,1:10]
dim(C4_insert_for_bgl)
bgl_hla[1:10,1:10]
dim(bgl_hla)

is.na(C4_insert_for_bgl) %>% table()

# NA -> 0으로 바꿔줌.
C4_insert_for_bgl[is.na(C4_insert_for_bgl) ] <- 0
dim(C4_insert_for_bgl)
dim(bgl_hla)

#############################################

#rm(list=ls())
#load("/kimlab_wd/yuo1996/C4_HLA_ref/try1/3.phasing_with_Beagle_C4_HLA/C4_bgl_insert.RData")

# 1. make marker file : insert C4 portion
output.markers <- read.table( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/2.MakeReference_with_HLA/", opt$snpver,"/integrated.markers_redefined"), header=F, sep=" ")
output.markers$V2 <- as.numeric(output.markers$V2)
head(output.markers); nrow(output.markers)

C4_markers$V2 <- c( rep( 31980581 ,12) , rep(31984684 ,  5) )
C4_markers$V3 <- "A"
C4_markers$V4 <- "T"
colnames(C4_markers)[1] <- "V1"
# C4 position 도 안겹치게 해줘야.
C4_markers$V2[1:12] <- 31980581:31980592
C4_markers$V2[13:17] <- 31984684:31984688

whole_markers <- rbind(C4_markers, output.markers)  %>% arrange( V2, V1) 
rownames(whole_markers) <- 1:nrow(whole_markers)
nrow(whole_markers)
unique(whole_markers$V2) %>% length()

write.table(whole_markers, paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "QUAL",opt$qual, "_", opt$haplonetver, "/whole.markers" ), col.names=F, row.names=F, quote=F, sep=" ")

whole_markers <- read.table(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "QUAL", opt$qual,  "_", opt$haplonetver,"/whole.markers"), header=F)
nrow(whole_markers) # 13431
head(whole_markers)
tail(whole_markers)


# 2. make bgl file : insert C4 portion
dim(bgl_hla)
dim(C4_insert_for_bgl)
bgl_hla[1:10,1:10]
C4_insert_for_bgl[1:18,1:10]
a <- bgl_hla[2,] %>% as.matrix() %>% as.vector() %>% .[-c(1:2) ]
b <- C4_insert_for_bgl[1,] %>%  as.matrix() %>% as.vector() %>% .[-1]
table( a==b )  # 446 *2 =892 / 441*2 =882

table( bgl_hla[ 6:nrow(bgl_hla) , 2] == output.markers$V1 ) #  원래 갯수 63777 -> 62573 ???
nrow(C4_markers)# 17 더해질 

C4_insert_for_bgl <- mutate( C4_insert_for_bgl, Vs="M" ) %>% relocate(Vs)
dim(C4_insert_for_bgl)
C4_insert_for_bgl[1:10,1:10]
dim(bgl_hla)
colnames(C4_insert_for_bgl) <- colnames(bgl_hla)


# whole_bgl start!!!
whole_bgl <- rbind(bgl_hla , C4_insert_for_bgl[-1,]) 
colnames(whole_markers)[1:2] <- c("V2", "pos")
whole_bgl[1:10,1:10]
nrow(whole_bgl)

whole_bgl <- merge( whole_bgl, whole_markers[1:2] , by="V2", all= TRUE  ) %>% arrange(  pos, V2)
whole_bgl <- whole_bgl  %>%  select(!pos) %>% relocate( V1) 

h  <- whole_bgl[c( which(whole_bgl$V2=="pedigree"), 
                                which(whole_bgl$V2=="id"),
                                which(whole_bgl$V2=="father"),
                                which(whole_bgl$V2=="mother"),
                                which(whole_bgl$V2=="gender"))  , ]
whole_bgl  <- whole_bgl[ -c( which(whole_bgl$V2=="pedigree"), 
                                which(whole_bgl$V2=="id"),
                                which(whole_bgl$V2=="father"),
                                which(whole_bgl$V2=="mother"),
                                which(whole_bgl$V2=="gender"))  , ]

whole_bgl <- rbind(h , whole_bgl)
rownames(whole_bgl) <- 1:nrow(whole_bgl)
nrow(whole_bgl)

table( whole_bgl$V2[-c(1:5)]==whole_markers$V2 ) #???????????????왜 안맞냐 127개??

whole_bgl[1:10,1:10]
dim(whole_bgl)


table( is.na(whole_bgl) )

######## sample QC step by C4/HERV/AB_QUAL QC -> JUST REMOVE UNDER opt$qual

total_table_uni <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/3rd_ver/total_table.txt", header=T )
total_table_kih <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/kih/total_table.txt", header=T )
total_table=rbind( total_table_uni, total_table_kih)
nrow(total_table)
filter(total_table, SAMPLE=="NIH21D7020955")

rmsams=filter(total_table, CNQ_C4 < opt$qual | AB_QUAL < opt$qual | CNQ_HERV < opt$qual)$SAMPLE 
length(rmsams)

head(total_table)
whole_bgl[1:10,1:10]
whole_bgl[1,] %>% as.matrix() %>% as.vector() %>% unique() %>% length()

whole_bgl = whole_bgl[ , -which(whole_bgl[1,] %in% rmsams)]
write.table( whole_bgl, paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/",  opt$region, "QUAL", opt$qual, "_", opt$haplonetver,"/whole.bgl"), col.names=F, row.names=F , sep=" ", quote=F)


# #####################################
# # dat, ped 파일에 C4 추가 (해서 plink로 바꾼다음 FRQ.frq 뽑는용도)

# setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_", opt$haplonetver) )
# path_whole.bgl <- "./whole.bgl" 
# system(paste0( "cat " , path_whole.bgl,  "| java -jar /kimlab_wd/yuo1996/tools/SNP2HLA_package_v1.0.3/MakeReference/beagle2linkage.jar whole"))

# ped <- read.table("./whole.ped", header=F)
# ped[1:10,1:10]
# ped$add1 <- -9
# ped <- ped %>% relocate(add1, .after=V5) 
# dim(ped)
# write.table(ped, "./whole.ped", col.names=F, row.names=F, sep=" ", quote=F)


# map <- read.table("./whole.markers", header=F)
# map <- data.frame( chr=6 , snp=map$V1 , c=0 , pos=map$V2) 
# write.table(map, "./whole.map", col.names=F, row.names=F ,quote=F, sep="\t")


# 직접 도려야댐.
# plink --file whole --make-bed 
# plink --bfile plink  --freq --out whole.FRQ # --keep-allele-order
