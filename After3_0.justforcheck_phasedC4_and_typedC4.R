# samples which are assinged as missing are below
# 1. phred scale score < 5 for C4, A/B, HERV 
# 2. disconcordant genotypes in diploid cn level between haplonet and typed C4

rm(list=ls())
library(dplyr);library(tidyr)
library(optparse)
library(stringr)

##################################  arguement setting ########################################
args_list <- list(
    make_option( "--snpver" , type="character" , default= "ver3" , help="snpversion", metavar="character"),
    make_option( "--region" , type="character" , default= "start31.7end32.2" , help="haplonet_region", metavar="character"),
    make_option( "--crossval" , type="character" , default= "crossval50" , help="snpversion", metavar="character"),
    make_option( "--haplonetver" , type="character" , default= "haplonet_RMSE_under0.9" , help="snpversion", metavar="character")

)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)

# opt=list(); opt$snpver="ver3"; opt$region="start31.7end32.2"; opt$haplonetver="ver1_best"; opt$crossval="crossval1"
# opt=list(); opt$snpver="ver3"; opt$region="start31.7end32.2QUAL30"; opt$haplonetver="ver2_best"; opt$crossval="crossval1"
# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2QUAL30"; opt$haplonetver="ver1_best"; opt$crossval="crossval1"
# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8QUAL20"; opt$haplonetver="ver1_best"; opt$crossval="crossval1"

###############################################################################################
setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_", opt$haplonetver) )
getwd()

#### before phasing: for 0-assigned sample extraction

whole_bgl.bfphasing <- read.table( "./whole.bgl" , header=F , sep=" " )

# # missing(0) 처리된 HLA 확인해범
# whole_bgl.bfphasing[1:10,1:10]
# whole_bgl.bfphasing$V2[-c(1:5)][ apply( whole_bgl.bfphasing[ -c(1:5), -c(1:2) ] , 1, function(x) length( which( x == 0) ) > 0  ) ] %>% grep("HLA", . , value=TRUE )
# filter(whole_bgl.bfphasing, V2=="HLA_DQB1*06:04:01G") %>% as.matrix() %>% as.vector() %>% table()

extract_C4 <- function(whole_bgl){
                phasedC4 <- whole_bgl[c(1, grep( "C4|HERV" ,whole_bgl$V2 ) ), -1]
#                colnames(phasedC4) <- phasedC4[1,]
#                phasedC4 <- phasedC4[-1,]
                rownames(phasedC4) <- phasedC4$V2
                phasedC4 <- phasedC4[,-1]
                return(phasedC4)
}

phasedC4_bf <- extract_C4(whole_bgl.bfphasing)
phasedC4_bf[1:10,1:10]
phasedC4_bf[,940:950]
phasedC4_bf[ which(phasedC4_bf=="0", arr.ind=TRUE)] 
phasedC4_bf[ which(phasedC4_bf=="0", arr.ind=TRUE)] 

missing_idx <- which(phasedC4_bf=="0", arr.ind=TRUE) %>% as.data.frame()
if( nrow(missing_idx) !=0){
missing_idx$type <- gsub( "\\.[0-9]+", "", rownames(missing_idx) ) 
missing_idx$SAMPLE <- as.vector( as.matrix( phasedC4_bf[1,] ))[ missing_idx$col ]
rownames(missing_idx) <- 1:nrow(missing_idx)
missing_idx <- select(missing_idx, type, SAMPLE) %>% unique 
missing_idx$gene <- str_extract( missing_idx$type , "[[:alnum:]]+_") %>% gsub("_", "", .)
missing_idx <- select(missing_idx, gene, SAMPLE) %>% unique 

}
head(missing_idx, n=20)
nrow(missing_idx)

# number of missing-assigned samples
missing_idx$SAMPLE %>% unique %>% length()
missing_idx$gene %>% table() # by the genes


#### after phasing: phased C4 (technically imputed C4 while phasing)
#whole_bgl.afphasing <- read.table( "./whole.bgl.phased.whole.bgl.phased" , header=F , sep=" " )
whole_bgl.afphasing <- read.table( "./whole.phased.sampleQC.bgl" , header=F , sep=" " )
whole_bgl.afphasing  <- rbind( whole_bgl.afphasing[1,] , whole_bgl.afphasing )
whole_bgl.afphasing  <- rbind(matrix(0, nrow=3, ncol=dim(whole_bgl.afphasing)[2]), whole_bgl.afphasing)
whole_bgl.afphasing  <- whole_bgl.afphasing[ c(4:5, 1:3, 6:dim(whole_bgl.afphasing)[1]  ) , ]
whole_bgl.afphasing[1,1] <- "P"; whole_bgl.afphasing[1,2] <- "pedigree"


make_imputed_C4_table <- function( output.test ) {

res_C4 <- output.test[c(1,2, grep("copy", output.test$V2) ),  ]
sam    <-  ( dim( res_C4 )[2]-2  )/2

C4_CN  <- data.frame( SAMPLE= NA, C4 = NA, C4A= NA, C4B  = NA, HERV =NA )
C4_CN_for2copy_check <- data.frame( SAMPLE= NA, C4 = NA, C4A= NA, C4B  = NA, HERV =NA )


for(i in 1:sam) {
    i=2*i+1
    # check for 2_1 case
    # tmp <- res_C4[-1, c(2,which(res_C4[1,] =="KOREA1K-178"))]

    tmp <- res_C4[ -1, c(2, i:(i+1) )]
    colnames(tmp) <- tmp[1,]
    tmp <- tmp[-1,]
    
    CN <- c()
    P <- c()

    for (k in c("C4_", "C4A_", "C4B_", "HERV")){

        indv <- tmp[grep(k , tmp$id),]        
        indv$cn <- as.numeric( substr( indv$id , nchar(indv$id) , nchar(indv$id) ) )


        chr1_P_num <- length(  indv[ grep("T", indv[,2]), "cn" ]  )
        if( chr1_P_num == 1 ){
            cn1 <- indv[ grep("T", indv[,2]), "cn" ]
        }else if(chr1_P_num == 0){
            cn1 <-0
        } else{
            cn1 <- sum( indv[ grep("T", indv[,2]), "cn" ] ) # 2_1 인 경우  sum으로 해놨음 !!!...
        }


        chr2_P_num <- length(  indv[ grep("T", indv[,3]), "cn" ]  )
        if( chr2_P_num == 1 ){
            cn2 <- indv[ grep("T", indv[,3]), "cn" ]
        }else if (chr2_P_num==0){
            cn2 <- 0
        }else{
            cn2 <- sum( indv[ grep("T", indv[,3]), "cn" ] )
        }
            
        
        total <- cn1+cn2
        CN <- c( CN, total)

        P_num <- paste0( c( chr1_P_num , chr2_P_num) , collapse="_" )
        P <- c( P , P_num)
 
    }

    t <- c( colnames(tmp)[2] , CN)  %>% data.frame() %>% t()
    colnames(t) <- c("SAMPLE", "C4", "C4A", "C4B", "HERV")
    C4_CN <- rbind(C4_CN, t)

    pa <- c( colnames(tmp)[2] , P)  %>% data.frame() %>% t()
    colnames(pa) <- c("SAMPLE", "C4", "C4A", "C4B", "HERV")
    C4_CN_for2copy_check<- rbind(C4_CN_for2copy_check, pa)
    
}


C4_CN <- C4_CN[-1,]
rownames(C4_CN) <- 1:nrow(C4_CN)
C4_CN$C4 <- as.numeric(C4_CN$C4)
C4_CN$C4A <- as.numeric(C4_CN$C4A)
C4_CN$C4B <- as.numeric(C4_CN$C4B)
C4_CN$HERV <- as.numeric(C4_CN$HERV)
C4_CN

C4_CN_for2copy_check <- C4_CN_for2copy_check[-1,]
rownames(C4_CN_for2copy_check) <- 1:nrow(C4_CN_for2copy_check)

return(list( C4_CN , C4_CN_for2copy_check) )
}


imputedC4 <- make_imputed_C4_table( whole_bgl.afphasing )[[1]]
nrow(imputedC4)

#### extract bf samples from af

head(missing_idx)
missing_idx$SAMPLE %>% unique %>% length() #304
head(imputedC4); dim(imputedC4)

C4_missing <- filter(missing_idx, gene=="C4")
C4A_missing <- filter(missing_idx, gene=="C4A")
C4B_missing <- filter(missing_idx, gene=="C4B")
HERV_missing <- filter(missing_idx, gene=="HERV")

C4_imputed <- imputedC4[ which(imputedC4$SAMPLE %in% C4_missing$SAMPLE), c("SAMPLE", "C4")] %>% rename(imputed=C4)
C4A_imputed <- imputedC4[ which(imputedC4$SAMPLE %in% C4A_missing$SAMPLE), c("SAMPLE", "C4A")] %>% rename(imputed=C4A)
C4B_imputed <- imputedC4[ which(imputedC4$SAMPLE %in% C4B_missing$SAMPLE), c("SAMPLE", "C4B")] %>% rename(imputed=C4B)
HERV_imputed <- imputedC4[ which(imputedC4$SAMPLE %in% HERV_missing$SAMPLE), c("SAMPLE", "HERV")] %>% rename(imputed=HERV)

imputed_C4AB <- merge(C4A_imputed, C4B_imputed, by="SAMPLE") %>% rename(C4A=imputed.x) %>% rename(C4B=imputed.y)
imputed_HERV <- HERV_imputed %>% rename( HERV=imputed)
nrow(imputed_HERV); nrow(imputed_C4AB)
write.table(imputed_C4AB, "./imputed_C4AB", col.names=T, row.names=F, sep="\t", quote=F)
write.table(imputed_HERV, "./imputed_HERV", col.names=T, row.names=F, sep="\t", quote=F)




# ##############################################################################################
# # typed C4 retrive
# pat <- "/kimlab_wd/yuo1996/C4_analysis/4.terra_result/3rd_ver/"
# bat <-paste0("batch", 1:9) 

# total_table <- data.frame()
# for( i in 1:length(bat)){
#     table <- read.table(  paste0( paste0(pat, bat)[i] , "/C4_hg38.C4_table.txt") , header=T)
#     table$batch <- bat[i]

#     total_table <- rbind( total_table, table)
# }


# C4_compare <- merge(total_table[c("SAMPLE", "C4", "CNQ_C4")], C4_imputed, by="SAMPLE" )  # 죄다 4,5,6
# C4A_compare <- merge(total_table[c("SAMPLE", "C4A", "AB_QUAL")], C4A_imputed, by="SAMPLE" ) 
# C4B_compare <- merge(total_table[c("SAMPLE", "C4B", "AB_QUAL")], C4B_imputed, by="SAMPLE" ) 
# HERV_compare <- merge(total_table[c("SAMPLE", "HERV", "CNQ_HERV")], HERV_imputed, by="SAMPLE" ) 



# table( C4_compare$C4 == C4_compare$imputed )
# table( C4A_compare$C4A == C4A_compare$imputed )
# table( C4B_compare$C4B == C4B_compare$imputed )
# table( HERV_compare$HERV == HERV_compare$imputed )

# c( filter( C4_compare, C4==imputed )$SAMPLE , 
# filter( C4A_compare, C4A==imputed )$SAMPLE,
# filter( C4B_compare, C4B==imputed )$SAMPLE,
# filter( HERV_compare, HERV==imputed )$SAMPLE ) %>% unique() %>% length() 
# # 304명 중에 typedC4랑 4개 유전자중 하나라도 일치하는 80명만데려갈까? 그럼 오버피팅 피할수잇나?
# # 근데그러면 phsing한번 더해야. 

# intersect( filter( C4_compare, C4==imputed )$SAMPLE, filter( C4A_compare, C4A==imputed )$SAMPLE ) %>% intersect(. ,filter( C4B_compare, C4B==imputed )$SAMPLE)


# #intersect( unique( missing_idx$SAMPLE ), multiple_aligned_sams[,1] )
# #multiple_aligned_sams <- read.table("/kimlab_wd/yuo1996/C4_korean1K/finish/batch11_98/batch11_98.list" , header=F)
# #head(multiple_aligned_sams)


# ###### 결과: typed C4랑 대부분 다름-> 혹시 haplonet_proximal 결과랑 일치하는지 함 볼까? 

# haplonet_path <- "/kimlab_wd/yuo1996/C4_HLA_ref/haplonet_proximal/4.haplonet/ver3/start31.7end32.2/genotpying_forcrossval/"
# C4_path <- paste0( haplonet_path, grep("C4.txt", list.files( haplonet_path ) ,value=T))
# C4A_path <- paste0( haplonet_path, grep("C4A.txt", list.files( haplonet_path ) ,value=T))
# C4B_path <- paste0( haplonet_path, grep("C4B.txt", list.files( haplonet_path ) ,value=T))
# HERV_path <- paste0( haplonet_path, grep("HERV.txt", list.files( haplonet_path ) ,value=T))

# #C4
# geno_prob_C4 <- read.table( C4_path , header=T)
# check_indiploid <- merge( C4_compare , geno_prob_C4[c("SAMPLE" , "estimated" )]  ) 
# table( check_indiploid$imputed==check_indiploid$estimated )

# #C4A
# geno_prob_C4A <- read.table( C4A_path , header=T)
# check_indiploid <- merge( C4A_compare , geno_prob_C4A[c("SAMPLE" , "estimated" )]  ) 
# table( check_indiploid$imputed==check_indiploid$estimated )

# #C4B
# geno_prob_C4B <- read.table( C4B_path , header=T)
# check_indiploid <- merge( C4B_compare , geno_prob_C4B[c("SAMPLE" , "estimated" )]  ) 
# table( check_indiploid$imputed==check_indiploid$estimated )

# #HERV
# geno_prob_HERV <- read.table( HERV_path , header=T)
# check_indiploid <- merge( HERV_compare , geno_prob_HERV[c("SAMPLE" , "estimated" )]  ) 
# table( check_indiploid$imputed==check_indiploid$estimated )

