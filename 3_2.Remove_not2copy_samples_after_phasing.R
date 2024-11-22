rm(list=ls())
library(dplyr);library(tidyr)
library(optparse)

##################################  arguement setting ########################################
args_list <- list(
    make_option( "--snpver" , type="character" , default= "ver3" , help="snpversion", metavar="character"),
    make_option( "--region" , type="character" , default= "start31.7end32.2" , help="haplonet_region", metavar="character"),
    make_option( "--haplonetver" , type="character" , default= "haplonet_RMSE_under0.9" , help="snpversion", metavar="character")
)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)

# opt=list(); opt$snpver="ver3"; opt$region="start31.7end32.2"; opt$haplonetver="ver1_best"
# opt=list(); opt$snpver="ver3"; opt$region="start31.7end32.2QUAL30"; opt$haplonetver="ver2_best"

###############################################################################################


setwd( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_", opt$haplonetver) )

whole.phased <- read.table("./whole.eagle.phased.bgl", header=F)
whole.phased[1:10,1:10]

output.test <- whole.phased

# make_imputed_C4_table <- function( output.test ) {

##################### 1. C4 accuracy ()

res_C4 <- output.test[c(1,1, grep("copy", output.test$V2) ), ]
sam    <- ( dim( res_C4 )[2]-2 )/2

C4_CN <- data.frame( SAMPLE= NA, C4 = NA, C4A= NA, C4B  = NA, HERV =NA)
C4_CN_for2copy_check <- data.frame( SAMPLE= NA, C4 = NA, C4A= NA, C4B  = NA, HERV =NA)


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
head(C4_CN)

table( C4_CN$C4 ==C4_CN$C4A + C4_CN$C4B )
table( C4_CN$C4 >= C4_CN$HERV )



C4_CN_for2copy_check <- C4_CN_for2copy_check[-1,]
rownames(C4_CN_for2copy_check) <- 1:nrow(C4_CN_for2copy_check)

#return(list( C4_CN , C4_CN_for2copy_check) )
#}


total_C4_CN <- C4_CN
total_C4_CN_for2copy_check <- C4_CN_for2copy_check
total_C4_CN_for2copy_check[2:5] %>% as.matrix() %>% as.vector() %>% table

total_C4_CN_for2copy_check [ which(total_C4_CN_for2copy_check[2:5] !="1_1" ,arr.ind=TRUE)[,1] , ]
not2copy_samples <- total_C4_CN_for2copy_check [ which(total_C4_CN_for2copy_check[2:5] !="1_1" ,arr.ind=TRUE)[,1] , ]$SAMPLE
not2copy_samples

total_table_uni <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/3rd_ver/total_table.txt", header=T )
total_table_kih <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/kih/total_table.txt", header=T )
total_table = rbind( total_table_uni, total_table_kih)

filter(total_table, SAMPLE %in% not2copy_samples)

# A+B =C4 & cnC4 < cnHERV?? -> 이런얘들도 다제거???
head(C4_CN)
table( C4_CN$C4==(C4_CN$C4A+C4_CN$C4B)  )
table( C4_CN$C4 >= (C4_CN$HERV)  )

not_matched_sum <- C4_CN[ which( C4_CN$C4!=(C4_CN$C4A+C4_CN$C4B)  ), ]$SAMPLE
not_matched_herv <- C4_CN[ which( C4_CN$C4 < C4_CN$HERV  ), ]$SAMPLE

length(not_matched_sum);length(not_matched_herv)

rmsamples <- union( union(not_matched_herv, not_matched_sum) , not2copy_samples) #%>% length()
length(rmsamples)
rmsamples

filter(total_table, SAMPLE %in% rmsamples)

if(length(rmsamples)==0){
    whole.phased.sampleQC=whole.phased
}else{
    whole.phased.sampleQC <- whole.phased[, -which(whole.phased[1, ] %in% rmsamples  )]

}

dim(whole.phased)
dim(whole.phased.sampleQC)


getwd()

rmsamples
survived_sams <- unique( as.vector(as.matrix( whole.phased.sampleQC[1,] ) )[-c(1:2)] ) 
head(survived_sams)

getwd()
write.table(rmsamples, "./rmsamples.txt", col.names=F, row.names=F, quote=F, sep="\n")
write.table(survived_sams, "./survived_sams.txt", col.names=F, row.names=F, quote=F, sep="\n")

# 나중에 70000명 imputation을 위해서 sampleQC버전도 저장해놔야.
write.table( whole.phased.sampleQC , "./whole.phased.sampleQC.bgl" , col.names=FALSE, row.names=FALSE,quote=F)


# HLA+C4 variants for masking (removing in test panel)
getwd()
markers <- read.table("./whole.markers", header=F)

masking <- markers[ grep("HLA|C4|HERV|SNP|AA", markers$V1 ), ] 
masked_alleles <- data.frame(chr=6, pos=masking$V2)

write.table(masked_alleles, "./masked_alleles.txt", col.names=F, row.names=F , quote=F, sep="\t")



# whole.phased.sampleQC[1:10,1:10]
# whole.phased.sampleQC$V2

# filter(C4_CN, SAMPLE %in% not_matched_herv)
# filter(total_table , SAMPLE %in% c("KOREA1K-15"))
# filter(total_table , SAMPLE %in% not_matched_herv)

# 845-87

# # 2copy call안된 얘들
# setwd( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try5/7.Accuracy/", opt$snpver ))
# total_C4_CN_for2copy_check[2:5] %>% as.matrix() %>% as.vector() %>% table
# non2copyimputed <- filter( total_C4_CN_for2copy_check , C4!="1_1" | C4A!="1_1"  | C4B!="1_1"  | HERV!="1_1" )
# write.table(non2copyimputed , "./non2copy_imputed", col.names=T, row.names=F, sep="\t", quote=F)
# # 얘들이 대부분 밑에 C4A+B=C4안되는 얘들임 . -> 걍 제거?? 





######################### not 2copy samples from HLA -> 있을수가 있나? d없음. 걍확인용
whole.phased[1:10,1:10]
setwd( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_", opt$haplonetver) )

whole.phased <- read.table("./whole.eagle.phased.bgl", header=F)
whole.phased[1:10,1:10]
grep("HLA", whole.phased$V2, value=T)

output.test <- whole.phased
f="2field"
f="Gfield"

make_HLA_imputed <- function( output.test , f) {

        res_HLA<- output.test[c(1,1, grep("HLA", output.test$V2) ),]

field_info= data.frame( V2=res_HLA$V2, field=NA)
field_info[, 2] = sapply( res_HLA$V2, function(x) strsplit(x, ":" )[[1]] %>% length() )  
field_info[which(field_info$V2=="HLA_DPA1*02:07:01G"), 2]="2"
field_info[1:2, 2] <- NA
field_info[ grep("gg", field_info$V2), "field" ] <- "G"
field_info$field = paste0(field_info$field , "field")
field_info[1:2, 2] <- NA
res_HLA_f = res_HLA[c(1:2, which(res_HLA$V2 %in% filter(field_info, field==f)$V2 ) ) , ] 
res_HLA <- res_HLA_f
res_HLA_f[1:10,1:10]

        sam <-  ( dim( res_HLA )[2]-2  )/2
        HLA_imputed <- data.frame()

        # take some minuntes
        i=1
        for(i in 1:sam) {
            i=2*i+1

            tmp <- res_HLA[ -1, c(2, i:(i+1) )]
            colnames(tmp) <- tmp[1,]
            tmp <- tmp[-1,]

            name <- colnames(tmp)[2]

            colnames(tmp)[2:3] <- c("GT1", "GT2")
            tmp$HLA <- sapply( tmp$id , function(x) strsplit(x , "*", fixed=TRUE )[[1]][1] )
            #tmp$HLA <- gsub("HLA_", "", tmp$id) %>% gsub( "_[a-zA-Z0-9]{1,}", "", .)  %>% paste0("HLA_", .)
            tmp_P <- tmp[which(tmp$GT1 =="T" |tmp$GT2 =="T"), ]


            for( h in unique(tmp_P$HLA) ) {

                #HLA_tmp <- NULL            
                #HLA_tmp$ID <- name
                #HLA_tmp <- rbind( data.frame(HLA_tmp) , data.frame(HLA_tmp) )
                
                HLA_tmp = data.frame(ID=name)
                #filter(tmp_P, HLA==h)[2:3]  어떻게  2copy 확인 가능?????
                # P를 3개이상잡거나한건 어떠헥되지? 확인되도록 코드 바꿔야할 듯..?
                HLA_tmp$imputed_allele_num =filter(tmp_P, HLA==h)[c("GT1", "GT2")] %>% as.matrix() %>% as.vector() %>% gsub( "T", 1, . ) %>% gsub( "A", 0, .) %>% as.numeric() %>% sum()
                
                HLA_tmp <- cbind(HLA_tmp, data.frame( GT= filter(tmp_P, HLA==h, GT1=="T"|GT2=="T")$id ) )

                HLA_tmp$type <- h     
                HLA_imputed <- rbind(HLA_imputed , HLA_tmp)
            }

            
        }

        return(HLA_imputed)
}



total_HLA_imputed_1f <- make_HLA_imputed( whole.phased , "1field")  
head(total_HLA_imputed_1f)
# 2copy check after imputation
total_HLA_imputed_1f %>% select( ID, imputed_allele_num , type) %>% unique() %>% filter(imputed_allele_num!=2)
write.table(total_HLA_imputed_1f , "./total_HLA_phased_1f",  col.names = T, row.names=F, sep="\t", quote=F)

total_HLA_imputed_2f <- make_HLA_imputed( whole.phased , "2field")  
head(total_HLA_imputed_2f)
filter(total_HLA_imputed_2f , type=="HLA_DPA1")$GT %>% unique()

total_HLA_imputed_2f %>% select( ID, imputed_allele_num , type) %>% unique() %>% filter(imputed_allele_num!=2) 
filter(total_HLA_imputed_2f , imputed_allele_num!=2)

write.table(total_HLA_imputed_2f , "./total_HLA_phased_2f",  col.names = T, row.names=F, sep="\t", quote=F)


total_HLA_imputed_Gf <- make_HLA_imputed( whole.phased , "Gfield")  

write.table(total_HLA_imputed_Gf , "./total_HLA_phased_Gf",  col.names = T, row.names=F, sep="\t", quote=F)
head(total_HLA_imputed_Gf)
total_HLA_imputed_Gf %>% select( ID, imputed_allele_num , type) %>% unique() %>% filter(imputed_allele_num!=2)  
total_HLA_imputed_Gf %>% select( ID, imputed_allele_num , type) %>% unique() %>% filter(imputed_allele_num > 2)  
total_HLA_imputed_Gf  %>% filter(imputed_allele_num != 2)  
total_HLA_imputed_Gf  %>% filter(imputed_allele_num > 2)  %>% select(GT, type) %>% unique()
filter(total_HLA_imputed_Gf, ID %in% filter(total_HLA_imputed_2f , imputed_allele_num!=2)$ID , type=="gg_HLA_DPA1")

# missing 하나도 없는거 맞는지 확인!!!!!!!!!!!!1





# rmsamples_hla = total_HLA_imputed_Gf %>% filter( !imputed_allele_num %in% 1:2)  %>% .$ID %>% unique()

# c4_rmsamples <- read.table("./rmsamples.txt", header=F)

# union(rmsamples, c4_rmsamples[,1])

# total_HLA_imputed_Gf %>% filter( imputed_allele_num !=2 )  %>% .$ID %>% unique()


# rmsamples = union(rmsample_c4, rmsamples_hla)

# whole.phased.sampleQC <- whole.phased[, -which(whole.phased[1, ] %in% rmsamples  )]
# whole.phased.sampleQC[1:10,1:10]

# dim(whole.phased)
# dim(whole.phased.sampleQC)
# rmsamples

# survived_sams <- unique( as.vector(as.matrix( whole.phased.sampleQC[1,] ) )[-c(1:2)] ) 
# head(survived_sams)

# write.table(rmsamples, "./rmsamples.txt", col.names=F, row.names=F, quote=F, sep="\n")
# write.table(survived_sams, "./survived_sams.txt", col.names=F, row.names=F, quote=F, sep="\n")

# # 나중에 70000명 imputation을 위해서 sampleQC버전도 저장해놔야.
# write.table( whole.phased.sampleQC , "./whole.phased.sampleQC.bgl" , col.names=FALSE, row.names=FALSE,quote=F)
