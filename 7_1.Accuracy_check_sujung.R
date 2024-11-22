rm(list=ls())
library(dplyr);library(tidyr)
library(optparse)
library(stringr)

##################################  arguement setting ########################################
args_list <- list(
    make_option( "--snpver" , type="character" , default= "ver3" , help="snpversion", metavar="character"),
    make_option( "--region" , type="character" , default= "start31.7end32.2" , help="haplonet_region", metavar="character"),
    make_option( "--crossval" , type="character" , default= "crossval50" , help="snpversion", metavar="character"),
    make_option( "--haplonetver" , type="character" , default= "crossval50" , help="snpversion", metavar="character")
)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)


# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8QUAL20"; opt$haplonetver="ver1_sujung"; opt$crossval="crossval1"

###############################################################################################


if(!dir.exists( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver) ) ){ dir.create(path=paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver) ) } 
if(!dir.exists( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver) ) ){ dir.create(path=paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver) ) } 
if(!dir.exists( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver,"/", opt$region, "_", opt$haplonetver, "/", opt$crossval) )){ dir.create(path=paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/", opt$crossval)) } 

setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval) )


# retrive all sets
imputed_list <- grep("vcf$", list.files("./") , value=TRUE)
#imputed_list <- grep("phased$", list.files("./") , value=TRUE)
head(imputed_list); length(imputed_list)

for(i in imputed_list){ 
    output.test <- read.delim( paste0( "./", i ) , header=T, sep="\t", skip=13) # impute5
#    output.test <- read.delim( paste0( "./", i ) , header=T, sep="\t", skip=11) # beagle
    num <- gsub("KOREA1K-","korea",i) %>% gsub("[^0-9]", "", . )
    assign( paste0("output.test", num), output.test , envir = .GlobalEnv )
}
ls()

# output.test <- output.test995

make_imputed_C4_table <- function( output.test ) {

##################### 1. C4 accuracy ()

#res_C4 <- output.test[c(1,2, grep("copy", output.test$V2) ),] # beagle? 
res_C4 <- output.test[ grep("copy", output.test$ID) ,] # impute5
id <- colnames( res_C4 )[dim(res_C4)[2]] %>% gsub( ".", "-", .,fixed=TRUE)
res_C4 <- rename(res_C4, sample=colnames( res_C4 )[dim(res_C4)[2]])
res_C4 <- separate( res_C4, sample, into=c("GT", "rest") , sep=":") %>% select( ID, GT ) 
res_C4$GT <- gsub("|", "abc", res_C4$GT ,fixed=TRUE)
res_C4 <- separate( res_C4, GT, into=c("V3","V4") , sep="abc" )
res_C4$V3 <- ifelse( res_C4$V3==1, "P", "A"); res_C4$V4 <- ifelse( res_C4$V4==1, "P", "A")
res_C4 <- res_C4 %>% mutate( V1="M") %>% relocate( V1, .before=ID) %>% rename( V2=ID)
res_C4 <- rbind( data.frame( V1=c("P","I"), V2=c("pedigree", "id"), V3=id, V4=id) , res_C4)

sam <-  ( dim( res_C4 )[2]-2  )/2

C4_CN <- data.frame( SAMPLE= NA, C4A= NA, C4B  = NA, HERV =NA)
C4_CN_for2copy_check <- data.frame( SAMPLE= NA, C4A= NA, C4B  = NA, HERV =NA)


for(i in 1:sam) {
    i=2*i+1
    # check for 2_1 case
    # tmp <- res_C4[-1, c(2,which(res_C4[1,] =="KOREA1K-178"))]

    tmp <- res_C4[ -1, c(2, i:(i+1) )]
    colnames(tmp) <- tmp[1,]
    tmp <- tmp[-1,]
    
    CN <- c()
    P <- c()

    for (k in c("C4A_", "C4B_", "HERV")){

        indv <- tmp[grep(k , tmp$id),]        
        indv$cn <- as.numeric( substr( indv$id , nchar(indv$id) , nchar(indv$id) ) )


        chr1_P_num <- length(  indv[ grep("P", indv[,2]), "cn" ]  )
        if( chr1_P_num == 1 ){
            cn1 <- indv[ grep("P", indv[,2]), "cn" ]
        }else if(chr1_P_num == 0){
            cn1 <-0
        } else{
            cn1 <- sum( indv[ grep("P", indv[,2]), "cn" ] ) # 2_1 인 경우  sum으로 해놨음 !!!...
        }


        chr2_P_num <- length(  indv[ grep("P", indv[,3]), "cn" ]  )
        if( chr2_P_num == 1 ){
            cn2 <- indv[ grep("P", indv[,3]), "cn" ]
        }else if (chr2_P_num==0){
            cn2 <- 0
        }else{
            cn2 <- sum( indv[ grep("P", indv[,3]), "cn" ] )
        }
            
        
        total <- cn1+cn2
        CN <- c( CN, total)

        P_num <- paste0( c( chr1_P_num , chr2_P_num) , collapse="_" )
        P <- c( P , P_num)
 
    }

    t <- c( colnames(tmp)[2] , CN)  %>% data.frame() %>% t()
    colnames(t) <- c("SAMPLE",  "C4A", "C4B", "HERV")
    C4_CN <- rbind(C4_CN, t)

    pa <- c( colnames(tmp)[2] , P)  %>% data.frame() %>% t()
    colnames(pa) <- c("SAMPLE", "C4A", "C4B", "HERV")
    C4_CN_for2copy_check<- rbind(C4_CN_for2copy_check, pa)
    
}


C4_CN <- C4_CN[-1,]
rownames(C4_CN) <- 1:nrow(C4_CN)
C4_CN$C4A <- as.numeric(C4_CN$C4A)
C4_CN$C4B <- as.numeric(C4_CN$C4B)
C4_CN$HERV <- as.numeric(C4_CN$HERV)
C4_CN

C4_CN_for2copy_check <- C4_CN_for2copy_check[-1,]
rownames(C4_CN_for2copy_check) <- 1:nrow(C4_CN_for2copy_check)

return(list( C4_CN , C4_CN_for2copy_check) )
}

#total_C4_CN <- make_imputed_C4_table(output.test2 )[[1]]  %>% mutate( set = paste0("set", 1)) 
#total_C4_CN_for2copy_check <- make_imputed_C4_table(output.test2 )[[2]]  %>% mutate( set = paste0("set", 1)) 

num <- gsub("KOREA1K-","korea", imputed_list) %>% gsub("[^0-9]", "", . )
total_C4_CN <- data.frame()
for(o in paste0("output.test", num) ){
    tmp <- make_imputed_C4_table( get(o) )[[1]]  # %>% mutate( set = paste0("set", a)) 
    total_C4_CN <- rbind(total_C4_CN, tmp)
}
head(total_C4_CN)
tail(total_C4_CN)
nrow(total_C4_CN)

total_C4_CN_for2copy_check <- data.frame()
for(o in paste0("output.test", num) ){
    tmp <- make_imputed_C4_table( get(o) )[[2]]  # %>% mutate( set = paste0("set", a)) 
    total_C4_CN_for2copy_check <- rbind(total_C4_CN_for2copy_check, tmp)
}
nrow(total_C4_CN_for2copy_check)
tail(total_C4_CN_for2copy_check)

total_C4_CN_for2copy_check %>% head
total_C4_CN_for2copy_check %>% nrow()

# 2copy call안된 얘들
total_C4_CN_for2copy_check[2:4] %>% as.matrix() %>% as.vector() %>% table
total_C4_CN_for2copy_check[which(total_C4_CN_for2copy_check[,-1]!="1_1", arr.ind=TRUE)[,1],]
not2copy_samples <- total_C4_CN_for2copy_check [ which(total_C4_CN_for2copy_check[2:4] !="1_1" ,arr.ind=TRUE)[,1] , ]$SAMPLE
not2copy_samples <- not2copy_samples %>% unique 
length(not2copy_samples)

not2copy <- total_C4_CN_for2copy_check[total_C4_CN_for2copy_check$SAMPLE %in% not2copy_samples, ] 

getwd()
write.table(not2copy, "../not2copy", col.names=T, row.names = F,quote=F, sep="\t")


#non2copyimputed <- filter( total_C4_CN_for2copy_check , C4!="1_1" | C4A!="1_1"  | C4B!="1_1"  | HERV!="1_1" )
#setwd( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try5/7.Accuracy/", opt$snpver ))
#write.table(non2copyimputed , "./non2copy_imputed", col.names=T, row.names=F, sep="\t", quote=F)
# 얘들이 대부분 밑에 C4A+B=C4안되는 얘들임 . -> 걍 제거?? 

# imputed cn of C4A + C4B == C4?? 
table(  ( total_C4_CN$C4A + total_C4_CN$C4B ) >= total_C4_CN$HERV  )
total_C4_CN[ ( total_C4_CN$C4A + total_C4_CN$C4B ) < total_C4_CN$HERV ,] 

not_matched_herv <- total_C4_CN[ which( ( total_C4_CN$C4A + total_C4_CN$C4B ) < total_C4_CN$HERV ), ]$SAMPLE

length(not2copy_samples);length(not_matched_herv)


rmsamples <- union( not_matched_herv , not2copy_samples) #%>% length()
length(rmsamples)
write.table(rmsamples, "../rmsamples", col.names=F, row.names = F,quote=F)

# total_C4_CN <- total_C4_CN[ -which(total_C4_CN$SAMPLE %in% rmsamples),] # 제거 원하느 경우에만..?
nrow(total_C4_CN)

###############  retrive C4 typing from WGS terra
total_table_uni <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/3rd_ver/total_table.txt", header=T )
total_table_kih <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/kih/total_table.txt", header=T )
total_table=rbind( total_table_uni, total_table_kih)

# lowqulity typing phased된 정보로 바꿔서 concordance rate 구하기위함..!
# sampleQC 버전에선 제외
if( !grepl("QUAL", opt$region) ){
imputed_C4AB <- read.table( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_" , opt$haplonetver, "/imputed_C4AB") , header=T)
imputed_HERV <- read.table( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/",  opt$region, "_",opt$haplonetver, "/imputed_HERV") , header=T)

total_table <- merge(total_table, imputed_C4AB, by="SAMPLE", all.x=TRUE) %>% merge(. , imputed_HERV, by="SAMPLE", all.x=TRUE) 
total_table[ !is.na(total_table$C4A.y), "C4A.x"]  <- total_table[ !is.na(total_table$C4A.y), "C4A.y"] 
total_table[ !is.na(total_table$C4B.y), "C4B.x"]  <- total_table[ !is.na(total_table$C4B.y), "C4B.y"] 
total_table[ !is.na(total_table$HERV.y), "HERV.x"]  <- total_table[ !is.na(total_table$HERV.y), "HERV.y"] 

total_table <- rename(total_table, C4A=C4A.x) %>% rename(C4B=C4B.x) %>% rename(HERV=HERV.x) 
}


round( ( table( total_table$C4 )/nrow(total_table) ) * 100 , 2)
round( ( table( total_table$C4A )/nrow(total_table) ) * 100 , 2)
round( ( table( total_table$C4B )/nrow(total_table) ) * 100 , 2)
round( ( table( total_table$HERV )/nrow(total_table) ) * 100 , 2)


# diploid freq. dist.

(  table( total_C4_CN$C4 )/nrow(total_table)  ) * 100
(  table( total_C4_CN$C4A )/nrow(total_table)  ) * 100
(  table( total_C4_CN$C4B )/nrow(total_table)  ) * 100
(  table( total_C4_CN$HERV )/nrow(total_table)  ) * 100



####
typedC4 <- select(total_table, SAMPLE, C4, CNQ_C4, HERV, CNQ_HERV , C4A, C4B, AB_QUAL, batch) 
colnames(total_C4_CN)[2:4] <- paste0("imputed_" ,colnames(total_C4_CN)[2:4] )
vali <- merge(typedC4, total_C4_CN, by="SAMPLE")
head(vali)
nrow(vali)

vali$imputed_C4=vali$imputed_C4A+ vali$imputed_C4B

vali$C4_mat <- ifelse(vali$C4==vali$imputed_C4, "yes", "no" )
vali$C4A_mat <- ifelse(vali$C4A==vali$imputed_C4A, "yes", "no" )
vali$C4B_mat <- ifelse(vali$C4B==vali$imputed_C4B, "yes", "no" )
vali$HERV_mat <- ifelse(vali$HERV==vali$imputed_HERV, "yes", "no" )

failed <- filter(vali, C4_mat =="no"  |C4A_mat =="no"| C4B_mat =="no"|HERV_mat =="no")
success <- filter(vali, C4_mat =="yes" , C4A_mat =="yes" ,  C4B_mat =="yes", HERV_mat =="yes")

nrow(success)/(nrow(success)+nrow(failed))*100

nrow( filter(vali, C4_mat=="yes") )/ nrow(vali)
nrow( filter(vali, C4A_mat=="yes") )/ nrow(vali)
nrow( filter(vali, C4B_mat=="yes") )/ nrow(vali)
nrow( filter(vali, HERV_mat=="yes") )/ nrow(vali)

# non 2copy 제거하고 해보기
#vali <- vali[-which( vali$SAMPLE %in% non2copyimputed$SAMPLE) , ]

# concordance rate

concor_cal <- function( con_table  , allele_num , sam_num ){
    con_table <- as.vector(con_table)
    ac <- (con_table[2])/(sam_num*allele_num)*100
    return(ac)
}

to <- c( concor_cal(table( vali$C4_mat ) , 1 , nrow(vali)),
        concor_cal(table( vali$C4A_mat ) , 1 , nrow(vali)),
        concor_cal(table( vali$C4B_mat ) , 1 , nrow(vali)),
        concor_cal(table( vali$HERV_mat ) , 1 , nrow(vali)),
        concor_cal( table( c( vali$C4_mat, vali$C4A_mat,  vali$C4B_mat, vali$HERV_mat ))  ,  4,  nrow(vali)) 
)

to # crossval1 이면 여기서 끝. 
# crossval1 인경우만!!
to <- t( as.data.frame(to) )
colnames(to) <- c("C4", "C4A", "C4B", "HERV", "total")
to

setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval) )
write.table(to, "./C4_diploid_contable.txt", col.names=T,row.names=F, sep="\t", quote=F)
ls()
rm( list=grep("output.test", ls(), value=T))
save.image( "./val_C4_result.RDdata" )


######################################## 2. HLA accurarcy
rm(list=ls())
library(stringr)
library(dplyr);library(tidyr)
library(optparse)
library(stringr)

##################################  arguement setting ########################################
args_list <- list(
    make_option( "--snpver" , type="character" , default= "ver3" , help="snpversion", metavar="character"),
    make_option( "--region" , type="character" , default= "start31.7end32.2" , help="haplonet_region", metavar="character"),
    make_option( "--crossval" , type="character" , default= "crossval50" , help="snpversion", metavar="character"),
    make_option( "--haplonetver" , type="character" , default= "crossval50" , help="snpversion", metavar="character")
)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)

# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8QUAL20"; opt$haplonetver="ver1_sujung"; opt$crossval="crossval1"
###############################################################################################


setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/", opt$crossval))

# retrive all sets
imputed_list <- grep("vcf$", list.files("./") , value=TRUE)
#imputed_list <- grep("phased$", list.files("./") , value=TRUE)
head(imputed_list); length(imputed_list)

for(i in imputed_list){ 
    output.test <- read.delim( paste0( "./", i ) , header=T, sep="\t", skip=13) # impute5
#    output.test <- read.delim( paste0( "./", i ) , header=T, sep="\t", skip=11) # beagle
    num <- gsub("KOREA1K-","korea",i) %>% gsub("[^0-9]", "", . )
    assign( paste0("output.test", num), output.test , envir = .GlobalEnv )
}
ls()

output.test <- output.test100
# take some minutes
# f"1field"/"2field"/"3field"
# f="Gfield"
make_HLA_imputed <- function( output.test , f) {
        #res_C4 <- output.test[c(1,2, grep("HLA", output.test$V2) ),] # beagle?
        res_C4 <- output.test[ grep("HLA", output.test$ID) ,] # impute5
        id <- colnames( res_C4 )[dim(res_C4)[2]] %>% gsub( ".", "-", . ,fixed=TRUE)
        res_C4 <- rename(res_C4, sample=colnames( res_C4 )[dim(res_C4)[2]])
        res_C4 <- separate( res_C4, sample, into=c("GT", "rest") , sep=":") %>% select( ID, GT ) 
        res_C4$GT <- gsub("|", "abc", res_C4$GT ,fixed=TRUE)
        res_C4 <- separate( res_C4, GT, into=c("V3","V4") , sep="abc" )
        res_C4$V3 <- ifelse( res_C4$V3==1, "P", "A"); res_C4$V4 <- ifelse( res_C4$V4==1, "P", "A")
        res_C4 <- res_C4 %>% mutate( V1="M") %>% relocate( V1, .before=ID) %>% rename( V2=ID)
        res_C4 <- rbind( data.frame( V1=c("P","I"), V2=c("pedigree", "id"), V3=id, V4=id) , res_C4)        
        res_HLA <- res_C4
#       res_HLA$V2 <- gsub( "*", "_", res_HLA$V2, fixed=TRUE) %>% gsub(":","", . , fixed=TRUE)

head(res_HLA)
field_info= data.frame( V2=res_HLA$V2, field=NA)
field_info[, 2] = sapply( res_HLA$V2, function(x) strsplit(x, ":" )[[1]] %>% length() )  
field_info[1:2, 2] <- NA
field_info[ grep("gg", field_info$V2), "field" ] <- "G"
field_info$field = paste0(field_info$field , "field")
field_info[1:2, 2] <- NA
res_HLA_f = res_HLA[c(1:2, which(res_HLA$V2 %in% filter(field_info, field==f)$V2 ) ) , ] 
res_HLA <- res_HLA_f
head(res_HLA)

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
            tmp_P <- tmp[which(tmp$GT1 =="P" |tmp$GT2 =="P"), ]


            for( h in unique(tmp_P$HLA) ) {

                #HLA_tmp <- NULL            
                #HLA_tmp$ID <- name
                #HLA_tmp <- rbind( data.frame(HLA_tmp) , data.frame(HLA_tmp) )
                
                HLA_tmp = data.frame(ID=name)
                #filter(tmp_P, HLA==h)[2:3]  어떻게  2copy 확인 가능?????
                # P를 3개이상잡거나한건 어떠헥되지? 확인되도록 코드 바꿔야할 듯..?
                HLA_tmp$imputed_allele_num =filter(tmp_P, HLA==h)[c("GT1", "GT2")] %>% as.matrix() %>% as.vector() %>% gsub( "P", 1, . ) %>% gsub( "A", 0, .) %>% as.numeric() %>% sum()
                
                HLA_tmp <- cbind(HLA_tmp, data.frame( GT= filter(tmp_P, HLA==h, GT1=="P"|GT2=="P")$id ) )

                HLA_tmp$type <- h     
                HLA_imputed <- rbind(HLA_imputed , HLA_tmp)
            }
        }
        return(HLA_imputed)
}



#total_HLA_imputed <- rbind( make_HLA_imputed(output.test1) %>% mutate( set = paste0("set", 1)) )


num <- gsub("KOREA1K-","korea", imputed_list) %>% gsub("[^0-9]", "", . )
total_HLA_imputed_1f <- data.frame()
for(o in paste0("output.test", sort(num)) ){
    print(o)
    tmp <- make_HLA_imputed( get(o) , "1field")  
    total_HLA_imputed_1f <- rbind(total_HLA_imputed_1f, tmp)
}
head(total_HLA_imputed_1f, n=50)
# 2copy check after imputation
total_HLA_imputed_1f %>% select( ID, imputed_allele_num , type) %>% unique() %>% filter(imputed_allele_num!=2)
total_HLA_imputed_1f %>% filter(imputed_allele_num!=2)


total_HLA_imputed_2f <- data.frame()
for(o in paste0("output.test", sort(num)) ){
    print(o)
    tmp <- make_HLA_imputed( get(o) , "2field")  
    total_HLA_imputed_2f <- rbind(total_HLA_imputed_2f, tmp)
}
head(total_HLA_imputed_2f)
total_HLA_imputed_2f %>% select( ID, imputed_allele_num , type) %>% unique() %>% filter(imputed_allele_num!=2) #%>% select( type) %>% table()
total_HLA_imputed_2f %>% filter(imputed_allele_num!=2)
total_HLA_imputed_2f %>% filter(ID=="KOREA1K-479")


total_HLA_imputed_Gf <- data.frame()
for(o in paste0("output.test", sort(num)) ){
    print(o)
    tmp <- make_HLA_imputed( get(o) , "Gfield")  
    total_HLA_imputed_Gf <- rbind(total_HLA_imputed_Gf, tmp)
}
head(total_HLA_imputed_Gf)
total_HLA_imputed_Gf %>% select( ID, imputed_allele_num , type) %>% unique() %>% filter(imputed_allele_num!=2)  #%>% select( type) %>% table()
total_HLA_imputed_Gf %>% filter(imputed_allele_num!=2)
total_HLA_imputed_Gf %>% filter(imputed_allele_num > 2)  


write.table(total_HLA_imputed_1f , paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_imputed_1f") , col.names=T, row.names=F, sep="\t", quote=F)
write.table(total_HLA_imputed_2f , paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_imputed_2f") , col.names=T, row.names=F, sep="\t", quote=F)
write.table(total_HLA_imputed_Gf , paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_imputed_Gf") , col.names=T, row.names=F, sep="\t", quote=F)

# 여까지 수정함. 

# # choose for field 1 or filed 2
# # bring all data base 2field HLA allleds
# hla2filed_list <- read.table("/kimlab_wd/yuo1996/tools/SNP2HLA_package_v1.0.3/MakeReference/HLA2field.list", header=F)
# hla2filed_list <- separate(hla2filed_list, V1, into= c("type", "f1", "f2"), sep=":") 
# hla2filed_list$type <- paste0("HLA_", hla2filed_list$type)
# hla2filed_list <- unite( hla2filed_list , allelle, c(f1, f2), sep="")  
# head(hla2filed_list)

total_HLA_imputed_1f <- read.table(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_imputed_1f") , header=T)
total_HLA_imputed_2f <- read.table(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_imputed_2f") , header=T)
total_HLA_imputed_Gf <- read.table(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_imputed_Gf") , header=T)


# retrive phased result (For groundtruth)
total_HLA_phased_1f <- read.table(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_phased_1f") , header=T)
total_HLA_phased_2f <- read.table(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_phased_2f") , header=T)
total_HLA_phased_Gf <- read.table(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_", opt$haplonetver,"/total_HLA_phased_Gf") , header=T)


total_HLA_imputed_Gf$type <- gsub("gg_", "", total_HLA_imputed_Gf$type)
total_HLA_phased_Gf$type <- gsub("gg_", "", total_HLA_phased_Gf$type)

# compare between phased vs imputed

head(total_HLA_imputed_Gf)
head(total_HLA_phased_1f)
head(total_HLA_phased_1f)

phased = total_HLA_phased_2f
imputed = total_HLA_imputed_2f

hla_performance <- function( phased, imputed){
    phased <- filter(phased, ID %in% unique( imputed$ID ))    
    res <- data.frame()
    for( i in unique( imputed$ID ) ){
        
        #i=unique( imputed$ID )[1]
        print(i)
        tmp = data.frame(ID=i)

        for(h in unique( imputed$type )){

            #h=unique( imputed$type )[1]
            
            l=length( intersect( filter(imputed, ID==i, type==h)$GT , filter(phased, ID==i, type==h)$GT ) )
            t=data.frame(h=l) 
            colnames(t)=h

            if(nrow( imputed[which(imputed$ID==i & imputed$type==h), ] ) !=0 ){ # 0 copy imputed된경우제외.. 
#            if( unique(imputed[which(imputed$ID==i & imputed$type==h), "imputed_allele_num"])==2 & l==1 ){ # 잘못된코드
            if( unique( imputed[which(imputed$ID==i & imputed$type==h), "imputed_allele_num"])==2 && length(unique( imputed[which(imputed$ID==i & imputed$type==h), "GT"]))==1 && l==1 ){
                t[1,h] <- 2 # homo 맞는경우
            }else if( unique(imputed[which(imputed$ID==i & imputed$type==h), "imputed_allele_num"])==3 & l==2 ){
                t[1,h] <- 1 
            }
}
            tmp = cbind(tmp, t)
            }
        res=rbind(res, tmp)
        }

    return(res)
}
# take minutes

res_1f= hla_performance(total_HLA_phased_1f, total_HLA_imputed_1f)
head(res_1f)
res_1f[-1] %>% as.matrix() %>% as.vector() %>% table()
res_2f= hla_performance(total_HLA_phased_2f, total_HLA_imputed_2f)
head(res_2f)
filter(res_2f, ID=="KOREA1K-102")
res_2f[-1] %>% as.matrix() %>% as.vector() %>% table()
res_Gf= hla_performance(total_HLA_phased_Gf, total_HLA_imputed_Gf)
res_Gf[-1] %>% as.matrix() %>% as.vector() %>% table()
head(res_Gf)

# hla_val = res~

concor_cal <- function( con_table  , type_num , sample_num){
    con_table <- as.vector(con_table)
    if(length(con_table)==3 ){
    ac <- (con_table[2]+con_table[3]*2)/(sample_num*type_num*2)*100
    }else if(length( con_table)==2){
    ac <- (con_table[1]+con_table[2]*2)/(sample_num*type_num*2)*100
    }else if(length( con_table)==1){
        ac=100
    }
    return(ac)
}

# 고민/////
#hla_val <- hla_val[-which(hla_val$ID %in% gsub("-", "", non2copyimputed$SAMPLE) ), ]

mk_concordance <- function( field ){
    if(field==1){
        res=res_1f
    }else if(field==2){
        res=res_2f
    }else if(field=="G"){
        res=res_Gf
    }
        con_A <- select( res, HLA_A)  %>% as.matrix() %>% as.vector() %>% table() 
        con_B <- select( res, HLA_B)  %>% as.matrix() %>% as.vector() %>% table() 
        con_C <- select( res, HLA_C)  %>% as.matrix() %>% as.vector() %>% table() 
        con_DPA1 <- select( res, HLA_DPA1)  %>% as.matrix() %>% as.vector() %>% table() 
        con_DPB1 <- select( res, HLA_DPB1)  %>% as.matrix() %>% as.vector() %>% table() 
        con_DQA1 <- select( res, HLA_DQA1)  %>% as.matrix() %>% as.vector() %>% table() 
        con_DQB1 <- select( res, HLA_DQB1)  %>% as.matrix() %>% as.vector() %>% table() 
        con_DRB1 <- select( res, HLA_DRB1)  %>% as.matrix() %>% as.vector() %>% table() 
        whole_table <- res[,-1] %>% as.matrix() %>% as.vector() %>% table()  

        tohla <- c( concor_cal(con_A, 1, nrow(res)),
                    concor_cal(con_B, 1, nrow(res)),
                    concor_cal(con_C, 1, nrow(res)),
                    concor_cal(con_DPA1, 1, nrow(res)),
                    concor_cal(con_DPB1, 1, nrow(res)),
                    concor_cal(con_DQA1, 1, nrow(res)),
                    concor_cal(con_DQB1, 1, nrow(res)),
                    concor_cal(con_DRB1, 1, nrow(res)),
                    concor_cal(whole_table, 8, nrow(res) )
                )
    tohla <- t( as.data.frame(tohla) )
    colnames(tohla) <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "total")

    return(tohla)
}

con1=mk_concordance(1)
con2=mk_concordance(2)
cong=mk_concordance("G")

con1
con2
cong

# only for crossval1
#round(tohla, digits=1)
path=paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/",opt$region, "_", opt$haplonetver, "/", opt$crossval) 
#setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/",opt$haplonetver, "/", opt$crossval) )
write.table(con1, paste0(path,"/HLA_diploid_contable1.txt"), col.names=T,row.names=F, sep="\t", quote=F)
write.table(con2, paste0(path,"/HLA_diploid_contable2.txt"), col.names=T,row.names=F, sep="\t", quote=F)
write.table(cong, paste0(path,"/HLA_diploid_contableg.txt"), col.names=T,row.names=F, sep="\t", quote=F)
#
rm( list=grep("output.test", ls(), value=T))
setwd(path)
save.image( "./val_HLA_result.RDdata") 
