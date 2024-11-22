rm(list=ls())
library(dplyr);library(tidyr)
library(optparse)
##################################  arguement setting ########################################
args_list <- list(
    make_option( "--snpver" , type="character" , default= "ver3" , help="snpversion", metavar="character"),
    make_option( "--region" , type="character" , default= "ver3" , help="snpversion", metavar="character"),
    make_option( "--crossval" , type="character" , default= "crossval50" , help="snpversion", metavar="character"),
    make_option( "--haplonetver" , type="character" , default= "crossval50" , help="snpversion", metavar="character")
)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)


# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8QUAL20"; opt$haplonetver="ver1_sujung"; opt$crossval="crossval1"

###############################################################################################

# tmp result file from imputation
setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval) )


##################  cross-val version ( several sets )

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


##################### 1. imputed C4 allele: 2개 이상 imputaion된거 고려해줘야??!!!
#output.test <- output.test942

make_C4_CN_allele <- function( output.test ){

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

    C4_CN_allele <- data.frame( SAMPLE= NA, C4A_1= NA, C4A_2= NA, C4B_1 = NA,C4B_2= NA, HERV_1 =NA, HERV_2 =NA)

    for(i in 1:sam) {
        i=2*i+1
        # check for 2_1 case
        # tmp <- res_C4[-1, c(2,which(res_C4[1,] =="KOREA1K-178"))]

        tmp <- res_C4[ -1, c(2, i:(i+1) )]
        colnames(tmp) <- tmp[1,]
        tmp <- tmp[-1,]
        
        allele <- c()
  

        for (k in c("C4A_", "C4B_", "HERV")){

            # k = "C4_"
            indv <- tmp[grep(k , tmp$id),]
            indv$cn <- as.numeric( substr( indv$id , nchar(indv$id) , nchar(indv$id) ) )

            chr1_P <-  indv[ grep("P", indv[,2]), 1 ]  
            chr2_P <-  indv[ grep("P", indv[,3]), 1 ] 

            # 한 copy이상 나오면 합쳐줌. 
            if( length(chr1_P) > 1 ){
                chr1_P <- paste0( k, "copy", sum( indv[ grep("P", indv[,2]), "cn" ] ) )
            }

            if( length(chr2_P) > 1 ){
                chr2_P <- paste0( k, "copy", sum( indv[ grep("P", indv[,3]), "cn" ] ) )
            }

            total <- c(chr1_P, chr2_P) 
            
            if(length(total) ==1 ){
               total <- c( total, "NULL") 
            }
            allele<- c( allele, total)

        }

        t <- c( colnames(tmp)[2] , allele)  %>% data.frame() %>% t()
        colnames(t) <- c("SAMPLE","C4A_1","C4A_2", "C4B_1","C4B_2", "HERV_1", "HERV_2")
        C4_CN_allele <- rbind(C4_CN_allele, t)
    }


    C4_CN_allele <- C4_CN_allele[-1,]
    rownames(C4_CN_allele) <- 1:nrow(C4_CN_allele)
    return(C4_CN_allele)


}

#total_C4_CN_allele <- make_C4_CN_allele( output.test1) %>% mutate( set = paste0("set", 1))
                    
num <- gsub("KOREA1K-","korea", imputed_list) %>% gsub("[^0-9]", "", . )

total_C4_CN_allele <- data.frame()
for(o in paste0("output.test", sort(num)) ){
    tmp <- make_C4_CN_allele( get(o) ) #%>% mutate( set = paste0("set", a)) 
    total_C4_CN_allele <- rbind(total_C4_CN_allele, tmp)
}


nrow(total_C4_CN_allele)
head(total_C4_CN_allele)
total_C4_CN_allele[-1] %>% as.matrix() %>% table() # NULL??


################## 2. retrive phased alleles : 2개 이상 imputaion된거 고려해줘야??!!!
# cross val set ver
library(stringr)

setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval) )

phased_list <- grep(".phased.vcf", list.files("./") , value=TRUE) 
length(phased_list)

i=phased_list[1]
for(i in phased_list){ 
    test_bgl <-read.delim( paste0( "./", i ) , header=T, sep="\t", skip=12)
#    test_bgl <-read.delim( paste0( "./", i ) , header=T, sep="\t", skip=8)
    num <- gsub("KOREA1K-","korea",i) %>% gsub("[^0-9]", "", . )
    assign( paste0("test_bgl", num), test_bgl, envir = .GlobalEnv )
    rm(test_bgl)
}
ls()


# test_bgl <- test_bgl901
make_C4_CN_masked <- function( test_bgl ){

        # test_bgl_C4 <- test_bgl[c(1,2, grep("copy", test_bgl$V2) ),] # beagle?
        res_C4 <- test_bgl[ grep("copy", test_bgl$ID) ,] # impute5
        id <- colnames( res_C4 )[dim(res_C4)[2]] %>% gsub( ".", "-", .,fixed=TRUE)
        res_C4 <- rename(res_C4, sample=colnames( res_C4 )[dim(res_C4)[2]])
        res_C4 <- select( res_C4, ID, sample ) %>% rename( GT=sample)
        res_C4$GT <- gsub("|", "abc", res_C4$GT ,fixed=TRUE)
        res_C4 <- separate( res_C4, GT, into=c("V3","V4") , sep="abc" )
        res_C4$V3 <- ifelse( res_C4$V3==1, "P", "A"); res_C4$V4 <- ifelse( res_C4$V4==1, "P", "A")
        res_C4 <- res_C4 %>% mutate( V1="M") %>% relocate( V1, .before=ID) %>% rename( V2=ID)
        res_C4 <- rbind( data.frame( V1=c("P","I"), V2=c("pedigree", "id"), V3=id, V4=id) , res_C4)
        test_bgl_C4 <- res_C4

        sam <-  ( dim( test_bgl_C4 )[2]-2  )/2

        C4_CN_masked <- data.frame( SAMPLE= NA, C4A_1= NA, C4A_2= NA, C4B_1 = NA,C4B_2= NA, HERV_1 =NA, HERV_2 =NA)

        for(i in 1:sam) {
            i=2*i+1
            tmp <- test_bgl_C4[ -1, c(2, i:(i+1) )]
            colnames(tmp) <- tmp[1,]
            tmp <- tmp[-1,]
            
            allele <- c()

            for (k in c( "C4A_", "C4B_", "HERV")){
                indv <- tmp[grep(k , tmp$id),]


                chr1_P <-  indv[ grep("P", indv[,2]), 1 ]  
                chr2_P <-  indv[ grep("P", indv[,3]), 1 ] 

                total <- c(chr1_P, chr2_P) 
                
                # phasing 후 2copy안나오는 경우땜에 추가. 아님 앞3_1에서제거하던가
                if(length(total)==1){
                    total <- c(total, "NULL")
                }
                #
                allele<- c( allele, total)

            }

            t <- c( colnames(tmp)[2] , allele)  %>% data.frame() %>% t()
            colnames(t) <- c("SAMPLE", "C4A_1","C4A_2", "C4B_1","C4B_2", "HERV_1", "HERV_2")
            C4_CN_masked <- rbind(C4_CN_masked, t)
        }


        C4_CN_masked  <- C4_CN_masked[-1,]
        rownames(C4_CN_masked ) <- 1:nrow(C4_CN_masked )
        return(C4_CN_masked)

}


#total_C4_CN_masked <- make_C4_CN_masked( test_bgl1 )  %>% mutate( set = paste0("set", 1) )

num <- gsub("KOREA1K-","korea", imputed_list) %>% gsub("[^0-9]", "", . )

total_C4_CN_masked <- data.frame()
for(o in paste0("test_bgl", sort(num)) ){
    tmp <- make_C4_CN_masked( get(o) ) #%>% mutate( set = paste0("set", a)) 
    total_C4_CN_masked <- rbind(total_C4_CN_masked, tmp)
}

#setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try5/7.Accuracy/", opt$snpver, "/", opt$crossval) )
#write.table( total_C4_CN_masked, "./total_C4_CN_masked.txt" ), col.names = TRUE, row.names=FALSE, sep="\t", quote=F)
head(total_C4_CN_masked)



######### compare

# # maf check
# head(total_C4_CN_masked)
# filter(total_C4_CN_masked, set=="set1") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set2") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set3") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set4") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set5") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set6") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set7") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set8") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set9") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
# filter(total_C4_CN_masked, set=="set10") %>% .[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()

total_C4_CN_masked[ grep("C4_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
total_C4_CN_masked[ grep("C4A_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
total_C4_CN_masked[ grep("C4B_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()
total_C4_CN_masked[ grep("HERV_", colnames(total_C4_CN_masked))] %>% as.matrix() %>% as.vector() %>% table()


#
head(total_C4_CN_allele)
head(total_C4_CN_masked)
#total_C4_CN_masked <- total_C4_CN_masked[which(total_C4_CN_masked$SAMPLE %in% total_C4_CN_allele$SAMPLE),]
nrow(total_C4_CN_allele)==nrow(total_C4_CN_masked)

is.na(total_C4_CN_masked) %>% table()
is.na(total_C4_CN_allele) %>% table()

if( names( table( total_C4_CN_allele$SAMPLE==total_C4_CN_masked$SAMPLE ) ) ){
    SAMPLE <- total_C4_CN_allele$SAMPLE
    print("can go!!")
}else{
    print("stop!!!!!!!!!!!!!!!!!!!!!")
}


C4_allele_val <- NULL

#id=SAMPLE[1]
for( id in SAMPLE) {

    #tmp <- data.frame( ID= id, set = total_C4_CN_allele[which( total_C4_CN_allele$SAMPLE==id) , "set" ])
    tmp <- data.frame( ID= id )

    for( c4type  in c( "C4A", "C4B", "HERV")  ){

        #c4type="HERV"
        
        imputed_GT <- filter(total_C4_CN_allele, SAMPLE==id ) %>% select( starts_with( paste0( c4type, "_"))) %>% as.matrix() %>% as.vector()
        masked_GT <- filter(total_C4_CN_masked, SAMPLE==id ) %>% select( starts_with( paste0( c4type, "_"))) %>% as.matrix() %>% as.vector()
        # imputed_GT;masked_GT


        if( length( unique(imputed_GT)) ==1  &  length( unique(masked_GT)) ==1 ){
            
            if( length( intersect(imputed_GT, masked_GT) )  == 1 ){
                tmp$t <-2  # homo이고 맞는 경우
            }else{
                tmp$t <-0  #homo이고 틀린경우
            }

        }else{
        tmp$t <- length( intersect(imputed_GT, masked_GT) ) ## homo 잘못 나왔을 수 있음 !!
        }
        colnames(tmp)[which(colnames(tmp)=="t")] <- c4type


    }

    C4_allele_val <- rbind(C4_allele_val, tmp)
}

# take 1 minutes under

head(C4_allele_val)
nrow(C4_allele_val)

setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/5.Imputation_C4HLA_using_Beagle/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval) )
rmsamples <- read.table("../rmsamples")
# C4_allele_val <- C4_allele_val[ -which(C4_allele_val$ID %in% rmsamples[,1]),]


table( C4_allele_val$C4 ) #%>% sum
table( C4_allele_val$C4A ) #%>% sum
table( C4_allele_val$C4B ) #%>% sum
table( C4_allele_val$HERV ) #%>% sum

concor_cal <- function( con_table  , allele_num , sample_num){
    con_table <- as.vector(con_table)
    if(length(con_table)==3 ){
    ac <- (con_table[2]+con_table[3]*2)/(sample_num*allele_num*2)*100
    }else if(length( con_table)==2){
    ac <- (con_table[1]+con_table[2]*2)/(sample_num*allele_num*2)*100
    }else if (length( con_table)==1){
    ac <- (con_table[1]*2)/(sample_num*allele_num*2)*100
    }
    return(ac)
}

#whole_table <- C4_allele_val[,-c(1:2)] %>% as.matrix() %>% as.vector() %>% table()  
# crossval1
whole_table <- C4_allele_val[,-1] %>% as.matrix() %>% as.vector() %>% table()  

concor_cal( whole_table , 4 , nrow(C4_allele_val))

to <- c( 
        concor_cal(table( C4_allele_val$C4A ) , 1 ,  nrow(C4_allele_val)) ,
        concor_cal(table( C4_allele_val$C4B ) , 1 ,  nrow(C4_allele_val)) ,
        concor_cal(table( C4_allele_val$HERV ) , 1 ,  nrow(C4_allele_val)),
        concor_cal( whole_table , 3 , nrow(C4_allele_val))
        )
to

# only for crossval1
to <- as.data.frame(to) %>% t()
colnames(to) <- c( "C4A", "C4B", "HERV", "total")
to
setwd(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval) )
write.table(to,  "./C4_allele_contable.txt" , col.names=T,row.names=F, sep="\t", quote=F)

rm( list=grep("test", ls(), value=T))
save.image("./val_C4_allele_result.RDdata")
load("./val_C4_allele_result.RDdata")




