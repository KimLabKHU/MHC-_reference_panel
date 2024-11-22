###################!!!!!!!! copy 수마다 acc 다른지 보기 ; cn2가 가장 freq이 높으니 가장 acc도 좋을것.
rm(list=ls())
library(dplyr);library(tidyr)
library(optparse)
library(stringr)

##################################  arguement setting ###########################################
args_list <- list(
    make_option( "--snpver" , type="character" , default= "ver3" , help="snpversion", metavar="character"),
    make_option( "--region" , type="character" , default= "start31.7end32.2" , help="haplonet_region", metavar="character"),
    make_option( "--crossval" , type="character" , default= "crossval50" , help="snpversion", metavar="character"),
    make_option( "--haplonetver" , type="character" , default= "crossval50" , help="snpversion", metavar="character")
)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)

# opt=list(); opt$snpver="ver3"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8"; opt$haplonetver="ver1_best"; opt$crossval="crossval1"
#################################################################################################
load(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval, "/val_C4_allele_result.RDdata") )
load("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/C4_genomstrip/ver3/C4_HERV_AB_GT_dist.RData")
ls()

head(total_C4_CN_allele)
head(total_C4_CN_masked)
nrow(total_C4_CN_allele)==nrow(total_C4_CN_masked)

# ref로 찢어서 고정한얘들
h1 <- "C4_copy2"; h2 <- "C4_copy3" # 1+4가 있는지 보려고
h1 <- "C4_copy2"; h2 <- "C4_copy2" # 1+3가 있는지 보려고
h1 <- "C4_copy2"; h2 <- "C4_copy4" # 3+3가 있는지 보려고 -> ver3 에서 좀 있어보임?
h1 <- "C4_copy1"; h2 <- "C4_copy1" # 2+0가 있는지 보려고
h1 <- "C4_copy1"; h2 <- "C4_copy2" # 만 가능
h1 <- "C4_copy3"; h2 <- "C4_copy4" # 만 가능
h1 <- "C4_copy4"; h2 <- "C4_copy4" # 만 가능
ref <- filter( total_C4_CN_masked , C4_1==h1 & C4_2==h2 | C4_2==h1 & C4_1==h2)  %>% select( SAMPLE , starts_with("C4_"), set) %>% .$SAMPLE
filter( total_C4_CN_allele , SAMPLE %in% ref) %>% select( SAMPLE , starts_with("C4_"), set) %>% filter( !(C4_1==h1 & C4_2==h2) & !(C4_2==h1 & C4_1==h2)) 

h1 <- "HERV_copy1"; h2 <- "HERV_copy1" # 2+0가 있는지 보려고
h1 <- "HERV_copy2"; h2 <- "HERV_copy1" # 0+3가 있는지 보려고
h1 <- "HERV_copy2"; h2 <- "HERV_copy2" # 1+3/0+4가 있는지 보려고 
h1 <- "HERV_copy2"; h2 <- "HERV_copy3" # 1+4가 있는지 보려고
h1 <- "HERV_copy2"; h2 <- "HERV_copy4" # 3+3가 있는지 보려고
ref <- filter( total_C4_CN_masked , HERV_1==h1 & HERV_2==h2 | HERV_2==h1 & HERV_1==h2)  %>% select( SAMPLE , starts_with("HERV_"), set) %>% .$SAMPLE
filter( total_C4_CN_allele , SAMPLE %in% ref) %>% select( SAMPLE , starts_with("HERV_"), set) %>% filter( !(HERV_1==h1 & HERV_2==h2) & !(HERV_2==h1 & HERV_1==h2)) 

h1 <- "C4A_copy1"; h2 <- "C4A_copy1" # 2+0가 있는지 보려고
h1 <- "C4A_copy1"; h2 <- "C4A_copy2" # 3+0가 있는지 보려고
h1 <- "C4A_copy1"; h2 <- "C4A_copy3" # 2+2가 있는지 보려고
ref <- filter( total_C4_CN_masked , C4A_1==h1 & C4A_2==h2 | C4A_2==h1 & C4A_1==h2)  %>% select( SAMPLE , starts_with("C4A_"), set) %>% .$SAMPLE
filter( total_C4_CN_allele , SAMPLE %in% ref) %>% select( SAMPLE , starts_with("C4A_"), set) %>% filter( !(C4A_1==h1 & C4A_2==h2) & !(C4A_2==h1 & C4A_1==h2))

h1 <- "C4B_copy1"; h2 <- "C4B_copy1" # 2+0가 있는지 보려고
h1 <- "C4B_copy1"; h2 <- "C4B_copy2" # 3+0가 있는지 보려고
h1 <- "C4B_copy1"; h2 <- "C4B_copy3" # 2+2가 있는지 보려고
ref <- filter( total_C4_CN_masked , C4B_1==h1 & C4B_2==h2 | C4B_2==h1 & C4B_1==h2)  %>% select( SAMPLE , starts_with("C4B_"), set) %>% .$SAMPLE
filter( total_C4_CN_allele , SAMPLE %in% ref) %>% select( SAMPLE , starts_with("C4B_"), set) %>% filter( !(C4B_1==h1 & C4B_2==h2) & !(C4B_2==h1 & C4B_1==h2))


# cp별 acc freq보는건 여기서부터

if( names(table( total_C4_CN_allele$SAMPLE==total_C4_CN_masked$SAMPLE ) )){
    SAMPLE <- total_C4_CN_allele$SAMPLE
    print("can go!!")
}else{
    print("stop!!!!!!!!!!!!!!!!!!!!!")
}


######### 1. C4 con/cor table according to diploid cn
###################!!!!!!!! copy 수마다 acc 다른지 보기 ; cn2가 가장 freq이 높으니 가장 acc도 좋을것.
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

# opt=list(); opt$snpver="ver3"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8"; opt$haplonetver="ver1_best"; opt$crossval="crossval1"
# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8QUAL20"; opt$haplonetver="ver1_best"; opt$crossval="crossval1"

###############################################################################################
load(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval, "/val_C4_result.RDdata") )
load("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/C4_genomestrip/C4_HERV_AB_GT_dist.RData")

head(vali)
vali$batch %>% unique %>% table
nrow(vali)

vali$batch
vali$SAMPLE

# diploid cn types
C4_alleles <- vali %>% select( C4 ) %>% as.matrix() %>% as.vector() %>% table() %>% names() %>% as.numeric()
C4A_alleles <- vali %>% select( C4A ) %>% as.matrix() %>% as.vector() %>% table() %>% names() %>% as.numeric()
C4B_alleles <- vali %>% select( C4B ) %>% as.matrix() %>% as.vector() %>% table() %>% names() %>% as.numeric()
HERV_alleles <- vali %>% select( HERV ) %>% as.matrix() %>% as.vector() %>% table() %>% names() %>% as.numeric()

C4_alleles;C4A_alleles;C4B_alleles;HERV_alleles

c4type="C4A"
make_dn_contigency_table <- function( c4type) {

    c4_type <- vali[ c( which( colnames(vali) ==c4type) , which( colnames(vali) == paste0( "imputed_", c4type)  ) )] 
    colnames(c4_type) = c("typed", "imputed")

    contable =table(c4_type)
contable
    sensitivity=diag(contable)/colSums(contable)
    sensitivity=diag(contable)/rowSums(contable)
    
    
    contable=rbind(contable, colSums(contable) )
    contable=cbind(contable,  rowSums(contable) )
    rownames(contable)[ length(rownames(contable) )]="sum"
    colnames(contable)[ length(colnames(contable) )]="sum"

contable[1:5, 1:5]


    
    dcn <- get(paste0(c4type, "_alleles") ) 

    contingency_table_set <- list()
    
    for(i in dcn){

        TP <- filter( c4_type, typed==i, imputed==i) %>% nrow 
        FN <- filter( c4_type, typed==i, imputed!=i) %>% nrow
        FP <- filter( c4_type, typed!=i , imputed==i) %>% nrow
        TN <- filter( c4_type, typed!=i , imputed!=i) %>% nrow

        t <- matrix( c(TP, FN, FP, TN), nrow=2)
        colnames(t) <- c("typed_P", "typed_N")
        rownames(t) <- c("imputed_P", "imputed_N")
        
        t <- list(t)
        names(t) <- i

        contingency_table_set <- append( contingency_table_set , t)
    }

    return( contingency_table_set )
}



# previous
# c4type="C4A"
# make_dn_contigency_table <- function( c4type) {

#     c4_type <- vali[ c( which( colnames(vali) ==c4type) , which( colnames(vali) == paste0( "imputed_", c4type)  ) )] 
#     colnames(c4_type) = c("typed", "imputed")
#     dcn <- get(paste0(c4type, "_alleles") ) 

#     contingency_table_set <- list()
#     for(i in dcn){

#         TP <- filter( c4_type, typed==i, imputed==i) %>% nrow 
#         FN <- filter( c4_type, typed==i, imputed!=i) %>% nrow
#         FP <- filter( c4_type, typed!=i , imputed==i) %>% nrow
#         TN <- filter( c4_type, typed!=i , imputed!=i) %>% nrow

#         t <- matrix( c(TP, FN, FP, TN), nrow=2)
#         colnames(t) <- c("typed_P", "typed_N")
#         rownames(t) <- c("imputed_P", "imputed_N")
        
#         t <- list(t)
#         names(t) <- i

#         contingency_table_set <- append( contingency_table_set , t)
#     }

#     return( contingency_table_set )
# }


make_dn_contigency_table( "C4")
make_dn_contigency_table( "C4A")
make_dn_contigency_table( "C4B")
make_dn_contigency_table( "HERV")

head(vali)
cor(vali$C4, vali$imputed_C4)#^2
cor(vali$C4A, vali$imputed_C4A)#^2
cor(vali$C4B, vali$imputed_C4B)#^2
cor(vali$HERV, vali$imputed_HERV)#^2


make_dn_cor_table <- function( c4type) {

    c4_type <- vali[ c( which( colnames(vali) ==c4type) , which( colnames(vali) == paste0( "imputed_", c4type)  ) )] 
    colnames(c4_type) = c("typed", "imputed")
    dcn <- get(paste0(c4type, "_alleles") ) 

    cor_table_set <- list()
    for(i in dcn){

        cor( filter( c4_type, typed==i)[,"typed"] , filter( c4_type, typed==i)[,"imputed"] )
# C4B에서 1을 2로 예측한 경우가 많다/




####### 2. C4 con/cor table according to haploid cn
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

# opt=list(); opt$snpver="ver3"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8"; opt$haplonetver="ver1_best"; opt$crossval="crossval1"

# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8QUAL20"; opt$haplonetver="ver1_best"; opt$crossval="crossval1"

###############################################################################################

load(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval, "/val_C4_result.RDdata") )
load(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval, "/val_C4_allele_result.RDdata") )
load("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/C4_genomestrip/C4_HERV_AB_GT_dist.RData")


if( names( table( total_C4_CN_allele$SAMPLE==total_C4_CN_masked$SAMPLE )  )){
    SAMPLE <- total_C4_CN_allele$SAMPLE
    print("can go!!")
}else{
    print("stop!!!!!!!!!!!!!!!!!!!!!")
}


C4_alleles <- total_C4_CN_masked %>% select( starts_with("C4_")) %>% as.matrix() %>% as.vector() %>% table() %>% names()
C4A_alleles <- total_C4_CN_masked %>% select( starts_with("C4A_")) %>% as.matrix() %>% as.vector() %>% table() %>% names()
C4B_alleles <- total_C4_CN_masked %>% select( starts_with("C4B_")) %>% as.matrix() %>% as.vector() %>% table() %>% names()
HERV_alleles <- total_C4_CN_masked %>% select( starts_with("HERV_")) %>% as.matrix() %>% as.vector() %>% table() %>% names()

all <- c(C4_alleles , C4A_alleles , C4B_alleles , HERV_alleles)
acc_for_alleles <- NULL

#ex 
data <- total_C4_CN_allele; c4type="C4"; alleles_set=C4_alleles
#
mkallele_count_table <- function( c4type, alleles_set, data){  # c4type = C4, C4A, C4B, HERV

tmp <- select(data, SAMPLE, starts_with( paste0( c4type, "_"))) 

for( a in alleles_set){
    tmp$tt <- 0
    colnames(tmp)[which( colnames(tmp) == "tt" )] <-a
    
}

for(i in 1:nrow(tmp)){

    ta <- as.vector(as.matrix( tmp[i, 2:3] )) %>% table()
    tmp[i, which(colnames(tmp) %in% names( ta ) ) ] <- ta[ which(names(ta) %in% colnames(tmp))]

}

return(tmp)

}


allele_count_masked_C4 <- mkallele_count_table( "C4" , C4_alleles, total_C4_CN_masked ) 
allele_count_masked_C4A <- mkallele_count_table( "C4A" , C4A_alleles, total_C4_CN_masked ) 
allele_count_masked_C4B <- mkallele_count_table( "C4B" , C4B_alleles, total_C4_CN_masked ) 
allele_count_masked_HERV <- mkallele_count_table( "HERV" , HERV_alleles, total_C4_CN_masked ) 


allele_count_imputed_C4 <- mkallele_count_table( "C4" , C4_alleles, total_C4_CN_allele ) 
allele_count_imputed_C4A <- mkallele_count_table( "C4A" , C4A_alleles, total_C4_CN_allele ) 
allele_count_imputed_C4B <- mkallele_count_table( "C4B" , C4B_alleles, total_C4_CN_allele ) 
allele_count_imputed_HERV <- mkallele_count_table( "HERV" , HERV_alleles, total_C4_CN_allele ) 

c4type="C4"; allele="C4_copy1"
calculate_per <- function( c4type, allele) {

masked <- get(paste0( "allele_count_masked_" , c4type))
imputed <- get(paste0( "allele_count_imputed_" , c4type))

head(masked)
head(imputed)

masked_not0=filter(masked, get(allele) != 0)
imputed_not0=imputed[which(imputed$SAMPLE %in% masked_not0$SAMPLE) , ]


if( any( masked_not0$SAMPLE==imputed_not0$SAMPLE ) ){
sum( masked_not0$C4_copy1)

masked <- masked_not0[which( colnames(masked_not0)==allele) ] %>% as.matrix() %>% as.vector()
imputed <- imputed_not0[which( colnames(imputed_not0)==allele) ] %>% as.matrix() %>% as.vector()

# 점수 계산
score <- sum(
  (masked == 2 & imputed == 2) * 2,  # masked가 2이고 imputed도 2이면 2점
  (masked == 2 & imputed == 1),      # masked가 2이고 imputed가 1이면 1점
  (masked == 1 & imputed == 1)       # masked가 1이고 imputed도 1이면 1점
)

}
return( score/sum(masked) )
}


calculate_per( "C4" , "C4_copy1")  %>% round(., digits = 4)*100
calculate_per( "C4" , "C4_copy2") %>% round(., digits = 4)*100
calculate_per( "C4" , "C4_copy3") %>% round(., digits = 4)*100# cor이 안좋다?
calculate_per( "C4" , "C4_copy4") %>% round(., digits = 4)*100# 베리굿!!!

calculate_per( "C4A" , "C4A_copy0") %>% round(., digits = 4) *100
calculate_per( "C4A" , "C4A_copy1") %>% round(., digits = 4) *100
calculate_per( "C4A" , "C4A_copy2") %>% round(., digits = 4) *100
calculate_per( "C4A" , "C4A_copy3") %>% round(., digits = 4) *100

calculate_per( "C4B" , "C4B_copy0") %>% round(., digits = 4) *100
calculate_per( "C4B" , "C4B_copy1") %>% round(., digits = 4)*100
calculate_per( "C4B" , "C4B_copy2") %>% round(., digits = 4)*100
calculate_per( "C4B" , "C4B_copy3") %>% round(., digits = 4)*100

calculate_per( "HERV" , "HERV_copy0") %>% round(., digits = 4)*100
calculate_per( "HERV" , "HERV_copy1") %>% round(., digits = 4)*100
calculate_per( "HERV" , "HERV_copy2") %>% round(., digits = 4)*100
calculate_per( "HERV" , "HERV_copy3") %>% round(., digits = 4)*100
calculate_per( "HERV" , "HERV_copy4") %>% round(., digits = 4)*100


c4type="C4"; allele="C4_copy1"

calculate_freq_and_acc <- function( c4type, allele) {

masked <- get(paste0( "allele_count_masked_" , c4type))
imputed <- get(paste0( "allele_count_imputed_" , c4type))

if( unique(masked$SAMPLE== imputed$SAMPLE ) ){

    masked_at <- masked[, which(colnames(masked) == allele )]  
    imputed_at <- imputed[, which(colnames(imputed) == allele )] 
# % ver    
#    masked_freq <-  round( table(masked_at)/nrow(masked)*100 , 2)
#    imputed_freq <- round( table(imputed_at)/nrow(masked)*100 , 2)
# indiv ver
    masked_freq <-  table(masked_at)
    imputed_freq <- table(imputed_at)

    con <- round( table( masked_at == imputed_at )["TRUE"]/nrow(masked)*100 , 2)
    cor <- round( cor( masked_at , imputed_at ) , 2)
}

    return( list(masked_freq , imputed_freq , con, cor) )

}

calculate_freq_and_acc( "C4" , "C4_copy1")  %>% .[[3]] # 베리굿!!!!! 
calculate_freq_and_acc( "C4" , "C4_copy2") %>% .[[3]]
calculate_freq_and_acc( "C4" , "C4_copy3") %>% .[[3]]# cor이 안좋다?
calculate_freq_and_acc( "C4" , "C4_copy4") %>% .[[3]]# 베리굿!!!

calculate_freq_and_acc( "C4A" , "C4A_copy0") %>% .[[3]]
calculate_freq_and_acc( "C4A" , "C4A_copy1") %>% .[[3]]
calculate_freq_and_acc( "C4A" , "C4A_copy2") %>% .[[3]]
calculate_freq_and_acc( "C4A" , "C4A_copy3") %>% .[[3]]

calculate_freq_and_acc( "C4B" , "C4B_copy0") %>% .[[3]]
calculate_freq_and_acc( "C4B" , "C4B_copy1") %>% .[[3]]
calculate_freq_and_acc( "C4B" , "C4B_copy2") %>% .[[3]]
calculate_freq_and_acc( "C4B" , "C4B_copy3") %>% .[[3]]

calculate_freq_and_acc( "HERV" , "HERV_copy0") %>% .[[3]]
calculate_freq_and_acc( "HERV" , "HERV_copy1") %>% .[[3]]
calculate_freq_and_acc( "HERV" , "HERV_copy2") %>% .[[3]]
calculate_freq_and_acc( "HERV" , "HERV_copy3") %>% .[[3]]
calculate_freq_and_acc( "HERV" , "HERV_copy4") %>% .[[3]]


nrow(allele_count_masked_C4)
freq = function (table){
a=table( table)
if(length(a)==3){
(a[2]+a[3]*2)/(1522*2) %>% return() 
}else if( length(a)==2){
(a[2])/(1522*2) %>% return() 
}
}

freq(allele_count_masked_C4$C4_copy1) *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4$C4_copy2)  *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4$C4_copy3)  *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4$C4_copy4) *100  %>%  round(. , digits=2)

freq(allele_count_masked_C4A$C4A_copy0) *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4A$C4A_copy1)  *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4A$C4A_copy2)  *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4A$C4A_copy3) *100  %>%  round(. , digits=2)

freq(allele_count_masked_C4B$C4B_copy0) *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4B$C4B_copy1)  *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4B$C4B_copy2)  *100  %>%  round(. , digits=2)
freq(allele_count_masked_C4B$C4B_copy3) *100  %>%  round(. , digits=2)


freq(allele_count_masked_HERV$HERV_copy0) *100  %>%  round(. , digits=2)
freq(allele_count_masked_HERV$HERV_copy1)  *100  %>%  round(. , digits=2)
freq(allele_count_masked_HERV$HERV_copy2)  *100  %>%  round(. , digits=2)
freq(allele_count_masked_HERV$HERV_copy3) *100  %>%  round(. , digits=2)
freq(allele_count_masked_HERV$HERV_copy4) *100  %>%  round(. , digits=2)


#  안고침
for( id in SAMPLE) {

    tmp <- data.frame( ID= id, set = total_C4_CN_allele[which( total_C4_CN_allele$SAMPLE==id) , "set" ])

    for( c4type  in c("C4", "C4A", "C4B", "HERV")  ){

        for( alleles in get(paste0(c4type, "_alleles")) ) {

            imputed_GT <- filter(total_C4_CN_allele, SAMPLE==id ) %>% select( starts_with( paste0( c4type, "_"))) %>% as.matrix() %>% as.vector()
            masked_GT  <- filter(total_C4_CN_masked, SAMPLE==id ) %>% select( starts_with( paste0( c4type, "_"))) %>% as.matrix() %>% as.vector()
        }

        
        imputed_GT <- filter(total_C4_CN_allele, SAMPLE==id ) %>% select( starts_with( paste0( c4type, "_"))) %>% as.matrix() %>% as.vector()
        masked_GT <- filter(total_C4_CN_masked, SAMPLE==id ) %>% select( starts_with( paste0( c4type, "_"))) %>% as.matrix() %>% as.vector()


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




###################

rm(list=ls())
library(dplyr); library(tidyr)
load("/kimlab_wd/yuo1996/C4_HLA_ref/try1/7.Accuracy/ver2_kkimsnps/val_C4_result.RDdata")

C4_typed <- read.table("/kimlab_wd/yuo1996/C4_analysis/4.terra_result/total_result.txt", header=T)
vali <- merge( vali, C4_typed[c("SAMPLE", "COUNTS" )], by="SAMPLE")


filter(vali, C4_mat=="no") %>% select( C4, imputed_C4 )# 주로 한개씩 차이-> C4B가 틀린경우가 대부분. 
filter(vali, C4A_mat=="no") %>% select( C4A, imputed_C4A ) 
filter(vali, C4B_mat=="no") %>% select( C4B, imputed_C4B ) # 1 을 2 로 잘못 imputation하는 경우가많음.
filter(vali, C4B_mat=="no", C4B==1)


###### 2copy 정상 call 안된 경우
rm(list=ls())
load("/kimlab_wd/yuo1996/C4_HLA_ref/try1/7.Accuracy/ver2_kkimsnps_crossval/val_C4_result.RDdata")
total_C4_CN_for2copy_check[, -c(1,6)] %>% as.matrix() %>% as.vector() %>% table()
head(total_C4_CN_for2copy_check)
filter(total_C4_CN_for2copy_check ,  C4 != "1_1" | C4A != "1_1" | C4B != "1_1" | HERV != "1_1" )


# C4 A/ B는 틀렸는데 C4는 맞는 경우?
filter(vali, C4_mat=="yes", C4A_mat=="no", , C4B_mat=="no")  # 10 %가 넘음. 



### allele (haploid) level에서 확인.

load("/kimlab_wd/yuo1996/C4_HLA_ref/try1/7.Accuracy/ver2_kkimsnps/val_C4_allele_result.RDdata")

head(C4_allele_val)
filter(C4_allele_val , ID %in%  filter(vali, C4B_mat=="no", C4B==1)$SAMPLE ) # 하나를 잘못잡음다

head(C4_CN_allele)
filter(C4_CN_allele , SAMPLE %in%  filter(vali, C4B_mat=="no", C4B==1)$SAMPLE ) # 하나를 잘못잡음다
filter(C4_CN_masked , SAMPLE %in%  filter(vali, C4B_mat=="no", C4B==1)$SAMPLE ) # 하나를 잘못잡음다


# allele freq checking
nrow(C4_CN_masked)
select( C4_CN_masked , starts_with("C4B")) %>% as.matrix() %>% as.vector() %>% table 
select( C4_CN_masked , starts_with("C4A")) %>% as.matrix() %>% as.vector() %>% table 


########### 정확도가 낮은 특정 set 확인
rm(list=ls())
library(tidyverse)
ver="ver3"
load(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try1/7.Accuracy/" , ver, "/val_C4_result.RDdata"))

vali %>% select( set, C4) %>% table
vali %>% select( set, C4, CNQ_C4) %>% group_by( set ) %>% summarise( cnqmean = mean( CNQ_C4))
filter(vali, set=="set1")$CNQ_C4 %>% summary()
filter(vali, set=="set2")$CNQ_C4 %>% summary()

t <- vali %>% select( set, HERV) %>% table
t <- vali %>% select( set, C4) %>% table # 2 1 3 4

apply( t, 1, function(x) sum( x[c(1,3,5,7)] )) %>% sort
apply( t, 1, function(x) sum( x[c(2,4,6)] )) %>% sort

# SNP check?????




######## haploid level에서 freq 확인
rm(list=ls())
load("/kimlab_wd/yuo1996/C4_HLA_ref/try1/7.Accuracy/ver2_kkimsnps_crossval/val_C4_allele_result.RDdata")
head(total_C4_CN_masked)

total_C4_CN_masked %>% select ( C4_1, C4_2) %>% as.matrix() %>% as.vector() %>% table()
total_C4_CN_masked %>% select ( C4A_1, C4A_2) %>% as.matrix() %>% as.vector() %>% table()
total_C4_CN_masked %>% select ( C4B_1, C4B_2) %>% as.matrix() %>% as.vector() %>% table()
total_C4_CN_masked %>% select ( HERV_1, HERV_2) %>% as.matrix() %>% as.vector() %>% table()

notset2 <- total_C4_CN_masked %>% filter( set!="set1") %>% select (  C4_1, C4_2) %>% as.matrix() %>% as.vector() %>% table()
set2 <- total_C4_CN_masked %>% filter( set=="set1") %>% select (  C4_1, C4_2) %>% as.matrix() %>% as.vector() %>% table()

notset2/sum(notset2)
set2/sum(set2)

total_C4_CN_masked %>% filter( set=="set1") %>% select (  C4_1, C4_2) %>% as.matrix() %>% as.vector() %>% table()

