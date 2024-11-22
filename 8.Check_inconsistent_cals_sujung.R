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

# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8QUAL20"; opt$haplonetver="ver1_sujung"; opt$crossval="crossval1"

###############################################################################################
load(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval, "/val_C4_result.RDdata") )
# load("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/C4_genomestrip/C4_HERV_AB_GT_dist.RData")

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
c4type="C4"
make_dn_contigency_table <- function( c4type) {

    c4_type <- vali[ c( which( colnames(vali) ==c4type) , which( colnames(vali) == paste0( "imputed_", c4type)  ) )] 
    colnames(c4_type) = c("typed", "imputed")

    contable = table(c4_type)    
#    contable=rbind(contable, colSums(contable) )
#    contable=cbind(contable,  rowSums(contable) )
#    rownames(contable)[ length(rownames(contable) )]="sum"
#    colnames(contable)[ length(colnames(contable) )]="sum"

if( dim(contable)[1]!=dim(contable)[2] ){

contable=cbind(contable, rep(0,8))

}

    contable = cbind( contable , diag(contable)/rowSums(contable) )
    contable = rbind( contable , diag(contable)/colSums(contable) )

contable

    return(contable)    
}




make_dn_contigency_table( "C4")
make_dn_contigency_table( "C4A")
make_dn_contigency_table( "C4B")
make_dn_contigency_table( "HERV")

head(vali)
cor(vali$C4, vali$imputed_C4)#^2
cor(vali$C4A, vali$imputed_C4A)#^2
cor(vali$C4B, vali$imputed_C4B)#^2
cor(vali$HERV, vali$imputed_HERV)#^2




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

# opt=list(); opt$snpver="ver4"; opt$region="start31.7end32.2_highLDfiltering_50_5_0.8QUAL20"; opt$haplonetver="ver1_sujung"; opt$crossval="crossval1"

###############################################################################################

load(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval, "/val_C4_result.RDdata") )
load(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/7.Accuracy/", opt$snpver, "/", opt$region, "_", opt$haplonetver, "/", opt$crossval, "/val_C4_allele_result.RDdata") )
# load("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/C4_genomestrip/C4_HERV_AB_GT_dist.RData")


if( names( table( total_C4_CN_allele$SAMPLE==total_C4_CN_masked$SAMPLE )  )){
    SAMPLE <- total_C4_CN_allele$SAMPLE
    print("can go!!")
}else{
    print("stop!!!!!!!!!!!!!!!!!!!!!")
}



#C4_alleles <- total_C4_CN_masked %>% select( starts_with("C4_")) %>% as.matrix() %>% as.vector() %>% table() %>% names()
C4A_alleles <- total_C4_CN_masked %>% select( starts_with("C4A_")) %>% as.matrix() %>% as.vector() %>% table() %>% names()
C4B_alleles <- total_C4_CN_masked %>% select( starts_with("C4B_")) %>% as.matrix() %>% as.vector() %>% table() %>% names()
HERV_alleles <- total_C4_CN_masked %>% select( starts_with("HERV_")) %>% as.matrix() %>% as.vector() %>% table() %>% names()

#all <- c(C4_alleles , C4A_alleles , C4B_alleles , HERV_alleles)
all <- c( C4A_alleles , C4B_alleles , HERV_alleles)
acc_for_alleles <- NULL

#ex 
data <- total_C4_CN_allele; c4type="C4A"; alleles_set=C4A_alleles
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


#allele_count_masked_C4 <- mkallele_count_table( "C4" , C4_alleles, total_C4_CN_masked ) 
allele_count_masked_C4A <- mkallele_count_table( "C4A" , C4A_alleles, total_C4_CN_masked ) 
allele_count_masked_C4B <- mkallele_count_table( "C4B" , C4B_alleles, total_C4_CN_masked ) 
allele_count_masked_HERV <- mkallele_count_table( "HERV" , HERV_alleles, total_C4_CN_masked ) 


#allele_count_imputed_C4 <- mkallele_count_table( "C4" , C4_alleles, total_C4_CN_allele ) 
allele_count_imputed_C4A <- mkallele_count_table( "C4A" , C4A_alleles, total_C4_CN_allele ) 
allele_count_imputed_C4B <- mkallele_count_table( "C4B" , C4B_alleles, total_C4_CN_allele ) 
allele_count_imputed_HERV <- mkallele_count_table( "HERV" , HERV_alleles, total_C4_CN_allele ) 

c4type="C4A"; allele="C4A_copy1"
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


#calculate_per( "C4" , "C4_copy1")  %>% round(., digits = 4)*100
#calculate_per( "C4" , "C4_copy2") %>% round(., digits = 4)*100
#calculate_per( "C4" , "C4_copy3") %>% round(., digits = 4)*100# cor이 안좋다?
#calculate_per( "C4" , "C4_copy4") %>% round(., digits = 4)*100# 베리굿!!!

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




nrow(allele_count_masked_C4)

nrow(allele_count_masked_C4A)

freq = function (table){
a=table( table)
if(length(a)==3){
(a[2]+a[3]*2)/(1537*2) %>% return() 
}else if( length(a)==2){
(a[2])/(1537*2) %>% return() 
}
}

#freq(allele_count_masked_C4$C4_copy1) *100  %>%  round(. , digits=2)
#freq(allele_count_masked_C4$C4_copy2)  *100  %>%  round(. , digits=2)
#freq(allele_count_masked_C4$C4_copy3)  *100  %>%  round(. , digits=2)
#freq(allele_count_masked_C4$C4_copy4) *100  %>%  round(. , digits=2)

freq(allele_count_masked_C4A$C4A_copy0) + freq(allele_count_masked_C4A$C4A_copy1) + freq(allele_count_masked_C4A$C4A_copy2) +freq(allele_count_masked_C4A$C4A_copy3)

round( freq(allele_count_masked_C4A$C4A_copy0) , digits=3) #  %>%  round(. , digits=2)
round(freq(allele_count_masked_C4A$C4A_copy1) , digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_C4A$C4A_copy2) , digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_C4A$C4A_copy3) , digits=3) #   %>%  round(. , digits=2)

round(freq(allele_count_masked_C4B$C4B_copy0), digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_C4B$C4B_copy1), digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_C4B$C4B_copy2), digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_C4B$C4B_copy3) , digits=3) #   %>%  round(. , digits=2)

round(freq(allele_count_masked_HERV$HERV_copy0), digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_HERV$HERV_copy1), digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_HERV$HERV_copy2) , digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_HERV$HERV_copy3) , digits=3) #   %>%  round(. , digits=2)
round(freq(allele_count_masked_HERV$HERV_copy4) , digits=3) #   %>%  round(. , digits=2)



