rm(list=ls())
library(dplyr);library(tidyr)
library(optparse)
library(stringr)

##################################  arguement setting ########################################
args_list <- list(
    make_option( "--snpver" , type="character" , default= "ver3" , help="snpversion", metavar="character")
)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)
# opt=list(); opt$snpver="ver3";
# opt=list(); opt$snpver="ver4";

###############################################################################################


ref=read.table(paste0("/kimlab_wd/yuo1996/C4_HLA_ref/try6/0.vcf_to_plink/", opt$snpver, "/plink.hwe"), header=T)
ref.bim=read.table(paste0("/kimlab_wd/yuo1996/C4_HLA_ref/try6/0.vcf_to_plink/", opt$snpver, "/2nd_plinkQC.bim" ), header=F)
ref_geno=cbind( ref.bim["V4"], ref[c("A1","A2", "GENO")])
head(ref_geno)
nrow(ref_geno) # ver3: 42141



kih=read.table( paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/", opt$snpver, "/plink.hwe" ), header=T)
kih.bim=read.table(paste0( "/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/", opt$snpver,"/2nd_plinkQC.bim" ), header=F)
kih_geno=cbind(kih.bim["V4"] , kih[c("A1","A2", "GENO")])
head(kih_geno)
nrow(kih_geno) # ver3 : 69800


mer = merge(ref_geno, kih_geno, by="V4") %>% filter(A1.x==A1.y )
head(mer)
nrow(mer) # ver3 38018
table(mer$A1.x==mer$A1.y)
table(mer$A2.x==mer$A2.y)
mer[mer$A1.x!=mer$A1.y , ] %>% head(n=100) # 대부분 ref/alt 가 걍 바뀐.




mer_fortest = separate( mer, GENO.x, into=c("AA", "AG", "GG") ,sep="/") %>% separate( GENO.y,into=c("BB", "BC", "CC") ,sep="/") 
mer_fortest$AA=as.numeric(mer_fortest$AA)
mer_fortest$AG=as.numeric(mer_fortest$AG)
mer_fortest$GG=as.numeric(mer_fortest$GG)
mer_fortest$BB=as.numeric(mer_fortest$BB)
mer_fortest$BC=as.numeric(mer_fortest$BC)
mer_fortest$CC=as.numeric(mer_fortest$CC)
head(mer_fortest)
fisherres= data.frame( pos=mer_fortest$V4, pval=NA)

a=1
for(i in 1:nrow( mer_fortest )){
    print(a)
  mat = rbind( mer_fortest[i, c("AA", "AG", "GG")] , mer_fortest[i, c("BB", "BC", "CC")] %>% rename(AA=BB, AG=BC, GG=CC) )  
  fisherres[i, "pval"] = fisher.test(mat)$p
   a=a+1 
}


head(fisherres)
nrow(fisherres)
arrange(fisherres, pval) %>% head(n=50)

filter(fisherres , pval < 0.005) %>% nrow()
filter(fisherres , pval < 0.001) %>% nrow()
filter(fisherres , pval < 10^(-1.5)) %>% nrow()
filter(fisherres , pval < 10^(-6)) %>% nrow()
filter(fisherres , pval < 0.05) %>% nrow()

p=c(1:nrow(fisherres))/nrow(fisherres)
plot( -log10(p), -log10(sort( fisherres$pval)))
abline(0,1, col="red")
hist(fisherres$pval)

filter(fisherres, pos=="32609531") #rs73410656
filter(fisherres, pos=="32433929") #rs117820904
filter(fisherres, pos=="32433423") #rs796666408
filter(fisherres, pos=="32728197") #rs144675976
filter(fisherres, pos=="32705844") #rs377300312
filter(fisherres, pos=="32689386") #rs188844922

filter(fisherres, pos=="32459228") # rs5875374
filter(fisherres, pos=="32459402") # rs1548306
filter(fisherres, pos=="32706948") # rs17206287


# masked_alleles= read.table(paste0("/kimlab_wd/yuo1996/C4_HLA_ref/try6/3.phasing_with_Beagle_C4_HLA/", opt$snpver, "/", opt$haplonetver, "/masked_alleles.txt" ), header=F)

newset = filter(fisherres , pval > 0.001) %>% mutate(chr="chr6") %>% .[c("chr", "pos")] %>% arrange(pos)
nrow(newset)

#intersect(masked_alleles$V2, newset$pos)
head(newset)

path=paste0("/kimlab_wd/yuo1996/C4_HLA_ref/try6_asnpqc_intronrm_with_Eagle_IMPUTE5/KIH_combine/0.vcf_to_plink/", opt$snpver, "/" )
write.table(newset, paste0(path, "/newsnps.txt"), quote=F, col.names=F, row.names=F, sep="\t" )


write.table(fisherres, paste0(path, "/fisherres.txt"), quote=F, col.names=T, row.names=F, sep="\t" )

# foreach(i = 1:nrow(mer_fortest), .combine = "rbind") %do% {
#   mat <- rbind(
#     mer_fortest[i, c("AA", "AG", "GG")],
#     mer_fortest[i, c("BB", "BC", "CC")] %>% rename(AA = BB, AG = BC, GG = CC)
#   )
#   p_val <- fisher.test(mat)$p.value
#   data.frame(pos = mer_fortest$V4[i], pval = p_val)
# } -> result_list

# fisherres <- do.call(rbind, result_list)