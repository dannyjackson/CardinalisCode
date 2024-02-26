bgzip /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf

bcftools index /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf.gz
bcftools view -r 'VYXE01005221.1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf.gz > /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/VYXE01005221.vcf

cd /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/all

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/VYXE01005221.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt


sed -i 's/VYXE//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fst$WEIR_AND_COCKERHAM_FST[fst$WEIR_AND_COCKERHAM_FST < 0] <- 0 
fst$POS <- (fst$POS - 2174500)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)
mydf2 <- subset(mydf, SNP>8)
mydf <- mydf2

pdf(file = "VYXE01005221_fst_all.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="POS",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()
q()
n

cd ../urban_urban

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/VYXE01005221.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt


sed -i 's/VYXE//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fst$WEIR_AND_COCKERHAM_FST[fst$WEIR_AND_COCKERHAM_FST < 0] <- 0 
fst$POS <- (fst$POS - 2174500)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)
mydf2 <- subset(mydf, SNP>8)
mydf <- mydf2

pdf(file = "VYXE01005221_fst_urban.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="POS",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()
q()
n



cd ../rural_rural

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/VYXE01005221.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt


sed -i 's/VYXE//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fst$WEIR_AND_COCKERHAM_FST[fst$WEIR_AND_COCKERHAM_FST < 0] <- 0 
fst$POS <- (fst$POS - 2174500)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)
mydf2 <- subset(mydf, SNP>8)
mydf <- mydf2



pdf(file = "VYXE01005221_fst_rural.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()
q()
n





cd ../pyrr

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/VYXE01005221.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt


sed -i 's/VYXE//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fst$WEIR_AND_COCKERHAM_FST[fst$WEIR_AND_COCKERHAM_FST < 0] <- 0 
fst$POS <- (fst$POS - 2174500)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)
mydf2 <- subset(mydf, POS>0)
mydf <- mydf2



pdf(file = "VYXE01005221_fst_pyrr.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()
q()
n



cd ../noca

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/VYXE01005221.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt


sed -i 's/VYXE//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fst$WEIR_AND_COCKERHAM_FST[fst$WEIR_AND_COCKERHAM_FST < 0] <- 0 
fst$POS <- (fst$POS - 2174500)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)
mydf2 <- subset(mydf, POS>0)
mydf <- mydf2



pdf(file = "VYXE01005221_fst_noca.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()
q()
n

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/VYXE01005221.vcf -o /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 