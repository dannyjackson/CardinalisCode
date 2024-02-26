# BayesScan
# Whole-genome SNP data unravel population structure andsignatures of selection for black plumage of indigenous chickenbreeds from Jiangxi province, China
# https://rpubs.com/lbenestan/outlier
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7768260/


interactive -t 1-00:00 --mem-per-cpu 300000
interactive -n 12 -t 3-00:00 #bayescan
interactive -n 12 -t 1-00:00 #bayescan
interactive -t 1-00:00 --mem-per-cpu 300000
interactive -t 7-00:00 --mem-per-cpu 500000

module purge
module load r-4.2.2-gcc-11.2.0
module load bayescan/2.1
module load htslib-1.16-gcc-11.2.0
module load bcftools-1.14-gcc-11.2.0

/scratch/dnjacks4/cardinalis/to_cardinalis/bayescan

module load openjdk-1.8.0_265-b01-gcc-12.1.0
module load sqlite-3.38.5-gcc-11.2.0
module load proj-8.2.1-gcc-11.2.0
pymodule load gdal-3.4.3-gcc-11.2.0
module load geos-3.9.1-gcc-11.2.0
module load gcc-12.1.0-gcc-11.2.0
module load libxc-5.1.7-gcc-12.1.0
module load libvterm-0.1.4-gcc-11.2.0

cd /scratch/dnjacks4/cardinalis/to_parus/bayescan

# filter vcf 

cp /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_snps_multiallelic.vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/parussnpsmultiallelic.vcf

sed -i 's/_//g' /scratch/dnjacks4/cardinalis/to_parus/vcfs/parussnpsmultiallelic.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/parussnpsmultiallelic.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out parusfiltered.geno25.maf1

grep 'CHROM' parusfiltered.geno25.maf1.vcf

Before main variant filters, 24 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.081039.
38271028 variants removed due to missing genotype data (--geno).
21084 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
14659 variants and 24 people pass filters and QC.
Note: No phenotypes present.
--recode vcf-iid to parus_filtered.geno25.maf1.vcf ... done.

4205 of 14659 variants removed.

10454 snps remaining

/scratch/dnjacks4/cardinalis/to_parus/bayescan/parus_filtered.geno25.maf1.vcf
# Genome-wide scan for selection signatures reveals novel insights into the adaptive capacity in local North African cattle

# First, we excluded rare SNPs with low minor allele frequencies (MAF) < 0.05. Then, the whole genotype dataset was subjected to linkage disequilibrium (LD) pruning using the default parameters of PLINK (SNP window size:50, step 5 SNPs, r2: 0.5). In total, 38,464 SNPs spread over all autosomal chromosomes were finally considered for population structure analyses.

# Bayescan uses a reversible-jump Markov Chain Monte Carlo to separate locus- specific effects of selection from population-specific effects of demography. Outliers are those loci that require the locus-specific component to explain observed genetic diversity. For the Markov chain Monte Carlo (MCMC) algorithm we used 20 pilot runs of 5,000 iterations, a burn-in of 50,000 iterations, a thinning interval of 10 (5,000 iterations were used for the estimation of posterior odds) with a resulting total number of 100,000 iterations. To control the number of false positives, significant SNPs were defined by applying a q-value threshold of 0.05.

# filter vcf by MAF < 0.10
# linkage disequilibrium (LD) pruning using the default parameters of PLINK (SNP window size:50, step 5 SNPs, r2: 0.5).

echo "INDIVIDUALS,STRATA" > /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt

echo 'MSB25201,PYRR_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'NOCA003,NOCA_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'NOCA004,NOCA_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'NOCA006,NOCA_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'NOCA008,NOCA_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'NOCA012,NOCA_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'NOCA013,NOCA_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'PYRR003,PYRR_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'PYRR004,PYRR_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'PYRR006,PYRR_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'PYRR007,PYRR_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'PYRR009,PYRR_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'PYRR011,PYRR_U' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM100619,NOCA_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM100620,NOCA_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM100621,NOCA_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM103345,NOCA_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM103346,PYRR_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM77548,PYRR_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM77718,PYRR_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM77780,PYRR_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM77781,PYRR_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM77856,NOCA_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt
echo 'UWBM77978,NOCA_R' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt


sed -i 's/,/\t/g' /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt

cd /scratch/dnjacks4/cardinalis/to_parus/bayescan



cd /scratch/dnjacks4/cardinalis/to_parus/bayescan


R


library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)


vcf <- read.vcfR("/scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf")

pop_map <- read.table("/scratch/dnjacks4/cardinalis/to_parus/reference_lists/bayescan_popfile.txt", header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn="cardinalis.filtered.bsc")


mkdir all_filtered
mkdir all_filtered_snp



bayescan -n 5000 -burn 50000 -pr_odds 10000 /scratch/dnjacks4/cardinalis/to_parus/bayescan/cardinalis.filtered.bsc -threads 12 -od /scratch/dnjacks4/cardinalis/to_parus/bayescan/all_filtered


grep -v "#" /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf | awk '{ print $1, $2 }' | tr ' ' '_' > /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.loci.txt

# pyrrhuloxia 
bcftools view -s 'MSB25201,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,UWBM103346,UWBM77548,UWBM77718,UWBM77780,UWBM77781' /scratch/dnjacks4/cardinalis/to_parus/vcfs/parussnpsmultiallelic.vcf > /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr_all.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr_all.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out pyrr.filtered.geno25.maf1


echo "INDIVIDUALS,STRATA" > /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt

echo 'MSB25201,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR003,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR004,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR006,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR007,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR009,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR011,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'UWBM103346,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'UWBM77548,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'UWBM77718,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'UWBM77780,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'UWBM77781,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt


sed -i 's/,/\t/g' /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt



/scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr_bayes_test.vcf

R

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

vcf <- read.vcfR("/scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr.filtered.geno25.maf1.vcf")

pop_map <- read.table("/scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt", header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn="pyrr.filtered.bsc")

mkdir filtered 

bayescan -n 5000 -burn 50000 -pr_odds 10000 /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr.filtered.bsc -od /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrrhuloxia/filtered/





grep -v "#" /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr.filtered.geno25.maf1.vcf | awk '{ print $1, $2 }' | tr ' ' '_' > /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr.filtered.geno25.maf1.loci.txt

R 

bayescan=read.table("pyrr.filtered_fst.txt") 


SNPb=read.table("/scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr.filtered.geno25.maf1.loci.txt",header=FALSE)

bayescan=cbind(SNPb, bayescan) 

colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

write.table(bayescan, "pyrr-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE) 

attach(bayescan)
class(bayescan$Q_VALUE)

bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 

positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

xtabs(data=bayescan, ~SELECTION) 

write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F) 




#northern cardinals 
bcftools view -s 'NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,UWBM100619,UWBM100620,UWBM100621,UWBM103345,UWBM77856,UWBM77978'  /scratch/dnjacks4/cardinalis/to_parus/vcfs/parussnpsmultiallelic.vcf > /scratch/dnjacks4/cardinalis/to_parus/bayescan/noca_all.vcf


plink --vcf /scratch/dnjacks4/cardinalis/to_parus/bayescan/noca_all.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out noca.filtered.geno25.maf1


echo "INDIVIDUALS,STRATA" > /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt

echo 'NOCA003,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'NOCA004,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'NOCA006,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'NOCA008,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'NOCA012,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'NOCA013,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'UWBM100619,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'UWBM100620,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'UWBM100621,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'UWBM103345,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'UWBM77856,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt
echo 'UWBM77978,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt


sed -i 's/,/\t/g' /scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt



R

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

vcf <- read.vcfR("/scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/noca.filtered.geno25.maf1.vcf")

pop_map <- read.table("/scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt", header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn="noca.filtered.bsc")




mkdir /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/filtered/

bayescan -n 5000 -burn 50000 -pr_odds 10000 /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/noca.filtered.bsc -threads 12 -od /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/filtered/






grep -v "#" /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/noca.filtered.geno25.maf1.vcf | awk '{ print $1, $2 }' | tr ' ' '_' > /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/noca.filtered.geno25.maf1.loci.txt

R 

bayescan=read.table("noca.filtered_fst.txt") 


SNPb=read.table("/scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/noca.filtered.geno25.maf1.loci.txt",header=FALSE)

bayescan=cbind(SNPb, bayescan) 

colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

write.table(bayescan, "noca-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE) 

attach(bayescan)
class(bayescan$Q_VALUE)

bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 

positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

xtabs(data=bayescan, ~SELECTION) 

write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F) 



# install_github("zhengxwen/gdsfmt")
# install_gi thub("zhengxwen/SeqArray")

# library("gdsfmt")
# load_all('/home/dnjacks4/R/x86_64-pc-linux-gnu-library/4.2/SeqArray/')

# library("SeqArray")
# setwd("/scratch/dnjacks4/cardinalis/to_parus/bayescan") 

# remotes::install_github("thierrygosselin/radiator")

# /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_snps_multiallelic.vcf




 awk '$2> 0.91 || NR==1' /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrrhuloxia/filtered_fullrun/pyrr.filtered_fst.txt

  awk '$2> 0.20 || NR==1' /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/filtered_fullrun/noca.filtered_fst.txt

 awk '{print $2}' /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrrhuloxia/filtered_fullrun/pyrr.filtered_fst.txt | sort -u

awk '{print $2}' /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/filtered_fullrun/noca.filtered_fst.txt | sort -u




# terminated at 45%
[dnjacks4@c035:/scratch/dnjacks4/cardinalis/to_parus/bayescan/all_filtered/test]$ bayescan -n 5000 -burn 50000 -pr_odds 10000 /scratch/dnjacks4/cardinalis/to_parus/bayescan/cardinalis.filtered.bsc -snp -od /scratch/dnjacks4/cardinalis/to_parus/bayescan/all_filtered/test  


# Terminated at 83% of pilot runs
[dnjacks4@c006:/scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrrhuloxia]$ bayescan -n 5 -burn 50 -pr_odds 10 /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr.filtered.bsc -od /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrrhuloxia/filtered/

# finished! but still lacking the q-value...
[dnjacks4@c038:/scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal]$ bayescan -n 5000 -burn 50000 -pr_odds 10000 /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/noca.filtered.bsc -threads 12 -od /scratch/dnjacks4/cardinalis/to_parus/bayescan/northerncardinal/filtered_fullrun/

# finished! but still lacking the q-value...
[dnjacks4@c053:/scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrrhuloxia/filtered_fullrun]$ bayescan -n 5000 -burn 50000 -pr_odds 10000 /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr.filtered.bsc -od /scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrrhuloxia/filtered_fullrun/ 








### repeating pyrr when aligned to scott edwards' shitty cardinal genome:

bcftools view -s 'MSB25201,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,UWBM103346,UWBM77548,UWBM77718,UWBM77780,UWBM77781' /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalissnpsmultiallelic.vcf > /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr_all.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr_all.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr.filtered.geno25.maf1


R

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

vcf <- read.vcfR("/scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr.filtered.geno25.maf1.vcf")

pop_map <- read.table("/scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt", header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn="pyrr.filtered.bsc")

mkdir filtered 

bayescan -n 5000 -burn 50000 -pr_odds 10000 /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr.filtered.bsc -od /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrrhuloxia/filtered/





grep -v "#" /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr.filtered.geno25.maf1.vcf | awk '{ print $1, $2 }' | tr ' ' '_' > /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr.filtered.geno25.maf1.loci.txt

R 

bayescan=read.table("pyrr.filtered_fst.txt") 


SNPb=read.table("/scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr.filtered.geno25.maf1.loci.txt",header=FALSE)

bayescan=cbind(SNPb, bayescan) 

colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

write.table(bayescan, "pyrr-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE) 

attach(bayescan)
class(bayescan$Q_VALUE)

bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 

# balancing   neutral 
# 10768     57013 
   
positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

xtabs(data=bayescan, ~SELECTION) 

write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F) 




## northern cardinals


bcftools view -s 'NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,UWBM100619,UWBM100620,UWBM100621,UWBM103345,UWBM77856,UWBM77978' /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalissnpsmultiallelic.vcf > /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca_all.vcf


plink --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca_all.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1


R

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

vcf <- read.vcfR("/scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1.vcf")

pop_map <- read.table("/scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt", header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn="noca.filtered.bsc")

mkdir northerncardinals 
cd northerncardinals
mkdir filtered 

bayescan -n 5000 -burn 50000 -pr_odds 10000 /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.bsc -od /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/northerncardinals/filtered/





grep -v "#" /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1.vcf | awk '{ print $1, $2 }' | tr ' ' '_' > /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1.txt

R 

bayescan=read.table("noca.filtered_fst.txt") 


SNPb=read.table("/scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1.txt",header=FALSE)

bayescan=cbind(SNPb, bayescan) 

colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

write.table(bayescan, "noca-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE) 

attach(bayescan)
class(bayescan$Q_VALUE)

bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 

# balancing   neutral 
# 10504     54688 

positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

xtabs(data=bayescan, ~SELECTION) 

write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F) 
