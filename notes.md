# sol
module purge
module load bcftools-1.14-gcc-11.2.0 
module load bedtools2-2.30.0-gcc-11.2.0
module load bwa-0.7.17-gcc-12.1.0 


# agave
module load java/8u241
module load trimmomatic/0.33

module purge

module load htstream/0.0.1  
module load admixture/1.3.0
module load angsd/0.921
module load angsd/936 
module load bamtools/2.5.1
module load bcftools/1.9 
module load bedtools2/2.24.0   
module load blast/2.9.0  
module load bowtie2/2.4.5 
module load bwa/0.7.17  
module load fastqc/0.11.7 
module load gatk/4.2.5.0 
module load ncbi-blast/2.6.0
module load numpy/python-3x
module load plink/1.9.0
module load samtools/1.9
module load r/4.1.0

module load picard/2.9.2

module purge
module load picard/2.18.3

picard AddOrReplaceReadGroups \
I=sorted_bam_files/MSB_25201_sorted.bam \
O=sorted_bam_files/MSB_25201__sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=MSB_25201

 picard/1.118
bcftools/1.3.1                                   htstream/0.0.1                                   picard/1.135
bcftools/1.4.0                                   hwloc/1.10.1                                     picard/2.18.3
bcftools/1.9                                     hwloc/1.11.13                                    picard/2.23.7
beagle/4.0                                       hwloc/1.11.2                                     picard/2.3.0
beagle/5.0                                       hwloc/1.11.3                                     picard/2.9.2
beagle-lib/2.1.2                                 hwloc/1.11.7                                     picard/latest
Sequences of adapter

5' Adapter:

5'-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT-3'

3' Adapter:

5'-GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG-3'

>fiveprime
TCTAGCCTTCTCGCAGCACATCCCTTTCTCACATCTAGAGCCACCAGCGGCATAGTAA
>Threeprime
CTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTGCCTACTGATAGAGCATACGGCAGAAGACGAAC

# COMPLIMENTS
>fiveprime
TCTAGCCTTCTCGCAGCACATCCCTTTCTCACATCTAGAGCCACCAGCGGCATAGTAA
>threeprime
CTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTGCCTACTGATAGAGCATACGGCAGAAGACGAAC





fastqc -t 6 \
    raw_fastas/*
    --outdir pre_trim_QC_files/

java -jar /packages/7x/trimmomatic/0.33/trimmomatic.jar PE -t 6 

./trim.sh -i /scratch/dnjacks4/cardinalis/filenames.txt -p /scratch/dnjacks4/cardinalis/raw_fastas/ -f 1.fq.gz -r 2.fq.gz -t 7
/scratch/dnjacks4/cardinalis/raw_fastas/

java -jar /packages/7x/trimmomatic/0.33/trimmomatic.jar PE   /scratch/dnjacks4/cardinalis/raw_fastas/UWBM_77978_CKDN230006124-1A_H5LY7DSX7_L2_1.fq.gz /scratch/dnjacks4/cardinalis/raw_fastas/UWBM_77978_CKDN230006124-1A_H5LY7DSX7_L2_2.fq.gz \
    -baseout /scratch/dnjacks4/cardinalis/raw_fastas/UWBM_77978_CKDN230006124-1A_H5LY7DSX7_L2_2.fq.gz_trimmed.fq.gz \
    ILLUMINACLIP:adapters.txt:1:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>trim_and_QC_log.txt


PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>
~/whole_genome_bioinformatics/trim-and-QC.sh -i /scratch/dnjacks4/cardinalis/filenames.txt -p /scratch/dnjacks4/cardinalis/raw_fastas/ -f L1_1.fq.gz -r L1_2.fq.gz -t 6

[-i] Sample list
[-p] Path to fastq files
[-f] Suffix of forward reads (e.g. R1_000.fastq.gz)
[-r] Suffix of reverse reads (e.g. R2_000.fastq.gz) \


fastqc -t 6 \
    trimmed_fastas/* --outdir post_trim_QC_files/

ls -d $PWD/trimmed_fastas/* > trimmedfilenames.txt

# i'm redoing this, because all individuals with multiple fastas came back with 100% missing snp frequencies
# i made a new directory with only the trimmed fasta with the largest file size for each individual, trimmed_fastas_select

ls -d $PWD/trimmed_fastas_select/* > trimmedfilenames.txt

./sort.sh -i /scratch/dnjacks4/cardinalis/sampleids.txt -p trimmed_fastas_select/ -r referencedata/cardinalreferencegenome.fa

ref="/scratch/dnjacks4/cardinalis/referencedata/cardinalreferencegenome.fa"
bamdir="/scratch/dnjacks4/cardinalis/sorted_bam_files/"
ID="cardinalis"
bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*_sorted_RGadded_dupmarked.bam | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf


plink --vcf /scratch/dnjacks4/cardinalis/cardinalis_snps_multiallelic.vcf --allow-extra-chr --missing --cluster-missing --freq

plink --vcf /scratch/dnjacks4/cardinalis/cardinalis_snps_multiallelic.vcf --allow-extra-chr --missing --cluster missing --within individualnames.txt --freq

plink --vcf /scratch/dnjacks4/cardinalis/cardinalis.vcf --allow-extra-chr --missing --cluster missing --within cluster_pop.txt --freq





#filters by quality
bcftools view -i 'QUAL>100' cardinalis_snps_multiallelic.vcf > cardinalis_qualitysort.vcf

#filters by depth and removes indels
vcftools --vcf cardinalis_qualitysort.vcf --min-meanDP 2 --max-meanDP 8 --remove-indels --recode --out cardinalis_filtered

# After filtering, kept 37653 out of a possible 585284 Sites

sed -i 's/\_//g' cardinalis_filtered.recode.vcf
sed -i 's/\UWBM/BM/g' cardinalis_filtered.recode.vcf


plink --vcf cardinalis_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid --out cardinalis_filtered_mind2

# 7973 variants and 21 people pass filters and QC.

plink --vcf cardinalis_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.1 --maf 0.01 --recode vcf-iid --out cardinalis_filtered_cleaned

# 22876 variants and 3 people pass filters and QC.


cp sula_flightless_cleaned.vcf sula_flightless_cleaned_zip.vcf
bgzip sula_flightless_cleaned_zip.vcf

bcftools index sula_flightless_cleaned_zip.vcf.gz


~/vcf2phylip/vcf2phylip.py -i cardinalis_filtered_mind2.vcf
mkdir pruned 

cd pruned 

python ~/sula/filter_invariants_all.py ../cardinalis_filtered_mind2.min4.phy
mv variantsites.phy ../variantsites_mind2.phy 
mv variantsites_kept.txt ../variantsites_mind2_kept.txt 
cd .. 
rm -r pruned

echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-4528130' > partitionfile.txt
~/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n cardinalis_1 -q /scratch/dnjacks4/cardinalis/mydata/partitionfile.txt -s /scratch/dnjacks4/cardinalis/mydata/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V


echo -e "PYRR" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "PYRR" >> pops.txt 
echo -e "NOCA" >> pops.txt 
echo -e "NOCA" >> pops.txt 


MSB25201    PYRR      
NOCA003 NOCA
NOCA004 NOCA
NOCA006 NOCA
NOCA008 NOCA
NOCA012 NOCA
NOCA013 NOCA
PYRR003 PYRR
PYRR004 PYRR
PYRR006 PYRR
PYRR007 PYRR
PYRR009 PYRR
PYRR011 PYRR
BM100619    NOCA
BM100620    NOCA
BM100621    NOCA
BM103345    NOCA
BM103346    PYRR
BM77548 PYRR
BM77718 PYRR
BM77780 PYRR
BM77781 PYRR
BM77856 NOCA
BM77978 NOCA


# required to run:
# library("devtools")
# install_github("zhengxwen/gdsfmt")
# install_github("zhengxwen/SNPRelate")

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf -o /scratch/dnjacks4/cardinalis/PCA/all/ -p /scratch/dnjacks4/cardinalis/PCA/all/pops.txt -n all -s y

MSB25201    PYRR      
PYRR003 PYRR
PYRR004 PYRR
PYRR006 PYRR
PYRR007 PYRR
PYRR009 PYRR
PYRR011 PYRR
BM103346    PYRR
BM77548 PYRR
BM77718 PYRR
BM77780 PYRR
BM77781 PYRR


bcftools view -s MSB25201,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,BM103346,BM77548,BM77718,BM77780,BM77781 /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --force-samples > /scratch/dnjacks4/cardinalis/mydata/pyrr_pca.vcf

cd /scratch/dnjacks4/cardinalis/PCA/pyrr
echo -e "rural" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/mydata/pyrr_pca.vcf -o /scratch/dnjacks4/cardinalis/PCA/pyrr/ -p /scratch/dnjacks4/cardinalis/PCA/pyrr/pops.txt -n pyrr -s y

  
NOCA003 NOCA
NOCA004 NOCA
NOCA006 NOCA
NOCA008 NOCA
NOCA012 NOCA
NOCA013 NOCA
BM100619    NOCA
BM100620    NOCA
BM100621    NOCA
BM103345    NOCA
BM77856 NOCA
BM77978 NOCA

bcftools view -s NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,BM100619,BM100620,BM100621,BM103345,BM77856,BM77978 /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --force-samples > /scratch/dnjacks4/cardinalis/mydata/noca_pca.vcf


cd /scratch/dnjacks4/cardinalis/PCA/noca

echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 


~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/mydata/noca_pca.vcf -o /scratch/dnjacks4/cardinalis/PCA/noca/ -p /scratch/dnjacks4/cardinalis/PCA/noca/pops.txt -n noca -s y



### pyrrhuloxia male 
cd /scratch/dnjacks4/cardinalis/PCA/pyrr_m

bcftools view -s MSB25201,PYRR004,PYRR006,PYRR007,BM77718,BM77781 /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --force-samples > /scratch/dnjacks4/cardinalis/mydata/pyrr_m_pca.vcf

echo -e "rural" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/mydata/pyrr_m_pca.vcf -o /scratch/dnjacks4/cardinalis/PCA/pyrr_m/ -p /scratch/dnjacks4/cardinalis/PCA/pyrr_m/pops.txt -n pyrr_m -s y


### pyrrhuloxia female
cd /scratch/dnjacks4/cardinalis/PCA/pyrr_f

bcftools view -s PYRR003,PYRR009,PYRR011,BM103346,BM77548,BM77780 /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --force-samples > /scratch/dnjacks4/cardinalis/mydata/pyrr_f_pca.vcf

echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "urban" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 
echo -e "rural" >> pops.txt 

PYRR_003        PYRR    female
PYRR_009        PYRR    female
PYRR_011        PYRR    female
UWBM_103346     PYRR    female
UWBM_77548      PYRR    female
UWBM_77780      PYRR    female

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/mydata/pyrr_f_pca.vcf -o /scratch/dnjacks4/cardinalis/PCA/pyrr_f -p /scratch/dnjacks4/cardinalis/PCA/pyrr_f/pops.txt -n pyrr_f -s y


### fst

#### species
cd /scratch/dnjacks4/cardinalis/fst/species


echo -e "MSB25201" >> pyrr_pops.txt 
echo -e "PYRR003" >> pyrr_pops.txt 
echo -e "PYRR004" >> pyrr_pops.txt 
echo -e "PYRR006" >> pyrr_pops.txt 
echo -e "PYRR007" >> pyrr_pops.txt 
echo -e "PYRR009" >> pyrr_pops.txt 
echo -e "PYRR011" >> pyrr_pops.txt 
echo -e "BM103346" >> pyrr_pops.txt 
echo -e "BM77548" >> pyrr_pops.txt 
echo -e "BM77718" >> pyrr_pops.txt 
echo -e "BM77780" >> pyrr_pops.txt 
echo -e "BM77781" >> pyrr_pops.txt 

echo -e "NOCA003" >> noca_pops.txt 
echo -e "NOCA004" >> noca_pops.txt 
echo -e "NOCA006" >> noca_pops.txt 
echo -e "NOCA008" >> noca_pops.txt 
echo -e "NOCA012" >> noca_pops.txt 
echo -e "NOCA013" >> noca_pops.txt 
echo -e "BM100619" >> noca_pops.txt 
echo -e "BM100620" >> noca_pops.txt 
echo -e "BM100621" >> noca_pops.txt 
echo -e "BM103345" >> noca_pops.txt 
echo -e "BM77856" >> noca_pops.txt 
echo -e "BM77978" >> noca_pops.txt 


vcftools --vcf /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --weir-fst-pop pyrr_pops.txt --weir-fst-pop noca_pops.txt --fst-window-size 25000 --fst-window-step 25000 

sed -i 's/JACDOX//g' out.windowed.weir.fst


library(qqman)
fst<-read.table("/scratch/dnjacks4/cardinalis/fst/species/out.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "species.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()


#### pyrr
cd /scratch/dnjacks4/cardinalis/fst/pyrr


echo -e "MSB25201" > rural.txt 
echo -e "PYRR003" > urban.txt 
echo -e "PYRR004" >> urban.txt 
echo -e "PYRR006" >> urban.txt 
echo -e "PYRR007" >> urban.txt 
echo -e "PYRR009" >> urban.txt 
echo -e "PYRR011" >> urban.txt 
echo -e "BM103346" >> rural.txt 
echo -e "BM77548" >> rural.txt 
echo -e "BM77718" >> rural.txt 
echo -e "BM77780" >> rural.txt 
echo -e "BM77781" >> rural.txt 

vcftools --vcf /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --weir-fst-pop rural.txt --weir-fst-pop urban.txt --fst-window-size 1000 --fst-window-step 1000 

sed -i 's/JACDOX//g' out.windowed.weir.fst

library(qqman)
fst<-read.table("/scratch/dnjacks4/cardinalis/fst/pyrr/out.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "pyrr_1000window.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()



vcftools --vcf /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --weir-fst-pop rural.txt --weir-fst-pop urban.txt 

sed -i 's/JACDOX//g' out.weir.fst

library(qqman)
fst<-read.table("/scratch/dnjacks4/cardinalis/fst/pyrr/out.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "pyrr_snps.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()

#### noca
cd /scratch/dnjacks4/cardinalis/fst/noca

echo -e "NOCA003" > urban.txt 
echo -e "NOCA004" >> urban.txt 
echo -e "NOCA006" >> urban.txt 
echo -e "NOCA008" >> urban.txt 
echo -e "NOCA012" >> urban.txt 
echo -e "NOCA013" >> urban.txt 
echo -e "BM100619" > rural.txt 
echo -e "BM100620" >> rural.txt 
echo -e "BM100621" >> rural.txt 
echo -e "BM103345" >> rural.txt 
echo -e "BM77856" >> rural.txt 
echo -e "BM77978" >> rural.txt 

vcftools --vcf /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --weir-fst-pop rural.txt --weir-fst-pop urban.txt --fst-window-size 1000 --fst-window-step 1000 

sed -i 's/JACDOX//g' out.windowed.weir.fst


library(qqman)
fst<-read.table("/scratch/dnjacks4/cardinalis/fst/noca/out.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "noca.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()

vcftools --vcf /scratch/dnjacks4/cardinalis/mydata/cardinalis_filtered.recode.vcf --weir-fst-pop rural.txt --weir-fst-pop urban.txt 

sed -i 's/JACDOX//g' out.weir.fst

library(qqman)
fst<-read.table("/scratch/dnjacks4/cardinalis/fst/noca/out.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "noca_snps.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()

### snpeff
# https://hackmd.io/@tlama/outlierFST#Evaluate-outliers-in-R-see-Rmarkdown-here

# downloaded snpeff

# download snpeff databases for chicken and great tit
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip



