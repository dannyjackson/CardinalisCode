interactive -t 2-00:00 

module purge
module load r-4.0.2-gcc-11.2.0 #r-4.2.2-gcc-11.2.0
module load bayescan/2.01
module load htslib-1.16-gcc-11.2.0
module load bcftools-1.14-gcc-11.2.0
module load vcftools-0.1.14-gcc-11.2.0

echo 'MSB25201,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'NOCA003,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'NOCA004,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'NOCA006,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'NOCA008,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'NOCA012,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'NOCA013,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'PYRR003,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'PYRR004,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'PYRR006,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'PYRR007,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'PYRR009,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'PYRR011,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM100619,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM100620,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM100621,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM103345,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM103346,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM77548,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM77718,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM77780,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM77781,PYRR' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM77856,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
echo 'UWBM77978,NOCA' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt
sed -i 's/,/\t/g' /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt




echo 'MSB25201' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt
echo 'UWBM103346' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt
echo 'UWBM77548' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt
echo 'UWBM77718' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt
echo 'UWBM77780' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt
echo 'UWBM77781' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt

echo 'PYRR003' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR004' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR006' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR007' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR009' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR011' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt

echo 'NOCA003' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA004' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA006' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA008' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA012' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA013' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt

echo 'UWBM100619' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt
echo 'UWBM100620' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt
echo 'UWBM100621' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt
echo 'UWBM103345' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt
echo 'UWBM77856' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt
echo 'UWBM77978' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt


echo 'MSB25201' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'UWBM103346' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'UWBM77548' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'UWBM77718' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'UWBM77780' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'UWBM77781' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt

echo 'PYRR003' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR004' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR006' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR007' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR009' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR011' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt

echo 'NOCA003' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'NOCA004' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'NOCA006' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'NOCA008' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'NOCA012' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'NOCA013' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt

echo 'UWBM100619' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'UWBM100620' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'UWBM100621' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'UWBM103345' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'UWBM77856' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt
echo 'UWBM77978' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt


vcftools --vcf /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt

sed -i 's/NC//g' out.weir.fst
sed -i 's/NW//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "pyrr_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()

vcftools --vcf /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt

awk ' gsub(".1","1",$4){for (i=1;i<=NF;i++){printf("%-9s",$i)}print "" }' out.weir.fst | head

awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

sed -i 's/NC//g' out.weir.fst
sed -i 's/NW//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "noca_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()

awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_parus/fst/pyrrhuloxia/out.weir.fst | awk '$3> 0.8 || NR==1' 

CHROM	POS	WEIR_AND_COCKERHAM_FST
# gene id 39089571
NC040875.1	482	0.868687
# gene id 39089587
NC040875.1	9649	1
NC040875.1	9696	1
NC040875.1	9702	1
NC040875.1	9736	1
# gene id 39089588
NC040875.1	10009	0.877637
# gene id 39089591
NC040875.1	12170	0.868687

awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_parus/fst/northerncardinal/out.weir.fst | awk '$3> 0.8 || NR==1' 

CHROM	POS	WEIR_AND_COCKERHAM_FST
# gene id 39089583
NC040875.1	5605	0.840426
NC040875.1	6247	1
NC040875.1	6556	0.840426
NC040875.1	6559	0.840426
# 39089565
NC040875.1	10692	0.840426
NC040875.1	11343	0.840426
NC040875.1	11359	0.840426



# this other paper said they searched 50kb on either side of every relevant snp, and this scaffold isn't that long
grep 'NC_040875.1' /scratch/dnjacks4/cardinalis/referencedata/ncbi_dataset_parus/data/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna.fai


grep 'NC_040875.1' /scratch/dnjacks4/cardinalis/referencedata/ncbi_dataset_parus/data/GCF_001522545.3/genomic.gff > relevant_items.txt

grep 'ID=gene' /scratch/dnjacks4/cardinalis/to_parus/fst/pyrrhuloxia/relevant_items.txt > relevant_genes.txt

awk '{print $9'} /scratch/dnjacks4/cardinalis/to_parus/fst/pyrrhuloxia/relevant_genes.txt  | awk 'BEGIN { FS = ";" } ; { print $2 }' | awk 'BEGIN { FS = ":" } ; { print $2 }' > genelist.txt

awk '$4 > 5000 || NR==1' relevant_genes.txt | awk '$4 < 13000 || NR==1'
 | awk '$3> 0.8 || NR==1' 


bgzip /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf
bcftools index /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf.gz 
bcftools view -r 'NC040875.1' /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf.gz > /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.NC040875.1.vcf
50,000


/scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.NC040875.1.vcf


mkdir /data5/sulidae/working

cd /data5/sulidae/working

echo -e "PYRR_rural" > pops.txt 
echo -e "NOCA_urban" >> pops.txt 
echo -e "NOCA_urban" >> pops.txt 
echo -e "NOCA_urban" >> pops.txt 
echo -e "NOCA_urban" >> pops.txt 
echo -e "NOCA_urban" >> pops.txt 
echo -e "NOCA_urban" >> pops.txt 
echo -e "PYRR_urban" >> pops.txt 
echo -e "PYRR_urban" >> pops.txt 
echo -e "PYRR_urban" >> pops.txt 
echo -e "PYRR_urban" >> pops.txt 
echo -e "PYRR_urban" >> pops.txt 
echo -e "PYRR_urban" >> pops.txt 
echo -e "NOCA_rural" >> pops.txt 
echo -e "NOCA_rural" >> pops.txt 
echo -e "NOCA_rural" >> pops.txt 
echo -e "NOCA_rural" >> pops.txt 
echo -e "PYRR_rural" >> pops.txt 
echo -e "PYRR_rural" >> pops.txt 
echo -e "PYRR_rural" >> pops.txt 
echo -e "PYRR_rural" >> pops.txt 
echo -e "PYRR_rural" >> pops.txt 
echo -e "NOCA_rural" >> pops.txt 
echo -e "NOCA_rural" >> pops.txt 











# running on chickadee because the R version on Sol can't install snprelate 
~/genomics/PCA_r.sh -v ~/parusfiltered.geno25.maf1.NC040875.1.vcf -o /data5/sulidae/working -p /data5/sulidae/working/pops.txt -n all -s y


module purge

module load picard/2.18.3
module load r/4.1.0

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








# RAxML Whole Genome

/data5/sulidae/working/wholegenome/parusfiltered.geno25.maf1.vcf
sed -i 's/UWBM/UW/g' /data5/sulidae/working/wholegenome/parusfiltered.geno25.maf1.vcf
~/vcf2phylip/vcf2phylip.py -i /data5/sulidae/working/wholegenome/parusfiltered.geno25.maf1.vcf

mkdir pruned 

cd pruned 

python ~/sula/filter_invariants_all.py /data5/sulidae/working/wholegenome/parusfiltered.geno25.maf1.min4.phy

mv variantsites.phy ../variantsites_mind2.phy 
mv variantsites_kept.txt ../variantsites_mind2_kept.txt 
cd .. 
rm -r pruned


echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-1205' > partitionfile.txt
raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n cardinalis_1 -q /data5/sulidae/working/wholegenome/raxml_wholegenome/partitionfile.txt -s /data5/sulidae/working/wholegenome/raxml_wholegenome/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V

# RAxML MT
~/parusfiltered.geno25.maf1.NC040875.1.vcf


sed -i 's/UWBM/UW/g' ~/parusfiltered.geno25.maf1.NC040875.1.vcf
~/vcf2phylip/vcf2phylip.py -i ~/parusfiltered.geno25.maf1.NC040875.1.vcf

mv ~/parusfiltered.geno25.maf1.NC040875.1.min4.phy .
mkdir pruned 

cd pruned 

python ~/sula/filter_invariants_all.py ../parusfiltered.geno25.maf1.NC040875.1.min4.phy

mv variantsites.phy ../variantsites_mind2.phy 
mv variantsites_kept.txt ../variantsites_mind2_kept.txt 
cd .. 
rm -r pruned


echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-75' > partitionfile.txt
raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n cardinalis_1 -q /data5/sulidae/working/partitionfile.txt -s /data5/sulidae/working/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V

vcftools --vcf /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt


cd /scratch/dnjacks4/cardinalis/to_parus/fst
vcftools --vcf /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.NC040875.1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt


/scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile.txt

# all all vs all all 
After filtering, kept 24 out of 24 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.074947
Weir and Cockerham weighted Fst estimate: 0.098933
After filtering, kept 14659 out of a possible 14659 Sites
Run Time = 1.00 seconds

# urban all v urban all
After filtering, kept 12 out of 24 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.062737
Weir and Cockerham weighted Fst estimate: 0.088914
After filtering, kept 14659 out of a possible 14659 Sites

# rural all v rural all 
Keeping individuals in 'keep' list
After filtering, kept 12 out of 24 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.069298
Weir and Cockerham weighted Fst estimate: 0.10738
After filtering, kept 14659 out of a possible 14659 Sites

# all mt v all mt 
After filtering, kept 24 out of 24 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.49523
Weir and Cockerham weighted Fst estimate: 0.63966
After filtering, kept 264 out of a possible 264 Sites

# urban mt v urban mt
After filtering, kept 12 out of 24 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.2494
Weir and Cockerham weighted Fst estimate: 0.33979
After filtering, kept 264 out of a possible 264 Sites

# rural mt v rural mt
After filtering, kept 12 out of 24 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.68622
Weir and Cockerham weighted Fst estimate: 0.86731
After filtering, kept 264 out of a possible 264 Sites



sed -i 's/NC//g' out.weir.fst
sed -i 's/NW//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "all_urban_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()






##### repeating on to_cardinalis edwards genome 
vcftools --vcf /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt


cp /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalis_snps_multiallelic.vcf /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalissnpsmultiallelic.vcf

sed -i 's/_//g' /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalissnpsmultiallelic.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalissnpsmultiallelic.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out cardinalisfiltered.geno25.maf1

## all 
# all v all 
vcftools --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalisfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr.txt

sed -i 's/JACDOX//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R 

library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "all_cardinalis_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()


# rural v rural
vcftools --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalisfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt

sed -i 's/JACDOX//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R 

library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "rural_cardinalis_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()


# urban v urban
vcftools --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalisfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt

sed -i 's/JACDOX//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R 

library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "urban_cardinalis_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()

# noca
vcftools --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalisfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_noca_rural.txt

sed -i 's/JACDOX//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R 

library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "noca_cardinalis_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()

# pyrr 
vcftools --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalisfiltered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_parus/reference_lists/fstpopfile_pyrr_rural.txt

sed -i 's/JACDOX//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R 

library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "pyrr_cardinalis_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()

## mt
# all v all 
# rural v rural
# urban v urban
# noca
# pyrr 


vcftools --vcf /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.NC040875.1.vcf --het --out mitochondrial.het
