module purge
module load r-4.0.2-gcc-11.2.0 #r-4.2.2-gcc-11.2.0
module load bayescan/2.1
module load htslib-1.16-gcc-11.2.0
module load bcftools-1.14-gcc-11.2.0
module load vcftools-0.1.14-gcc-11.2.0

# module load r-4.2.2-gcc-11.2.0
module load sqlite-3.38.5-gcc-11.2.0
module load proj-8.2.1-gcc-11.2.0
module load gdal-3.4.3-gcc-11.2.0
module load geos-3.9.1-gcc-11.2.0
module load  gcc-12.1.0-gcc-11.2.0
module load libxc-5.1.7-gcc-12.1.0
module load libvterm-0.1.4-gcc-11.2.0

bwa index /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna


cd /scratch/dnjacks4/cardinalis/to_parus

/scratch/dnjacks4/cardinalis/sort.sh -i /scratch/dnjacks4/cardinalis/sampleids.txt -p /scratch/dnjacks4/cardinalis/genomicdata/trimmed_fastas_select/ -r  /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna

/scratch/dnjacks4/cardinalis/to_parus/sort_postbam.sh -i /scratch/dnjacks4/cardinalis/sampleids.txt -p /scratch/dnjacks4/cardinalis/genomicdata/trimmed_fastas_select/ -r  /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna


ref="/scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna"
bamdir="/scratch/dnjacks4/cardinalis/to_b10k/sorted_bam_files/"
ID="b10k"
bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*_sorted_RGadded_dupmarked.bam | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf






bcftools stats -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_snps_multiallelic.vcf > b10k_snps_multiallelic.stats.txt

# bcftools stats /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf > parusfiltered.geno25.maf1.stats.txt

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/b10k_snps_multiallelic.vcf --allow-extra-chr --missing --cluster-missing --freq

# plink --vcf /scratch/dnjacks4/cardinalis/to_parus/parus_snps_multiallelic.vcf --allow-extra-chr --missing --cluster missing --within individualnames.txt --freq

# plink --vcf /scratch/dnjacks4/cardinalis/to_parus/parus_snps_multiallelic.vcf --allow-extra-chr --missing --cluster missing --within cluster_pop.txt --freq


#filters by quality
bcftools view -i 'QUAL>100' /scratch/dnjacks4/cardinalis/to_b10k/b10k_snps_multiallelic.vcf  > b10k_qualitysort.vcf

vcftools --vcf b10k_qualitysort.vcf --min-meanDP 2 --remove-indels --recode --out b10k_filtered

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out b10k_filtered.geno25.maf1

bcftools stats -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf > b10k_filtered.recode.stats.txt
bcftools stats -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > b10k_filtered.geno25.maf1.stats.txt

# with a max-meanDP 8:
# After filtering, kept 16120 out of a possible 3246716 Sites
# After filtering, kept 20825 out of a possible 3246716 Sites

sed -i 's/\_//g' b10k_filtered.recode.vcf
sed -i 's/\UWBM/BM/g' b10k_filtered.recode.vcf

plink --vcf b10k_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02  --recode vcf-iid --out b10k_filtered_geno02

plink --vcf b10k_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.2  --recode vcf-iid --out b10k_filtered_geno2

plink --vcf b10k_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.15  --recode vcf-iid --out b10k_filtered_geno15

plink --vcf b10k_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.1  --recode vcf-iid --out b10k_filtered_geno1


plink --vcf b10k_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid --out b10k_filtered_mind2

plink --vcf b10k_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.1 --maf 0.01 --recode vcf-iid --out b10k_filtered_cleaned


cp b10k_filtered_cleaned.vcf b10k_filtered_cleaned_zip.vcf
bgzip b10k_filtered_cleaned_zip.vcf

bcftools index b10k_filtered_cleaned_zip.vcf.gz

~/vcf2phylip/vcf2phylip.py -i /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf

python ~/sula/filter_invariants_all.py b10k_filtered.geno25.maf1.min4.phy

mkdir raxml
cd raxml
echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-2906' > partitionfile.txt
~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n cardinalis_1 -q /scratch/dnjacks4/cardinalis/to_b10k/raxml/partitionfile.txt -s/scratch/dnjacks4/cardinalis/to_b10k/variantsites.phy -T 6 -p 12345 -N 10 ­-b 12345 -V







## pca

### all
cd /scratch/dnjacks4/cardinalis/to_b10k/PCA/all

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

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf -o /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 

### pyrr urban vs rural
bcftools view -s MSB25201,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,BM103346,BM77548,BM77718,BM77780,BM77781 /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf --force-samples > /scratch/dnjacks4/cardinalis/to_b10k/pyrr_pca.vcf

cd /scratch/dnjacks4/cardinalis/to_b10k/PCA/pyrr
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

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf -o /scratch/dnjacks4/cardinalis/to_b10k/PCA/pyrr/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/pyrr/pops.txt -n pyrr 

### noca urban vs rural 

bcftools view -s NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,BM100619,BM100620,BM100621,BM103345,BM77856,BM77978 /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf --force-samples > /scratch/dnjacks4/cardinalis/to_b10k/noca_pca.vcf


cd /scratch/dnjacks4/cardinalis/to_b10k/PCA/noca

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


~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/noca_pca.vcf -o /scratch/dnjacks4/cardinalis/to_b10k/PCA/noca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/noca/pops.txt -n noca 
















# FST and tajima's D

echo 'MSB25201,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'NOCA003,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'NOCA004,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'NOCA006,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'NOCA008,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'NOCA012,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'NOCA013,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'PYRR003,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'PYRR004,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'PYRR006,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'PYRR007,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'PYRR009,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'PYRR011,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM100619,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM100620,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM100621,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM103345,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM103346,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM77548,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM77718,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM77780,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM77781,PYRR' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM77856,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
echo 'BM77978,NOCA' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt
sed -i 's/,/\t/g' /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile.txt




echo 'MSB25201' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt
echo 'BM103346' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt
echo 'BM77548' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt
echo 'BM77718' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt
echo 'BM77780' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt
echo 'BM77781' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt

echo 'PYRR003' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR004' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR006' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR007' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR009' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt
echo 'PYRR011' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt

echo 'NOCA003' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA004' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA006' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA008' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA012' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt
echo 'NOCA013' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt

echo 'BM100619' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt
echo 'BM100620' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt
echo 'BM100621' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt
echo 'BM103345' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt
echo 'BM77856' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt
echo 'BM77978' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt


echo 'MSB25201' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'BM103346' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'BM77548' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'BM77718' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'BM77780' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'BM77781' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt

echo 'PYRR003' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR004' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR006' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR007' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR009' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt
echo 'PYRR011' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt

echo 'NOCA003' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'NOCA004' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'NOCA006' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'NOCA008' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'NOCA012' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'NOCA013' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt

echo 'BM100619' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'BM100620' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'BM100621' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'BM103345' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'BM77856' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt
echo 'BM77978' >>  /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt




cd /scratch/dnjacks4/cardinalis/to_b10k/fst/pyrr
vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt

sed -i 's/VYXE//g' out.weir.fst
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

awk '$3 != "-nan" || NR==1' out.weir.fst | awk '$3 == 1 || NR==1' | wc -l

awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_b10k/fst/noca/out.weir.fst | awk '$3> 0.75 || NR==1' 


awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_b10k/fst/pyrr/out.weir.fst | awk '$3> 0.75 || NR==1' | awk '{print $1}' | sort -u

01005221.1
01006626.1
01008328.1
01008382.1
01020886.1
01021115.1
01021741.1

CHROM	POS	WEIR_AND_COCKERHAM_FST
01005221.1	2174717	0.78166
01005221.1	2174752	0.752871
01005221.1	2174761	0.78166
01006626.1	142206	0.771763
01006626.1	142332	0.840426
01006626.1	142337	0.840426
01006626.1	142340	0.840426
01006626.1	142346	0.840426
01008328.1	194860	0.777778
01008328.1	194864	0.777778
01008328.1	194882	0.777778
01008382.1	741075	0.840426
01020886.1	118926	0.777778
01020886.1	118927	0.777778
01020886.1	118929	0.777778
01020886.1	118935	0.777778
01020886.1	118945	0.777778
01020886.1	118948	0.777778
01020886.1	118950	0.777778
01020886.1	119022	0.777778
01020886.1	119043	0.777778
01020886.1	119046	0.777778
01020886.1	119055	0.777778
01020886.1	119058	0.840426
01020886.1	119067	0.840426
01020886.1	119094	0.777778
01021115.1	168762	0.777778
01021115.1	168765	0.777778
01021741.1	238957	0.771763

# unique chroms
01005221.1
01006626.1
01008328.1
01008382.1
01020886.1
01021115.1
01021741.1

grep '01005221.1\|01006626.1\|01008328.1\|01008382.1\|01020886.1\|01021115.1\|01021741.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff > relevant_items.txt



grep 'ID=gene' relevant_items.txt > relevant_genes.txt

awk '{print $9'} relevant_genes.txt  | awk 'BEGIN { FS = ";" } ; { print $2 }' | awk 'BEGIN { FS = "=" } ; { print $2 }' > genelist.txt






cd /scratch/dnjacks4/cardinalis/to_b10k/tajimad/pyrr/urban

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_pyrr.filtered.geno25.maf1.vcf --TajimaD 20000 


sed -i 's/VYXE//g' out.Tajima.D

awk '{sub(/\./,"",$1)}1' out.Tajima.D | column -t > out.Tajima.D.formanhattan


# Tajima's D
R

taj.all <- read.table("out.Tajima.D",header=T)
taj.subset<-taj.all[complete.cases(taj.all),]

SNP<-c(1: (nrow(taj.subset)))

lower = min(taj.subset$TajimaD)
upper = max(taj.subset$TajimaD)
lower_cutoff = lower + ((upper-lower)*0.025)
upper_cutoff = upper - ((upper-lower)*0.025)

MoreThanLower <- taj.subset$TajimaD > lower_cutoff
LessThanUpper <- taj.subset$TajimaD < upper_cutoff
significant <- MoreThanLower & LessThanUpper 

myBg <- !significant


mydf<-data.frame(SNP,myBg,taj.subset)

sigdf <-  mydf[which(mydf$myBg),]

write.table(sigdf, file = "sigtd.tsv")

pdf(file = "pyrr_urban_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(taj.subset$TajimaD,br=20)
dev.off()

pdf(file = "pyrr_urban_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(TajimaD ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()


cd /scratch/dnjacks4/cardinalis/to_b10k/tajimad/pyrr/rural

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_pyrr.filtered.geno25.maf1.vcf --TajimaD 20000 


sed -i 's/VYXE//g' out.Tajima.D

awk '{sub(/\./,"",$1)}1' out.Tajima.D | column -t > out.Tajima.D.formanhattan


# Tajima's D
R

taj.all <- read.table("out.Tajima.D",header=T)
taj.subset<-taj.all[complete.cases(taj.all),]

SNP<-c(1: (nrow(taj.subset)))

lower = min(taj.subset$TajimaD)
upper = max(taj.subset$TajimaD)
lower_cutoff = lower + ((upper-lower)*0.025)
upper_cutoff = upper - ((upper-lower)*0.025)

MoreThanLower <- taj.subset$TajimaD > lower_cutoff
LessThanUpper <- taj.subset$TajimaD < upper_cutoff
significant <- MoreThanLower & LessThanUpper 

myBg <- !significant



mydf<-data.frame(SNP,myBg,taj.subset)

pdf(file = "pyrr_rural_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(taj.subset$TajimaD,br=20)
dev.off()

pdf(file = "pyrr_rural_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(TajimaD ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()

cd /scratch/dnjacks4/cardinalis/to_b10k/tajimad/noca/urban

vcftools --gzvcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz --TajimaD 20000 


sed -i 's/VYXE//g' out.Tajima.D

awk '{sub(/\./,"",$1)}1' out.Tajima.D | column -t > out.Tajima.D.formanhattan


# Tajima's D
R

taj.all <- read.table("out.Tajima.D",header=T)
taj.subset<-taj.all[complete.cases(taj.all),]

SNP<-c(1: (nrow(taj.subset)))

lower = min(taj.subset$TajimaD)
upper = max(taj.subset$TajimaD)
lower_cutoff = lower + ((upper-lower)*0.025)
upper_cutoff = upper - ((upper-lower)*0.025)

MoreThanLower <- taj.subset$TajimaD > lower_cutoff
LessThanUpper <- taj.subset$TajimaD < upper_cutoff
significant <- MoreThanLower & LessThanUpper 

myBg <- !significant



mydf<-data.frame(SNP,myBg,taj.subset)

pdf(file = "noca_urban_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(taj.subset$TajimaD,br=20)
dev.off()

pdf(file = "noca_urban_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(TajimaD ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()



cd /scratch/dnjacks4/cardinalis/to_b10k/tajimad/noca/rural

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.filtered.geno25.maf1.vcf --TajimaD 20000 


sed -i 's/VYXE//g' out.Tajima.D

awk '{sub(/\./,"",$1)}1' out.Tajima.D | column -t > out.Tajima.D.formanhattan


# Tajima's D
R

taj.all <- read.table("out.Tajima.D",header=T)
taj.subset<-taj.all[complete.cases(taj.all),]

SNP<-c(1: (nrow(taj.subset)))

lower = min(taj.subset$TajimaD)
upper = max(taj.subset$TajimaD)
lower_cutoff = lower + ((upper-lower)*0.025)
upper_cutoff = upper - ((upper-lower)*0.025)

MoreThanLower <- taj.subset$TajimaD > lower_cutoff
LessThanUpper <- taj.subset$TajimaD < upper_cutoff
significant <- MoreThanLower & LessThanUpper 

myBg <- !significant



mydf<-data.frame(SNP,myBg,taj.subset)

pdf(file = "noca_rural_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(taj.subset$TajimaD,br=20)
dev.off()

pdf(file = "noca_rural_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(TajimaD ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()

cd /scratch/dnjacks4/cardinalis/to_b10k/nucleotidediversity/pyrr/urban
vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_pyrr.vcf --window-pi 20000


sed -i 's/VYXE//g' out.windowed.pi

awk '{sub(/\./,"",$1)}1' out.windowed.pi | column -t > out.windowed.pi.formanhattan

R

pi.all <- read.table("out.windowed.pi",header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$PI)
upper = max(pi.subset$PI)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$PI < cutoff

myBg <- !LessThanCutoff

mydf<-data.frame(SNP,myBg,pi.subset)

sigdf <-  mydf[which(mydf$myBg),]

write.table(sigdf, file = "sigpi.tsv")


pdf(file = "pyrr_urban_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(pi.subset$PI,br=20)
dev.off()

pdf(file = "pyrr_urban_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(PI ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()

cd /scratch/dnjacks4/cardinalis/to_b10k/nucleotidediversity/pyrr/rural

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_pyrr.vcf --window-pi 20000


sed -i 's/VYXE//g' out.windowed.pi

awk '{sub(/\./,"",$1)}1' out.windowed.pi | column -t > out.windowed.pi.formanhattan

R

pi.all <- read.table("out.windowed.pi",header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$PI)
upper = max(pi.subset$PI)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$PI < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,pi.subset)

pdf(file = "pyrr_rural_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(pi.subset$PI,br=20)
dev.off()

pdf(file = "pyrr_rural_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(PI ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()

cd /scratch/dnjacks4/cardinalis/to_b10k/nucleotidediversity/noca/urban


vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.vcf --window-pi 20000


sed -i 's/VYXE//g' out.windowed.pi

awk '{sub(/\./,"",$1)}1' out.windowed.pi | column -t > out.windowed.pi.formanhattan

R

pi.all <- read.table("out.windowed.pi",header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$PI)
upper = max(pi.subset$PI)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$PI < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,pi.subset)

pdf(file = "noca_urban_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(pi.subset$PI,br=20)
dev.off()

pdf(file = "noca_urban_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(PI ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()

cd /scratch/dnjacks4/cardinalis/to_b10k/nucleotidediversity/noca/rural

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.vcf --window-pi 20000


sed -i 's/VYXE//g' out.windowed.pi

awk '{sub(/\./,"",$1)}1' out.windowed.pi | column -t > out.windowed.pi.formanhattan

R

pi.all <- read.table("out.windowed.pi",header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$PI)
upper = max(pi.subset$PI)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$PI < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,pi.subset)

pdf(file = "noca_rural_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(pi.subset$PI,br=20)
dev.off()

pdf(file = "noca_rural_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(PI ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()























cd /scratch/dnjacks4/cardinalis/to_b10k/fst/noca
vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt

sed -i 's/VYXE//g' out.weir.fst
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



awk '$3 != "-nan" || NR==1' out.weir.fst | awk '$3 == 1 || NR==1' | wc -l

awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_b10k/fst/noca/out.weir.fst | awk '$3> 0.75 || NR==1' 

List of genes of interest, both species:
Plekhf2,Ch037,Ctnna3_0,Lrrtm3,Mrps5,Ift122,Rho,CARCAR_R15577,CARCAR_R00463,CARCAR_R00463,Tmcc1,Fxr1_0,Fxr1_1,Pcbp3_0,Col6a1,Chga,Itpk1

# noca 
List of genes of interest:
Plekhf2,Ch037,Ctnna3_0,Lrrtm3,Mrps5

Scaffold: 01005221.1 
Total SNPs: 3
True range: 2174717-2174761
Search range: 2124717-2224761

grep 'VYXE01005221.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: Plekhf2,Ch037

Scaffold: 01008768.1
Total SNPs: 11
True range: 38618-38684
Search range: 0 - 88684

grep 'VYXE01008768.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: None

Scaffold: 01009024.1 
Total SNPs: 1
True range: 321029
Search range: 271029-371029

grep 'VYXE01009024.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: Ctnna3_0, Lrrtm3

Scaffold: 01011506.1  
Total SNPs: 2
True range: 97212-97217
Search range: 47212-147217

grep 'VYXE01011506.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: Mrps5


Scaffold: 01018132.1
Total SNPs: 2
True range: 45204-45210
Search range: 0-95210

grep 'VYXE01018132.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: None (one unidentified)


Scaffold: 01019412.1
Total SNPs: 1
True range: 27631
Search range: 0-77631

grep 'VYXE01019412.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: None (one unidentified)

# pyrr
awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_b10k/fst/pyrr/out.weir.fst | awk '$3> 0.75 || NR==1' 

List of genes of interest:
Plekhf2, Ch037,Ift122,Rho,CARCAR_R15577,CARCAR_R00463,CARCAR_R00463,Tmcc1,Fxr1_0,Fxr1_1,Pcbp3_0,Col6a1,Chga,Itpk1

Scaffold: 01005221.1 
Total SNPs: 3
True range: 2174717-2174761
Search range: 2124717-2224761

grep 'VYXE01005221.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'


Genes in range: Plekhf2, Ch037

Scaffold: 01006626.1 
Total SNPs: 5
True range: 142206-142346
Search range: 92206-192346

grep 'VYXE01006626.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: Ift122,Rho,CARCAR_R15577,CARCAR_R00463


Scaffold: 01008328.1 
Total SNPs: 3
True range: 194860-194882
Search range: 144860-244882

grep 'VYXE01006626.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'


Genes in range: CARCAR_R00463,Tmcc1


Scaffold: 01008382.1
Total SNPs: 1
True range: 741075
Search range: 691075-791075

grep 'VYXE01008382.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'


Genes in range: Fxr1_0,Fxr1_1


Scaffold: 01020886.1 
Total SNPs: 14
True range: 118926-119094
Search range: 68926-169094

grep 'VYXE01020886.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: Pcbp3_0,Col6a1


Scaffold: 01021115.1 
Total SNPs: 3
True range: 168762-238957
Search range: 118762-288957

grep 'VYXE01021115.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene'

Genes in range: Chga,Itpk1

awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_b10k/fst/noca/out.weir.fst | awk '$3> 0.75 || NR==1' | awk '{print $1}' | sort -u

awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_b10k/fst/pyrr/out.weir.fst | awk '$3> 0.75 || NR==1' | awk '{print $1}' | sort -u

awk '$3 != "-nan" || NR==1' /scratch/dnjacks4/cardinalis/to_b10k/fst/noca/out.weir.fst | awk '$3> 0.75 || NR==1' | awk '{print $1}' | sort -u

Scaffold VYXE01005221.1 has high fst in both species

grep 'VYXE01005221.1' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff


01005221.1      2174717 0.78166
01005221.1      2174752 0.752871
01005221.1      2174761 0.78166
01005221.1      2174717 0.78166
01005221.1      2174752 0.752871
01005221.1      2174761 0.78166

(range of snps with buffer of 50k on either side)
2124717 - 2224761

VYXE01005221.1     Genbank gene    2141860 2142606 .       +       .       ID=gene-CARCAR_R01774;Name=Plekhf2;end_range=2142606,.;gbkey=Gene;gene=Plekhf2;gene_biotype=protein_coding;locus_tag=CARCAR_R01774;partial=true;start_range=.,2141860

VYXE01005221.1     Genbank gene    2176060 2178234 .       -       .       ID=gene-CARCAR_R12099;Name=Ch037;end_range=2178234,.;gbkey=Gene;gene=Ch037;gene_biotype=protein_coding;locus_tag=CARCAR_R12099;partial=true;start_range=.,2176060


cd /scratch/dnjacks4/cardinalis/to_b10k/fst/all
vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr.txt

sed -i 's/VYXE//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "all_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()




mkdir urban_urban

cd /scratch/dnjacks4/cardinalis/to_b10k/fst/urban_urban
vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_urban.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_urban.txt

sed -i 's/VYXE//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "urban_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()



cd /scratch/dnjacks4/cardinalis/to_b10k/fst/rural_rural
vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_noca_rural.txt --weir-fst-pop /scratch/dnjacks4/cardinalis/to_b10k/reference_lists/fstpopfile_pyrr_rural.txt

sed -i 's/VYXE//g' out.weir.fst
awk '{sub(/\./,"",$1)}1' out.weir.fst | column -t > out.weir.formanhattan.fst

R
library(qqman)
fst<-read.table("out.weir.formanhattan.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "rural_fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()







# bayescan
interactive -n 12 -t 3-00:00 #bayescan
module purge
module load r-4.2.2-gcc-11.2.0
module load bayescan/2.1
module load htslib-1.16-gcc-11.2.0
module load bcftools-1.14-gcc-11.2.0


cd /scratch/dnjacks4/cardinalis/to_b10k/bayescan

cp /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf.gz /scratch/dnjacks4/cardinalis/to_b10k/b10kfiltered.recode.vcf.gz

gunzip /scratch/dnjacks4/cardinalis/to_b10k/b10kfiltered.recode.vcf.gz

sed -i 's/_//g' /scratch/dnjacks4/cardinalis/to_b10k/b10kfiltered.recode.vcf

# pyrrhuloxia 
bcftools view -s 'MSB25201,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,BM103346,BM77548,BM77718,BM77780,BM77781' /scratch/dnjacks4/cardinalis/to_b10k/b10kfiltered.recode.vcf > /scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr_all.vcf


plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr_all.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out pyrr.filtered.geno25.maf1


echo "INDIVIDUALS,STRATA" > /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt

echo 'MSB25201,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR003,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR004,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR006,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR007,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR009,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'PYRR011,Urban' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'BM103346,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'BM77548,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'BM77718,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'BM77780,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt
echo 'BM77781,Rural' >>  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt


sed -i 's/,/\t/g' /scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt



/scratch/dnjacks4/cardinalis/to_parus/bayescan/pyrr_bayes_test.vcf

R

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

vcf <- read.vcfR("/scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr.filtered.geno25.maf1.vcf")

pop_map <- read.table("/scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt", header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn="pyrr.filtered.bsc")

mkdir filtered 

bayescan -n 5000 -burn 50000 -pr_odds 10 /scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr/pyrr.filtered.bsc -od /scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr/filtered/




grep -v "#" /scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr.filtered.geno25.maf1.vcf | awk '{ print $1, $2 }' | tr ' ' '_' > /scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr.filtered.geno25.maf1.loci.txt

R 

bayescan=read.table("pyrr.filtered_fst.txt") 


SNPb=read.table("/scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr.filtered.geno25.maf1.loci.txt",header=FALSE)

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
bcftools view -s 'NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,BM100619,BM100620,BM100621,BM103345,BM77856,BM77978'   /scratch/dnjacks4/cardinalis/to_b10k/b10kfiltered.recode.vcf > /scratch/dnjacks4/cardinalis/to_b10k/bayescan/noca_all.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/bayescan/noca_all.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out noca.filtered.geno25.maf1


R

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

vcf <- read.vcfR("/scratch/dnjacks4/cardinalis/to_b10k/bayescan/noca.filtered.geno25.maf1.vcf")

pop_map <- read.table("/scratch/dnjacks4/cardinalis/to_parus/reference_lists/noca_bayescan_popfile.txt", header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn="noca.filtered.bsc")

cd noca
mkdir filtered 

bayescan -n 5000 -burn 50000 -pr_odds 10 /scratch/dnjacks4/cardinalis/to_b10k/bayescan/noca/noca.filtered.bsc -od /scratch/dnjacks4/cardinalis/to_b10k/bayescan/noca/filtered/





grep -v "#" /scratch/dnjacks4/cardinalis/to_b10k/bayescan/noca.filtered.geno25.maf1.vcf | awk '{ print $1, $2 }' | tr ' ' '_' > /scratch/dnjacks4/cardinalis/to_b10k/bayescan/noca.filtered.geno25.maf1.loci.txt

R 

bayescan=read.table("noca.filtered_fst.txt") 


SNPb=read.table("/scratch/dnjacks4/cardinalis/to_b10k/bayescan/noca.filtered.geno25.maf1.loci.txt",header=FALSE)

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






































# pcadapt
module purge
module load r-4.2.2-gcc-11.2.0
module load bayescan/2.1
module load htslib-1.16-gcc-11.2.0
module load bcftools-1.14-gcc-11.2.0


# PCAdapt

# pyrrhuloxia

bcftools view -s 'MSB25201,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,BM103346,BM77548,BM77718,BM77780,BM77781' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > /scratch/dnjacks4/cardinalis/to_b10k/pcadapt/pyrr_all.vcf


plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/pcadapt/pyrr_all.vcf --allow-extra-chr --out pyrr.filtered --make-bed

cd /scratch/dnjacks4/cardinalis/to_b10k/pcadapt/pyrr

R 

library(pcadapt)
library(qvalue)

path_to_file <- "/scratch/dnjacks4/cardinalis/to_b10k/pcadapt/pyrr.filtered.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
x <- pcadapt(input = filename, K = 4)

pdf(file = "scree.pdf", width = 10, height = 10, useDingbats=FALSE)
    plot(x, option = "screeplot")
    dev.off()

poplist.names <- c(rep("Rural", 1), rep("Urban", 6), rep("Rural", 5))

pdf(file = "pca.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "scores", pop = poplist.names)
    dev.off()
    
pdf(file = "hist_pvalues.pdf", width = 10, height = 10, useDingbats=FALSE)
  hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
    dev.off()



# E.1. q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)
write.csv(outliers, file = "qvalue_outliers.csv")

# E.2. Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
write.csv(outliers, file = "BenjaminiHochberg_outliers.csv")


# E.3. Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
write.csv(outliers, file = "Bonferroni_outliers.csv")


plot(x, option = "screeplot", K = 2)

# With integers

# create a scree plot to determine the actual number of K)

pdf(file = "scree.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "screeplot", pop = poplist.int)
    dev.off()

pdf(file = "pc12.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "scores", pop = poplist.int)
    dev.off()

pdf(file = "pc34.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "scores", i = 3, j = 4, pop = poplist.int)
    dev.off()

pdf(file = "manhattan.pdf", width = 40, height = 10, useDingbats=FALSE)
  plot(x, option = "manhattan")
    dev.off()

pdf(file = "qqplot.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "qqplot")
    dev.off()

pdf(file = "hist.pdf", width = 20, height = 10, useDingbats=FALSE)
  hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "blue")
    dev.off()

sumx <- summary(x)
cat(paste(sumx),file=paste0("pcadapt_summary.txt"),sep="\n",append=TRUE)

#q-values
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers_q <- which(qval < alpha)
lo_q <- length(outliers_q)


cat(paste("Length of q-value outliers= ", lo_q),file=paste0("pcadapt_stats.txt"),sep="\n",append=TRUE)
cat(paste(outliers_q),file=paste0("outliers_q.txt"),sep="\n",append=TRUE)


#Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
outliers_bh <- which(padj < alpha)
lo_bh <- length(outliers_bh)

cat(paste("Length of Benjamini-Hochberg outliers= ", lo_bh),file=paste0("pcadapt_stats.txt"),sep="\n",append=TRUE)
cat(paste(outliers_bh),file=paste0("outliers_bh.txt"),sep="\n",append=TRUE)

#Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
outliers_bf <- which(padj < alpha)
lo_bf <- length(outliers_bf)

cat(paste("Length of Bonferroni outliers= ", lo_bf),file=paste0("pcadapt_stats.txt"),sep="\n",append=TRUE)
cat(paste(outliers_bf),file=paste0("outliers_bf.txt"),sep="\n",append=TRUE)




























# northern cardinals

bcftools view -s 'NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,UWBM100619,UWBM100620,UWBM100621,UWBM103345,UWBM77856,UWBM77978'  /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > /scratch/dnjacks4/cardinalis/to_b10k/pcadapt/noca_all.vcf


plink --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1.vcf --allow-extra-chr --out noca.filtered --make-bed

cd /scratch/dnjacks4/cardinalis/to_b10k/pcadapt/noca

R 

library(pcadapt)
library(qvalue)

path_to_file <- "/scratch/dnjacks4/cardinalis/to_b10k/pcadapt/noca.filtered.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
x <- pcadapt(input = filename, K = 10)

pdf(file = "scree.pdf", width = 10, height = 10, useDingbats=FALSE)
    plot(x, option = "screeplot")
    dev.off()

poplist.names <- c(rep("Urban", 6), rep("Rural", 6))
pdf(file = "pca.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "scores") #, pop = poplist.names)
    dev.off()
    
pdf(file = "hist_pvalues.pdf", width = 10, height = 10, useDingbats=FALSE)
  hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
    dev.off()



# E.1. q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)
write.csv(outliers, file = "qvalue_outliers.csv")

# E.2. Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
write.csv(outliers, file = "BenjaminiHochberg_outliers.csv")


# E.3. Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
write.csv(outliers, file = "Bonferroni_outliers.csv")


plot(x, option = "screeplot", K = 2)

# With integers

# create a scree plot to determine the actual number of K)

pdf(file = "scree.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "screeplot", pop = poplist.int)
    dev.off()

pdf(file = "pc12.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "scores", pop = poplist.int)
    dev.off()

pdf(file = "pc34.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "scores", i = 3, j = 4, pop = poplist.int)
    dev.off()

pdf(file = "manhattan.pdf", width = 40, height = 10, useDingbats=FALSE)
  plot(x, option = "manhattan")
    dev.off()

pdf(file = "qqplot.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "qqplot")
    dev.off()

pdf(file = "hist.pdf", width = 20, height = 10, useDingbats=FALSE)
  hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "blue")
    dev.off()

sumx <- summary(x)
cat(paste(sumx),file=paste0("pcadapt_summary.txt"),sep="\n",append=TRUE)

#q-values
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers_q <- which(qval < alpha)
lo_q <- length(outliers_q)


cat(paste("Length of q-value outliers= ", lo_q),file=paste0("pcadapt_stats.txt"),sep="\n",append=TRUE)
cat(paste(outliers_q),file=paste0("outliers_q.txt"),sep="\n",append=TRUE)


#Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
outliers_bh <- which(padj < alpha)
lo_bh <- length(outliers_bh)

cat(paste("Length of Benjamini-Hochberg outliers= ", lo_bh),file=paste0("pcadapt_stats.txt"),sep="\n",append=TRUE)
cat(paste(outliers_bh),file=paste0("outliers_bh.txt"),sep="\n",append=TRUE)

#Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
outliers_bf <- which(padj < alpha)
lo_bf <- length(outliers_bf)

cat(paste("Length of Bonferroni outliers= ", lo_bf),file=paste0("pcadapt_stats.txt"),sep="\n",append=TRUE)
cat(paste(outliers_bf),file=paste0("outliers_bf.txt"),sep="\n",append=TRUE)





















# SweeD 
# Tests for selective sweeps

bcftools view -s 'PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011' /scratch/dnjacks4/cardinalis/to_b10k/b10ksnpsmultiallelic.vcf > /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_pyrr.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_pyrr.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid  --out urban_pyrr.filtered.geno25.maf1

bcftools view -s 'MSB25201,UWBM103346,UWBM77548,UWBM77718,UWBM77780,UWBM77781' /scratch/dnjacks4/cardinalis/to_b10k/b10ksnpsmultiallelic.vcf > /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_pyrr.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_pyrr.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid  --out rural_pyrr.filtered.geno25.maf1

 SweeD -name ms1 -grid 1000 -length 100000 -input ms1.out

~/programs/sweed/SweeD -name pyrr_urban -input /scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/urban_pyrr.filtered.geno25.maf1.vcf -grid 100 -length 100000

~/programs/omegaplus/OmegaPlus -input /data5/sulidae/working/urban_noca.filtered.geno25.maf1.vcf -name urban_noca -grid 1000 -length 100000 -minwin 100 -maxwin 10000 -seed  9382


bgzip /scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/urban_pyrr.filtered.geno25.maf1.vcf
tabix /scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/urban_pyrr.filtered.geno25.maf1.vcf.gz
zcat /scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/urban_pyrr.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names.txt



grep 'Position' SweeD_Report.pyrr_urban | head -n 1 > SweeD_Report.header
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.pyrr_urban_manhattan


while IFS=$'\t' read -r col1 col2 col3 col4 col5
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 $col3 $col4 $col5" >> SweeD_Report.pyrr_urban_manhattan
    fi
done < SweeD_Report.pyrr_urban

grep 'Position' SweeD_Report.pyrr_urban | head -n 1 > SweeD_Report.pyrr_urban_manhattan_scaffnames
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.pyrr_urban_manhattan

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> SweeD_Report.pyrr_urban_manhattan_scaffnames
  fi
done < SweeD_Report.pyrr_urban_manhattan

sed -i 's/VYXE//g' SweeD_Report.pyrr_urban_manhattan_scaffnames
# sed -i 's/\(.\{1\}\).1/\1/' SweeD_Report.pyrr_urban_manhattan_scaffnames
awk '{sub(/\./,"",$1)}1' SweeD_Report.pyrr_urban_manhattan_scaffnames | column -t > SweeD_Report.pyrr_urban_manhattan_scaffnames.formanhattan

mv SweeD_Report.pyrr_urban_manhattan SweeD_Report.pyrr_urban_manhattan_scaffnames.formanhattan

R
library(qqman)
sweed<-read.table("SweeD_Report.pyrr_urban_manhattan_scaffnames.formanhattan", header=TRUE)
sweedsubset<-sweed[complete.cases(sweed),]
CLR<-c(1: (nrow(sweedsubset)))
mydf<-data.frame(CLR,sweedsubset)

pdf(file = "pyrr_urban_sweed.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()


### rural pyrrhuloxia 
~/programs/sweed/SweeD -name pyrr_rural -input /scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/rural_pyrr.filtered.geno25.maf1.vcf -grid 100


bgzip /scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/rural_pyrr.filtered.geno25.maf1.vcf
tabix /scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/rural_pyrr.filtered.geno25.maf1.vcf.gz
zcat /scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/rural_pyrr.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names.txt



grep 'Position' SweeD_Report.pyrr_rural | head -n 1 > SweeD_Report.header
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.pyrr_rural_manhattan


while IFS=$'\t' read -r col1 col2 col3 col4 col5
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 $col3 $col4 $col5" >> SweeD_Report.pyrr_rural_manhattan
    fi
done < SweeD_Report.pyrr_rural

grep 'Position' SweeD_Report.pyrr_rural | head -n 1 > SweeD_Report.noca_pyrr_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> SweeD_Report.pyrr_rural_manhattan_scaffnames
  fi
done < SweeD_Report.pyrr_rural_manhattan

sed -i 's/VYXE//g' SweeD_Report.pyrr_rural_manhattan_scaffnames
# sed -i 's/\(.\{1\}\).1/\1/' SweeD_Report.pyrr_rural_manhattan_scaffnames
awk '{sub(/\./,"",$1)}1' SweeD_Report.pyrr_rural_manhattan_scaffnames | column -t > SweeD_Report.pyrr_rural_manhattan_scaffnames.formanhattan

R
library(qqman)
sweed<-read.table("SweeD_Report.pyrr_rural_manhattan_scaffnames", header=TRUE)
sweedsubset<-sweed[complete.cases(sweed),]
CLR<-c(1: (nrow(sweedsubset)))
mydf<-data.frame(CLR,sweedsubset)

pdf(file = "pyrr_rural_sweed.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()



















### northern cardinals 
bcftools view -s 'NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013' /scratch/dnjacks4/cardinalis/to_b10k/b10ksnpsmultiallelic.vcf > /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid  --out urban_noca.filtered.geno25.maf1

bcftools view -s 'UWBM100619,UWBM100620,UWBM100621,UWBM103345,UWBM77856,UWBM77978'   /scratch/dnjacks4/cardinalis/to_b10k/b10ksnpsmultiallelic.vcf > /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid  --out rural_noca.filtered.geno25.maf1

cd northerncardinals

~/programs/sweed/SweeD -name noca_urban -input /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf -grid 100

bgzip /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf
tabix /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz
zcat /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names.txt

cat -n scaff_names.txt > scaff_names_numbers.txt



grep 'Position' SweeD_Report.noca_urban | head -n 1 > SweeD_Report.header
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.noca_urban_manhattan


while IFS=$'\t' read -r col1 col2 col3 col4 col5
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 $col3 $col4 $col5" >> SweeD_Report.noca_urban_manhattan
    fi
done < SweeD_Report.noca_urban

grep 'Position' SweeD_Report.noca_urban | head -n 1 > SweeD_Report.noca_urban_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> SweeD_Report.noca_urban_manhattan_scaffnames
  fi
done < SweeD_Report.noca_urban_manhattan

sed -i 's/VYXE//g' SweeD_Report.noca_urban_manhattan_scaffnames
sed -i 's/\(.\{1\}\).1/\1/' SweeD_Report.noca_urban_manhattan_scaffnames
awk '{sub(/\./,"",$1)}1' SweeD_Report.noca_urban_manhattan_scaffnames | column -t > SweeD_Report.noca_urban_manhattan_scaffnames.formanhattan

R
library(qqman)
sweed<-read.table("SweeD_Report.noca_urban_manhattan_scaffnames", header=TRUE)
sweedsubset<-sweed[complete.cases(sweed),]
CLR<-c(1: (nrow(sweedsubset)))
mydf<-data.frame(CLR,sweedsubset)

pdf(file = "noca_sweed.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()


#### noca rural

~/programs/sweed/SweeD -name noca_rural -input /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.filtered.geno25.maf1.vcf -grid 100

bgzip /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.filtered.geno25.maf1.vcf
tabix /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.filtered.geno25.maf1.vcf.gz
zcat /scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names.txt

rm scaff_names_numbers.txt



grep 'Position' SweeD_Report.noca_rural | head -n 1 > SweeD_Report.header
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.noca_rural_manhattan


while IFS=$'\t' read -r col1 col2 col3 col4 col5
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 $col3 $col4 $col5" >> SweeD_Report.noca_rural_manhattan
    fi
done < SweeD_Report.noca_rural

grep 'Position' SweeD_Report.noca_rural | head -n 1 > SweeD_Report.noca_rural_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> SweeD_Report.noca_rural_manhattan_scaffnames
  fi
done < SweeD_Report.noca_rural_manhattan

sed -i 's/VYXE//g' SweeD_Report.noca_rural_manhattan_scaffnames
# sed -i 's/\(.\{1\}\).1/\1/' SweeD_Report.noca_rural_manhattan_scaffnames
awk '{sub(/\./,"",$1)}1' SweeD_Report.noca_rural_manhattan_scaffnames | column -t > SweeD_Report.noca_rural_manhattan_scaffnames.formanhattan

R
library(qqman)
sweed<-read.table("SweeD_Report.noca_rural_manhattan_scaffnames", header=TRUE)
sweedsubset<-sweed[complete.cases(sweed),]
CLR<-c(1: (nrow(sweedsubset)))
mydf<-data.frame(CLR,sweedsubset)

pdf(file = "noca_rural_sweed.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()


module purge
module load r-4.2.2-gcc-11.2.0
module load sqlite-3.38.5-gcc-11.2.0
module load proj-8.2.1-gcc-11.2.0
module load gdal-3.4.3-gcc-11.2.0
module load geos-3.9.1-gcc-11.2.0
module load  gcc-12.1.0-gcc-11.2.0
module load libxc-5.1.7-gcc-12.1.0
module load libvterm-0.1.4-gcc-11.2.0

# couldn't install on Sol, running on chickadee
scp dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz .

scp dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/cardinalis/to_b10k/sweed/rural_noca.filtered.geno25.maf1.vcf.gz .

scp dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/urban_pyrr.filtered.geno25.maf1.vcf.gz .

scp dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/cardinalis/to_b10k/sweed/pyrrhuloxia/rural_pyrr.filtered.geno25.maf1.vcf.gz .

scp urban_noca.filtered.geno25.maf1.vcf.gz daja5529@chickadee.colorado.edu:/data5/sulidae/working/

scp rural_noca.filtered.geno25.maf1.vcf.gz daja5529@chickadee.colorado.edu:/data5/sulidae/working/

scp urban_pyrr.filtered.geno25.maf1.vcf.gz daja5529@chickadee.colorado.edu:/data5/sulidae/working/

scp rural_pyrr.filtered.geno25.maf1.vcf.gz daja5529@chickadee.colorado.edu:/data5/sulidae/working/

# Northern Cardinals
# Urban
~/programs/omegaplus/OmegaPlus -input /data5/sulidae/working/urban_noca.filtered.geno25.maf1.vcf -name urban_noca -grid 1000 -length 100000 -minwin 100 -maxwin 10000 -seed  9382


echo -e 'Scaffold\tPosition\tLikelihood' > OmegaPlus_Report.urban_noca_manhattan
while IFS=$'\t' read -r col1 col2 
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 " >> OmegaPlus_Report.urban_noca_manhattan
    fi
done < OmegaPlus_Report.urban_noca

bgzip ../urban_noca.filtered.geno25.maf1.vcf
tabix ../urban_noca.filtered.geno25.maf1.vcf.gz
zcat ../urban_noca.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names_noca_urban.txt

echo -e 'Scaffold\tPosition\tLikelihood' > OmegaPlus_Report.noca_urban_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names_noca_urban.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> OmegaPlus_Report.noca_urban_manhattan_scaffnames
  fi
done < OmegaPlus_Report.urban_noca_manhattan

sed -i 's/VYXE//g' OmegaPlus_Report.noca_urban_manhattan_scaffnames

R
library(qqman)
omega<-read.table("OmegaPlus_Report.noca_urban_manhattan_scaffnames", header=TRUE)
omega.subset<-omega[complete.cases(omega),]
SNP<-c(1: (nrow(omega.subset)))

lower = min(omega.subset$Likelihood)
upper = max(omega.subset$Likelihood)
cutoff = upper - ((upper-lower)*0.05)
LessThanCutoff <- omega.subset$Likelihood < cutoff
myBg <- !LessThanCutoff
mydf<-data.frame(SNP,myBg,omega.subset)
sigdf <-  mydf[which(mydf$myBg),]
write.table(sigdf, file = "sigomega.tsv")

mydf<-data.frame(SNP,omega.subset)

pdf(file = "noca_urban_omega.pdf", width = 20, height = 7, useDingbats=FALSE)
    print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()


# Rural

~/programs/omegaplus/OmegaPlus -input /data5/sulidae/working/rural_noca.filtered.geno25.maf1.vcf -name rural_noca -grid 1000 -length 100000 -minwin 100 -maxwin 10000 -seed  9382


echo -e 'Scaffold\tPosition\tLikelihood' > OmegaPlus_Report.rural_noca_manhattan
while IFS=$'\t' read -r col1 col2 
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 " >> OmegaPlus_Report.rural_noca_manhattan
    fi
done < OmegaPlus_Report.rural_noca

bgzip  /data5/sulidae/working/rural_noca.filtered.geno25.maf1.vcf
tabix /data5/sulidae/working/rural_noca.filtered.geno25.maf1.vcf.gz
zcat /data5/sulidae/working/rural_noca.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names_noca_rural.txt

echo -e 'Scaffold\tPosition\tLikelihood' > OmegaPlus_Report.noca_rural_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names_noca_rural.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> OmegaPlus_Report.noca_rural_manhattan_scaffnames
  fi
done < OmegaPlus_Report.rural_noca_manhattan

sed -i 's/VYXE//g' OmegaPlus_Report.noca_rural_manhattan_scaffnames

R
library(qqman)
omega<-read.table("OmegaPlus_Report.noca_rural_manhattan_scaffnames", header=TRUE)
omega.subset<-omega[complete.cases(omega),]
SNP<-c(1: (nrow(omega.subset)))

lower = min(omega.subset$Likelihood)
upper = max(omega.subset$Likelihood)
cutoff = upper - ((upper-lower)*0.05)
LessThanCutoff <- omega.subset$Likelihood < cutoff
myBg <- !LessThanCutoff
mydf<-data.frame(SNP,myBg,omega.subset)
sigdf <-  mydf[which(mydf$myBg),]
write.table(sigdf, file = "sigomega.tsv")

mydf<-data.frame(SNP,omega.subset)

pdf(file = "noca_rural_omega.pdf", width = 20, height = 7, useDingbats=FALSE)
    print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()




# Pyrrhuloxia
# Urban
~/programs/omegaplus/OmegaPlus -input /data5/sulidae/working/urban_pyrr.filtered.geno25.maf1.vcf -name urban_pyrr -grid 1000 -length 100000 -minwin 100 -maxwin 10000 -seed  9382


echo -e 'Scaffold\tPosition\tLikelihood' > OmegaPlus_Report.urban_pyrr_manhattan
while IFS=$'\t' read -r col1 col2 
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 " >> OmegaPlus_Report.urban_pyrr_manhattan
    fi
done < OmegaPlus_Report.urban_pyrr

bgzip /data5/sulidae/working/urban_pyrr.filtered.geno25.maf1.vcf
tabix /data5/sulidae/working/urban_pyrr.filtered.geno25.maf1.vcf.gz
zcat /data5/sulidae/working/urban_pyrr.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names_pyrr_urban.txt

echo -e 'Scaffold\tPosition\tLikelihood' > OmegaPlus_Report.pyrr_urban_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names_pyrr_urban.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> OmegaPlus_Report.pyrr_urban_manhattan_scaffnames
  fi
done < OmegaPlus_Report.urban_pyrr_manhattan

sed -i 's/VYXE//g' OmegaPlus_Report.pyrr_urban_manhattan_scaffnames

R
library(qqman)
omega<-read.table("OmegaPlus_Report.pyrr_urban_manhattan_scaffnames", header=TRUE)
omega.subset<-omega[complete.cases(omega),]
SNP<-c(1: (nrow(omega.subset)))

lower = min(omega.subset$Likelihood)
upper = max(omega.subset$Likelihood)
cutoff = upper - ((upper-lower)*0.05)
LessThanCutoff <- omega.subset$Likelihood < cutoff
myBg <- !LessThanCutoff
mydf<-data.frame(SNP,myBg,omega.subset)
sigdf <-  mydf[which(mydf$myBg),]
write.table(sigdf, file = "sigomega.tsv")

mydf<-data.frame(SNP,omega.subset)

pdf(file = "pyrr_urban_omega.pdf", width = 20, height = 7, useDingbats=FALSE)
    print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()


# Rural

~/programs/omegaplus/OmegaPlus -input /data5/sulidae/working/rural_pyrr.filtered.geno25.maf1.vcf -name rural_pyrr -grid 1000 -length 100000 -minwin 100 -maxwin 10000 -seed  9382


echo -e 'Scaffold\tPosition\tLikelihood' > OmegaPlus_Report.rural_pyrr_manhattan
while IFS=$'\t' read -r col1 col2 
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 " >> OmegaPlus_Report.rural_pyrr_manhattan
    fi
done < OmegaPlus_Report.rural_pyrr

bgzip /data5/sulidae/working/rural_pyrr.filtered.geno25.maf1.vcf
tabix /data5/sulidae/working/rural_pyrr.filtered.geno25.maf1.vcf.gz
zcat /data5/sulidae/working/rural_pyrr.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names_pyrr_rural.txt

echo -e 'Scaffold\tPosition\tLikelihood' > OmegaPlus_Report.pyrr_rural_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names_pyrr_rural.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> OmegaPlus_Report.pyrr_rural_manhattan_scaffnames
  fi
done < OmegaPlus_Report.rural_pyrr_manhattan

sed -i 's/VYXE//g' OmegaPlus_Report.pyrr_rural_manhattan_scaffnames

R
library(qqman)
omega<-read.table("OmegaPlus_Report.pyrr_rural_manhattan_scaffnames", header=TRUE)
omega.subset<-omega[complete.cases(omega),]
SNP<-c(1: (nrow(omega.subset)))

lower = min(omega.subset$Likelihood)
upper = max(omega.subset$Likelihood)
cutoff = upper - ((upper-lower)*0.05)
LessThanCutoff <- omega.subset$Likelihood < cutoff
myBg <- !LessThanCutoff
mydf<-data.frame(SNP,myBg,omega.subset)
sigdf <-  mydf[which(mydf$myBg),]
write.table(sigdf, file = "sigomega.tsv")

mydf<-data.frame(SNP,omega.subset)

pdf(file = "pyrr_rural_omega.pdf", width = 20, height = 7, useDingbats=FALSE)
    print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()
















# hapflk
# agave

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

plink --vcf ../b10k_filtered.geno25.maf1.vcf --allow-extra-chr --out b10k_filtered --recode


awk '{sub(/MSB25201/,"PYRR_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/NOCA003/,"NOCA_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/NOCA004/,"NOCA_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/NOCA006/,"NOCA_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/NOCA008/,"NOCA_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/NOCA012/,"NOCA_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/NOCA013/,"NOCA_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/PYRR003/,"PYRR_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/PYRR004/,"PYRR_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/PYRR006/,"PYRR_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/PYRR007/,"PYRR_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/PYRR009/,"PYRR_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/PYRR011/,"PYRR_urban",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM100619/,"NOCA_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM100620/,"NOCA_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM100621/,"NOCA_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM103345/,"NOCA_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM103346/,"PYRR_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM77548/,"PYRR_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM77718/,"PYRR_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM77780/,"PYRR_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM77781/,"PYRR_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM77856/,"NOCA_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped
awk '{sub(/BM77978/,"NOCA_rural",$1)}1' b10k_filtered.ped > tmp && mv tmp b10k_filtered.ped


../hapflk1.1_linux64 --file b10k_filtered

python -m pip install hapflk -e /scratch/dnjacks4/cardinalis/to_b10k/hapflk-1.4/hapflk

practical/data/NEU --outgroup=Soay -p practical/NEU
/scratch/dnjacks4/cardinalis/to_b10k/pcadapt/pyrr.filtered.bed



Annotate SNPs in SnpEff (Cingolani et al. 2012)
Genes within islands of differentiation were identified using BEDTools v2.27.1 (Quinlan et al. 2010). Gene Ontology (GO) term enrichment analyses were performed using the R package clusterProfiler

library(clusterProfiler)


,"Plekhf2","Ctnna3_0","Lrrtm3","Mrps5","Dlg2","Rexo1","Ccdc27","Cep104","Dffb","Kcna6","Lrrc47","Mtco2_0","Slc26a5","Tp73","Tp73as1")

# gene <- names(geneList)[abs(geneList) > 2]


genelist <- c('157657','427538','420234','423648','423649','416706','419024','428327','101751184','419379','374158','100858067','419380','4513','417715','419382','57212','423420','396000','424968','416123','423421','423421','100858906','751791','418379','420327','420325','420326','420286')

nocaurban_geneList <- c('157657','427538','420234','423648','423649','416706','419024','428327')

nocarural_geneList <- c('157657','427538','420234','423648','423649','416706','101751184','419379','374158','100858067','419380','4513','417715','419382','57212')

pyrrurban_geneList <- c('157657','427538','420234','419024','423420','396000','424968','416123','423421','423421','100858906','751791','418379')

pyrrrural_geneList <- c('157657','427538','420234','423420','396000','424968','416123','423421','423421','100858906','751791','420327','420325','420326','420286')

Change <- clusterProfiler:: bitr("ENSG00000172765",
fromType = 'ENSEMBL',
toType = c('ENTREZID','GENENAME'),
OrgDb = "org.Hs.eg.db",
drop = TRUE)

library(GOSemSim)
ggGO <- godata('org.Gg.eg.db', ont="MF")
hsGO <- godata('org.Hs.eg.db', ont="MF")

GOSemSim::geneSim('157657','427538','420234','423648','423649','416706','419024','428327', semData=ggGO, measure="Wang", combine="BMA")

mgeneSim(genes=c("157657","427538","420234","423648"), semData=ggGO, measure="Wang",verbose=FALSE)

clusterSim(nocaurban_geneList, nocarural_geneList, semData=ggGO, measure="Wang", combine="BMA")


nocaurban_ggo <- groupGO(gene     = nocaurban_geneList,
               OrgDb    = "org.Gg.eg.db",
               ont      = "CC",
               level    = 3,
               readable = TRUE)

nocarural_ggo <- groupGO(gene     = nocarural_geneList,
               OrgDb    = "org.Gg.eg.db",
               ont      = "CC",
               level    = 3,
               readable = TRUE)

pyrrurban_ggo <- groupGO(gene     = pyrrurban_geneList,
               OrgDb    = "org.Gg.eg.db",
               ont      = "CC",
               level    = 3,
               readable = TRUE)

pyrrrural_ggo <- groupGO(gene     = pyrrrural_geneList,
               OrgDb    = "org.Gg.eg.db",
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)
nrow(ggo)

ego <- enrichGO(gene          = pyrrrural_geneList,
                universe      = names(geneList),
                OrgDb         = "org.Gg.eg.db",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
head(ego)
nrow(ego)

pyrr.geneList <- data.frame(gene.pyrr.interest=c("Chga","Col6a1","Fxr1_0","Fxr1_1","Ift122","Itpk1","Pcbp3_0","Rho","Tmcc1","Dcbld2","Efr3a","Kcnq3","Oc90","Rad21"))


 java -jar snpEff.jar download GRCg6a.105

java -Xmx20g -jar snpEff.jar build -v b10k 2>&1 | tee b10k.build 

java -jar snpEff.jar build -gtf22 -v b10k -noCheckProtein
java -jar snpEff.jar build -gtf22 -v b10k


java -Xmx8g -jar snpEff.jar -v b10k VYXE01005221.vcf > VYXE01005221.ann.vcf
java -Xmx8g -jar snpEff.jar -v b10k noca.filtered.geno25.maf1.vcf > noca.filtered.geno25.maf1.ann.vcf
java -Xmx8g -jar snpEff.jar -v b10k pyrr.filtered.geno25.maf1.vcf > pyrr.filtered.geno25.maf1.ann.vcf


java -jar snpEff.jar build -gff3 -v ncbi_cardinalis 

java -Xmx4g -jar /path/to/snpEff/snpEff.jar -c /path/to/snpEff/snpEff.config -v b10k VYXE01005221.vcf > VYXE01005221.ann.vcf

https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_b10k.zip

awk '!/CARCAR/' b10k_original.gff


	CDS check:	b10k	OK: 0	Warnings: 0	Not found: 17322	Errors: 0	Error percentage: NaN%
FATAL ERROR: No CDS checked. This is might be caused by differences in FASTA file transcript IDs respect to database's transcript's IDs.
Transcript IDs from database (sample):
	'gene-CARCAR_R09680'
	'rna-gnl|WGS:VYXE|CARCAR_R15298_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R15299_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R15296_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R15297_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05925_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05926_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05927_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05928_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05921_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05922_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05923_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05924_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05920_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05929_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05914_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05915_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05916_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05917_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05910_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05911_mrna'
	'rna-gnl|WGS:VYXE|CARCAR_R05912_mrna'
Transcript IDs from database (fasta file):
	'lcl|VYXE01019281'
	'lcl|VYXE01019280'
	'lcl|VYXE01019285'
	'lcl|VYXE01005875.1_cds_NWT24985.1_3657'
	'lcl|VYXE01019282'
	'lcl|VYXE01000495.1_cds_NWT22072.1_401'
	'lcl|VYXE01009046.1_cds_5814'
	'lcl|VYXE01012859.1_cds_8286'
	'lcl|VYXE01005915.1_cds_NWT25076.1_3765'
	'lcl|VYXE01019288'
	'lcl|VYXE01019296'
	'lcl|VYXE01009593.1_cds_NWT27093.1_6071'
	'lcl|VYXE01019295'
	'lcl|VYXE01023648.1_cds_NWT35018.1_14961'
	'lcl|VYXE01023692.1_cds_15104'
	'lcl|VYXE01012799.1_cds_NWT29024.1_8240'
	'lcl|VYXE01017982.1_cds_NWT31844.1_11413'
	'lcl|VYXE01007327.1_cds_NWT25878.1_4683'
	'lcl|VYXE01019297'
	'lcl|VYXE01019263'
	'lcl|VYXE01019262'
	'lcl|VYXE01019261'

In order to fix this, i need to convert the transcript IDs in the GTF file to reflect the Transcript ID in the CDS.fa file




for line in gtf file,
if it contains Genbank gene AND DOES NOT contain pseudogene,
search the cds file for a match to $16 (without the quotes and semi colon), take the first column of the cds file without th >, and replace the gtf file's $12 with it.
done


# actual script



awk '/gene/ {print "\"" substr($1,2) "\";", substr($2,7, length($2)-7);}' cds.fa > dict.txt


head -n 3 old_genes.gtf > header.txt
tail -n +4 old_genes.gtf > working_genes.txt 
cat header.txt > genes.gtf 

awk '{OFS = "\t"} NR==FNR { genes[$2]=$1; genes[$3]=$2; next }  {if ($0 ~ /gbkey/) {$12=genes[substr($16,2,length($16)-3)]} else {$12=genes[substr($14,2,length($14)-3)]}}1' dict.txt working_genes.txt > working_genes_temp.gtf

awk -v OFS='\t' -v col=9 '{FS = "\t"}
    {$(col)=$(col)" "$(col+1)" "$(col+2)" "$(col+3)" "$(col+4)" "$(col+5)" "$(col+6)" "$(col+7)" "$(col+8)" "$(col+9)" "$(col+10)" "$(col+11)" "$(col+12)" "$(col+13)" "$(col+14)" "$(col+15)" "$(col+16)" "$(col+17)" "$(col+18)" "$(col+19)" "$(col+20);              # merge col and col+1
    print $1,$2,$3,$4,$5,$6,$7,$8,$9;                # remove the last field
}1' working_genes_temp.gtf > genes.gtf





WARNING_TRANSCRIPT_NOT_FOUND: Cannot find transcript 'TRANSCRIPT_CDS_VYXE01000004.1_19854_19961'. Created transcript 'TRANSCRIPT_CDS_VYXE01000004.1_19854_19961' and gene 'GENE_CDS_VYXE01000004.1_19854_19961' for VYXE01000004.1	Genbank	CDS	19853	19960	-
. File '/Users/danjack/snpEff/./data/b10k/genes.gtf' line 33	'VYXE01000004.1	Genbank	CDS	19854	19961	.	-	0	gene_id	"CARCAR_R04070";	transcript_id	"lcl|VYXE01000004.1_cds_NWT21712.1_1";	gbkey	"CDS";	gene	"Atp2b4";	locus_tag	"CARCAR_R04070";	orig_transcript_id	"gnl|WGS:VYXE|CARCAR_R04070_mrna";	partial	"true";	product	"AT2B4	ATPase";	protein_id	"NWT21712.1";	exon_number	"11";'

# to fix the protein file, i need to take the $2 from the protein.fa, use it to search for and find the 

awk '/>/ {print $2}' protein.fa

awk '{OFS = "\t"} NR==FNR { genes[$2]=$1; genes[$3]=$2; next }  {if ($0 ~ /gbkey/) {$12=genes[substr($16,2,length($16)-3)]} else {$12=genes[substr($14,2,length($14)-3)]}}1' dict.txt working_genes.txt > working_genes_temp.gtf

awk 'NR==FNR { genes[$2]=$1; genes[$3]=$2; next }  {if ($0 ~ />/) {print genes[$2]}}1' dict.txt temp.txt 

awk '/>/' protein.fa > temp.txt
> working_genes_temp.gtf




### didn't really work. it worked but it's impossible to read now
# i need to instead replace all copies of "transcript ID" with the value for "gene"
awk '/gene/ {print "\"" substr($1,2) "\";", substr($2,7, length($2)-7);}' cds.fa > dict.txt


awk '{OFS = "\t"} {if ($0 ~ /gbkey/) {$12="\""substr($16,2,length($16)-3)"\";"} else {$12="\""substr($14,2,length($14)-3)"\";"}}1' working_genes.txt > working_genes_temp.gtf

awk -v OFS='\t' -v col=9 '{FS = "\t"}
    {$(col)=$(col)" "$(col+1)" "$(col+2)" "$(col+3)" "$(col+4)" "$(col+5)" "$(col+6)" "$(col+7)" "$(col+8)" "$(col+9)" "$(col+10)" "$(col+11)" "$(col+12)" "$(col+13)" "$(col+14)" "$(col+15)" "$(col+16)" "$(col+17)" "$(col+18)" "$(col+19)" "$(col+20);              # merge col and col+1
    print $1,$2,$3,$4,$5,$6,$7,$8,$9;                # remove the last field
}1' working_genes_temp.gtf > genes.gtf

mv cds.fa old_cds.fa 

awk '{if ($0 ~ />/) {$1=">"substr($2,7,length($2)-7)}}1' old_cds.fa > cds.fa

mv protein.fa old_protein.fa

awk '{FS = "\t"} {print $24}' genes.gtf | head
awk '{OFS = "\t"} {if ($0 ~ /gbkey/) {print $24}}' genes.gtf | head
genes.gtf


awk '{OFS = "\t"} NR==FNR { genes[$2]=$1; genes[$3]=$2; next }  {if ($0 ~ /gbkey/) {$12=genes[substr($16,2,length($16)-3)]} else {$12=genes[substr($14,2,length($14)-3)]}}1' dict.txt working_genes.txt > working_genes_temp.gtf



java -jar snpEff.jar build -gtf22 -v b10k


java -jar snpEff.jar build -gtf22 -v b10k -noCheckProtein

grep 'Rad21' noca.filtered.geno25.maf1.ann.vcf | wc -l
grep 'Rad21' pyrr.filtered.geno25.maf1.ann.vcf | wc -l

# phylogenetics on just Ch037
VYXE01005221.1 2176060 2178234

awk '$0 ~ /\##/' VYXE01005221.vcf > Ch037.vcf

awk '{$2 > 2176060 && $2 < 2178234}1' VYXE01005221.vcf | awk '$0 !~ /\##/' >> Ch037.vcf

bgzip Ch037.vcf

bcftools index Ch037.vcf.gz

~/vcf2phylip/vcf2phylip.py -i Ch037.vcf.gz

python ~/sula/filter_invariants_all.py Ch037.min4.phy

MSB25201   RYWCRARNNAGCTCCTGTTCACCTAGCACACACYYTCYYAYCKMMYCAGCCATCCCAC
NOCA003    ACRCANNNNNNNNNNNNNNNNCTARGCAYGCRYCTTTCCACYKYATCARYCRNNNNNN
NOCA004    AYRCRNNNNNNNNNNNNNNNNCTARGCAYGCRYYYTYCCACCKCNNNNNNNNNNNNNN
NOCA006    AYGCRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCACCGCATCANNNNNCCYAM
NOCA008    AYGCGNNNNNNNNNNNNCCCACTAGGCATGCRYYYTYCCRSCGCATCANNNNNNNNNN
NOCA012    AYGMRNNNNNNNNNNNNNNNNNNNNNNNNNNACYYTCCCACCKCATCNNNNNNNNNNN
NOCA013    ACRMRAANNNNNNNNNNNNNNNNNNNNNNNNACYYTCNCACCKCATCANNNNNNNNNN
PYRR003    ACACAAAYCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
PYRR004    RYWCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYCCACCKYATCARYNNNNNNNN
PYRR006    AYRYANNNNNNNNYYYSCCYMYTMARYMCGYACYYTCCCACCGCATCANNNNNNNNNN
PYRR007    ACRCRNNNNRRYWNNNNNNNNNNNNNNNNNNNNNNNNCCACCGCATCANNNNNNNNNN
PYRR009    NNNCANNNNNNNNNNNNNNNNNNNAGCACGCACYYTCYCACCKCATCANNNNNNNNNN
PYRR011    AYWYRNNNNAGCTCCTGTTCACCTAGCACGCACCTTCYYATCKMMYCAGCCATCCCAC
BM100619   AYRCRRACCARYWYYYSCYYMYTMRRYATGCGTCTTTCCGGCTAMYYRMYYACACCGC
BM100620   AYACRGACCNNNNNNNNCCCACTAGRYMYGCRYCTTTCCRSCTMCCTGCCTACACCGC
BM100621   AYRMRGACCNNNNNNNNCCCACTAGGCATGCGTCTTTCCGGCTACCYGCCTACACCGC
BM103345   RTGYRNNNNNNNNNNNNCCCACTAGGCATGCRYYYTYCCRSCKMMYYRCCTNYMNCRC
BM103346   RYWYRAGCCAGCTCCTGTTCACCTAGCACGCACCTTCTTATCTACTCAGCCATCCCAC
BM77548    ACACAAGCCAGCTCCTGTTCACCTAGCACGCACCTTCTTATCTACCCAGCCATCCCAC
BM77718    ACRCAAGCCAGCTCCTGTTCACCTAGCACGCACCTTCTTATCTACCCAGCCATCCCAC
BM77780    RYWYRAGCCAGCTCCTGTTCATCTAGCACGCACCTTCTTATCTACCCAGCCATCCCAC
BM77781    NNNNNAGCCAGCTCCTGTTCACCTAGCACGCACCTTCTTATCTACCCAGCCATCCCAC
BM77856    ACRYRGACCARYTCCTGCCCACTAGGCATGCGTCTTTCCGGCTACCYRCCTACACCGC
BM77978    NNNNNGACCNNNNNNNNCCCACTAGGCATGCGTCTTTCCGGCTACCTGCCTACACCGC

# not enough variant sites for RAxML to be effective

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/Ch037.vcf -o /scratch/dnjacks4/cardinalis/to_b10k/fst/shared_region/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 





# phylogenetics on just RHO
awk '$0 ~ /\"Rho\"/ {print $1,$4,$5}' genes.gtf
VYXE01006626.1 115846 118575

awk '$0 ~ /\#/' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > rho.vcf

awk '{$2 > 115846 && $2 < 118575}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01006626/' >> rho.vcf

bgzip rho.vcf

bcftools index rho.vcf.gz

~/vcf2phylip/vcf2phylip.py -i rho.vcf.gz

python ~/sula/filter_invariants_all.py rho.min4.phy
24 18
MSB25201   CCTGGCCCTCTCCTTTAA
NOCA003    NNNNNNNNNNTYCCCCCA
NOCA004    NNNNNNNNNNNNNNNNNN
NOCA006    TYCAATTYYYYCYCCCCR
NOCA008    TYCAATTYYTCCTCCCCG
NOCA012    NNNNNNNYYCYCYCCCCR
NOCA013    TYCAATTNNNNNNNNNNN
PYRR003    NNNNNNNNNNNNNNNNNN
PYRR004    YCYRRYYCTCTYCNNNNA
PYRR006    TNCAATTCTCTYNNNNNN
PYRR007    YCYRRNYCTCTYCCCCCA
PYRR009    YNNNNNNNNNNNNNNNNN
PYRR011    TCCRRYYCTCTCCYYYMA
BM100619   TTCAATTTCTCCTCCCCG
BM100620   TTCAATTTCTCCTCCCCG
BM100621   TTCAATTTCTCCTCCCCG
BM103345   TTCAATTTCTCCTCCCCG
BM103346   CCTGGCCCTCTCCTTTAA
BM77548    CCTGGCCCTCTCCTTTAA
BM77718    CCTGGCCCTCTCCTTTAA
BM77780    CCTGACCCTCTCCTTTAA
BM77781    CCTGGCCCTCTCCTTTAA
BM77856    TTCAATTTCTCCTCCCCG
BM77978    TTCAATTTCTCCTCCCCG


mkdir raxml
cd raxml
echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-10' > partitionfile.txt
~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n rho -q /scratch/dnjacks4/cardinalis/to_b10k/fst/rho/raxml/partitionfile.txt -s/scratch/dnjacks4/cardinalis/to_b10k/fst/rho/variantsites.phy -T 6 -p 12345 -N 10 ­-b 12345 -V

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/fst/rho/rho.vcf.gz -o /scratch/dnjacks4/cardinalis/to_b10k/fst/rho/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 




# phylogenetics on just Hydin
awk '$0 ~ /\"Hydin_0\"/ {print $1,$3,$4,$5}' genes.gtf
VYXE01006187.1 gene 5864 13079

cd 

awk '$0 ~ /\#/' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > hydin.vcf

awk '{$2 > 5864 && $2 < 13079}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01006187/' >> hydin.vcf

bgzip hydin.vcf

bcftools index hydin.vcf.gz

~/vcf2phylip/vcf2phylip.py -i hydin.vcf.gz

python ~/sula/filter_invariants_all.py hydin.min4.phy
24 602
MSB25201   CYKRRSRMRCCYYCRRSTCYGYWAYGKGAMRRYGCGCTMYAASRRGACRMTTMKRCCYMWYMRTTGSCRCGYGMGSTGCYAGRWKYTRMTCYSCCCGAGCSRCCTKGTAYARASRCYKKACRSCYGMACCCAGGCRYSKRYKMKYKMGCATCGAARMYKASSMRAMTRMTYGCGKATTMTCAYYCGAAGCAAAYSKMRSTARSGYKCSWYMAARSRYTWASMCYMMTYTGRTTSRTYMSTMWGRRCYYMKMRKMTWYMKYCRKRSRYSSSRMSMSWMYMGRYKSRRRSSYCTKRRCKMSRYTKCSASRKKYRRMRYSTWKKSGATTMKKSYMMCKRCGYSMGTRSGCGCTGTARRMRARAAAYGGYAATCRCRAGCGKACTCGTGTCGGCGGCGGGCCCASYSTGCSTCSGATGYYSRTGGGKCAKRTGWWWSCAYYRSYYGCMGYCYACCAGKRTMKSRCASMYKYYCACAGKGCCCWRMMRGSSSSYMKMWYSRGTGCYGGYYACCCRYRCKYCSGKKYTNYRYARWMYYKMCATCRRSKSRMKKSTYRCTRSRSKSMASSSSGTGYGRWRNSKMRSSYASGYCYRCMRYCYYNSRYSKY
NOCA003    CCGGGCGMNCCTNCAAGTSNAYNRYRGGNNNRYGYGCNNCARCRGSACAMTYAGRYSTMACCATTSSSRMSNNAGSTGNTMTATGCTRCKCCSCCCGRGCGGCYKKGYACRARCACCTKRCGCMCGAAAMSRGKCRYSKRYKAKYKAGSRKSKARRMCKMSSMARMYGCYYGCGTATKMYYRCCCGWGGCMMRTGGCGGTRGCGYKCSWYMARAGGCKAACMMYMMTYTSAYYCGYNNNTCAGGGCTTCTCRTCTNNMGTYRNRSGTSGSGMCMSWMTCRAYTCAAASCTCTTGAYGMSRYYKCSRGGKTCGACNYGTWKGGSMTKMKTCTMMCKGSGTSARTGGKMKCYRYWARMARGRARYRSYRMYMRYRRSSRKAYTMRYNNNNNNNNNNNGCSSRSTGTGCSKSGGRYSNTCGYRGSKYWGRTRNTTCCACTRCCTGYMGTMYWSSAGTRTMNGGCANNNNTCCACAGGNCNNNNCCGGSGGSYMTMWYSANNNNNNNNCAYNYAYRCKYNNNTKCTTTNYARWATYKCCRTCRRSKSRMKKSKYAYKGGGGTCMRSCNSGYGTGRTGGGKCGCSYASSTYYANMRYMYYWSRYSTC
NOCA004    CCGGGCRCACCTTMAASYSYRTTRYGKGWARACSCRCYMYARCGRSRYRMTYMGGCSYMACCAWYSGCGMSCGMNSYRYTMTATKNNNNKCYCCMCGRRCSRYCKKKYAYAARSRSCTTACRSMCGMMMMSARGCRYSKGYKMGTTMGCAKSKMRAMYGMSSMARMYGCYYGYGKRTTCYCACYCGWAGCAAATGGCGCKARCGTGCSWYMARRSRYKWWSMCYMMWYWCACCCGYTACKCAKGGYTTMKMATCTWYMGTCATGGGTGGSGMCCGANYCGGCGGGRGGSYSKKRRCKASGCYGCGRSRKTYGACRTSYWGGGSMTKMKKCTMAMKRSGYSAGTRGKCGCNNNNAGMAARRARYRSYRMYCRCGAGSRKAYKMRYRWMSRTRSYGGGCCSAGTSKGSSTSGGRYSTTCGYRGSKYWGRYRAWWCCACYRSYYRYMGTCTWSSRKTNNNNNNNNSCCKCCMMSMKNKCCCWRMCNSSGGGTCKMWCNRNNNNNNNNNNCSYRYRYKTYCSKKCTKYRYAAAMYYKMSAYYRAGGCGMKKSKYRYKGGGSKSMRSSCSSYRYKAWRKSTMACGYASSYYTACARYCCCTCGTCTC
NOCA006    CCGGGCRMRCCTYCAASTCCATTRYRKKAAARCGCGSTCYRACGRGACRAWTMGRTSYMACCATTSNCGCGCGMNSTGCTMKRTKCTACKCYCCCCGRGCGGCCKTGTRYRRRCRSCTKRCGCMYGMMMCCAGGCRCSKGYKMNTKMGCAKSGARRMCKMSSMAAMYRCYYSCGKATTCYYRNYMSWRGCMMRTSKMRSKRGSGYTCCACAARAGGCKAWSCMTCCTYWSRYYCGYYMSKMWKRRYYYMKMRKCTWYMKYYRNRNRYSSSGMCANNNNNNAYTCAAANCTCTTGAYGAGRYYKYSRGGGKTRRCGTSTWTGSGANNNNNNNNMMKRSGTSAGTRGGCGCTGYAARMRRGRARYRSYRMYMAYRRSSRKAYKMRYRWMSRTRSYGKGCSSAGTGKKSSKSGSRYSTYCGYRGSKYWKATRNWWCCACYACYTGCMGTMTASSRKKRYATSNMMSMYTTCCNNNGKGNNNNNNNNGCGGGTMKMWYCARYKYCSSYYACSYRYRYKYYCGGKYTKYRYARWAYYKMSAYYAACGCRCKKCKYRYKRSRSKSCAGGCSSTRYKRWRKSKMRSSYASGYCYRYMGYMYYWSRYCKC
NOCA008    CCGRGCGCACCTTMRAGYCCATTATGKGAARRYSCGCYMNAANNRSRCRMTYAGGCCTMACCATTSNNNNNCGASSYGCTAGRTKCYRCTCTSCMMKRRSGGCYTKGYRYRRRCASCTTRCRCCCGMMACCAGGCRYSKGTGMGTTAGSRKSKARAMCKMSSCARMYGCYYSYGKRKKCTCACYCGARKSAAAYSKMRSKNRSRTGYSWTMWAASRYKWWGCCTCCTCWSAYYSGYYMSTNWGRRYTTAGANNCTWYCGTYANGGGTGGGGCCAGACYCGRYKSRRRSSYSKKRRCGMSRYYKYSAGRKTYRACAYGYWTGSSAYTMKKCYAMCKGCGYSANTRGGCGCYRYAAACARGARRYRSYRMYMAYRRSSGKMYKAAYGTMSRTGCYGNKSSSASTGKGCCKSSGRYSYTSGYRRSKYAGATRAWWSCACYRSYYRYMGTMTWSSRKTRYMKSGCNCCCTNNCACAGGGSYSWAMCRNNNGNTCTAACNRRYKTCGGCTACSYRYRYKYYSSKKYKKYNCRAAATYKMSRYYGASGCGMKKSKYRYKGGGSKSMRSSCSGYGCGRWRNSTCASSTACGTCYRNMGYMYYWSRYCKC
NOCA012    CCGRGCRCACCTYMRASYSCATTNYGGGWARRCGCNNTCNANCRRSAYAMTYMGGCSTMAYMATTSGCGCGYNANSTGCNAKRTNCTACKCCSCCCGRRSGNCCNGGYRYRRRCACCKKRNASMYGMMNMCARGCRYSKRTGMKTTANNNKSKMRRMCKMSSCRAMYGCYYSYGKRKKCYCRCCCGARGCMMRYSKMRSKRRCGTTCGWYMAAASRYKWWGMCYCCTCWSRYYCGYYMSTMWGRRCYYMGAAGCYTTAGCTGTACGCCCSRCCMSWMCCGRCKNRRRNSYSKKNNYGMSRYYKYSASGGKTRRMGTGTWKGGGAYTCKKCYMACKGSGYSMRYRGGCGCYRYAARMARRRRRYRSYRMYMRYRRSSRKMTNAAYNNACATGCCRGKSSSRSTGNGCSTCGGATGYTGGYGNNNNANRTGAWWSCACYACYTGCMSYMTACCAGTRYATSRMMSMCKYCCACMKGGSYSWAMCGGSSSSYMTMAYNRRYKYYSGCCACSYRYRCKTYSSKKYKTTNCARAAYYKMSRYYGASKSRMTTSTYRCTRGGCGGCAGGNNNNNNKAWAKSTMASSTACSTYTACAGCCCYTSRYSTY
NOCA013    CCGRRSRMRCCTTCRASYSCRTTRCGKKAMRACGCRSTCCAACGRGAYRMTYMGRCSTAAYMRTTGGCNCNCGAGSTGCTMTRWGCTACNNYSCCCGRGCGRYYKKKYACRRASRSCTKRCRSMYGMAMMSRRKMRYSKRYKMGTTCGCAKSGARRMCKMSSMRRMYGCYYSCGKRKKMYCRCYCGWRGCAARTGGCGSKRRSGYKCSAYMARASRYKWWSMCYCCWCTSRYYSGYYMSTMWGRRCTYMKMRKMTTCAKTCGKACRTGGSGMSMGAMYCGRYKSRRRNSYSKKRRYKMSRYTKCSAGRKKYRACRYSYWKKSSAYTMKKCTCMCKGSGTSARTRGGCGCTRYWARMRRGAARYRSYRMYMAYRRSSGKAYKMRYRWMSRYRSCRNNNSCRGTGNKCSKSSGRYSYTSGYRGSKYWGAYRAWACSACYACYTRYASTMYACNAGTATMKSRCACCYKYYCACMNKGNNNNRCCNNNNSSYMKMWYSARYKYNNSYNNCCYRYRCKYCCGNNTTGCACAAAAYCKCCRTCARCKSGMKKSKYRYKGGGSKSMRSSCSGYGTGATGGGTCGCGCAGGYCTRYMGCCCCTCGTCGC
PYRR003    NNNGGCRMRSYTYCGASTCCRTTATRNGAMGGNGCGCYMCNASRRGNNAMTTANRCCTAACCATTGGCGCGYKMGSYRYTATAWKCTRCTMTCYMMKAGCSGCYTKGTACNNNNRCCKKRCASCYGMAMCNNGGCACGTGTTAKYKAGSATCGAAAACTAGSMRAMTRMTTGCRTATTCTCAYYCGAAGCAAATGGCGCNNNNNNNNNACAAARGGCTAAGACYMMWCTGRTTCGTYAGTMWGRRCYYAKAAKCYNTMGYYGTRSGYSSSRMCMSWMYCGRYKSRRRSSYSKKRRYGMSRYTKYSANNNNNNNNNNSTWKKGGATTMKKCCAMCGGCGTSMGYRSGCKMTRTWARMAARAAAYGGYAATMRYRRSSGKMYNMRYGTACRTGCCRGKCSSRSTSTGCCTCGGATGYYSRYGGGGCAKRTGATTCCMYTRSCTGCNGTCTACCAGKRTATNRCACNNNNNCASAGGGCCCTRCMRSSSSGTCKMACCAGTGYYSGCYMCSYANNNNNNNNKKYKTTRYRRWAYNKNNNNNNNGGSRMTTGTYRCTGGGCGGCAGGSSSTGYGATNNNNNRCGYASGTCTNNARNCCYWCGTCTN
PYRR004    CYGRRCRMACCTYMRAGYSYATWAYRKKWAARYNCRSTCYAASRRGACRMTYMKGCCTAWYMRTTGSCGCGCGMSSYGCTMTRWGCYACTCYSCCCGAGCNGCYTKGTAYARASRSCTKAYRSMYGMAMCSRGKMACSTGTKMGTKAGCATCGAARACTACSCRAMTRCYTGCRTATKMTCAYYCGAAGCAAATGGCRSTAGCNNNYSWYMWAASRYTWASMMTCCTTTNNTTSRTYMSTMWGRRCYYMKMAKNYWYMGYYRTRSGYSSGGCSCSAMYCGRCKCRRRGSYCTKRRYKMSRYTKYSASRTTYGRCRYSTWKGSSMTKMKKNYAACKRCGTSMGYRGKMKMTRTAARMRARAAAYGGYAATMRCRAGSGKAYKMRYRWMSRYGSCGGGCSSRSTSTGCCKCGGATGYTSGYGGGGCAKRTGNWWSCAYTGCCYGCMSYCTWSSRKTAYMKSRCACCYKYCCACMGKGCCSWRMCRSSSSSYMKMWCCAGTGYYGGCYMCSYRYRYKYYSSKKYKKYRYRRWMYYKMCRTCRRSKSRMTKSTYRCTRSRSKSMASSSCGTGNGRTRGSKMGSGYRSSYYYACAGYCCCWCGTCTC
PYRR006    SCGGGCGMANNTTNAAGTNNNYTRNGNGWARGCSYGCTCCARSGRSACAMTYMKGCCTMANNNWYGSCGMSYKMGSYRYYMKAWKCYACTCYCCCCGAGCGRYCTKKTAYARACACYKKRCASCCGCAMCCAGGCACSNNTKAGTKAGSATCGMRRAYTACSCRAMTGCTTGCGTATKMTCANYMSWAKCAAAYGGNGNNNNCGTGCNWYMAARSRYTWNGACYCCTCTGRTTCRTYASTMWGRRCNNNNNNNMTTTMKTCRKRSRTSSGGCCMSWMTCGAYTCAAASCTCTTGANNMSGCTGCNAGRNTCGACACGTAKKGGAYTMKKNYMAMKRCSTSMGTRGKCGCTGTAAACGARAAAYGGYAATMACRAGSGKMYTMRYRWMSRYRSCGGGCCCRSYGTGCCTSSSATGYTSGTGRGGCAKRTGATTCCMYTRSCCGCMGYCTACCAGTRTMTGGCASCYKTCCACMGKGCCCWRCCGGSSSGTMKMAYSRGTGCYGGCYACSYAYRCGCCCGKKYTKYRYRRWMYCGCCANNRASNSRMKKSKYRYKGSRSKSMASSCGSYRTKRWGGGKCGCGCRGGYCTACMRTCCYWCGTCTY
PYRR007    SYKRRSRMACCTYCGASYCCRYTAYGKGANNAYGNGCTCCNNSRRGACAATTAGGCCYMWYCATTGSCGCGYKAGSTGNTATATGCYANTCYSCCCGAGCNGCCTGGTAYAAASGCYKGAYACMCGMAMCSRGKMAYSKGTGMKYKASCATCGAARAYTACSCAAMYGCYTGCRTATKCTCACYCGARGCAARYSKMRSKARSGTGYGWYMAAASRYTWASMMTCCTYTGRYTSRTYMSTMWGRRCYYMKMAKCYTTMGYYRTRSGYSSSRMCMSWMYMRRYKSRRRSSYSKKRRYGANRYTGCGAGRTTYRAMRYGTAKKGGATTCKTCCMACKGCSYSMGYRSGMKCTRYAARMRARAAAYGGYAATMAYRRSCGKMYNMRYGTMSRYGSCRGGCCCASYSTGSSTSSGATGTTCGTGGSKYAGRTRATTSCAYTGCCYGCMSTMTACCAKTRTCGGGCACCYTTCMMSAKKKSYSWRMCRGSSGGTCGAWYCAGTGCCGGCCAYSCAYRYKYYNGKKYTNYRYARWMYYKMCRTCRRSKSRMTKSTYRCTGSRSKSMASSCCGTGTGGWRGGTCRCGYRSGYCTACARYCCYWCGTCTC
PYRR009    CCGGGCGNNSYYYCRRSTCCAYWAYGNGWMRRNNCGSTCCAACRRSRCRMTYMGGYCTNAYCANNNGCGCGCGASCTGCTMTNNNCTRCKCYCCMMKARSSRCCTKGTACARACASCTTRCACCCGCAMCNNNGCAYSKRTKMKYTAGCATCGARRAYTACSCRAMTGCYTGYGKATKMTYAYYMSWAGCAAAYGGCGSTAGCGNNYGWYMAAAGGCTAAGACTMMTYTGGTTSRTYASTMWGRRCYYMKMAKCYTTAGCCRTRSGYSSSRMSMSWMYMRRYKSRRRSSYSKKGAYGMSRYTKYSAGRTTYRRMRTSTATGSGATTCKKNYMMCKRCGTSMGYGSGCGCTGTAARMRARAAAYGGYAATMRYRRSSGKMYNMRYGTMGGTGCYGKKSCCRSYSTGSCKCSSATGTYSRYGGGGCAKRTGATTSCAYTRCYYRYMGYCTACCAGTRTATGGCACCCTYNCACMKKGNYSARMMRGSCGGNNNNWYSAGTGCYGGCYACNNNNNCKYCSSKKYTKCACAAAMYCKCCRTCNNNNNNMTKSTYRCTGSRSKSMASSSSGYRYKRWRKSKMGCGCAGGTCNACAGYCCYWCGTCKC
PYRR011    SYKRRCAMACCTYCGRSYSCATTRYRKKWAARYNCGCTCYAACGRSRCRMTYAGRCCYAWYMRTYGSSGCGCGMSSYGCTMKRTKCTRCTMYSYMMGRGCNGYYKTKTACRRRCRCCTGACACCYGCAMCCAGGCACGKGTKAGTKAGCATCGMARAYTASSCRACTGCTTGCGTATKCTCACYCGARGCAAAYSKMRSTAGSGTGYGNNNNNNGGCTAASCMTCNTYTGRTTSRTYAGTMWGRRCYYMKMRKMYWYMKYYRKRSRYSSSGMSMSWMYMRRYKSRRRGSYSKTRRYGMSRYTKCSAGRTTCGACRYGTAKGGGATTCKTNYMACKRCGTSAGTRSGCGCTRTARRCRAGAAAYGGYAATMAYRRSSGKMYKMRYGTMSRYGSYGKGCCCRSYSTGCCTSSGATGYYSRYGGGGCAGRTGWWWSSMCYRSYYRCMGYMTWSSRKKRYATSGMMCCYKYCCMSMKKKCCCTRCCGGSCGGTCKMACCARYKYYSSYYACCCRCRCGCCCGGGTTNCACAAAMCCGCCATNRRSKSRMKKSKYRYKRSRSKSMASSSSSYGYKRWRKSKMRSGCAGGTCYACAGYCCYWSRYSKY
BM100619   CCGRRCRMRSYYYARRSTCYRTTRTRGGWMRANSYGNYCCAACGRSAYRMWYAGGYSTMAYMRTTSGCGCGCGAGCTGCTAKRTKYTRMKMYSYMMKARSSRCCKTGTRCRRACASYKKRCRSMCGMMNMCARGCRYSKRTKMGTTMGSRKSGMRRMCKMSSMRRMTRMYYGYGKRKTCYYRCYCGWRGSMMRYSGMRSTRGSRYKCCACAARRSRYKWWSMMYMMTYWSRYYSRYYASKCAKGGYTYMKMRTMYWYMGTCATGGGYSGSRMSMSWMYMRRYKSRRRSSYCTTGAYGMSGCYGCGRGRTTYRAMRYSTWKKGGATTCKKCCAMCKRSSTSMRTRSGCGMYRYWARMRRGARRYRSYRMYCACRASSRKAYKMRYGTMGGYGSYGKKSSSRSTGTKCCKSGSATSYTSGTRRGGYWGAYRAWWSSACYRCCYRCMSYMTACCAKTATAKSRMMSCCKYCMMSMKKKCYSWRMMRNNGSSYMKMNNNRGYKTCSSYTACCNNCGCGCCCGKKYKKYRYARWMYYKMSRYYRRSKSRMKKCKYRYKRSRSKSMASSSCGTGYGRTGGGTCRSGYASSYYYRYMGNCCCTCGTCKC
BM100620   CCGGGCGMRCCTYCAASTCCATTRTGKNWAARYSCRCTCCRRCRRSAYRCWYMKGCSTAACCATYSSSRMSYKMSCTGCYAKRTKYYAMKCYCYCCGAGCGRCCKKGTACRRASACCNTACRSCYGMMMMCARGCRYGTGYKMGTTMGSRKSGMRAMCGMGSCARMYGCYYSYGTAKKCYYRCYCGWRGCAAAYSKMRSKRRSGTGNNACAARRGGCKANSMMYCCWYTSAYYCGYTAGTNWGRRCTYMKMRKNTWYMGYCNNNNNTGGSGMSMSAMYMRACKSRRRSSYSKKRRNNMSGNTGCGRGGKTYRNCRTSYWKGGGAYTCKKCTCMCKGCGYSMRTRGGNNCYGYAAACARRARRYRSYRMYMAYRRSSRKMYTMRYNNNCRTGCCRGKSCCRSYGKGCCTSSSATGYTSGTGRSGCAGATGWWWSSACYRCYTGCCGTCTACCAGKATNNNNNNNNNNNNCASMGGGSYSTAMCRGSSGGTCKMWYSAGTGYCNNNNNNNYAYRCGYCNNKGYKTYNCAAAMYCGCCATCRRSKSRMTTSTYRCTRSRSKSMASSCSGYRYKRWRKSTMNNNNNNNNNNACARCCCYWSRTCKY
BM100621   CYGGGCGMACCTTCNRGTCYATTRYGKKWAAAYSCGSTCYAACGGSRYRMWYMKRYNTAACMRTTSGCGCGYGANGTRYTMKRTKYNANTCYSCCCGRGCSGYCTKKYACRRRCRCCTKRCRSMCRMMNCCAGGCRYSKGYGMNTKMGSRKGKAGAMCKMSSCRRCTGCYYSYGKAKKNTCACTMSWRKSAARYSKMRSTRRSRYKCCACAARRGGCKAWSMCYMMWNTCAYYCGYYACKCAKGGYTTCKCRTMYTYAKYYGKACRYCCCRACACTATMRATTCAAANCTCTTGACGMSRYYKYSRGGKTTRRMRYGTWKKGGATTMKKCYMMCKGCGYSAGTRGGCGCYGYAAACARGRARYRSYRMYMAYRRSCGKMYNMRYNNNCRTGCCRGGCCCRSYGKGNCKSGGRYSYTSGYRRSKYWGAYRWTTCCMCYRCYYGYMGTMYWSSRKTNNATSNMMSMCKTYMMCMGGGCYSWANNNGSNSSYMKAACCAGTGCYGGCCACSYRYRYKYYSSKKYTKYNCRAAAYYKMSAYYRRSKSRMKKSKYACTRSRSGGCAGGNNNNRTKAWRNSTMRSSYASSTYYACAGCCCCTCNNNTC
BM103345   CCGGGCRCACCTTCAAGTCCATTRYGTKWARNNSYGCTCCAACGGGACAMWTMGRYSTAAYCATYSSSRMSYGAGSTGCYMKRWGCTANKCCSCCCGAGCNGCCKKGYRCRARCASCTTRYRCMCRMMMCCAGGCRYSKGYGMGTKANNNKSKARAMCKMSSCRAMYRCYYGYGKRKTCTCAYYCGARGCAAAYGGCGSKARSGTGCGWYMAAAGGCTAWSCCTMCTNWCAYYCGYTASKMWKRRYTTMKMATMTTYAKYYRTRSGYSSSRMSMSWMYMRRYTCARRSSYSKKRRCKAGRYTTYCAGATTCRAMRYSYWKGGGAYTCKKCTCAMKRCSTSMRTRSKCGMYRYWARMRRGARRYRSYRMYMACRAGSRKAYTMRTRWMSRTRSCRGGCCCRGTGKGCCTSSSATSYTSGYRRSKYWKAYRAWWCSACYRCYYGCMGTMTACCAGGATCGSRCACCYTYYMANNNKGNYSWRMCRSSGGGTCGMAYSARYKYYSSYNNCSYAYRCKTCSCKKYTKYRYARWMYYKASRYYRASKSRMTTSTYRCTRSRSKSMASSNNNNNNGRTNGGNNGCGCAGGYCTACMGCMYYWSRYSKY
BM103346   CCGGGCGNANNYTMRRSTCNANNATGNKNNNACGCGCNNNAACGASACAMNTMGNCCTAWYMRWYGGSRMSYKAGCNNNTATNTKYTAMTNTNCCCGNNNSGCCTGGTAYAAACNNNTKAYASMCRMAMCSRGGCACSTRTKMKYKASCAKSGMRRACKMSSCRAMTGMTTGYRKATTCTCACNNNNAGCAAAYSKMRSTARSGTGNNNNNNNNNNNNNNNNNNNNNNNGATTSRTNNNNNNNNNNNNAGANKCYTTNGYNGTACGYSSSRMSMSAMYCGGCKSGRRSGCGGKAGCKMGGYTKYSAGNNNNNNNNNNNNKGGSMTKMKTCCAAMKRCSTSMGTAGGMKMTGTWARMAAGRAAYGGYAATMRYRASSGKMYNMRYRWMSRYGSCRGKCCCASTGTGCSTCGSATGTTGGYGGGGCANATGWTTSCMCTRCCYRYMSYCTACCAGKNTATGGCACCCKYYMMSAKTKCYSWAMCRGSSSGTCKMAYSAGTGYYGGCYMYSYACGCGYCSSKKYTKYRYARAATCKCCRTCARSNGRMKKSKTRYKRGRSKSMASSSSGYRTKRWRKSKMGCGCAGGYCTACAGYCCCTCGYNNN
BM77548    CCGGGCRMACCYYMRASTCCNYWATGKKNNNRCGCGCTCNAACRRGANAMTYNNRCCYAAYMRWYNGCGCGCGAGSYGCTMKAWKCTRCTCTCCCCGAGCSGYYTKKTACARACACCKGAYASMYRMAMCSRGKMAYSTGTKMGTTAGCAKCGMARMCNACSCRACTGCTTGYRTATKCTCAYYCGAAGCAAAYSKMRSKARSGTGCSWYMWARSRYTWNSMCYMMWYTGRYYSGTNNNNNNNNNCTYAKAAKMYWYMGYCATGGGTGGGGCSMGACCCGRCKSRGGGGCGGGAGNNNNGCTGCGASRNTCRAMRYSTATGSGMTKMTKCYAMCKRCGTSMRTRSGMKCTRYAARMRARAAAYGGYAATMRCRAGSGKATKMRTRWMNNNNNNNNGCCCRSYGTGCSKCGGATGTYSRYGGGGCAGATGWTTSCACTACCYRYMGYMTACCAGKRTCGNGCACCTKTCCASMGKGCYSTNMCGSSSGSYMTMWYCAGTGCTGGCYACCCRCGCGCCNNGGTTKCACAAAMNNNNNNNNRRSKSRMTKSTYRCTGGGCGGCAGGSCGTGYGRTGGGKCGCGCAGGYCTRCMRYMYYWCGTCTC
BM77718    CNGRRCGMRSYNYCRRSTCCANTAYGKGAAARYGCGCTCYAASGAGACRATTMKGCCTAWNNANNGSCGCGYGASSTNCYAKRWGCTAMTMYSCCCGAGCSRYCTGKTACAAACRSCTKACRCCCGANNNNNNNNACGTGTKAGTTASCATSGARRMYNACSCRAMYRMTTGCRKATTCTYAYYCGWAKSAAATGGCGSKARCGTGCSWYMWAASRYTWASMMYCCTYTNNTTSNTYMSTMWGRRCYTNKAAKCYTTAKYYGKACRYCSSRMSMSWMYCGRYKSRRRGSYNNNNNNNCCGCTGCGNNNNNCGACGTGTWKGSGMTKMKKCTMMCKRCGTSAGYGGKCGCTRYARRMRARAAAYGGYAATMRYRASSGKMYTMRYGTMSRTGSYRKGCCCRSYSTGCSTCGGATGTYSRTGGGGCNNNTGATTSCAYTRCYYRYMSYMYWSSRKKRTMKSRCASCYKTYMMSMKNKCCCTRCCGNSCGGTCKMAYCARYKYCSSYYAYSYRCNCKYCSCKKYTKYACAAAMYYKMCATCARSKSRATTGTYRCTRSRSKSCASSSCGTGYGGTGGSKMRSGYAGGCCTACARNMYYWSRYCTY
BM77780    CCGRGCRMRCCNTMARSYSCRNNNNGKGAMGRCGYGCTMYRACRRGACACTYMKGCCTAANNNNNNGCGCNYKMGSYRYTMTNNNNNNNTCTSCCCGAGCSRYYTKKTACRAAGNCCKGAYASMCRMAACSRGKMACSTGTGMGTKASCAKSGMARACTACSCRAMTRMTTGCGKATKCTCACTCGAAKSAAAYSKMRSKNRSRTKCGTTCWAASRYKWWGMCTCCTCWGGTTSGTYMGTMWGRRCYYMKMRKCCTTAGYYGTACGYSSSRMSMSWMYCGRYKSRRRSSTSKKGACKMGRYTKYCASRKTCGRCRYGTAKGGGATTCTKCCAACKRCGYSMGYGGGNNMTGTWARMRARRAAYGGYAATMACRAGSGKMYKMRYRWASRTGCYRKGCCCRSTSTKSCTSSGATGTYSRTGGSKYAKRTGWWWCCACTRCYYRCMSYMTACCAGTRTATSRCACMYKTCCASAKKGNCCTRMCRGSSGGTCTMACCANNGYCGGCTAYSYRCGCGCCCGGGTTKCACAAAMCCKCCRTCARSKSRCTTSTTRCTRSRSKSMASSSCGTGNGRTGGGKCRCGCRSSTYTRCMGYMYYTNGTCTN
BM77781    SYKRRSRMACCYYCARSYCCNTWRTGKKWMARYNCGCYMYRASRRSACRMTYMKGCCYAWYMRNTGSCGMSYKMGSYRYYMKRWKYYAMTMYSCMMKARSCGCCTKGTRYRRRSRCYTGACASMCRMMMCSRGKMAYSKRTKMGTKASCATCGAARACNASSMRAMYGCYTGCGKATTMTYACTMSWAKCAAAYSKMRSKNRSRTKYGWYMAARSRYTWASMCYMMTYTGRTTCGTTAGTMWGARCYYMKMRKCYWYMGYYRKRSRYSSSRMSMSWMYMRRYKSRRRSSYSKKRRYGMSRYTKYSASRKKYRRMRYSTWKGGSMTKMKKSYMMCKRCSTSMGYRSGCGMTRTARRMRARAAAYGGYAATMAYRRSCGKMYKMRYGTMSRTGCYRKGCCCRSTSTGCSTCSGATGYYSRYGGGGCAGRYRAWWSSACYRSYYRYMGYMYACSAKKATCKSRMASMTKTYCASMKGGNYSWRMMRGSSGGTMKMWCNRGTGCYGGCNNYSCRYRYKYYSSKKYTKYRYAAWMYCKCCATCRASKSRMTTSTYRYKRSRSKSMRSSCSSTGYGRWRNSKCRSSYRSGYCYACARYCYYWSRYSTY
BM77856    SYKGRSRCANNNTCRNGTCCATTRYRKGWARNCGYRSNNNNRCGGGACAMTYMGRYSYAACCAWYSNNNNNYGAGCNNNNNNRWGCTRCKCYSYCCGAGCSGCCKKGTRCGARSASCTTRCGCMYGAANMCARGCRYGKGYGAGTTMNNNKSGMRAMYTMSSMRRMYGCYYSYGKRKTCYCRCCCGARGCMMRTGGCGSKRRCGTKCSWCMARRSRCKAWGACYCCTNNNNNNNNNYMSTCAGGNCTTMKMAKCYTTAKYCGKANRTSSSGMCMSAMYMGRYKSARRSSTCTTGACKAGGYYKCSRGGKTYGRCRYSTAGGGGATTCKTSCAMMKRSGTGNRNGGNMKMYGYWARMARGRARYRSYRMYMAYRRGSRKMYKAAYRWMGGTGSYGKGCCCRGTGTKSSKSGGATGTTNGTGGSGYWGAYRAWTCCANTRCNYGYMSTMYNCCRKTAYMKSNMMSMYGYYMMSMKNKCYSWAMMNGNGGGTCKMAYCARYKYCSGCYMNCCAYRCKYCSSKKYTKYRYARWMYYKMCRTCRRSKSRCTTSTNACTNNNNNNNNNNCSGYGTGATANSTNRSGYASGYCYGYMGCCCCTCGNCTC
BM77978    CCGGGCGMANNNTNRRSTSCANTRYRKGWARACNNGCNNNRACGRSACAMTYMGGYCTAACCATTGSCGCGNNNNNNNNNNNATKYTAMTCYSCNNNNNNGRCCKKGYRYRRRSASNTTRCRSMYRMMNCSRGKMRYGNNNNAGTKASCAKSGARRMCGMSSCARCTRCTYSYGKAKKCTYANYMSWRKSMMRYSKMRGTRGSGTKCCACAARAGGCKAWSMCTCCWYWSATYCGTYMCKCAKGGYTTMNNNTNTTYMKYTRNRSGYSSSRMCMSWNTCRAYTCAAANNNNNNNNNNNNGCTGCSRGRKNNNACGTGTATGSSMYKMTTSYMMMKRSSTSAGTRSKCGCTGTWRRMRARARAYGGTAATMRCRAGSGKAYTCGTGTCSGYGSCGGKSSSRGYSTGSSTCGGATGYTGGTGGGGCAGATGWTTCCACTRCYTNYMSTMTACNAGKATMNSRCACNYKYCCASMKKGCCCTGCCGNSGGSYMTANYSARYKYCSSYCNCNCAYRCGYCCGGGYTKYNCRAAAYYKMSAYCRASKSRMTTSTYGCTGGGCGSCAGGNNNNNNKRWRGGKCRCGYASGTCTNNNNNNNNNCGTCKN


mkdir raxml
cd raxml
echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-28' > partitionfile.txt
~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n hydin -q /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/hydin/raxml/partitionfile.txt -s/scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/hydin/variantsites.phy -T 6 -p 12345 -N 10 ­-b 12345 -V

cd ..
mkdir pca

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/hydin/hydin.vcf.gz -o /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/hydin/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 






# phylogenetics on just Dlg2
awk '$0 ~ /\"Dlg2\"/ {print $1,$3,$4,$5}' genes.gtf
VYXE01009321.1 gene 525674 1108545


awk '$0 ~ /\#/' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > dlg2.vcf

awk '{$2 > 525674 && $2 < 1108545}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01009321/' >> dlg2.vcf

bgzip dlg2.vcf

bcftools index dlg2.vcf.gz

~/vcf2phylip/vcf2phylip.py -i dlg2.vcf.gz

python ~/sula/filter_invariants_all.py dlg2.min4.phy

24 49
MSB25201   TAWARTYSSATTTCKTYWSGGGKMRYSTYYYRGGASTKSMCGTAGRTAA
NOCA003    KRTMAKTSCMYKGSTYTACRRKTAGTCWCTCASCMCKKSASSYRSNNNN
NOCA004    KRTMAKTSCMTKKSTYTACRRKTAGTCWCTCASSMCKKSASSYRSNNGG
NOCA006    KRTMAKTSCMTKKSTYTACRRKTAGTCWCTCASSMCKKSASSYRSRYRR
NOCA008    KRTMAKTSCMYKGSTYTACRRKKAGTCWCTCASSMCKKSASSYRSRYAA
NOCA012    KRTMAKTSCMYKKSTYTACRRKTAGTCWCTCASSMCKKSASSYRSRCRR
NOCA013    KRTMAKTSCMYKKSTYTACRRKTAGTCWCTCASSMCKKSASSYRSGYRR
PYRR003    TAWARTYSSATTTCKTYWSGGGKMRYSTYYCASGASTTSMCGTAGNNNN
PYRR004    TAWARTYSSATTTCKTYWSGGGKMRYSTYYYRGGASTKSMCGTAGGYGG
PYRR006    TAWARTYSSATTTCKTYWCGGGKMRYSTYYYRGGASTKSMCGTAGGTGG
PYRR007    TAWARTTSSATTTCKTYWSGGGTMRYSTYYCASGASTKSMCGTAGRCGG
PYRR009    TAWARTYSSATTTCKTYWSGGGKMRYSTYYYRSGASTKSMCGTAGRYRR
PYRR011    TAWARTTSSATTTCKTYWCGGGKMRYSTYYYRSGASTKSMCGTAGRYRR
BM100619   KRTMAKTSCMYKGSTYTACRRKTAGTCWCTCASSMCKKSASSYRSGCGG
BM100620   KRTMAKTSCMYKGSTYTACRRKTAGTCWCTCASSMCKKSASSYRGGCGG
BM100621   KRTMAKTSCMTKGSTYTACRRKTAGTCWCTCASSMCKKSASSYRSRYAA
BM103345   KRTMAKTSCMTGKSTYTACRRKTAGTCWCTCASSMCKKSASSYRSACAA
BM103346   TAWARTYSSATTTCKTYWSGGGKMRYSTYYCASGASTKSMCGTAGRYRR
BM77548    TAWARTTSSATTTCKTYWCGGGKMRYSTYYCRSGASTKSMCGTAGNYNR
BM77718    TAWARTYSSATTTCKTYWSRGGTMRYSTYYYRGGASTKSMCGTAGRYRR
BM77780    TAWARTYSSATTTCKTYWCRGGTMRYSTYYCASGASTTSMCGTAGRYRR
BM77781    TAWARTYSSATTTCKTYWSGGGTMRYSTYYCRSGASTKSMCGTAGRTRG
BM77856    KRTMAKTSCMTKGSTYTACRRKTAGTCWCTCASSMCKKSASSYRSAYRR
BM77978    KRTMAKTSCMTGGSTYTACRRKTAGTCWCTCASSMCKKSASSYRGRCRR

mkdir raxml
cd raxml
echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-4' > partitionfile.txt
~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n dlg2 -q /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/dlg2/raxml/partitionfile.txt -s/scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/dlg2/variantsites.phy -T 6 -p 12345 -N 10 ­-b 12345 -V

cd ..
mkdir pca

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/dlg2/dlg2.vcf.gz -o /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/dlg2/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 





# phylogenetics on just Fxr1_0
awk '$0 ~ /\"Fxr1_0\"/ {print $1,$3,$4,$5}' genes.gtf
VYXE01008382.1 gene 732481 733607

awk '$0 ~ /\#/' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > fxr1.vcf

awk '{$2 > 732481 && $2 < 733607}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01008382/' >> fxr1.vcf

bgzip fxr1.vcf

bcftools index fxr1.vcf.gz

~/vcf2phylip/vcf2phylip.py -i fxr1.vcf.gz

python ~/sula/filter_invariants_all.py fxr1.min4.phy

cat fxr1.min4.phy
24 6
MSB25201   AGAAGG
NOCA003    RGNNNN
NOCA004    AGRRGA
NOCA006    RRAGRA
NOCA008    RRRGAR
NOCA012    NNNNNN
NOCA013    NRAGRA
PYRR003    NNNNNN
PYRR004    RGNNNN
PYRR006    NNNNNN
PYRR007    RGRGGR
PYRR009    NNNNNN
PYRR011    AGARGR
BM100619   AGRRRA
BM100620   AGRRRA
BM100621   AAAGAA
BM103345   RRAGAA
BM103346   AGAAGG
BM77548    AGAAGG
BM77718    AGAAGG
BM77780    AGAAGG
BM77781    AGAAGG
BM77856    AAAGAA
BM77978    AAAGAA

mkdir raxml
cd raxml
echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-1' > partitionfile.txt
~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n fxr1 -q /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/fxr1/raxml/partitionfile.txt -s/scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/fxr1/variantsites.phy -T 6 -p 12345 -N 10 ­-b 12345 -V

cd ..
mkdir pca

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/fxr1/fxr1.vcf.gz -o /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/fxr1/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 





# phylogenetics on just Rexo1
awk '$0 ~ /\"Rexo1\"/ {print $1,$3,$4,$5}' genes.gtf
VYXE01016584.1 gene 2987 31838

mkdir Rexo1
cd Rexo1
awk '$0 ~ /\#/' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > Rexo1.vcf

awk '{$2 > 732481 && $2 < 733607}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01008382/' >> Rexo1.vcf

bgzip Rexo1.vcf

bcftools index Rexo1.vcf.gz

~/vcf2phylip/vcf2phylip.py -i Rexo1.vcf.gz

python ~/sula/filter_invariants_all.py Rexo1.min4.phy

cat Rexo1.min4.phy
24 6
MSB25201   AGAAGG
NOCA003    RGNNNN
NOCA004    AGRRGA
NOCA006    RRAGRA
NOCA008    RRRGAR
NOCA012    NNNNNN
NOCA013    NRAGRA
PYRR003    NNNNNN
PYRR004    RGNNNN
PYRR006    NNNNNN
PYRR007    RGRGGR
PYRR009    NNNNNN
PYRR011    AGARGR
BM100619   AGRRRA
BM100620   AGRRRA
BM100621   AAAGAA
BM103345   RRAGAA
BM103346   AGAAGG
BM77548    AGAAGG
BM77718    AGAAGG
BM77780    AGAAGG
BM77781    AGAAGG
BM77856    AAAGAA
BM77978    AAAGAA


mkdir pca

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Rexo1/Rexo1.vcf.gz -o /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Rexo1/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 








# phylogenetics on just Ctnna3_0
awk '$0 ~ /\"Ctnna3_0\"/ {print $1,$3,$4,$5}' genes.gtf
VYXE01009024.1 gene 127289 302091

mkdir Ctnna3
cd Ctnna3

awk '$0 ~ /\#/' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > Ctnna3_0.vcf

awk '{$2 > 127289 && $2 < 302091}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01009024/' >> Ctnna3_0.vcf

bgzip Ctnna3_0.vcf

bcftools index Ctnna3_0.vcf.gz

~/vcf2phylip/vcf2phylip.py -i Ctnna3_0.vcf.gz

python ~/sula/filter_invariants_all.py Ctnna3_0.min4.phy

cat Ctnna3_0.min4.phy
24 64
MSB25201   GGCGGYSGAYRGGYYRGARRRYKAYRGRYRRRKRYYSARASYRRYAYTAGTTCAACCAATCCGC
NOCA003    SNNNNTGCRYGRRYYRGWRRRYKAYRRRYGRRTRYYSARASYGRYAYTTANNNNNTTAAWYCNN
NOCA004    SRCSGYSSRYGGGYYRGWRRRCTRCNRRYGRRKAYYSRRRSYGRTAYTTACYTARNNNNNNANN
NOCA006    GNSSSTGCRYGRGYYRGWRRRCTRCNRGCGRRTRYYSRRRGCGRTAYTTANNTCGTTAANNNNN
NOCA008    GRNNNYSSRYGGGYYRGWRRRCTRCARGCGRRTAYYSRRRSCGRYMYTTACYTMGTTAANNATY
NOCA012    GASSSNNNRYGRRYYRGWRRRCTRYRRRCGRRTRYYSRGRSYGRYAYTTACYTCGTTARWYNNN
NOCA013    CASCGCCSRYGGGYYRGARRRYKACARRCGRRTRYYSARRSCGRYAYYNNNNTCGTTARWYMKC
PYRR003    GASCSYSGAYRGRYYRGARRRYKACRRRCRRRKRYYSARAGYRRYAYTAGTTNNNNNRAWYCKC
PYRR004    SRCCSYSSAYRGRYYRGARRRYKACAGGCRRRKRYYSARASYRRYAYTWRYTNNNNNNNNNCGC
PYRR006    NNSCSTGCAYRGGYCRGARRRYKACAGGCRRAKRYCSARASYRRYAYTWRYTYARYYRATCNNN
PYRR007    SNSCGYSSAYRGGYCRRARRRYKACARGCRRRKRYYSARASYRRYAYTNNNNNNNNNNNNNNNN
PYRR009    GACGGTGSAYRGGYYRGARRRYKACRGRYRRRKRYCSARASYRRYAYTNNNNYARNNNNNNCGC
PYRR011    SRCCGYSGAYRGGYYRGARRRYKAYRRRCRRRKRYYSAGAGYRRYAYTWRNTCAACCGAWYCKC
BM100619   SASSSNNNRYGGRYYRGWRRRCKRCNRRCGRRTRYYSRRRSYGRTAYTTACCTCGTTAAATATT
BM100620   GNNNNYSSRYGRRYYRGWRRRYKRYNRRCGRRTRYYSARRSYGRYMYTTACCTCGTTAAATATT
BM100621   SRSCGCCSRYGGGYYRGARRRCTACARGCGRRTRYYSRGASYRRYMYCTACCTCGTTAAATATT
BM103345   SASCGTGSRYGGRYYRRWRRRCTRYNRRYGRRTRYYSRRRSYGRTAYTTACCTCGTTAGATATT
BM103346   SNCSSTGCAYRGGYYRRARRRYKAYRGRYRRRKRYCSARASYRRYAYTAGTTCAACCGATCCGC
BM77548    SANNNTGCAYRGRYYRGARRRYKAYRRRYRRRKRYYSARAGYRRYAYNAGTTCAACCGATCCGC
BM77718    SACSGTGSAYRGGYCRGARRRYKACARRYRRRKRYYSARASYRRYAYYAGTTCAACCGATCCGC
BM77780    SRSCGYSSAYRGGYCRRARRRYKAYRGRYGRRKRYCSARASYRRYAYTAGTTCAACCGGTCCGC
BM77781    SGSSGYSSAYRGGYCAGARRRYKACAGRCRRRKRYCSARASCGRTAYTAGTTCAACCGATCCGC
BM77856    GRGCGYSSRYGRGYYRRWRRRCKACRRRYGRRTRYYSARASYGRYMYTTACCTCGTTAAATATT
BM77978    SRSSSNNCRYGGRYYRGWRRRCTRCNRRCGRRTRYYSARRSYGRTMYYTACCTCGTTAAATATT

mkdir raxml
cd raxml
echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-10' > partitionfile.txt
~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n Ctnna3_0 -q /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Ctnna3/raxml/partitionfile.txt -s/scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Ctnna3/variantsites.phy -p 12345 -N 10 ­-b 12345 -V

cd ..
mkdir pca

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Ctnna3/Ctnna3_0.vcf.gz -o /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Ctnna3/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 





# phylogenetics on just Dcbld2
awk '$0 ~ /\"Dcbld2\"/ {print $1,$3,$4,$5}' genes.gtf
VYXE01002912.1 gene 17908 40432

mkdir Dcbld2
cd Dcbld2

awk '$0 ~ /\#/' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > Dcbld2.vcf

awk '{$2 > 17908 && $2 < 40432}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01002912/' >> Dcbld2.vcf

bgzip Dcbld2.vcf

bcftools index Dcbld2.vcf.gz

~/vcf2phylip/vcf2phylip.py -i Dcbld2.vcf.gz

python ~/sula/filter_invariants_all.py Dcbld2.min4.phy

cat Dcbld2.min4.phy

24 1
MSB25201   G
NOCA003    N
NOCA004    G
NOCA006    G
NOCA008    G
NOCA012    G
NOCA013    G
PYRR003    R
PYRR004    R
PYRR006    N
PYRR007    G
PYRR009    R
PYRR011    N
BM100619   G
BM100620   R
BM100621   G
BM103345   G
BM103346   G
BM77548    G
BM77718    G
BM77780    R
BM77781    G
BM77856    G
BM77978    G


mkdir pca

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Dcbld2/Dcbld2.vcf.gz -o /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Dcbld2/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 





# phylogenetics on just Col6a1
awk '$0 ~ /\"Col6a1\"/ {print $1,$3,$4,$5}' genes.gtf
VYXE01020886.1 gene 140268 142418

mkdir Col6a1
cd Col6a1

awk '$0 ~ /\#/' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > Col6a1.vcf

awk '{$2 > 140268 && $2 < 142418}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01020886/' >> Col6a1.vcf

bgzip Col6a1.vcf

bcftools index Col6a1.vcf.gz

~/vcf2phylip/vcf2phylip.py -i Col6a1.vcf.gz

python ~/sula/filter_invariants_all.py Col6a1.min4.phy

cat Col6a1.min4.phy

24 134
MSB25201   GYSSKCCGGGGCCCGTCYTKWCCSRKYYRRMTSMRYRRYRRYASWSSWKTGRARRCCCYRNASYYCYRCASRRRGSAYCGYYACCRYMYSKCSTRMSRSWRWWRCRSSSTCCGCCCTGTGTTAAATCTTGACCG
NOCA003    GYGSTSSRRSSSSYSTCCYGTYCSRGCTAGCTSNGTRAYRRYASASSWTTGRARRCCSYRSRSYYCYRCWGRRRSSATNNYYNCCRYMYSKYSKRASRSWGWWRCRSGGKTTATTACACAACCGRCTCAAMACN
NOCA004    KCSSKCCGGSSNCYGYCCTKTCCSAGCTAGCYSCGTRAYAATAGASSWTTGRARRCSSYRSRSYYYYRSWSRRRSSMCCGCCRYYRCMYSTCSKRMSGSWGAWRYGSGGGTTATTACACAACCRGCTCAACMCG
NOCA006    KYGSTCSRRSSSSYSTCCYKTCCSAKYYRRCYSMRYRAYAAYASASCATTGRARRSSCYRSRSYYYYRSWSRRNSSACCGYYNYYRCAYSTCSTRASRSWGAWRYGSGSKNNATTACACAACCGGCTCAACAMA
NOCA008    KCSGTSSRRSSCSYSTCCYGTYCSAGCTAGCYSMRYRRCRRCRSAGSWTTGRARRSSSYRSRSCYYCRSWSARRGSACCRCCACCRYMYSKYSKRMSRGWRWWRCGSGSKTTATTACACAACCGGCTCAACAMR
NOCA012    KYSSKCSGRSSSCYSTYCTKTCCSAGCTAGCTGCGTRRYRAYRGASCATTGGAGRCCSYRNASYYCCRCWGRRRGSMCCRYYNCCRCMYSTCSKRASRSWRWWRYGSSSTTTATTACACAACCGGCTCAACAAN
NOCA013    KYSSKCCGGSSCCYSTCCTKTCSSAKCTAGMYSCGYGRYRRYASWSCAKTGGAGRCCSYRSAGCYCYRCWSRRRGSAYNGCCACYRCMYSTCSKRASGSWGAWRYGSGSKTTATTACACAACCGGCTCNNNNNG
PYRR003    KYSSKCCGGSSNCCGYCCTKTCCSAGYYRGCTSCGYRAYARYASWGSWTTGRRARCCSYRSRGCYYYRCAGRRRGSMYNGYYACCRYMYSTCSTGAGRSWRWWRCRSGSGYYRYYNYRYRWYMRRNNNNNNNNG
PYRR004    KYSSKCCGGGGCCYGTYYTGTCSSAGYYRRMTSMGYRRYRRYAGAGSWTTGRARRSCSYRGRGCYYYRSASRARSSACCRCCACCRCMYSTCSKRASRSWRWWRCRSSSKNNNNNNNNNNNNNNNNNNNACMCG
PYRR006    GYSSKCCGGGSCCYGTYYTKWCSSAGCTAGCTSCGTGACAAYRSASSWKTGRARRSSCYRSRGCYYCRSASRRRSSAYCRCCACYRCMYSTCSKRASRSWRWAGCRSSNNNNNNNNNNNNNNNNNNNNNNNNNG
PYRR007    KCSSKCCGGGSSCYGYCCTKWCCSAKYYRRMYCMRYRRCANYAGAGSWTTGRAARCCSYRNRSYYCCRCASARRGGAYCGCCACCRCMYSTCCKRASRSWRWWRYGSGGGNNNNNNNNNNNNNNNNNNNNNNNG
PYRR009    KYSSKCCGGSSCCYGTYYTKTCSCRKCTAGMYCCRTRACARYASASCAKTGRARRSSCYRSRGYYYCRSASARNSSMCCGCCAYCRCMYSTCSKRASRSWRWWRCRSSSKNNNNNNNNNNNNNNNNNNNNNNNR
PYRR011    GYSSKCCGGSSCCCGTCCTGWCCCRGYYRRCTCNRYGRCRRYAGASSAKKSRRARCCSYRSAGCYYYRSASRRRGSMYCGYYRCCRCMYSTYSTRAGRSWRWWRYRSGGTYNNNNNNNNRWYMRRYYYWRMNNR
BM100619   KYSSTSSRRSSSSYSTYCYKTYCSRKYYRRCYSMRYRAYRRYASAGSATTGRARRSSCYRSRSYYYCGSASRARSSACCGYYRCCRCMYSKYSKRASGSWRWWRCGSSSTTTATTACACAACCGGCTCAACAAR
BM100620   KYSSKCCGRSSSCYSYCYTKWYCSRKYYRRMTSMRYGRYRRYRSASCATTGGAGRCCCYRSASYYCCRSWSRRRGGACCGCCACYRCACSTCSKRASGSWRWWRCGSSGGTTATTACACAACCGGCTCAACAAG
BM100621   KCGSTSSRRSSSSYSTCYYGTYCSANNTAGMNCCGYGRYRACASWSSWTTGRARRCSCYRSRSYYCCRCASRANGCANCGNNACCRYMYSKCSKGAGGSWRWWRCGSSGGTTATTACACAACCGGCTCAACAAR
BM103345   GYGSTCSRRSSSSYSTCCYKTYCSRGCTAGCYSCGTRAYARCRCWSCAKTGGARRSSCYRSRGCYYYGSASRANSSMCMRYYRYYRCMYSTCSKRMSRSWRWWRCGSGGKTTATTACACAACCGGCTCAACAAN
BM103346   KTSSKSSGRSSCCYGTYYTKWCCSAGCTAGCYCCGTRACAAYAGASCATKSGRGRCCSYRNAGCYYCRCAGRRGGSACNGCCRYYRYMYSKYSKRMSGSWRWWRCRSGNGCCGCCCTGTGTTAAATCTTGACCG
BM77548    GYSSKCCGRSSCSYGTCCTGWCCCAGYYRRCTSCRYGRYRRYRSWSSWTKSRRRRCCSYRNAGYYCYRCAGARRGSAYNGCCRCCRCMYSTCSTRASRSWRWWRCRSSSKCCGCCCTGTGTTAAATCTTGACCG
BM77718    GYSSKCCGGGSCCYGTYYTKTCCCRGCTAGMTSNRYGRYRRCANNGSWTKSRRRRSSCYRGRGCYYCRSASRRRSSANCGYYRCCRYMYSTYSTRAGRSWRWWRCRSSSKCCGCCCTGTGTTAAATCTTGACCR
BM77780    KYSSKCCGGGGCCCGTYYTKTCSCRKCTAGMTCNRCRRCRRCRSWNNNNNNNNNRCCCYRNAGYYYYRSASRANSCACMGYYACCRCMYSKYSTRASRSWRWWRYRSSSTCCGCCCTGTGTTAAATCTTGACCG
BM77781    KYSSKCSGRSSNCYGYCCTKWCCSRKYYRRMYSMRYRRYARYASWSSWTKSRRRRSSCYRGRSCYYCRCASRRRGSMYMRYYACCRCMYSTCSKRMSRSWRWWRCRSSSKCCGCCCTGTGTTAAATCTTGACCG
BM77856    KCSSKCCGRSSSCCGTCCTKTCSSRGCTAGCTGCGTGAYARCASWGCATTGGARRCCSYRSASYTCYGCWSRAGGCMYNGYYNCYRCMYSTCSKRMSRGTRWWGCGSSGGTTATTACACAACCGGCTCAACAAN
BM77978    GYSSKSCGRSSSCYSTYCYKTYCSRGYYRRCTSMRYGAYAAYAGNGCATTGGAGRSCCYRSRGCYYCGSASRAGGSMCMGYYNYYRYACSKCSKRMSRSWGAWRCGGGNNTTATTACACAACCGGCTCAACAAN

mkdir raxml
cd raxml
echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-23' > partitionfile.txt
~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n Col6a1 -q /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Col6a1/raxml/partitionfile.txt -s/scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Col6a1/variantsites.phy -p 12345 -N 10 ­-b 12345 -V

cd ..
mkdir pca

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Col6a1/Col6a1.vcf.gz -o /scratch/dnjacks4/cardinalis/to_b10k/genesofinterest/Col6a1/pca/ -p /scratch/dnjacks4/cardinalis/to_b10k/PCA/all/pops.txt -n all 






# Structure
# 1 = NOCA urban, 2 = NOCA rural, 3 = PYRR urban, 4 = PYRR rural
plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf --recode structure --out b10k_filtered.geno25.maf1 --allow-extra-chr 

/scratch/dnjacks4/cardinalis/to_b10k/structure/b10k_filtered.geno25.maf1.recode.strct_in


head -1 /scratch/dnjacks4/cardinalis/to_b10k/structure/b10k_filtered.geno25.maf1.recode.strct_in  > /scratch/dnjacks4/cardinalis/to_b10k/structure/b10k_filtered.geno25.maf1.recode.strct_in.nomap
tail -n +3 /scratch/dnjacks4/cardinalis/to_b10k/structure/b10k_filtered.geno25.maf1.recode.strct_in >> /scratch/dnjacks4/cardinalis/to_b10k/structure/b10k_filtered.geno25.maf1.recode.strct_in.nomap

awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' popfile.txt /scratch/dnjacks4/cardinalis/to_b10k/structure/b10k_filtered.geno25.maf1.recode.strct_in.nomap > /scratch/dnjacks4/cardinalis/to_b10k/structure/b10k_filtered.geno25.maf1.recode.strct_in.txt

~/programs/structure_kernel_src/structure 
~/programs/distruct1.1/distructLinux1.1

