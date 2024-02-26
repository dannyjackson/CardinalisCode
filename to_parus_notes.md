# interactive -t 15 -p debug -q wildfire
# interactive -n 7 -t 1-00:00
#   Long wait times? See if there are available resources in different partitions with the command, `showparts`!

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



zcat NOCA_004_CKDN230006109-1A_H5HWMDSX7_L2_1.fq.gz > NOCA_004_1.fq
zcat NOCA_004_CKDN230006109-1A_H5HWMDSX7_L2_2.fq.gz > NOCA_004_2.fq

zcat NOCA_004_CKDN230006109-1A_H5KJTDSX7_L2_1.fq.gz >> NOCA_004_1.fq
zcat NOCA_004_CKDN230006109-1A_H5KJTDSX7_L2_2.fq.gz >> NOCA_004_2.fq

zcat NOCA_004_CKDN230006109-1A_H5MVYDSX7_L1_1.fq.gz >> NOCA_004_1.fq
zcat NOCA_004_CKDN230006109-1A_H5MVYDSX7_L1_2.fq.gz >> NOCA_004_2.fq

bgzip NOCA_004_1.fq 
bgzip NOCA_004_2.fq 



zcat PYRR_003_CKDN230006114-1A_H5MT2DSX7_L1_1.fq.gz > PYRR_003_1.fq
zcat PYRR_003_CKDN230006114-1A_H5MT2DSX7_L1_2.fq.gz > PYRR_003_2.fq
zcat PYRR_003_CKDN230006114-1A_H5MVYDSX7_L1_1.fq.gz >> PYRR_003_1.fq
zcat PYRR_003_CKDN230006114-1A_H5MVYDSX7_L1_2.fq.gz >> PYRR_003_2.fq

bgzip PYRR_003_1.fq
bgzip PYRR_003_2.fq

zcat PYRR_007_CKDN230006117-1A_H5HNVDSX7_L1_1.fq.gz > PYRR_007_1.fq
zcat PYRR_007_CKDN230006117-1A_H5HNVDSX7_L1_2.fq.gz > PYRR_007_2.fq
zcat PYRR_007_CKDN230006117-1A_H5MT2DSX7_L1_1.fq.gz >> PYRR_007_1.fq
zcat PYRR_007_CKDN230006117-1A_H5MT2DSX7_L1_2.fq.gz >> PYRR_007_2.fq

bgzip PYRR_007_1.fq
bgzip PYRR_007_2.fq

# gotta redo these 

zcat UWBM_103346_CKDN230006128-1A_H5MTCDSX7_L2_1.fq.gz > UWBM_103346_1.fq
zcat UWBM_103346_CKDN230006128-1A_H5MTCDSX7_L2_2.fq.gz > UWBM_103346_2.fq
zcat UWBM_103346_CKDN230006128-1A_H5MTGDSX7_L3_1.fq.gz >> UWBM_103346_1.fq
zcat UWBM_103346_CKDN230006128-1A_H5MTGDSX7_L3_2.fq.gz >> UWBM_103346_2.fq

bgzip UWBM_103346_1.fq
bgzip UWBM_103346_2.fq

zcat UWBM_77856_CKDN230006121-1A_H5LY7DSX7_L2_1.fq.gz > UWBM_77856_1.fq
zcat UWBM_77856_CKDN230006121-1A_H5LY7DSX7_L2_2.fq.gz > UWBM_77856_2.fq
zcat UWBM_77856_CKDN230006121-1A_H5MNLDSX7_L2_1.fq.gz >> UWBM_77856_1.fq
zcat UWBM_77856_CKDN230006121-1A_H5MNLDSX7_L2_2.fq.gz >> UWBM_77856_2.fq

bgzip UWBM_77856_1.fq
bgzip UWBM_77856_2.fq


NOCA_004_1.fq.gz
NOCA_004_2.fq.gz
PYRR_003_1.fq.gz
PYRR_003_2.fq.gz
PYRR_007_1.fq.gz
PYRR_007_2.fq.gz
UWBM_103346_1.fq.gz
UWBM_103346_2.fq.gz
UWBM_77856_1.fq.gz
UWBM_77856_2.fq.gz

sed -i 's/\_1\.fq\.gz//g' combined_filenames.txt
sed -i 's/\_2\.fq\.gz//g' combined_filenames.txt

./trim.sh -i /scratch/dnjacks4/cardinalis/combined_filenames.txt -p /scratch/dnjacks4/cardinalis/genomicdata/raw_fastas/ -f 1.fq.gz -r 2.fq.gz -t 7

java -jar /packages/7x/trimmomatic/0.33/trimmomatic.jar PE   /scratch/dnjacks4/cardinalis/genomicdata/raw_fastas/UWBM_77856_1.fq.gz /scratch/dnjacks4/cardinalis/genomicdata/raw_fastas/UWBM_77856_2.fq.gz \
    -baseout /scratch/dnjacks4/cardinalis/genomicdata/raw_fastas/UWBM_77856_trimmed.fq.gz \
    ILLUMINACLIP:adapters.txt:1:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>trim_and_QC_log.txt


bwa index /scratch/dnjacks4/cardinalis/referencedata/ncbi_dataset_parus/data/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna

cd /scratch/dnjacks4/cardinalis/to_parus

/scratch/dnjacks4/cardinalis/sort.sh -i /scratch/dnjacks4/cardinalis/sampleids.txt -p /scratch/dnjacks4/cardinalis/genomicdata/trimmed_fastas_select/ -r /scratch/dnjacks4/cardinalis/referencedata/ncbi_dataset_parus/data/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna

/scratch/dnjacks4/cardinalis/to_parus/sort_postbam.sh -i /scratch/dnjacks4/cardinalis/sampleids.txt -p /scratch/dnjacks4/cardinalis/genomicdata/trimmed_fastas_select/ -r /scratch/dnjacks4/cardinalis/referencedata/ncbi_dataset_parus/data/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna


ref="/scratch/dnjacks4/cardinalis/referencedata/ncbi_dataset_parus/data/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna"
bamdir="/scratch/dnjacks4/cardinalis/sorted_bam_files/"
ID="parus"
bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*_sorted_RGadded_dupmarked.bam | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf


bcftools stats -v /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_snps_multiallelic.vcf > parus_snps_multiallelic.stats.txt

bcftools stats /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf > parusfiltered.geno25.maf1.stats.txt

plink --vcf /scratch/dnjacks4/cardinalis/to_parus/parus_snps_multiallelic.vcf --allow-extra-chr --missing --cluster-missing --freq

plink --vcf /scratch/dnjacks4/cardinalis/to_parus/parus_snps_multiallelic.vcf --allow-extra-chr --missing --cluster missing --within individualnames.txt --freq

plink --vcf /scratch/dnjacks4/cardinalis/to_parus/parus_snps_multiallelic.vcf --allow-extra-chr --missing --cluster missing --within cluster_pop.txt --freq


#filters by quality
bcftools view -i 'QUAL>100' parus_snps_multiallelic.vcf > parus_qualitysort.vcf

vcftools --vcf parus_qualitysort.vcf --min-meanDP 2 --remove-indels --recode --out parus_filtered
# with a max-meanDP 8:
# After filtering, kept 16120 out of a possible 3246716 Sites
# After filtering, kept 20825 out of a possible 3246716 Sites

sed -i 's/\_//g' parus_filtered.recode.vcf
sed -i 's/\UWBM/BM/g' parus_filtered.recode.vcf

plink --vcf parus_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02  --recode vcf-iid --out parus_filtered_geno02

plink --vcf parus_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.2  --recode vcf-iid --out parus_filtered_geno2

plink --vcf parus_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.15  --recode vcf-iid --out parus_filtered_geno15

plink --vcf parus_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.1  --recode vcf-iid --out parus_filtered_geno1


plink --vcf parus_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid --out parus_filtered_mind2

plink --vcf parus_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.1 --maf 0.01 --recode vcf-iid --out parus_filtered_cleaned


cp parus_filtered_cleaned.vcf parus_filtered_cleaned_zip.vcf
bgzip parus_filtered_cleaned_zip.vcf

bcftools index parus_filtered_cleaned_zip.vcf.gz

~/vcf2phylip/vcf2phylip.py -i /scratch/dnjacks4/cardinalis/to_parus/parus_snps_multiallelic.vcf

mkdir pruned 

cd pruned 

python ~/sula/filter_invariants_all.py ../parus_filtered_mind2.min4.phy
mv variantsites.phy ../variantsites_mind2.phy 
mv variantsites_kept.txt ../variantsites_mind2_kept.txt 
cd .. 
rm -r pruned


echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-93' > partitionfile.txt
~/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n cardinalis_1 -q /scratch/dnjacks4/cardinalis/to_parus/partitionfile.txt -s /scratch/dnjacks4/cardinalis/to_parus/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V




# Filtering for phylonet

# on Sol
module purge
module load bcftools-1.14-gcc-11.2.0 
module load htslib-1.16-gcc-11.2.0
module load samtools-1.9-gcc-12.1.0
module load bedtools2-2.30.0-gcc-11.2.0

bcftools view -m2 -M2 -v snps /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_cleaned.vcf -o /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf

bgzip /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf

bcftools index /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf.gz 

zcat /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf.gz | cut -f 1 | grep -v '#' > /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames.txt

zcat /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf.gz | grep '##contig' > /scratch/dnjacks4/cardinalis/to_parus/vcfs/contiglengths.txt 

awk -F '[=,]' '{if(substr($5, 1, length($5)-1)*1 > 5000)print $3} ' /scratch/dnjacks4/cardinalis/to_parus/vcfs/contiglengths.txt  > /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt


parallel -j 30 -a /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt 'C={}; bcftools view -O z -r ${C} -o /scratch/dnjacks4/cardinalis/to_parus/vcfs/split_vcf_parallel/split.${C}.vcf.gz /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf.gz' 

parallel -j 30 -a /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt 'C={}; bcftools index /scratch/dnjacks4/cardinalis/to_parus/vcfs/split_vcf_parallel/split.${C}.vcf.gz' 


### there are a couple of i think unnecessary lines of code here


mkdir /scratch/dnjacks4/cardinalis/to_parus/reference_lists/

mkdir /scratch/dnjacks4/cardinalis/to_parus/reference_lists/myPIRsLists/

mkdir /scratch/dnjacks4/cardinalis/to_parus/haplotypes/

ls -d /scratch/dnjacks4/cardinalis/sorted_bam_files/*bam | awk 'BEGIN { FS = "/" } ; {print $7} ' | awk 'BEGIN { FS = "_" } ; {print $1, '\t', "/scratch/dnjacks4/cardinalis/sorted_bam_files/", $0 '\t', "replacethis"} ' | sed 's/files\/ /files\//g' > /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bamlist_reference.txt

 


while read -r line
do 
  sed 's/replacethis/'"$line"'/g' /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bamlist_reference.txt > /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bamlists/bamlist_"$line".txt 
done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt


while read -r line
do 
  ~/programs/extractPIRs.v1.r68.x86_64/extractPIRs --bam /scratch/dnjacks4/cardinalis/to_parus/reference_lists/bamlists/bamlist_"$line".txt  \
              --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/split_vcf_parallel/split."$line".vcf.gz \
              --out /scratch/dnjacks4/cardinalis/to_parus/reference_lists/myPIRsLists/myPIRsList."$line".txt

done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt

while read -r line
do 
  ~/programs/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -assemble \
        --input-vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/split_vcf_parallel/split."$line".vcf.gz \
        --input-pir /scratch/dnjacks4/cardinalis/to_parus/reference_lists/myPIRsLists/myPIRsList."$line".txt \
        -O /scratch/dnjacks4/cardinalis/to_parus/haplotypes/myHaplotypeData."$line".txt 
done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt

mkdir /scratch/dnjacks4/cardinalis/to_parus/phylonet/

mkdir /scratch/dnjacks4/cardinalis/to_parus/phylonet/stats


module load gsl-2.7.1-gcc-11.2.0
module load mamba/latest
mamba create -n genomics -c conda-forge scipy numpy cython=0.26.1 python=2.7

source activate genomics

mamba repoquery search plink
mamba install

while read -r line
do 
  # mkdir /scratch/dnjacks4/cardinalis/to_parus/phylonet/stats/"$line"

  cd /scratch/dnjacks4/cardinalis/to_parus/phylonet/stats/"$line"

  ~/programs/plink --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/split_vcf_parallel/split."$line".vcf.gz --missing --allow-extra-chr

done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt


## I skipped a lot of these ... come back to them later if needed

while read -r chrom
do 
  LIST=$(awk -F ' ' '{print $6;}' /scratch/dnjacks4/cardinalis/to_parus/phylonet/stats/"$chrom"/plink.imiss | sed 1d | tr '\n' ' ')

  PRESENT=$(cat /scratch/dnjacks4/cardinalis/to_parus/phylonet/stats/"$chrom"/plink.imiss | wc -l)

  MISSCOUNT=0
  if (( $(echo "$PRESENT < 5" | bc -l) )); then
        MISSCOUNT=$(($MISSCOUNT + 1))
  fi

  
  for item in ${LIST}; do
    if (( $(echo "$item > 0.07" | bc -l) )); then
        MISSCOUNT=$(($MISSCOUNT + 1))
    fi
  done

  if (( $(echo "$MISSCOUNT < 1" | bc -l) )); then 
    echo $chrom >> /scratch/dnjacks4/cardinalis/to_parus/reference_lists/goodchroms.07.txt 
    cat /scratch/dnjacks4/cardinalis/to_parus/reference_lists/goodchroms.07.txt | wc -l
  fi
done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt







while read -r line
do 
  ~/programs/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit  -convert \
        --input-haps /scratch/dnjacks4/cardinalis/to_parus/haplotypes/myHaplotypeData."$line".txt \
        --output-vcf /scratch/dnjacks4/cardinalis/to_parus/shapeit_vcf/myHaplotypeData."$line".vcf
done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt


while read -r line
do
  grep '##' /scratch/dnjacks4/cardinalis/to_parus/shapeit_vcf/myHaplotypeData."$line".vcf > header.tmp
  zcat /data5/sulidae/revisions/sula_flightless_filtered_biallelic.vcf.gz | grep $line | head -1 >> header.tmp
  grep -v '##' /data5/sulidae/revisions/shapeit_vcf/myHaplotypeData."$line".vcf >> header.tmp
  mv header.tmp /data5/sulidae/revisions/shapeit_vcf_new/myHaplotypeData."$line".vcf
done < /data5/sulidae/revisions/chromosomenames_5000min.txt

while read -r line
do
  bgzip /data5/sulidae/revisions/shapeit_vcf/myHaplotypeData."$line".vcf
done < /data5/sulidae/revisions/chromosomenames_5000min.txt

while read -r line
do
  bcftools index /data5/sulidae/revisions/shapeit_vcf/myHaplotypeData."$line".vcf.gz
done < /data5/sulidae/revisions/chromosomenames_5000min.txt




# picking it back up here
# i need to convert these into fasta files

MSB25201,NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,BM100619,BM100620,BM100621,BM103345,BM103346,BM77548,BM77718,BM77780,BM77781,BM77856,BM77978

echo 'MSB25201,NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,BM100619,BM100620,BM100621,BM103345,BM103346,BM77548,BM77718,BM77780,BM77781,BM77856,BM77978' >  /scratch/dnjacks4/cardinalis/to_parus/reference_lists/phylonet_samplenames.txt
sed -i 's/,/\n/g' /scratch/dnjacks4/cardinalis/to_parus/reference_lists/phylonet_samplenames.txt


sed 's/NC\_/NC/g' /scratch/dnjacks4/cardinalis/referencedata/ncbi_dataset_parus/data/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna > /scratch/dnjacks4/cardinalis/referencedata/parus_reference_nounderscore.fna

sed -i 's/NW\_/NW/g' /scratch/dnjacks4/cardinalis/referencedata/parus_reference_nounderscore.fna

cd /scratch/dnjacks4/cardinalis/referencedata/


samtools faidx /scratch/dnjacks4/cardinalis/referencedata/parus_reference_nounderscore.fna

cd /scratch/dnjacks4/cardinalis/to_parus/shapeit_vcf/

ls > ../shapeitvcffilenames.txt

while read file; do 
    bgzip $file
done < ../shapeitvcffilenames.txt

ls > ../shapeitvcffilenames.txt

while read file; do 
    bcftools index $file
done < ../shapeitvcffilenames.txt




while read -r file; do
  
   if [[ $(awk -F '\t' '{print $2}' /data5/sulidae/revisions/phylonet/5000kb_masked/"$file"  | sed 's/[^N]//g' | awk '{ print length }' | uniq -c | awk '{print $2}') -le 350 ]]; then
    echo $file >> scaffolds.07missing.txt
  fi
  
done < /data5/sulidae/revisions/phylonet/5000kb_masked_filenames.txt 


while read -r chrom; do
  while read -r bird; do

    if [[ $(zcat /scratch/dnjacks4/cardinalis/to_parus/shapeit_vcf/myHaplotypeData."$chrom".vcf.gz | wc -l) -ge 7 ]]; then

        echo '>' $bird >> /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/"$chrom".fa

        samtools faidx /scratch/dnjacks4/cardinalis/referencedata/parus_reference_nounderscore.fna $chrom| bcftools consensus /scratch/dnjacks4/cardinalis/to_parus/shapeit_vcf/myHaplotypeData."$chrom".vcf.gz -s $bird | sed 1d >> /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/"$chrom".fa
    fi

  done < /scratch/dnjacks4/cardinalis/to_parus/reference_lists/phylonet_samplenames.txt
done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt



while read -r chrom; do

    if [[ $(zcat /scratch/dnjacks4/cardinalis/to_parus/shapeit_vcf/myHaplotypeData."$chrom".vcf.gz | wc -l) -ge 7 ]]; then

        echo "$chrom" >> /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomes_fromhaplotypes.txt
    fi
done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomenames_5000min.txt



while read -r chrom; do

 bedtools maskfasta -fi /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/"$chrom".fa -bed /scratch/dnjacks4/cardinalis/to_parus/vcfs/split_vcf_parallel/split."$chrom".vcf.gz -fo /scratch/dnjacks4/cardinalis/to_parus/phylonet/masked/"$chrom".fa

done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomes_fromhaplotypes.txt



# and i need to split these fasta files into 6kb chunks

while read -r chrom; do
  echo "starting $chrom"
  cat /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/"$chrom".fa| awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}' > /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/temp.fa

  VAR1=$(zcat /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf.gz | grep $chrom | head -1 | awk -F'=' '{print $4;}' | sed 's/>//g')

  VARCOUNT=$"0"
  echo "looping $chrom"

  until [ $VAR1 -le 5000 ]; do 
    echo "round $VARCOUNT"

    awk -F'\t' 'BEGIN {OFS = '\t'} {ind=$1;} {fasta=substr($2,0,5000);} {print ind, "\t", fasta;}' /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/temp.fa > /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/5000kb/"$chrom"_"$VARCOUNT".fa 

    awk -F'\t' 'BEGIN {OFS = '\t'} {ind=$1;} {fasta=substr($2,15001);} {print ind, "\t", fasta;}' /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/temp.fa > /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/temp2.fa

    mv /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/temp2.fa /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/temp.fa

    VARCOUNT=$(($VARCOUNT + 1))
    VAR1=$(($VAR1 - 15000))
  done
  rm /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/temp.fa
done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/chromosomes_fromhaplotypes.txt


ls /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/5000kb/ | sed 's/.fa//g' > /scratch/dnjacks4/cardinalis/to_parus/vcfs/5000kbhaplotypes.txt

while read -r haplotype; do

  cat /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/5000kb/"$haplotype".fa | tr "\t" "\n" |  fold -w 60 > /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/5000kb/"$haplotype".folded.fa

  mkdir /scratch/dnjacks4/cardinalis/to_parus/raxml_trees/"$haplotype"/

  ~/programs/standard-RAxML/raxmlHPC -s /scratch/dnjacks4/cardinalis/to_parus/phylonet_fastas/5000kb/"$haplotype".folded.fa -m GTRCAT -p 12345 -n test -w /scratch/dnjacks4/cardinalis/to_parus/raxml_trees/"$haplotype"/

  cat /scratch/dnjacks4/cardinalis/to_parus/raxml_trees/"$haplotype"/RAxML_bestTree.test >> /scratch/dnjacks4/cardinalis/to_parus/raxml_trees//BestTrees.tre

done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/5000kbhaplotypes.txt


mkdir /scratch/dnjacks4/cardinalis]$ cd /scratch/dnjacks4/cardinalis/to_parus/phylonet/all/
## four reticulations 
echo -e "#NEXUS\n\nBEGIN TREES;\n" > /scratch/dnjacks4/cardinalis/to_parus/phylonet/all/script_4ret.nex


TREECOUNT="0"
while read -r haplotype; do

  if [[ $(cat /scratch/dnjacks4/cardinalis/to_parus/raxml_trees/"$haplotype"/RAxML_bestTree.test | wc -l) -ge 1 ]]; then 
    echo -n "Tree gt$TREECOUNT=" |   cat - /scratch/dnjacks4/cardinalis/to_parus/raxml_trees/"$haplotype"/RAxML_bestTree.test >> /scratch/dnjacks4/cardinalis/to_parus/phylonet/all/script_4ret.nex

    TREECOUNT=$(($TREECOUNT + 1))
  fi

done < /scratch/dnjacks4/cardinalis/to_parus/vcfs/5000kbhaplotypes.txt



echo -e "\n\nEND;\n\nBEGIN PHYLONET;\n\nInferNetwork_MPL (all) 1 -x 100 resultOutputFile phylonet_allspecies_4ret.txt -di phylonet_allspecies_4ret.tre;\n\nEND;" >> /scratch/dnjacks4/cardinalis/to_parus/phylonet/all/script_4ret.nex


java -jar ~/programs/Phylonetv3_8_2.jar /scratch/dnjacks4/cardinalis/to_parus/phylonet/all/script_4ret.nex > /scratch/dnjacks4/cardinalis/to_parus/phylonet/all/phylonet_allspecies_4ret.txt 


java -jar ~/programs/Phylonetv3_8_2.jar script_4ret.nex > phylonet_allspecies_4ret_cardinalis.txt 


echo -e "\n\nEND;\n\nBEGIN PHYLONET;\n\nInferNetwork_MPL (all) 1 -a <PYRR_urban:PYRR_003,PYRR_004,PYRR_006,PYRR_007,PYRR_009,PYRR_011;PYRR_rural:MSB_25201,UWBM_103346,UWBM_77548,UWBM_77718,UWBM_77780,UWBM_77781; NOCA_urban:NOCA_003,NOCA_004,NOCA_006,NOCA_008,NOCA_012,NOCA_013; NOCA_rural:UWBM_100619,UWBM_100620,UWBM_100621,UWBM_103345,UWBM_77856,UWBM_77978> -x 100 resultOutputFile phylonet_allspecies_4ret.txt -di phylonet_allspecies_4ret.tre;\n\nEND;" >> /scratch/dnjacks4/cardinalis/to_parus/phylonet/all/script_4ret.nex


















## pca

### all
cd /scratch/dnjacks4/cardinalis/to_parus/PCA/all

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

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf -o /scratch/dnjacks4/cardinalis/to_parus/PCA/all/ -p /scratch/dnjacks4/cardinalis/to_parus/PCA/all/pops.txt -n all -s y

### pyrr urban vs rural
bcftools view -s MSB25201,PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011,UWBM103346,UWBM77548,UWBM77718,UWBM77780,UWBM77781 /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf --force-samples > /scratch/dnjacks4/cardinalis/to_parus/pyrr_pca.vcf

cd /scratch/dnjacks4/cardinalis/to_parus/PCA/pyrr
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

~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_parus/pyrr_pca.vcf -o /scratch/dnjacks4/cardinalis/to_parus/PCA/pyrr/ -p /scratch/dnjacks4/cardinalis/to_parus/PCA/pyrr/pops.txt -n pyrr -s y

### noca urban vs rural 

bcftools view -s NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013,UWBM100619,UWBM100620,UWBM100621,UWBM103345,UWBM77856,UWBM77978 /scratch/dnjacks4/cardinalis/to_parus/bayescan/parusfiltered.geno25.maf1.vcf --force-samples > /scratch/dnjacks4/cardinalis/to_parus/noca_pca.vcf


cd /scratch/dnjacks4/cardinalis/to_parus/PCA/noca

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


~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_parus/noca_pca.vcf -o /scratch/dnjacks4/cardinalis/to_parus/PCA/noca/ -p /scratch/dnjacks4/cardinalis/to_parus/PCA/noca/pops.txt -n noca -s y




# run structure on it instead of pca
# i'm doing this on sol instead because agave is down

module load gsl-2.7.1-gcc-11.2.0
module load mamba/latest
# mamba create -n fastStructrueEnv_2 -c conda-forge scipy numpy cython=0.26.1 python=2.7
source activate fastStructrueEnv_2

cd ~/programs/

# git clone https://github.com/rajanil/fastStructure
cd fastStructure

# We now have an environment with the source code for fastStructure Our gsl is in a non-standard path, so we have to define some environment variables (LD_LIBRRAY_PATH is set for us by the module load)

export CFLAGS="-l/packages/apps/spack/18/opt/spack/gcc-11.2.0/gsl-2.7.1-cid/include/gsl/"
export LDFLAGS="-L/packages/apps/spack/18/opt/spack/gcc-11.2.0/gsl-2.7.1-cid/lib/"

Now we can build fastStructure
cd vars
python setup.py build_ext --inplace
cd ../
python setup.py build_ext --inplace


# Ran these on agave then copied the bed file over



# Filtering for Structure (one SNP per 10k)

module purge
module load r-4.2.2-gcc-11.2.0
module load sqlite-3.38.5-gcc-11.2.0
module load proj-8.2.1-gcc-11.2.0
module load gdal-3.4.3-gcc-11.2.0
module load geos-3.9.1-gcc-11.2.0
module load  gcc-12.1.0-gcc-11.2.0
module load libxc-5.1.7-gcc-12.1.0
module load libvterm-0.1.4-gcc-11.2.0

R
install.packages("SNPfiltR")

library(SNPfiltR)
library(vcfR)
vcf <- read.vcfR( "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_namedforplink.vcf", verbose = FALSE )

vcf

 
> vcf
***** Object of Class vcfR *****
24 samples
1047 CHROMs
3,246,716 variants
Object size: 1224.8 Mb
0 percent missing data
*****        *****         *****

vcf2 <- distance_thin(vcf, min.distance = 10000)

> distance_thin(vcf, min.distance = 10000)
  |======================================================================| 100%

84572 out of 3246716 input SNPs were not located within 10000 base-pairs of another SNP and were retained despite filtering
***** Object of Class vcfR *****
24 samples
1047 CHROMs
84,572 variants
Object size: 35.2 Mb
0 percent missing data

write.vcf(vcf2, "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_10kfilter.structure.vcf")

### repeat the above with the biallelic file 

bcftools view -m2 -M2 -v snps /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_geno02.vcf -o /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf

bcftools view -m2 -M2 -v snps /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_geno1.vcf -o /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.1.vcf

bcftools view -m2 -M2 -v snps /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_geno15.vcf -o /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.15.vcf

bcftools view -m2 -M2 -v snps /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_geno2.vcf -o /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.2.vcf


sed 's/_//g' /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic.vcf > /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic_namedforplink.vcf

sed 's/_//g' /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.1.vcf > /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.1.plink.vcf

sed 's/_//g' /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.15.vcf > /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.15.plink.vcf

sed 's/_//g' /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.2.vcf > /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.2.plink.vcf

R
library(SNPfiltR)
library(vcfR)

## Biallelic 
vcf <- read.vcfR( "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered_biallelic_namedforplink.vcf", verbose = FALSE )

vcf2 <- distance_thin(vcf, min.distance = 10000)

write.vcf(vcf2, "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_biallelic_10kfilter.structure.vcf")

> vcf2
***** Object of Class vcfR *****
24 samples
179 CHROMs
486 variants
Object size: 0.2 Mb
0 percent missing data

# bigger geno filters 

vcf <- read.vcfR( "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.1.plink.vcf", verbose = FALSE )

vcf2 <- distance_thin(vcf, min.distance = 10000)

# 787 out of 12668 input SNPs were not located within 10000 base-pairs of another SNP and were retained despite filtering

write.vcf(vcf2, "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.1.plink.thinned.vcf")


vcf <- read.vcfR( "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.15.plink.vcf", verbose = FALSE )

vcf2 <- distance_thin(vcf, min.distance = 10000)

# 898 out of 14711 input SNPs were not located within 10000 base-pairs of another SNP and were retained despite filtering

write.vcf(vcf2, "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.15.plink.thinned.vcf")



vcf <- read.vcfR( "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.2.plink.vcf", verbose = FALSE )

vcf2 <- distance_thin(vcf, min.distance = 10000)

# 989 out of 16375 input SNPs were not located within 10000 base-pairs of another SNP and were retained despite filtering


write.vcf(vcf2, "/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.2.plink.thinned.vcf")

## these seem to be the same as the file without the bcftools filter...

#  biallelic 0.02
plink --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_biallelic_10kfilter.structure.vcf --recode structure --out forStructure_biallelic_10kfilter --allow-extra-chr 

mv forStructure_biallelic_10kfilter.* /scratch/dnjacks4/cardinalis/to_parus/vcfs/

cd /scratch/dnjacks4/cardinalis/to_parus/vcfs/

head -1 forStructure_biallelic_10kfilter.recode.strct_in  > forStructure_biallelic_10kfilter.recode.strct_in.nomap
tail -n +3 forStructure_biallelic_10kfilter.recode.strct_in >> forStructure_biallelic_10kfilter.recode.strct_in.nomap

awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' popfile.txt /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure_biallelic_10kfilter.recode.strct_in.nomap > /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure_biallelic_10kfilter.txt

# biallelic 0.1
cd /scratch/dnjacks4/cardinalis/to_parus/vcfs/

plink --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.1.plink.thinned.vcf --recode structure --out forStructure.1.thinned --allow-extra-chr 

head -1 forStructure.1.thinned.recode.strct_in  > forStructure.1.thinned.nomap
tail -n +3 forStructure.1.thinned.recode.strct_in >> forStructure.1.thinned.nomap

awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' popfile.txt /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.1.thinned.nomap > /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.1.thinned.txt

# biallelic 0.15
plink --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.15.plink.thinned.vcf --recode structure --out forStructure.15.thinned --allow-extra-chr 

head -1 forStructure.15.thinned.recode.strct_in  > forStructure.15.thinned.nomap
tail -n +3 forStructure.15.thinned.recode.strct_in >> forStructure.15.thinned.nomap

awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' popfile.txt /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.15.thinned.nomap > /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.15.thinned.txt

# biallelic 0.2
plink --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_filtered.2.plink.thinned.vcf --recode structure --out forStructure.2.thinned --allow-extra-chr 

head -1 forStructure.2.thinned.recode.strct_in  > forStructure.2.thinned.nomap
tail -n +3 forStructure.2.thinned.recode.strct_in >> forStructure.2.thinned.nomap

awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' popfile.txt /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.2.thinned.nomap > /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.2.thinned.txt



echo -e ".\n4\n1\n1\n1\n1\n1\n1\n3\n3\n3\n3\n3\n3\n2\n2\n2\n2\n4\n4\n4\n4\n4\n2\n2" > popfile.txt


awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' popfile.txt /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure_biallelic_10kfilter.recode.strct_in.nomap > /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure_biallelic_10kfilter.txt


# geno 1
cd /scratch/dnjacks4/cardinalis/to_parus/june1_structure/geno1

sed -i 's/forStructure_biallelic_10kfilter/forStructure.1.thinned.txt/g' mainparams

sed -i 's/486/787/g' mainparams

sed -i 's/biallelic_structure/june1_structure/geno1/g' mainparams


~/programs/structure_kernel_src/structure 

## edit and revise mainparams for 2k, 3k, and 4k then rerun 

~/programs/distruct1.1/distructLinux1.1

# geno 15
cd /scratch/dnjacks4/cardinalis/to_parus/june1_structure/geno15

sed -i 's/forStructure_biallelic_10kfilter/forStructure.15.thinned/g' mainparams

sed -i 's/486/898/g' mainparams

~/programs/structure_kernel_src/structure 

# geno 2
cd /scratch/dnjacks4/cardinalis/to_parus/june1_structure/geno2

sed -i 's/forStructure_biallelic_10kfilter/forStructure.2.thinned/g' mainparams

sed -i 's/486/989/g' mainparams

~/programs/structure_kernel_src/structure 



## former code, too big and not filtered by linkage 

plink --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_toplink.vcf.gz --make-bed --out /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort.vcf --allow-extra-chr 

sed 's/_//g' /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort.vcf > /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_namedforplink.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_namedforplink.vcf --recode structure --out forStructure --allow-extra-chr 




# plink exports a header with the first row of maker names and the second of the difference between the current variant's base-pair coordinate and the previous variant's bp coordinate (or -1 when the current variant starts a new chromosome)

head -1 forStructure.recode.strct_in  > forStructure.recode.strct_in.nomap
tail -n +3 forStructure.recode.strct_in >> forStructure.recode.strct_in.nomap
# https://www.biostars.org/p/313943/
# http://rajanil.github.io/fastStructure/


# I need to recode population assignments using awk. It's currently assigning each individual a random number.

.
MSB25201    4   PYRR
NOCA003 1
NOCA004 1
NOCA006 1
NOCA008 1
NOCA012 1
NOCA013 1
PYRR003 3
PYRR004 3
PYRR006 3
PYRR007 3
PYRR009 3
PYRR011 3
UWBM100619 2 NOCA
UWBM100620  2   NOCA
UWBM100621  2   NOCA
UWBM103345  2   NOCA
UWBM103346  4   PYRR
UWBM77548   4   PYRR
UWBM77718   4   PYRR
UWBM77780   4   PYRR
UWBM77781   4   PYRR
UWBM77856   2   NOCA
UWBM77978   2   NOCA


echo -e ".\n4\n1\n1\n1\n1\n1\n1\n3\n3\n3\n3\n3\n3\n2\n2\n2\n2\n4\n4\n4\n4\n4\n2\n2" > popfile.txt

# awk {'print $1'} /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.recode.strct_in.nomap

# rewrite the first row to read 'Label'

# awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' labels.txt /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.recode.strct_in.nomap > temp.txt


awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' popfile.txt /scratch/dnjacks4/cardinalis/to_parus/vcfs/forStructure.recode.strct_in.nomap > /scratch/dnjacks4/cardinalis/to_parus/vcfs/structurefile.txt

# rm temp.txt
awk {'print $1'} /scratch/dnjacks4/cardinalis/to_parus/vcfs/structurefile.txt
awk {'print $2'} /scratch/dnjacks4/cardinalis/to_parus/vcfs/structurefile.txt

~/programs/structure_kernel_src/structure 

~/programs/distruct1.1/distructLinux1.1



nohup python ~/programs/fastStructure/structure.py -K 1 --input=/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_bed --output=/scratch/dnjacks4/cardinalis/to_parus/structure/output --full --seed=100 &

nohup python ~/programs/fastStructure/structure.py -K 2 --input=/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_bed --output=/scratch/dnjacks4/cardinalis/to_parus/structure/output --full --seed=100 &

nohup python ~/programs/fastStructure/structure.py -K 3 --input=/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_bed --output=/scratch/dnjacks4/cardinalis/to_parus/structure/output --full --seed=100 &

nohup python ~/programs/fastStructure/structure.py -K 4 --input=/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_bed --output=/scratch/dnjacks4/cardinalis/to_parus/structure/output --full --seed=100 &

nohup python ~/programs/fastStructure/structure.py -K 5 --input=/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_bed --output=/scratch/dnjacks4/cardinalis/to_parus/structure/output --full --seed=100 &

python ~/programs/fastStructure/structure.py -K 6 --input=/scratch/dnjacks4/cardinalis/to_parus/vcfs/parus_qualitysort_bed --output=/scratch/dnjacks4/cardinalis/to_parus/structure/output --full --seed=100 


python ~/programs/fastStructure/chooseK.py --input=/scratch/dnjacks4/cardinalis/to_parus/structure/output


python ~/programs/fastStructure/distruct.py -K 6 --input=/scratch/dnjacks4/cardinalis/to_parus/structure/output --output=/scratch/dnjacks4/cardinalis/to_parus/structure/output_k6.svg



python ~/programs/fastStructure/distruct.py -K 6 --input=/scratch/dnjacks4/cardinalis/to_parus/structure/output --output=/scratch/dnjacks4/cardinalis/to_parus/structure/output_k6.svg



Visualizing it, becasue distruct sucks
R
library(pophelper)
x <- readQ("output.6.meanQ")