# SweeD 
# Tests for selective sweeps
# urban northern cardinals 

module purge
module load r-4.0.2-gcc-11.2.0 #r-4.2.2-gcc-11.2.0
module load htslib-1.16-gcc-11.2.0
module load bcftools-1.14-gcc-11.2.0
module load vcftools-0.1.14-gcc-11.2.0

cd /scratch/dnjacks4/cardinalis/to_b10k/sweed/

bcftools view -s 'NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf.gz > /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid  --out urban_noca.filtered.geno25.maf1

cd northerncardinals

~/programs/sweed/SweeD -name urban_noca -input /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf -grid 100 -length 100000


bgzip /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf
tabix /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz
zcat /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names_noca_urban.txt

scp -r to_b10k/omega dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/cardinalis/to_b10k/sweed

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
done < SweeD_Report.urban_noca

grep 'Position' SweeD_Report.urban_noca | head -n 1 > SweeD_Report.header
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.noca_urban_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names_noca_urban.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> SweeD_Report.noca_urban_manhattan_scaffnames
  fi
done < SweeD_Report.noca_urban_manhattan

sed -i 's/VYXE//g' SweeD_Report.noca_urban_manhattan_scaffnames
# sed -i 's/\(.\{1\}\).1/\1/' SweeD_Report.noca_urban_manhattan_scaffnames
# awk '{sub(/\./,"",$1)}1' SweeD_Report.noca_urban_manhattan_scaffnames | column -t > SweeD_Report.noca_urban_manhattan_scaffnames.formanhattan

R
library(qqman)
sweed<-read.table("SweeD_Report.noca_urban_manhattan_scaffnames", header=TRUE)
sweed.subset<-sweed[complete.cases(sweed),]
SNP<-c(1: (nrow(sweed.subset)))

lower = min(sweed.subset$Likelihood)
upper = max(sweed.subset$Likelihood)
cutoff = upper - ((upper-lower)*0.05)
LessThanCutoff <- sweed.subset$Likelihood < cutoff
myBg <- !LessThanCutoff
mydf<-data.frame(SNP,myBg,sweed.subset)
sigdf <-  mydf[which(mydf$myBg),]
write.table(sigdf, file = "sigsweed.tsv")

mydf<-data.frame(SNP,sweedsubset)

pdf(file = "noca_urban_sweed.pdf", width = 20, height = 7, useDingbats=FALSE)
    print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()

scp dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/cardinalis/to_b10k/sweed/northerncardinals/noca_urban_sweed.pdf .