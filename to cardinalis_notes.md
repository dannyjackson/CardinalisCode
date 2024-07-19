#to cardinalis
# scott edwards' cardinal genome

# PCA
cd /scratch/dnjacks4/cardinalis/to_cardinalis/
mkdir pca
cd pca
mkdir all 
cd all 
/scratch/dnjacks4/cardinalis/to_cardinalis/pca/all
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


~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_cardinalis/mydata/cardinalisfiltered.geno25.maf1.vcf -o /scratch/dnjacks4/cardinalis/to_cardinalis/pca/all -p /scratch/dnjacks4/cardinalis/to_cardinalis/pca/all/pops.txt -n all -s y


cd /scratch/dnjacks4/cardinalis/to_cardinalis/pca
mkdir noca 
cd noca 

echo -e "Urban" > pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 


~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1.vcf -o /scratch/dnjacks4/cardinalis/to_cardinalis/pca/noca -p /scratch/dnjacks4/cardinalis/to_cardinalis/pca/noca/pops.txt -n all -s y



cd /scratch/dnjacks4/cardinalis/to_cardinalis/pca
mkdir pyrr 
cd pyrr 

echo -e "Rural" > pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Urban" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 
echo -e "Rural" >> pops.txt 


~/genomics/PCA_r.sh -v /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/pyrr.filtered.geno25.maf1.vcf -o /scratch/dnjacks4/cardinalis/to_cardinalis/pca/pyrr -p /scratch/dnjacks4/cardinalis/to_cardinalis/pca/pyrr/pops.txt -n all -s y
