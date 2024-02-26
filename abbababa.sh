interactive -t 2-00:00 

#Convert vcf to geno
grep 'CHROM' /data5/sulidae/working/mt/raxml_MT/parusfiltered.geno25.maf1.NC040875.1.vcf

python ~/programs/cardinalis/genomics_general/VCF_processing/parseVCF.py -i /data5/sulidae/working/mt/raxml_MT/parusfiltered.geno25.maf1.NC040875.1.vcf -o mitogenome.geno.gz --addRefTrack

echo -e "MSB25201\tP1" > pops.txt 
echo -e "UW103346\tP1" >> pops.txt 
echo -e "UW77548\tP1" >> pops.txt 
echo -e "UW77718\tP1" >> pops.txt 
echo -e "UW77780\tP1" >> pops.txt 
echo -e "UW77781\tP1" >> pops.txt 

echo -e "PYRR003\tP2" >> pops.txt 
echo -e "PYRR004\tP2" >> pops.txt 
echo -e "PYRR006\tP2" >> pops.txt 
echo -e "PYRR007\tP2" >> pops.txt 
echo -e "PYRR009\tP2" >> pops.txt 
echo -e "PYRR011\tP2" >> pops.txt 

echo -e "NOCA003\tP3" >> pops.txt 
echo -e "NOCA004\tP3" >> pops.txt 
echo -e "NOCA006\tP3" >> pops.txt 
echo -e "NOCA008\tP3" >> pops.txt 
echo -e "NOCA012\tP3" >> pops.txt 
echo -e "NOCA013\tP3" >> pops.txt 

echo -e "REF\tP4" >> pops.txt 

bcftools view -s NOCA_003,NOCA_004,NOCA_006,NOCA_008,NOCA_012,NOCA_013,UWBM_100619,UWBM_100620,UWBM_100621,UWBM_103345,UWBM_77856,UWBM_77978 /scratch/dnjacks4/cardinalis/to_parus/parus_qualitysort.vcf --force-samples > /scratch/dnjacks4/cardinalis/to_parus/noca_pca.vcf



~/sula/abba.sh -a mitogenomeabba -b /data5/sulidae/working/mt/abba/output -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/working/mt/abba/mitogenome.geno.gz -f /data5/sulidae/working/mt/abba/pops.txt -g P1 -h P2 -i P3 -j P4 -k y -n 2

~/sula/abba.sh -a mitogenomeabba -b [PATH TO OUTPUT] -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/NAMEOFRUN/pops.txt -g P1 -h P2 -i P3 -j P4 -k y -n 2


(rural PYRR) (all Urban PYRR) (all NOCAs) (Outgroup)