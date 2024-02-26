# Significant Regions
# FST 

# noca 
List of genes of interest:
Plekhf2,Ch037,Ctnna3_0,Lrrtm3,Mrps5

grep 'Plekhf2' /scratch/dnjacks4/cardinalis/to_b10k/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' 

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


Tajima's D:
Northern Cardinals
Urban
"SNP" "myBg" "CHROM" "BIN_START" "N_SNPS" "TajimaD"
"1693" 1039 TRUE 1009321.1 580000 27 2.82764 Genes in range: Dlg2
"2951" 1847 TRUE 1016584.1 0 5 -1.50034 Genes in range: Rexo1

Rural
"SNP" "myBg" "CHROM" "BIN_START" "N_SNPS" "TajimaD"
"1510" 884 TRUE 1007098.1 0 8 -1.6253 Genes in range: Mtco2_0,Kcna6
"2361" 1383 TRUE 1011537.1 0 6 -1.55221 Genes in range: Tp73,Tp73as1,Ccdc27,Lrrc47,Cep104,Dffb
"3736" 2342 TRUE 1019264.1 160000 9 -1.65206 Genes in range: Slc26a5
"4843" 3087 TRUE 1026539.1 0 18 2.82101 Genes in range: None 

Pyrrhuloxia
Urban
"SNP" "myBg" "CHROM" "BIN_START" "N_SNPS" "TajimaD"
"26" 26 TRUE 1000233.1 0 7 -1.59273 Genes in range: None 
"193" 134 TRUE 1001023.1 0 18 2.63265 Genes in range: None
"556" 351 TRUE 1002912.1 20000 5 -1.50034 Genes in range: Dcbld2
"1824" 1106 TRUE 1009321.1 580000 23 2.67515 Genes in range: Dlg2

Rural 
"SNP" "myBg" "CHROM" "BIN_START" "N_SNPS" "TajimaD"
"218" 171 TRUE 1001302.1 0 13 2.52497 Genes in range: None
"1663" 950 TRUE 1007741.1 620000 4 -1.43145 Genes in range: Rad21
"2081" 1227 TRUE 1010057.1 0 4 -1.43145 Genes in range: None 
"2139" 1269 TRUE 1010419.1 20000 12 -1.433 Genes in range: None
"2213" 1296 TRUE 1010703.1 0 144 2.49775 Genes in range: None
"3883" 2420 TRUE 1019758.1 0 5 -1.50034 Genes in range: None
"3888" 2425 TRUE 1019804.1 280000 7 2.53271 Genes in range: Kcnq3,Oc90,Efr3a


# Nucleotide Diversity 
Northern Cardinals
Urban
"SNP" "myBg" "CHROM" "BIN_START" "BIN_END" "N_VARIANTS" "PI"
"921" 921 TRUE 1006187.1 1 20000 666 0.0123613 Genes in range: Hydin_0
Rural
"SNP" "myBg" "CHROM" "BIN_START" "BIN_END" "N_VARIANTS" "PI"
"936" 936 TRUE 1006187.1 1 20000 708 0.0121765 Genes in range: Hydin_0

Pyrrhuloxia
Urban
"SNP" "myBg" "CHROM" "BIN_START" "BIN_END" "N_VARIANTS" "PI"
"904" 904 TRUE 1006187.1 1 20000 652 0.0110707 Genes in range: Hydin_0

Rural
"SNP" "myBg" "CHROM" "BIN_START" "BIN_END" "N_VARIANTS" "PI"
"906" 906 TRUE 1006187.1 1 20000 688 0.0117187 Genes in range: Hydin_0


# Sweed
Northern cardinals
Urban
"SNP" "myBg" "Scaffold" "Position" "Likelihood" "Alpha" "StartPos" "EndPos"
"39018" 39018 TRUE 1006187.1 2650.1213 6.195905 0.003933158 198 5561 Genes in range: Hydin_0
Rural
"SNP" "myBg" "Scaffold" "Position" "Likelihood" "Alpha" "StartPos" "EndPos"
"163110" 163110 TRUE 1023961.1 2405.5455 10.30593 0.005772233 1073 2405 Genes in range: None

Pyrrhuloxia
Urban
"SNP" "myBg" "Scaffold" "Position" "Likelihood" "Alpha" "StartPos" "EndPos"
"40708" 40708 TRUE 1006187.1 1076.4444 12.6683 0.009598233 56 2175 Genes in range: Hydin_0
"40709" 40709 TRUE 1006187.1 1222.2222 12.10249 0.003522958 56 4457 Genes in range: Hydin_0

Rural
"SNP" "myBg" "Scaffold" "Position" "Likelihood" "Alpha" "StartPos" "EndPos"
"155790" 155790 TRUE 1022801.1 5319.6971 8.272938 0.0370888 5198 5593 Genes in range: None



# Omega
Northern cardinals
Urban
1013494.1 364 - 442 Genes in range: None

Rural
1017116.1 (116 regions) 2124.3555 - 2138.897 Genes in range: None
1023961.1 (542 regions) 4638.4604 - 12327.0342 Genes in range: None

Pyrrhuloxia
Urban
1008477.1 966.0082 - 1218.4209 Genes in range: None


Rural
1015323.1 152.8487 - 1532.1931


# Bayescan
None 
# Hapflk
???
# 