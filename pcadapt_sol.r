module purge
module load r-4.2.2-gcc-11.2.0
module load bayescan/2.1
module load htslib-1.16-gcc-11.2.0
module load bcftools-1.14-gcc-11.2.0


path_to_file <- "/scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1.vcf"






# PCAdapt

# northern cardinals

plink --vcf /scratch/dnjacks4/cardinalis/to_cardinalis/bayescan/noca.filtered.geno25.maf1.vcf --allow-extra-chr --out noca.filtered.geno25.maf1 --make-bed


R 

library(pcadapt)
library(qvalue)

path_to_file <- "/scratch/dnjacks4/cardinalis/to_cardinalis/pcadapt/noca.filtered.geno25.maf1.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
x <- pcadapt(input = filename, K = 1)

pdf(file = "hist_pvalues.pdf", width = 10, height = 10, useDingbats=FALSE)
  hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
    dev.off()

poplist.names <- c(rep("Urban", 6), rep("Rural", 6))
pdf(file = "pca.pdf", width = 10, height = 10, useDingbats=FALSE)
  plot(x, option = "scores") #, pop = poplist.names)
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
