# pDC-project
#install limma package
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
#read files
#you must create first a targets /t .txt file. First column is sample number, 2nd filename for single array(txt) and 3rd=condition
targets <- readTargets("targets.txt")
x <- read.maimages(targets, path="E:/R/pDC project", source="agilent",green.only=TRUE)
#correct background
y <- backgroundCorrect(x, method="normexp", offset=25)
#normalize between arrays
y <- normalizeBetweenArrays(y, method="quantile")
#average replicates
y.ave <- avereps(y, ID=y$genes$ProbeName)
#Build the design matrix for the linear modelling function.
f <- factor(targets$Condition, levels = unique(targets$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
#Apply the intensity values to lmFit.
fit <- lmFit(y.ave, design)
#Create a contrast matrix. In this example, all combinations of contrasts can be set up as below.
contrast.matrix <- makeContrasts("prediab-diab", "prediab-aCD3", "diab-aCD3", levels=design)
#Apply this contrast matrix to the modeled data and compute statistics for the data.
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
#Output the statistics for the dataset and write them to disk for further analysis.
output <- topTable(fit2, adjust="BH", coef="prediab-diab", genelist=y.ave$genes, number=Inf)
write.table(output, file="prediab_vs_diab.txt", sep="\t", quote=FALSE)
output <- topTable(fit2, adjust="BH", coef="prediab-aCD3", genelist=y.ave$genes, number=Inf)
write.table(output, file="prediab_vs_aCD3.txt", sep="\t", quote=FALSE)
output <- topTable(fit2, adjust="BH", coef="diab-aCD3", genelist=y.ave$genes, number=Inf)
write.table(output, file="diab_vs_aCD3.txt", sep="\t", quote=FALSE)
#Output the normalized log2 intensity values to analyze with qlucore
pDC_qlucore <- write.table(y.ave, file = "df.csv", sep=";")
