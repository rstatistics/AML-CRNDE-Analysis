library(affy)
library(affyPLM)
library(hgu133plus2hsensgcdf)
setwd("~/CRNDE/Analyze/04.APL.vs.Pros.microarray")
setwd("~/CRNDE/Analyze/04.APL.vs.Pros.microarray")
celFiles <- list.celfiles(path = "./GSE12662", full.names=TRUE)
data.affy <- ReadAffy(filenames = celFiles, cdfname = "HGU133Plus2_Hs_ENSG")
eset <- rma(data.affy)

MAplot(eset, pairs = TRUE, which=c(1,3,5), plot.method = "smoothScatter")
hist(eset)
boxplot(data.frame(exprs(eset)), col=c(1,3,5))

eset@protocolData@data
pheno <- read.table("./type.txt",sep = "\t",header = TRUE)
dist_mat <- dist(t(exprs(eset)))
clustering <- hclust(dist_mat, method = "complete")
plot(clustering, labels = pheno$type)

library(hgu133plus2.db)

Expr.raw <- exprs(eset)
rownames(Expr.raw) <- gsub("\\_at", "", rownames(Expr.raw))

# prepare gene_biotype database
library(EnsDb.Hsapiens.v79)
## Making a "short cut"
edb <- EnsDb.Hsapiens.v79
## Get the ENSEMBL ID
keys <- rownames(Expr.raw)
## Retrieve all gene IDs of all lncRNAs

lncRNA <- keys(edb, filter=
                     list(GeneBiotypeFilter(c(
                     "processed_transcript",
                     "3prime_overlapping_ncrna",
                     "antisense",
                     "non_coding",
                     "sense_intronic",
                     "sense_overlapping",
                     "lincRNA")
                     )),
               keytype="GENEID")
lncRNA <- lncRNA[!duplicated(lncRNA)]
lncRNA <- select(edb, keys=lncRNA, columns=c("SYMBOL", "GENEBIOTYPE"), keytype = "GENEID")
# de-duplication SYMBOL
lncRNA <- lncRNA[!duplicated(lncRNA$SYMBOL),]

ENSEMBL2DB <- data.frame(ID=keys(edb, keytype="GENEID"))

Expr.raw <- merge(data.frame(ID=ENSEMBL2DB$ID), Expr.raw, by.x=1, by.y=0, all = FALSE)

rownames(Expr.raw) <- Expr.raw$ID
Expr.raw <- Expr.raw[, -1]

## Making a "short cut"
Expr <- Expr.raw

design <- c("M3","M3","M3","M3","M3","M3","M3","M3","M3","M3","M3","M3","M3", "M3", "CD34","PMNs","PMNs","CD34","CD34", "Pros","PMNs", "Pros",
            "Pros","PMNs","PMNs","CD34","CD34","Pros","Pros")
colnames(Expr) <- design
design <- factor(design, levels = c("M3","CD34","Pros","PMNs")) 
design.mat <- model.matrix(~0+design)
colnames(design.mat) <- levels(design)
fit <- lmFit(Expr, design.mat)
cont.matrix <- makeContrasts("M3-CD34","M3-Pros","M3-PMNs",levels = design.mat)
fit2 <- contrasts.fit(fit, cont.matrix) 
fit2 <- eBayes(fit2)

# APL-vs-CD34+
allDiff_CD34 <- topTable(fit2,coef=1,adjust.method="BH",number = Inf)
# add Expression matrix
allDiff_CD34 <- merge(allDiff_CD34, Expr.raw[,which(colnames(Expr) %in% c("M3","CD34"))], by.x=0, by.y=0, all.x=TRUE)
# annotation
allDiff_CD34 <- merge(data.frame(ID=lncRNA$GENEID, SYMBOL=lncRNA$SYMBOL), allDiff_CD34, by.x=1, by.y =1, all = FALSE)
allDiff_CD34 <- allDiff_CD34[!duplicated(allDiff_CD34$SYMBOL),]

allDiff_CD34 <- allDiff_CD34[order(allDiff_CD34$logFC, decreasing = TRUE), ]

write.table(allDiff_CD34, file="APL-vs-CD34.microarray.all.txt", row.names = FALSE, col.names = TRUE, append = FALSE, sep="\t", quote=FALSE)
# filter DGE entries
allDiff_CD34 <- allDiff_CD34[with(allDiff_CD34,abs(allDiff_CD34$logFC)>log2(1.5) & allDiff_CD34$adj.P.Val < 0.05),]

write.table(allDiff_CD34, file="APL-vs-CD34.microarray.txt", row.names = FALSE, col.names = TRUE, append = FALSE, sep="\t", quote=FALSE)

# APL-vs-Promyelocytes
allDiff_Pros <- topTable(fit2,coef=2,adjust.method="BH",number = Inf)
# add Expression matrix
allDiff_Pros <- merge(allDiff_Pros, Expr.raw[,which(colnames(Expr) %in% c("M3","Pros"))], by.x=0, by.y=0, all.x=TRUE)
# annotation
allDiff_Pros <- merge(data.frame(ID=lncRNA$GENEID, SYMBOL=lncRNA$SYMBOL), allDiff_Pros, by.x=1, by.y =1, all = FALSE)
allDiff_Pros <- allDiff_Pros[!duplicated(allDiff_Pros$SYMBOL),]

allDiff_Pros <- allDiff_Pros[order(allDiff_Pros$logFC, decreasing = TRUE), ]

write.table(allDiff_Pros, file="APL-vs-Pros.microarray.all.txt", row.names = FALSE, col.names = TRUE, append = FALSE, sep="\t", quote=FALSE)
# filter DGE entries
allDiff_Pros <- allDiff_Pros[with(allDiff_Pros,abs(allDiff_Pros$logFC)>log2(1.5) & allDiff_Pros$adj.P.Val < 0.05),]

write.table(allDiff_Pros, file="APL-vs-Pros.microarray.txt", row.names = FALSE, col.names = TRUE, append = FALSE, sep="\t", quote=FALSE)
