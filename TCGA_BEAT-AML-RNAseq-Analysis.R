# Joint-RNAseq-Analysis

library(limma)
library(sva)
setwd("~/CRNDE/Analyze")

# read TCGA gene expression matrix
tcgaAML_expr <- read.table("tcgaAML.expr.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)
# remove genes with expression(FPKM) less than 0.1 in less than 25% samples.
tcgaAML_expr <- tcgaAML_expr[rowSums(tcgaAML_expr>0.1)>floor(0.25 * ncol(tcgaAML_expr)),]
# log2 transform
tcgaAML_expr <- log2(tcgaAML_expr + 1)

# read beatAML gene expression matrix
beatAML_expr <- read.table("beatAML.expr.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)
# remove genes with expression(FPKM) less than 0.1 in less than 25% samples.
beatAML_expr <- beatAML_expr[rowSums(beatAML_expr>0.1)>floor(0.25 * ncol(beatAML_expr)),]
beatAML_expr <- log2(beatAML_expr + 1)

# combine two cohorts
expr <- merge(tcgaAML_expr, beatAML_expr, by.x = 0, by.y = 0, all = FALSE, sort = FALSE)
rownames(expr) <- expr$Row.names
expr <- expr[, -1]
samples <- sort(colnames(expr))
expr <- expr[, samples]

# remove batch effect
pheno <- data.frame(batch = as.numeric(grepl('TCGA', colnames(expr)))+1)
rownames(pheno) <- colnames(expr)
batch <- pheno$batch
modcombat<-model.matrix(~1, data=pheno)
expr <- ComBat(dat = as.matrix(expr), batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
expr <- data.frame(expr)
write.table(expr, file="Joint_RNAseq_expression.txt", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

# diffrential analysis
design <- c(rep("AML", 350), rep("APL", 29), rep("Normal", 21))
#colnames(expr) <- design
design <- factor(design, levels = c("AML", "Normal", "APL")) 
design.mat <- model.matrix(~0+design)
colnames(design.mat) <- levels(design)
fit <- lmFit(expr,design.mat)
cont.matrix <- makeContrasts("APL-Normal","APL-AML",levels = design.mat)
fit2 <- contrasts.fit(fit, cont.matrix) 
fit2 <- eBayes(fit2)

# DGE of APL.vs.Normal
APL.vs.Normal <- topTable(fit2, coef = 1, number = Inf)
anno <- read.table("EnsemblToSymbol.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)
biotype <- read.table("EnsemblToBiotype.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)
APL.vs.Normal <- merge(biotype, APL.vs.Normal, by.x=0, by.y=0, all.x=FALSE, all.y=TRUE)
lncRNA <- c("processed_transcript", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "lincRNA")
APL.vs.Normal <- APL.vs.Normal[APL.vs.Normal$biotype %in% lncRNA, ][, -2]
APL.vs.Normal <- merge(anno, APL.vs.Normal, by.x=0, by.y=1, all.x=FALSE, all.y=TRUE)
APL.vs.Normal <- merge(APL.vs.Normal, expr, by.x=1, by.y=0, all.x=TRUE, all.y=FALSE)
colnames(APL.vs.Normal)[1] <- "ID"
write.table(APL.vs.Normal, file="APL.vs.Normal_lncRNA_limma.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

APL.vs.Normal.diff <- APL.vs.Normal[rownames(APL.vs.Normal[!is.na(APL.vs.Normal$adj.P.Val) & APL.vs.Normal$adj.P.Val<0.05 & abs(APL.vs.Normal$logFC)>=log2(1.5),]), ]
write.table(APL.vs.Normal.diff, file="APL.vs.Normal_lncRNA_limma.diff.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# DGE of APL.vs.AML
APL.vs.AML <- topTable(fit2, coef = 2, number = Inf)
anno <- read.table("EnsemblToSymbol.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)
biotype <- read.table("EnsemblToBiotype.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)
APL.vs.AML <- merge(biotype, APL.vs.AML, by.x=0, by.y=0, all.x=FALSE, all.y=TRUE)
lncRNA <- c("processed_transcript", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "lincRNA")
APL.vs.AML <- APL.vs.AML[APL.vs.AML$biotype %in% lncRNA, ][, -2]
APL.vs.AML <- merge(anno, APL.vs.AML, by.x=0, by.y=1, all.x=FALSE, all.y=TRUE)
APL.vs.AML <- merge(APL.vs.AML, expr, by.x=1, by.y=0, all.x=TRUE, all.y=FALSE)
colnames(APL.vs.AML)[1] <- "ID"
write.table(APL.vs.AML, file="APL.vs.AML_lncRNA_limma.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

APL.vs.AML.diff <- APL.vs.AML[rownames(APL.vs.AML[!is.na(APL.vs.AML$adj.P.Val) & APL.vs.AML$adj.P.Val<0.05 & abs(APL.vs.AML$logFC)>=log2(1.5),]), ]
write.table(APL.vs.AML.diff, file="APL.vs.AML_lncRNA_limma.diff.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# Heatmap of APL.vs.AML
Expr <- read.table("03.APL.vs.AML_heatmap/Joint_expr_addLabel.txt", header = TRUE, row.names = 2, stringsAsFactors = FALSE, sep = "\t")[, -1]
sampleNames <- colnames(Expr)
order_samples <- c(sampleNames[grepl("APL", sampleNames)], sampleNames[grepl("AML", sampleNames)])
group <- gsub(".*.TCGA.*", "TCGA", order_samples)
group <- gsub(".*.BA.*", "BEAT", group)
Expr <- Expr[order_samples]
Expr <- as.matrix(Expr)
#Expr <- log2(Expr + 1)
Expr <- Expr - apply(Expr, 1, mean)
Expr[which(Expr < -1.5)] <- -1.5
Expr[which(Expr > 1.5)] <- 1.5

nsamples <- ncol(Expr)
ngenes <- nrow(Expr)
width <- 180/nsamples
height <- 240/ngenes

library(gplots)
library(pheatmap)
annotation_col = data.frame(
  Group=factor(group),
  Type=factor(c(rep("APL", 29), rep("Non-APL", 350)))
)
rownames(annotation_col) <- colnames(Expr)
annotation_colors = list(
  Group = c("BEAT" = "red", "TCGA" = "blue"),
  Type = c("APL" = "#219418", "Non-APL" = "#FF6A00")
)

p <- pheatmap(as.matrix(Expr), annotation_colors=annotation_colors, 
         scale="none", show_rownames=TRUE, show_colnames=FALSE, 
         col = colorRampPalette(c( "blue", "white", "red"))( 100 ), 
         annotation_col=annotation_col, main="", annotation_names_col=FALSE, 
         annotation_legend=TRUE, cellheight=height, cellwidth=width, fontsize=6, 
         cluster_cols = TRUE, treeheight_col=10, treeheight_row=10, cutree_cols = 2,
         border = NA,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "canberra")
pdf(file="~/CRNDE/Analyze/03.APL.vs.AML_heatmap/APL.vs.Non-APL.heatmap.addLabel.pdf", width = 5.7, height = 5)
p
dev.off()

# Heatmap of APL.vs.Normal
Expr <- read.table("04.APL.vs.Normal_heatmap/Joint_expr_addLabel.txt", header = TRUE, row.names = 2, stringsAsFactors = FALSE, sep = "\t")[, -1]
sampleNames <- colnames(Expr)
order_samples <- c(sampleNames[grepl("Normal", sampleNames)], sampleNames[grepl("APL", sampleNames)])
group <- gsub(".*.TCGA.*", "TCGA", order_samples)
group <- gsub(".*.BA.*", "BEAT", group)
Expr <- Expr[order_samples]
Expr <- as.matrix(Expr)
#Expr <- log2(Expr + 1)
Expr <- Expr - apply(Expr, 1, mean)
Expr[which(Expr < -2)] <- -2
Expr[which(Expr > 2)] <- 2

nsamples <- ncol(Expr)
ngenes <- nrow(Expr)
width <- 180/nsamples
height <- 240/ngenes

library(gplots)
library(pheatmap)
annotation_col = data.frame(
  Group=factor(group),
  Type=factor(c(rep("Normal", 21), rep("APL", 29)))
)
rownames(annotation_col) <- colnames(Expr)
annotation_colors = list(
  Group = c("BEAT" = "red", "TCGA" = "blue"),
  Type = c("APL" = "#FF6A00", "Normal" = "#219418")
)

p <- pheatmap(as.matrix(Expr), annotation_colors=annotation_colors, 
         scale="none", show_rownames=TRUE, show_colnames=FALSE, 
         col = colorRampPalette(c( "blue", "white", "red"))( 100 ), 
         annotation_col=annotation_col, main="", annotation_names_col=FALSE, 
         annotation_legend=TRUE, cellheight=height, cellwidth=width, fontsize=6, 
         cluster_cols = TRUE, treeheight_col=10, treeheight_row=10, cutree_cols = 2,
         border = NA, cutree_rows = 1,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "manhattan")
pdf(file="~/CRNDE/Analyze/04.APL.vs.Normal_heatmap/APL.vs.Normal.heatmap.addLabel.pdf", width = 5.7, height = 5)
p
dev.off()

# Venn Plot
library(venn)
options(stringsAsFactors = FALSE)
APL.vs.AML_Up <- read.table("~/CRNDE/Analyze/03.Specific_and_Oncogenic.venn/Up.AML")[,1]
APL.vs.AML_Dn <- read.table("~/CRNDE/Analyze/03.Specific_and_Oncogenic.venn/Dn.AML")[,1]
APL.vs.Normal_Up <- read.table("~/CRNDE/Analyze/03.Specific_and_Oncogenic.venn/Up.Normal")[,1]
APL.vs.Normal_Dn <- read.table("~/CRNDE/Analyze/03.Specific_and_Oncogenic.venn/Dn.Normal")[,1]

d <- list(
  APL.vs.AML_Up = APL.vs.AML_Up,
  APL.vs.Normal_Up = APL.vs.Normal_Up,
  APL.vs.AML_Dn = APL.vs.AML_Dn,
  APL.vs.Normal_Dn = APL.vs.Normal_Dn
)

pdf(file = "~/CRNDE/Analyze/03.Specific_and_Oncogenic.venn/Venn.pdf", width = 5, height = 5)
venn(d, box = FALSE, zcolor = c("red", "orange"), opacity = 0.4)
dev.off()
