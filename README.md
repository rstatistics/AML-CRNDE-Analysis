# AML-CRNDE-Analysis
Data analysis for LncRNA CRNDE in AML

## Analysis of The Cancer Genome Atlas and beatAML RNA-seq data
The available RNA-seq data of 171 AML patient samples, which included 16 APL and 155 non-APL AML patient samples, were downloaded from The Cancer Genome Atlas ([TCGA](https://tcga-data.nci.nih.gov/tcga/)) Data Portal. RNA-seq reads were aligned to GRCh38 whole genome with Ensembl v84 using Hisat2 2.0.4. StringTie 1.2.3 and Ballgown 1.0.1 were used to assemble the alignments into full and partial transcripts and estimate the expression levels of all lncRNAs. Limma algorithm [53] was used to analyze the differential expression of the lncRNAs. The utilization of patient data was performed in accordance with the TCGA Human Subjects Protection and Data Access Policies.

## workflow
这里绘制流程图

## Database
The Human genome and annotation files were retrieved from [Genome](ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna) 
and [GTF](ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz).
The hisat2-index files were retrieved from [genome_snp_tran](https://cloud.biohpc.swmed.edu/index.php/s/grch38_snp_tran/download).  

## Interpretation biotype in Ensembl
#### Defination of [LncRNA](https://m.ensembl.org/info/genome/genebuild/biotypes.html) in Ensembl GTF file:
Processed transcript  
3' overlapping ncRNA  
Antisense  
Non coding  
Retained intron  
Sense intronic  
Sense overlapping  
lincRNA  

## Align to the reference genome
```bash
hisat2 -x hisat2_index/genome_snp_tran -1 Input_R1.fq.gz -2 Input_R2.fq.gz \
       --fr --threads 10 --dta 2> Input.summary.txt | samtools view -F 4 - -b | \
       samtools sort -o Input.bam > Input.stdout.log 2> Input.stderr.log

samtools index Input.bam
```
## Transcriptome assembly
```bash
stringtie Input.bam -p 10 -G Homo_sapiens.GRCh38.84.gtf -eB -o Input.gtf
```

## Gene quantification
```r
library(ballgown)
bg <- ballgown(dataDir="~/CRNDE/Analyze/", samplePattern="BA", meas='all')
gene_expression = gexpr(bg)
write.table(gene_expression,"AML_Expr.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
```

## Differential Gene Expression
```sh
Rscript AML-CRNDE-Analysis.R
```
