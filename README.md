# AML-CRNDE-Analysis
Data analysis for LncRNA CRNDE in AML

## Analysis of The Cancer Genome Atlas and beatAML RNA-seq data
The available RNA-seq data of 171 AML patient samples, which included 16 APL and 155 non-APL AML patient samples, were downloaded from The Cancer Genome Atlas ([TCGA](https://tcga-data.nci.nih.gov/tcga/)) Data Portal. RNA-seq reads were aligned to [GRCh38 whole genome](http://asia.ensembl.org/Homo_sapiens/Info/Index) with Ensembl v84 using [Hisat2 2.2.0](https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download), [StringTie 2.1.2](http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.2.Linux_x86_64.tar.gz) and [Ballgown 2.18.0](https://bioconductor.org/packages/release/bioc/html/ballgown.html) were used to assemble the alignments into full and partial transcripts and estimate the expression levels of all lncRNAs. [Limma](https://bioconductor.org/packages/release/bioc/html/limma.html) algorithm was used to analyze the differential expression of the lncRNAs. The utilization of patient data was performed in accordance with the TCGA Human Subjects Protection and Data Access Policies.

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

## Session Information
```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 10 (buster)

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] venn_1.9                    sva_3.34.0                  genefilter_1.68.0          
 [4] mgcv_1.8-31                 nlme_3.1-144                limma_3.42.2               
 [7] RColorBrewer_1.1-2          gplots_3.0.4                lattice_0.20-38            
[10] tximport_1.14.2             pheatmap_1.0.12             DESeq2_1.26.0              
[13] SummarizedExperiment_1.16.1 DelayedArray_0.12.3         BiocParallel_1.20.1        
[16] matrixStats_0.56.0          Biobase_2.46.0              GenomicRanges_1.38.0       
[19] GenomeInfoDb_1.22.1         IRanges_2.20.2              S4Vectors_0.24.4           
[22] SRAdb_1.48.2                RCurl_1.98-1.2              graph_1.64.0               
[25] BiocGenerics_0.32.0         RSQLite_2.2.0              

loaded via a namespace (and not attached):
 [1] bitops_1.0-6           bit64_4.0.5            tools_3.6.3            backports_1.1.9       
 [5] R6_2.4.1               KernSmooth_2.23-16     rpart_4.1-15           Hmisc_4.4-1           
 [9] DBI_1.1.0              colorspace_1.4-1       nnet_7.3-12            tidyselect_1.1.0      
[13] gridExtra_2.3          bit_4.0.4              compiler_3.6.3         htmlTable_2.0.1       
[17] xml2_1.3.2             caTools_1.18.0         scales_1.1.1           checkmate_2.0.0       
[21] readr_1.3.1            stringr_1.4.0          digest_0.6.25          foreign_0.8-75        
[25] GEOquery_2.54.1        XVector_0.26.0         base64enc_0.1-3        jpeg_0.1-8.1          
[29] pkgconfig_2.0.3        htmltools_0.5.0        htmlwidgets_1.5.1      rlang_0.4.7           
[33] rstudioapi_0.11        farver_2.0.3           generics_0.0.2         gtools_3.8.2          
[37] dplyr_1.0.2            magrittr_1.5           GenomeInfoDbData_1.2.2 Formula_1.2-3         
[41] Matrix_1.2-18          Rcpp_1.0.5             munsell_0.5.0          lifecycle_0.2.0       
[45] stringi_1.4.6          zlibbioc_1.32.0        grid_3.6.3             blob_1.2.1            
[49] gdata_2.18.0           crayon_1.3.4           splines_3.6.3          annotate_1.64.0       
[53] hms_0.5.3              locfit_1.5-9.4         knitr_1.29             pillar_1.4.6          
[57] admisc_0.8             geneplotter_1.64.0     XML_3.99-0.3           glue_1.4.2            
[61] latticeExtra_0.6-29    data.table_1.13.0      png_0.1-7              vctrs_0.3.4           
[65] gtable_0.3.0           purrr_0.3.4            tidyr_1.1.2            ggplot2_3.3.2         
[69] xfun_0.16              xtable_1.8-4           survival_3.2-3         tibble_3.0.3          
[73] AnnotationDbi_1.48.0   memoise_1.1.0          cluster_2.1.0          ellipsis_0.3.1        
```
