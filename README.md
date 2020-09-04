# AML-CRNDE-Analysis
Data analysis for LncRNA CRNDE in AML

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

stringtie Input.bam -p 10 -G Homo_sapiens.GRCh38.84.gtf -eB -o Input.gtf

```

