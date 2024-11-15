---
title: "Coverage Calculation"
author: "QiongJia"
date: "2024-11-15"
output: 
  html_document:
    highlight: pygments
    toc: true
    keep_md: true
---



## Design  

To calculate the average coverage(sequencing depth) for the given sample BAM file, two approach are presented below.  

1. The first one is based on the coverage calculation equation:   

$$C = \frac{LN}{G}$$  

where    
 - $C$ is coverage.  
 - $G$ is the haploid genome length.  
 - $L$ is the read length in the sequencing.  
 - $N$ is the number of reads.  
 
   The haploid reference genome (GRCh38) size is $3,099,734,149$ bases ([Reference](https://www.ncbi.nlm.nih.gov/grc/human/data)). The read length and number can be obtained from ***.bam.bas*** file. 

2. The second one is more precise. The ***samtool depth*** can be used to calculate the coverage at each genomic position and the average coverage of the given BAM file.  Same as ***bedtools genomecov***. 

3. At the end, one another tools ***mosdepth*** is presented briefly in coverage calculation.   

## Download BAM file 


``` bash
# make working directory
cd ~
mkdir TakeHomeFulgent
cd TakeHomeFulgent

# download cram file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram

# download cram index file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram.crai

# download .bam.bas file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bam.bas

```


``` bash
# define work variable
cram=NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram
bas=NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bam.bas
```

## Quick estimate 


``` bash
# reference genome size
G=3099734149
# read length equals total bases number / total reads number
L=$(awk 'NR==2 {print $9 / $11}' "$bas")
# reads number equals total reads number - duplicate reads number
N=$(awk 'NR==2 {print $11-$21}' "$bas")
# coverage calculation
C=$(echo "$N*$L/$G" | bc -l)
printf "The average coverage is: %.2f\n", "$C"
```

The average coverage is: **5.77**  

## Coverage calculation for each position 

### samtools  

**Step1**: calculate the coverage at each genomic position. 


``` bash
samtools depth -a $cram > NA12878_cram_coverage.txt
```

Overlook the coverage output.   


``` bash
head -n 5 NA12878_cram_coverage.txt
```
  

```
##     V1 V2 V3
## 1 chr1  1  0
## 2 chr1  2  0
## 3 chr1  3  0
## 4 chr1  4  0
## 5 chr1  5  0
```
  
Each line represents a genomic position. Three columns are included int he coverage output:   

- Chromosome;  
- Position;  
- Reads covered this position.   

**Step2**: calculate the average coverage. 


``` bash
awk '{sum+=$3} END { print "Average coverage from samtools = ",sum/NR}' NA12878_cram_coverage.txt
```

Average coverage from samtools = **5.14261**

### bedtools  

***bedtools genomecov -d*** also reports the genome coverage per base as below:  


``` bash
# To use -ibam flag in bedtools genomecov, the bam file is needed to be sorted by position
samtools sort $cram | bedtools genomecov -ibam stdin -d > NA12878_cram_genomecov.txt
```

Then the average coverage would be: 


``` bash
awk '{sum+=$3} END { print "Average coverage from bedtools = ",sum/NR}' NA12878_cram_genomecov.txt
```

Average coverage from bedtools = **5.59374**


## Other avaialble tools 

***mosdepth*** can report coverage for both per-base and summary result at the same time. It can take **BAM** file directly, but it needs reference genome file with **CRAM** file as input. 


``` bash
mosdepth NA12878 $cram
```

The file ended with ***.mosdepth.summary.txt*** contain the average coverage result.

Besides, It can also report coverage based on the user defined region by using ***--by <bed|window>***. 

## Runing time comparison  

- *samtools depth* :  1176s;
- *bedtools genomecov* : 32812s;


```
## R version 4.2.3 (2023-03-15)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] digest_0.6.37     R6_2.5.1          lifecycle_1.0.4   jsonlite_1.8.8   
##  [5] evaluate_0.24.0   cachem_1.1.0      rlang_1.1.4       cli_3.6.3        
##  [9] rstudioapi_0.15.0 jquerylib_0.1.4   bslib_0.8.0       rmarkdown_2.28   
## [13] tools_4.2.3       xfun_0.47         yaml_2.3.10       fastmap_1.2.0    
## [17] compiler_4.2.3    htmltools_0.5.8.1 knitr_1.48        sass_0.4.9
```
