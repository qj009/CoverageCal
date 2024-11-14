# CoverageCal

## Design  

To calculate the average coverage(sequencing depth) for the given sample BAM file, two approach are presented below.  

1. The first one is based on the coverage calculation equation:   
$$ C = \frac{LN}{G} $$
where 
 - $C$ is coverage.  
 - $G$ is the haploid genome length. 
 - $L$ is the read length in the sequencing. 
 - $N$ is the number of reads. 

    The *samtool idxstats* can be used to get chromosome lengths and number of mapped reads.

2. The second one is more precise. The *samtool depth* can be used to calculate the coverage at each genomic position and the average coverage of the given BAM file.  

3. At the end, two other tools *bedtools genomecov* and *mosdepth* are presented briefly in coverage calculation.

The full report can be found here[]
