# Coverage Calculation

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
  

## Example file

The **CRAM** and **.bam.bas** file of NA12878 from 1000 Genomes on GRCh38 is used for demonstration in [Coverage_Calulation.md](https://github.com/qj009/CoverageCal/blob/main/Coverage_Calulation.md). The source code of coverage calculation is provided in [RunningTimeCram.sh](https://github.com/qj009/CoverageCal/blob/main/RunningTimeCram.sh). 
