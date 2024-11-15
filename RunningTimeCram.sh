#!/bin/bash

#SBATCH --job-name=CoverageCalculation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=CoverageCalculation_run_%j.log

cd ~/bigdata/TakeHomeFulgent
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram.crai
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bam.bas
cram=NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram

module load samtools
module load bedtools
module load mosdepth

# quick calculation

# reference genome size
G=3099734149
# read length equals total bases number / total reads number
L=$(awk 'NR==2 {print $9 / $11}' "$bas")
# reads number equals total reads number - duplicate reads number
N=$(awk 'NR==2 {print $11-$21}' "$bas")
# coverage calculation
C=$(echo "$N*$L/$G" | bc -l)
printf "The average coverage is: %.2f\n", "$C"

# samtools
start_time=$(date +%s)

samtools depth -a $cram > NA12878_cram_coverage.txt
awk '{sum+=$3} END { print "Average coverage from samtools = ",sum/NR}' NA12878_cram_coverage.txt

end_time=$(date +%s)

runtime=$((end_time - start_time))

echo "The samtools took $runtime seconds to execute."

# bedtools
start_time=$(date +%s)

samtools sort $cram | bedtools genomecov -ibam stdin -d > NA12878_cram_genomecov.txt
awk '{sum+=$3} END { print "Average coverage from bedtools = ",sum/NR}' NA12878_cram_genomecov.txt

end_time=$(date +%s)

runtime=$((end_time - start_time))

echo "The bedtools took $runtime seconds to execute."

# mosdepth
start_time=$(date +%s)

mosdepth NA12878_cram $cram # error: need reference genome file.
awk 'NR==1 {print} {last=$0} END {print last}' NA12878_cram.mosdepth.summary.txt
awk '{last=$4} END {print "Average coverage from mosdepth = ",last}' NA12878_cram.mosdepth.summary.txt

end_time=$(date +%s)

runtime=$((end_time - start_time))

echo "The mosdepth took $runtime seconds to execute."
