#!/bin/bash

#SBATCH --job-name=CoverageCalculation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=CoverageCalculation_run_%j.log

cd ~/bigdata/TakeHomeFulgent
bam=NA12878.mapped.illumina.mosaik.CEU.exome.20110411.bam

module load samtools
module load bedtools
module load mosdepth

# samtools
start_time=$(date +%s)

samtools depth -a $bam > NA12878_coverage.txt
awk '{sum+=$3} END { print "Average coverage from samtools = ",sum/NR}' NA12878_coverage.txt

end_time=$(date +%s)

runtime=$((end_time - start_time))

echo "The samtools took $runtime seconds to execute."

# bedtools
start_time=$(date +%s)

samtools sort $bam | bedtools genomecov -ibam stdin -d > NA12878_genomecov.txt
awk '{sum+=$3} END { print "Average coverage from bedtools = ",sum/NR}' NA12878_genomecov.txt

end_time=$(date +%s)

runtime=$((end_time - start_time))

echo "The bedtools took $runtime seconds to execute."

# mosdepth
start_time=$(date +%s)

mosdepth NA12878 $bam
awk 'NR==1 {print} {last=$0} END {print last}' NA12878.mosdepth.summary.txt
awk '{last=$4} END {print "Average coverage from mosdepth = ",last}' NA12878.mosdepth.summary.txt

end_time=$(date +%s)

runtime=$((end_time - start_time))

echo "The mosdepth took $runtime seconds to execute."
