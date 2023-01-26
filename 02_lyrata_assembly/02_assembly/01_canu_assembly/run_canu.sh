#!/bin/sh

module load Java/1.8.0_212 
module load gnuplot/5.2.6-GCCcore-8.2.0

# CANU will restart where it left off.
# This can be confusing. You fix an option and restart, and the same option gets used.
# If you change this script, consider running in a new directory, possibly with -d option.

# cd /cluster/home/jasonrm/work/GenomeAssembly/A.arenosa/canu_asm1
export CANU=/cluster/home/jasonrm/Source/CANU/canu-2.0/Linux-amd64/bin/canu
export PACBIO=/cluster/projects/nn9525k/lyrata_genome/01_data/lyrata_pacbio.fastq.gz

echo $PACBIO
echo $NANOPORE

#$CANU -correct   genomeSize=250m \
# gridOptions="--account=${ACCOUNT} --job-name=LYRA --time=96:00:00 --mem-per-cpu=16G --cpus-per-task=4" \
#  -p lyrata -d ver1  -pacbio  $PACBIO 

#$CANU -trim   genomeSize=250m \
# gridOptions="--account=${ACCOUNT} --job-name=LYRA --time=96:00:00 --mem-per-cpu=16G --cpus-per-task=4" \
#  -p lyrata -d ver1  -pacbio  $PACBIO 

$CANU -assemble   genomeSize=250m \
 gridOptions="--account=${ACCOUNT} --job-name=LYRA --time=96:00:00 --mem-per-cpu=16G --cpus-per-task=4" \
  -p lyrata -d ver1  -pacbio  $PACBIO 

