#!/bin/bash
#SBATCH --job-name=pdstep123
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
##SBATCH --qos=devel
##SBATCH --time=00:30:00
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END

module --quiet purge
module load minimap2/2.17-GCC-8.2.0-2.31.1
module load Python/3.7.2-GCCcore-8.2.0
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load BWA/0.7.17-GCC-8.2.0-2.31.1
module list

#purge_dups Pipeline
#https://github.com/dfguan/purge_dups#pplg

#echo "Step 1. Run minimap2 to align pacbio data and generate paf files, then calculate read depth histogram and base-level read depth"
echo "Q: Can purge_dups work with Illumina short reads?"

#bwa mem $pri_asm $sr_1.fq $sr_2.fq | samtools view -b -o - > $sr.bam
#ngscstat $sr.bam... # The program will generate two/three outputs, TX.stat and TX.base.cov which functions the same way as PB.stat and PB.base.cov respectively.

echo "make dir: purge_work"
mkdir purge_work
cd purge_work

ASSEMBLY=../arenosa_pilon12.fasta
READ1=../Sample_06-F-2Az1-9_R1_rep.fastq
READ2=../Sample_06-F-2Az1-9_R2_rep.fastq
PREFIX=$(basename "$ASSEMBLY" .fasta )

echo "
bwa index "$ASSEMBLY""
bwa index "$ASSEMBLY"

echo "bwa mem -t "$SLURM_CPUS_PER_TASK" "$ASSEMBLY" "$READ1" "$READ2" | samtools view -b -o - > "${PREFIX}_bwa.bam""
bwa mem -t "$SLURM_CPUS_PER_TASK" "$ASSEMBLY" "$READ1" "$READ2" | samtools view -b -o - > "${PREFIX}_bwa.bam"

echo "The program will generate two/three outputs, TX.stat and TX.base.cov"
/cluster/projects/nn9525k/Programs/purge_dups/bin/ngscstat "${PREFIX}_bwa.bam"

/cluster/projects/nn9525k/Programs/purge_dups/bin/calcuts TX.stat > cutoffs 2>calcults.log

echo "Step 1. Split an assembly and do a self-self alignment."
/cluster/projects/nn9525k/Programs/purge_dups/bin/split_fa  "$ASSEMBLY" > pri_asm.split
minimap2 -xasm5 -t "$SLURM_CPUS_PER_TASK" -DP pri_asm.split pri_asm.split | gzip -c - > pri_asm.split.self.paf.gz

echo "make the histogram plot"
/cluster/projects/nn9525k/Programs/purge_dups/scripts/hist_plot.py -c cutoffs TX.stat TX_hist

echo "Step 2. Purge haplotigs and overlaps"

/cluster/projects/nn9525k/Programs/purge_dups/bin/purge_dups -2 -T cutoffs -c TX.base.cov pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

echo "Step 3. Get purged primary and haplotig sequences from draft assembly."

/cluster/projects/nn9525k/Programs/purge_dups/bin/get_seqs dups.bed "${ASSEMBLY}"
