#!/bin/bash
#SBATCH --job-name=pdstep1
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
#SBATCH --qos=devel
#SBATCH --time=00:30:00
##SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END

module --quiet purge
module load minimap2/2.17-GCC-8.2.0-2.31.1
module load Python/3.7.2-GCCcore-8.2.0
#module load SAMtools/1.9-GCC-8.2.0-2.31.1
module list 

#purge_dups Pipeline
#https://github.com/dfguan/purge_dups#pplg

echo "Step 1. Run minimap2 to align pacbio data and generate paf files, then calculate read depth histogram and base-level read depth"

echo "make dir: purge_work"
mkdir purge_work
cd purge_work

REF=../petraea_pilon9/petraea_pilon9.fasta
READS=../petraea_pilon9/pb.fofn

echo "minimap2 -xmap-pb -t "$SLURM_CPUS_PER_TASK" $REF $READS | gzip -c - > pb.paf.gz"

minimap2 -xmap-pb -t "$SLURM_CPUS_PER_TASK" $REF $READS | gzip -c - > pb.paf.gz

echo "produces PB.base.cov and PB.stat files"
/cluster/projects/nn9525k/Programs/purge_dups/bin/pbcstat *.paf.gz 

/cluster/projects/nn9525k/Programs/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log

echo "Step 1. Split an assembly and do a self-self alignment."
/cluster/projects/nn9525k/Programs/purge_dups/bin/split_fa $REF > pri_asm.split
minimap2 -xasm5 -t "$SLURM_CPUS_PER_TASK" -DP pri_asm.split pri_asm.split | gzip -c - > pri_asm.split.self.paf.gz

echo "make the histogram plot"
/cluster/projects/nn9525k/Programs/purge_dups/scripts/hist_plot.py -c cutoffs PB.stat PB_hist

echo "Step 2. Purge haplotigs and overlaps"

/cluster/projects/nn9525k/Programs/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

echo "Step 3. Get purged primary and haplotig sequences from draft assembly."

/cluster/projects/nn9525k/Programs/purge_dups/bin/get_seqs dups.bed "$REF"