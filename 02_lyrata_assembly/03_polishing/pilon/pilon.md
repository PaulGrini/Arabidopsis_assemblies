# pilon_petraea_canu

``` bash
mkdir pilon_petrea_canu
cd pilon_petrea_canu
mkdir data
cd data

cp /cluster/projects/nn9525k/hybrids/jasonrm/GenomeAsssembly/A.lyrata.ver1/lyrata.contigs.fasta petrea.contigs.fasta

ln -s /cluster/projects/nn9525k/hybrids/molbar_illumina_DNA/trimmed_reads/Sample_02-B-2Lz3-4/Sample_02-B-2Lz3-4_R?_rep.fastq.gz .

cd ..

sbatch pilon_loop.slurm #job 810515 1-19:51:49   06:18:01 59 617 972K
```

## input: pilon_round3.fasta

mapping changed to bowtie2 default without --score-min 'C,0,-1'

``` bash
sbatch pilon_loop2.slurm #job 811487 3-04:23:01   10:32:39   67 088 900K
```

## input: pilon_round6.fasta

mapping changed to `bowtie2 --threads ${THREADS} --no-unal --sensitive --local -x ${PREFIX}_round${ROUND} -1 ../${R1} -2 ../${R2} -S ${PREFIX}_round${ROUND}.sam`

``` bash
sbatch pilon_loop3.slurm #job 815787
```

``` bash
mkdir summary
cd summary/
cp ../bowtie_round?/*.flagstats .
cp ../bowtie_round?/*.flagstat .
cp ../bowtie_round?/*.stats.SN .
cp ../pilon_round?/*.changes .
cd ..

ls pilon_round3/
mkdir quast
cd quast/
ln -s ../pilon_round?/*_round?.fasta .
ln -s ../data/petrea.contigs.fasta .
cp ../../QUAST/quast.slurm .
sbatch quast.slurm
```
