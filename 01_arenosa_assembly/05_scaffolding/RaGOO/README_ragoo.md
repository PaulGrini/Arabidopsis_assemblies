RaGOO for Assemblytics
================
2020-09-08

https://github.com/malonge/RaGOO/wiki/Calling-Structural-Variants

## setup dir

```{bash, eval=FALSE}
cd /cluster/work/users/jonathbr/arenosa_assembly
mkdir -p 05_scaffolding/ragoo
cd 05_scaffolding/ragoo
```

## links

```{bash, eval=FALSE}
## Reference
mkdir -p links

ln -s /cluster/projects/nn9525k/lyrata_genome_phytomzome/Alyrata_384_v1.fa links/Alyrata_384_v1.fasta

## Reads 
ln -s /cluster/projects/nn9525k/arenosa_genome_pb_uncorr/canu_out/arenosa_pp_uncorr_assembly.correctedReads.fasta.gz links

ln -s /cluster/projects/nn9525k/arenosa_genome/Illumina_data/Sample_06-F-2Az1-9_R?_rep.fastq.gz links

## Assembly
zcat /cluster/work/users/jonathbr/arenosa_assembly/04_phasing/purge_results/arenosa_purged_sr.fasta.gz > links/arenosa_purged_sr.fasta
```

## RaGOO
```{bash, eval=FALSE}
sbatch ragoo.slurm
```
## tar copy and save

see [README_ragtag.md](../RagTag/README_ragtag.md)



