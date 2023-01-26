salsa juicer visualization v3
================
2020-11-27

  - [Index the draft fasta sequence using
    BWA](#index-the-draft-fasta-sequence-using-bwa)
  - [copy merged\_nodups.txt](#copy-merged_nodups.txt)
  - [Visualize candidate assembly](#visualize-candidate-assembly)
  - [copy and save](#copy-and-save)

using
/cluster/projects/nn9525k/arenosa\_genome/02\_scaffolds/05\_juicer/12\_arenosa\_purged\_sr\_salsa\_juicer\_NEW

"just need to run the visualization module. See page 5 …"

https://github.com/aidenlab/Juicebox

``` bash
srun --cpus-per-task=4 --mem-per-cpu=14G --partition=bigmem --time=03:00:00 --account=nn9525k --pty bash -i

cd /cluster/work/users/jonathbr/salsa/

mkdir -p juicer_v3
cd juicer_v3

ml BWA/0.7.17-GCC-8.2.0-2.31.1
ml SAMtools/1.9-GCC-8.2.0-2.31.1
ml Java/11.0.2
ml parallel/20190622-GCCcore-8.2.0
```

## Index the draft fasta sequence using BWA

``` bash
mkdir references
cd references

cp /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_salsa_juicer_NEW/reference/* .

#bwa index #done by Anders
cd ..
```

``` bash
cp -r /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_salsa_juicer_NEW/restriction_sites/ .
#awk '{print $1,$NF}' restriction_sites/*.txt > reference/chrom.size #done by Anders
```

## copy merged\_nodups.txt

``` bash
cp /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_salsa_juicer_NEW/merged_nodups.txt .
```

## Visualize candidate assembly

``` bash
awk -f /cluster/projects/nn9525k/Programs/3d-dna/utils/generate-assembly-file-from-fasta.awk references/scaffolds_FINAL.fasta > scaffolds_FINAL.fasta.assembly 

bash /cluster/projects/nn9525k/Programs/3d-dna/visualize/run-assembly-visualizer.sh scaffolds_FINAL.fasta.assembly merged_nodups.txt
```

![12\_arenosa\_purged\_sr\_salsa\_v3\_juciebox.png](12_arenosa_purged_sr_salsa_v3_juciebox.png)

## copy and save

``` bash
cp temp.*mnd.txt /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_salsa_juicer_NEW/assembly_hic
cp *superscaf_track.txt /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_salsa_juicer_NEW/assembly_hic
cp *scaffold_track.txt /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_salsa_juicer_NEW/assembly_hic

cp scaffolds_FINAL.fasta.assembly /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_salsa_juicer_NEW/assembly_hic
cp scaffolds_FINAL.fasta.hic /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_salsa_juicer_NEW/assembly_hic
```
