salsa juicer visualization v2
================
2020-11-24

  - [Index the draft fasta sequence using
    BWA](#index-the-draft-fasta-sequence-using-bwa)
  - [set up a working directory](#set-up-a-working-directory)
  - [Visualize candidate assembly](#visualize-candidate-assembly)
  - [copy and save](#copy-and-save)

"just need to run the visualization module. See page 5 …"

https://github.com/aidenlab/Juicebox

``` bash
srun --cpus-per-task=4 --mem-per-cpu=14G --partition=bigmem --time=03:00:00 --account=nn9525k --pty bash -i

cd /cluster/work/users/jonathbr/salsa/12_arenosa_purged_sr_sala

mkdir -p juicer_v2
cd juicer_v2

ml BWA/0.7.17-GCC-8.2.0-2.31.1
ml SAMtools/1.9-GCC-8.2.0-2.31.1
ml Java/11.0.2
ml parallel/20190622-GCCcore-8.2.0
```

## Index the draft fasta sequence using BWA

``` bash
mkdir references
cd references

cp /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/04_scaffolds_purge_dups/12_arenosa_purged_sr_salsa_purge_dups_purged.fasta .

bwa index 12_arenosa_purged_sr_salsa_purge_dups_purged.fasta
cd ..
```

``` bash
cp /cluster/projects/nn9525k/Programs/juicer/misc/generate_site_positions.py .
chmod u+x generate_site_positions.py

#./generate_site_positions.py <restriction enzyme> <genome> [location]
./generate_site_positions.py Arima arenosa_purged_sr references/12_arenosa_purged_sr_salsa_purge_dups_purged.fasta  #creates arenosa_purged_sr_Arima.txt 

mkdir restriction_sites
mv arenosa_purged_sr_Arima.txt restriction_sites/
```

## set up a working directory

``` bash
cp /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_sala_juicer/merged_nodups.txt .
```

``` bash
awk '{print $1,$NF}' restriction_sites/*.txt > chrom.size
```

## Visualize candidate assembly

``` bash
wget https://github.com/aidenlab/3d-dna/raw/master/utils/generate-assembly-file-from-fasta.awk

awk -f generate-assembly-file-from-fasta.awk references/12_arenosa_purged_sr_salsa_purge_dups_purged.fasta > 12_arenosa_purged_sr_salsa_purge_dups_purged.fasta.assembly 

bash /cluster/projects/nn9525k/Programs/3d-dna/visualize/run-assembly-visualizer.sh 12_arenosa_purged_sr_salsa_purge_dups_purged.fasta.assembly merged_nodups.txt
```

![12\_arenosa\_purged\_sr\_salsa\_purge\_dups\_purged\_juciebox.png](12_arenosa_purged_sr_salsa_purge_dups_purged_juciebox.png)

## copy and save

``` bash
cp temp.12_arenosa_purged_sr_salsa_purge_dups_purged.fasta.asm_mnd.txt /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_sala_juicer/assembly_hic
cp 12_arenosa_purged_sr_salsa_purge_dups_purged.fasta_asm.scaffold_track.txt /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_sala_juicer/assembly_hic
cp 12_arenosa_purged_sr_salsa_purge_dups_purged.fasta_asm.superscaf_track.txt /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_sala_juicer/assembly_hic

cp 12_arenosa_purged_sr_salsa_purge_dups_purged.fasta.assembly /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_sala_juicer/assembly_hic
cp 12_arenosa_purged_sr_salsa_purge_dups_purged.fasta.hic /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/05_juicer/12_arenosa_purged_sr_sala_juicer/assembly_hic
```
