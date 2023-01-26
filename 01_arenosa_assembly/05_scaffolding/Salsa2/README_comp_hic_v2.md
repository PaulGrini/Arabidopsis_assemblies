
## compare_hic

``` bash
cd /cluster/work/users/jonathbr/
mkdir compare_hic
```

## add 00_ref_assemblies for 11 and 12_arenosa_purged_sr

``` bash
cd /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/00_ref_assemblies
mkdir 11_arenosa_pilon12
cp /cluster/projects/nn9525k/jonathan/arenosa_assembly/03_polishing_assembly/arenosa_pilon12.fasta.gz 11_arenosa_pilon12/
mkdir 12_arenosa_purged_sr
cp /cluster/projects/nn9525k/jonathan/arenosa_assembly/04_phasing/purge_results/arenosa_purged_sr.fasta.gz 12_arenosa_purged_sr/

cd /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/00_ref_assemblies/
for DIR in $(ls -d */ | tr -d "//")
do
cd $DIR
FILE=$(ls *fasta.gz)
echo "cp $FILE /cluster/work/users/jonathbr/compare_hic/${DIR}.fasta.gz"
cp $FILE /cluster/work/users/jonathbr/compare_hic/${DIR}.fasta.gz
cd .. 
done

cd /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/01_Salsa_scaffolds/

for DIR in $(ls -d */ | tr -d "//")
do
 cd $DIR
 for SUBDIR in $(ls -d */ | tr -d "//")
 do 
  cd $SUBDIR
  for SUBSUB in $(ls -d */ | tr -d "//")
  do
   cd $SUBSUB
   FILE=$(ls *fasta.gz)
   echo "cp -n $FILE /cluster/work/users/jonathbr/compare_hic/${DIR}_${SUBDIR}_${SUBSUB}.fasta.gz"
   cp -n $FILE /cluster/work/users/jonathbr/compare_hic/${DIR}_${SUBDIR}_${SUBSUB}.fasta.gz
   cd ..
  done
  cd ..
 done
 cd ..
done

cd /cluster/work/users/jonathbr/compare_hic
ls -l

for g in *.gz
do 
echo $g
gunzip $g
done 

sbatch quast-lg.slurm

```

## compare_hic_v2 02. November 2020

## arenosa_genome/02_scaffolds/03_purged_scaffolds

``` bash
cd /cluster/work/users/jonathbr/compare_hic
mkdir compare_purged
cd mkdir compare_purged
```

## add 02_Salsa_scaffolds_with_matlock

``` bash
cd /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/02_Salsa_scaffolds_with_matlock
for FILE in $(find . -maxdepth 2 -type f -name "*.fasta")
do 
 #echo $FILE 
 FASTA=$(basename $FILE)
 #echo $FASTA
 echo "
 cp -n $FILE /cluster/work/users/jonathbr/compare_hic/compare_purged/$FASTA"
 cp -n $FILE /cluster/work/users/jonathbr/compare_hic/compare_purged/$FASTA
done

cd /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/03_purged_scaffolds
for FILE in $(find . -maxdepth 2 -type f -name "*.fasta")
do 
 #echo $FILE 
 FASTA=$(basename $FILE)
 #echo $FASTA
 echo "
 cp -n $FILE /cluster/work/users/jonathbr/compare_hic/compare_purged/$FASTA"
 cp -n $FILE /cluster/work/users/jonathbr/compare_hic/compare_purged/$FASTA
done

cd /cluster/work/users/jonathbr/compare_hic/compare_purged
ll

cp ../quast-lg.slurm .
sbatch quast-lg.slurm
```
