# PBJelly_v3
```{bash, eval=FALSE}
## Reads 
zcat /cluster/projects/nn9525k/arenosa_genome_pb_uncorr/canu_out/arenosa_pp_uncorr_assembly.correctedReads.fasta.gz > pbReads.fasta

## Assembly
zcat /cluster/work/users/jonathbr/arenosa_assembly/05_scaffolding/ragtag/arenosa_ragtag.fasta.gz > arenosa.fasta.bak

../PBJelly_v2/simplifyFastaHeaders.pl .
./simplifyFastaHeaders.pl arenosa.fasta.bak arenosa arenosa.fasta header.map

mkdir jellyout

```

```{bash, eval=FALSE}
split -l 30000 pbReads.fasta pbReads.split.

for f in pbReads.split.* ; do  mv $f ${f}.fasta ; done
for f in pbReads.split.* ; do  echo "<job>${f}</job>" ; done

```

```
------- Protocol.xml
<jellyProtocol>
    <reference>/cluster/work/users/jonathbr/arenosa_assembly/06_gap_filling/PBJelly_v3/arenosa.fasta</reference>
    <outputDir>/cluster/work/users/jonathbr/arenosa_assembly/06_gap_filling/PBJelly_v3/jellyout/</outputDir>
    <blasr>--minMatch 8 --sdpTupleSize 8 --minPctSimilarity 75 --bestn 1 --nCandidates 10 --maxScore -500 --noSplitSubreads</blasr>
    <input baseDir="/cluster/work/users/jonathbr/arenosa_assembly/06_gap_filling/PBJelly_v3/">
        <job>pbReads.split.aa.fasta</job>
		<job>pbReads.split.ab.fasta</job>
		<job>pbReads.split.ac.fasta</job>
		<job>pbReads.split.ad.fasta</job>
		<job>pbReads.split.ae.fasta</job>
		<job>pbReads.split.af.fasta</job>
		<job>pbReads.split.ag.fasta</job>
		<job>pbReads.split.ah.fasta</job>
		<job>pbReads.split.ai.fasta</job>
		<job>pbReads.split.aj.fasta</job>
		<job>pbReads.split.ak.fasta</job>
		<job>pbReads.split.al.fasta</job>
		<job>pbReads.split.am.fasta</job>
		<job>pbReads.split.an.fasta</job>
		<job>pbReads.split.ao.fasta</job>
		<job>pbReads.split.ap.fasta</job>
		<job>pbReads.split.aq.fasta</job>
		<job>pbReads.split.ar.fasta</job>
		<job>pbReads.split.as.fasta</job>
		<job>pbReads.split.at.fasta</job>
		<job>pbReads.split.au.fasta</job>
		<job>pbReads.split.av.fasta</job>
    </input>
</jellyProtocol>
----------
```

## add fake Phred score to fasta

```{bash, eval=FALSE}
conda activate /cluster/projects/nn9525k/Programs/conda_envs/PBJelly

##fakeQuals.py data/reference/lambda.fasta data/reference/lambda.qual
fakeQuals.py arenosa.fasta arenosa.qual

# fix read names # Sequence names CANNOT have spaces.
#cat arenosa_pb_uncorr_assembly.correctedReads.fasta.original | sed 's/ id=.*$//' > arenosa_pb_uncorr_assembly.correctedReads.fasta

for FASTA in pbReads.split.*.fasta;
do 
 PREFIX=$(basename $FASTA .fasta)
 mv $FASTA ${FASTA}.original
 echo -e "\n
 cat ${FASTA}.original | sed 's/ id=.*$//' > ${FASTA}"
 cat ${FASTA}.original | sed 's/ id=.*$//' > ${FASTA}
 echo -e "\n
 fakeQuals.py $FASTA ${PREFIX}.qual"
 fakeQuals.py $FASTA ${PREFIX}.qual
done

```

```{bash, eval=FALSE}
summarizeAssembly.py arenosa.fasta > sum_arenosa.txt

```
## Jelly.py setup
```{bash, eval=FALSE}
Jelly.py setup Protocol.xml
```
## Jelly.py mapping
```{bash, eval=FALSE}
sbatch jelly_mapping.slurm #job 1239739 25-19:36:+ 1-09:37:31 12817260
```

## Support The Gaps & Extract Useful Reads as slurm script

```{bash, eval=FALSE}
# Jelly.py support Protocol.xml #and 
# Jelly.py extraction Protocol.xml #included 
sbatch jelly_support.slurm #8G job 1248294 00:06:46   00:09:11 3452984K
```

## Assemble The Gaps
```{bash, eval=FALSE}
#Jelly.py assembly Protocol.xml -x "--nproc=16"
sbatch jelly_assembly.slurm  #job 1249214 3-05:23:14   11:51:58 41478248K
```

## Output Your Results

```{bash, eval=FALSE}

Jelly.py output Protocol.xml

grep -c "^>" jellyout/jelly.out.fasta
# 462

grep -Ho N jellyout/jelly.out.fasta | uniq -c
# 12409 jellyout/jelly.out.fasta:N

summarizeAssembly.py jellyout/jelly.out.fasta > sum_arenosa_jelly.txt
```
## tar copy and save

```{bash, eval=FALSE}
du -sh .
#46G

#remove input 
rm *.fasta*
rm *.qual

rm core.*

tar --remove-files -czvf slurm.base.tar.gz slurm-*

cp jellyout/jelly.out.fasta arenosa_jelly.fasta
gzip arenosa_jelly.fasta

tar --remove-files -czvf jelly_out.tar.gz jellyout


du -sh . #171M

cd ..
rm -r PBJelly/
rm -r PBJelly_v2/
rm -r PBJelly_test/

cp -rv /cluster/work/users/jonathbr/arenosa_assembly/06_gap_filling /cluster/projects/nn9525k/jonathan/arenosa_assembly/
```
## continue with
```{bash, eval=FALSE}
wc /cluster/projects/nn9525k/jonathan/arenosa_assembly/06_gap_filling/PBJelly_v3/arenosa_jelly.fasta.gz
#160696   875776 41798295 /cluster/projects/nn9525k/jonathan/arenosa_assembly/06_gap_filling/PBJelly_v3/arenosa_jelly.fasta.gz
```

## Extras

cp (https://github.com/esrice/PBJelly) .

```{bash, eval=FALSE}
blasrToBed.py
This script will convert blasr's .m4 or .m5 format into a
BED Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format1 )
If you would like to visualize the alignments, I
reccommend using IGB ( http://bioviz.org/igb/index.html ).

bedToCoverageWig.py
Turn a bed file with alignments into a depth of coverage
WIG Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format6 ).
```
