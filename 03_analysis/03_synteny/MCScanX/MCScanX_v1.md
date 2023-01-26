MCScanX Arabidopsis
================
2021-03-03

  - [dowload and install MCScanX](#dowload-and-install-mcscanx)
-   [setup directory](#setup-directory)
    -   [copy GFFs and proteins](#copy-gffs-and-proteins)
-   [pepare files - change to .bed
    format](#pepare-files---change-to-.bed-format)
    -   [simplify gff](#simplify-gff)
        -   [read acronym.csv file](#read-acronym.csv-file)
        -   [walk input](#walk-input)
    -   [simplify .peps](#simplify-.peps)
        -   [rename .peps](#rename-.peps)
        -   [single line & shorten fasta
            header](#single-line-shorten-fasta-header)
    -   [compare fasta to bed header](#compare-fasta-to-bed-header)
        -   [small fixes](#small-fixes)
        -   [compare gene count](#compare-gene-count)
    -   [rename files](#rename-files)
        -   [sort\_bed function](#sort_bed-function)
        -   [detect splice variants in
            .pep](#detect-splice-variants-in-.pep)
-   [test plots](#test-plots)

``` r
library(tidyverse)
#library(seqinr)
#library(biostrings)
#library(rtracklayer)
```

## dowload and install MCScanX

``` bash
cd /cluster/projects/nn9525k/Programs/
mkdir MCScanX
wget http://chibba.pgml.uga.edu/mcscan2/MCScanX.zip
#wget http://chibba.pgml.uga.edu/mcscan2/transposed/MCScanX-transposed.zip
unzip MCScanX.zip
rm MCScanX.zip
module load Java/11.0.2
make 

MCS=/cluster/projects/nn9525k/Programs/MCScanX/
$MCS/MCScanX
```

used as guidance 

<https://github.com/zhaotao1987/SynNet-Pipeline>

<https://github.com/Jnthnoaa/MADS_synteny_network/blob/main/README_SynNet_v2.md>

# setup directory

-   Whole genome protein files in fasta format.
-   GFF/BED file for each genome

### copy GFFs and proteins

``` bash
cd /cluster/work/users/jonathbr/
mkdir -p mcscanx/gff
mkdir -p mcscanx/pep
cd mcscanx

# petraea 
cp /cluster/projects/nn9525k/jonathan/petraea_curated/braker_curated_sort.gtf.gz gff/petraea_curated.gtf.gz
cp /cluster/projects/nn9525k/jonathan/petraea_curated/braker_curated_sort.aa.gz pep/petraea_curated.aa.gz

# arenosa
cp /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/13_braker_arenosa_remasked/braker_prot/augustus.hints.gtf gff/arenosa_remasked.gtf
cp /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/13_braker_arenosa_remasked/braker_prot/augustus.hints.aa pep/arenosa_remasked.aa

# thaliana
cp /cluster/projects/nn9525k/mads_paper/01_genomes/Araport11/annotation/Athaliana_447_Araport11.gene.gff3.gz gff/thaliana.gff3.gz
cp /cluster/projects/nn9525k/mads_paper/01_genomes/Araport11/annotation/Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa.gz pep/thaliana.aa.gz

# lyrata reference 
cp /cluster/projects/nn9525k/mads_paper/01_genomes/01_included/GCF_000004255.2_v.1.0_Alyrata/GCF_000004255.2_v.1.0_genomic.gff.gz gff/lyrata.gff.gz
cp /cluster/projects/nn9525k/mads_paper/01_genomes/01_included/GCF_000004255.2_v.1.0_Alyrata/GCF_000004255.2_v.1.0_protein.faa.gz pep/lyrata.aa.gz

for GZ in */*.gz 
do
 echo $GZ
 gunzip $GZ &
done 
wait 
```

# pepare files - change to .bed format

Index your genome files, using 3-5 letters, for example for “Arabidopis
thaliana”, rename genome file and gff file as “Ath.pep” and “Ath.bed”

``` bash
module purge
module load R/3.6.2-foss-2019b
R
```

## simplify gff

### read acronym.csv file

``` r
library(tidyverse)

input <- read_csv2("acronym_short_v2.csv")
```

#### get\_cds\_name function

``` r
#gff_t1 %>% print(n=40)
#mrna_id = "XM_013814943.2"
#str_detect(mrna_id, "^XM_")
get_cds_name  <- function(gff, mrna_id) {
    gff %>% 
      filter(type == "CDS", parent == str_c("rna-", mrna_id)) %>%
      pull(name) %>% first()
}
#get_cds_name(gff = gff_t1, mrna_id = mrna_id)
```

#### read\_gff3 function

``` r
read_gff3 <- function(gff_file) {
  read_tsv(file = gff_file,
           comment = "#",
           col_names = c("seqid", "source", "type", "start", "end",
                         "score", "strand", "phase", "attributes")) %>%
    mutate(id = str_extract(attributes, "ID=.*?;") %>% str_remove_all("ID=|;")) %>%
    mutate(name = str_extract(attributes, "Name=.*?;") %>% str_remove_all("Name=|;")) %>%
    mutate(longest = str_extract(attributes, "longest=.*?;") %>% str_remove_all("longest=|;")) %>%
    mutate(parent = str_extract(attributes, "Parent=.*?;") %>% str_remove_all("Parent=|;")) %>%
    mutate(pacid = str_extract(attributes, "pacid=.*?;") %>% str_remove_all("pacid=|;"))
}
```

#### read\_gff\_simple function

``` r
read_gff_simple <- function(gff_file) {
  read_tsv(file = gff_file,
           comment = "#",
           col_names = c("seqid", "source", "type", "start", "end",
                         "score", "strand", "phase", "attributes")) %>%
    mutate(id = str_extract(attributes, "ID=.*?;") %>% str_remove_all("ID=|;")) %>%
    mutate(name = str_extract(attributes, "Name=.*?;") %>% str_remove_all("Name=|;")) %>%
    mutate(longest = str_extract(attributes, "longest=.*?;") %>% str_remove_all("longest=|;")) %>%
    mutate(parent = str_extract(attributes, "Parent=.*?;") %>% str_remove_all("Parent=|;")) %>%
    mutate(pacid = str_extract(attributes, "pacid=.*?;") %>% str_remove_all("pacid=|;"))
}
```

#### read\_gff function

``` r
read_gff <- function(gff_file) {
  read_tsv(file = gff_file,
           comment = "#",
           col_names = c("seqid", "source", "type", "start", "end",
                         "score", "strand", "phase", "attributes")) %>%
    mutate(attributes_2 = str_remove_all(attributes, "\"")) %>% 
    mutate(transcript_id = str_extract(attributes_2, "transcript_id .*?; ") %>% str_remove_all("transcript_id |; ")) %>%
    mutate(gene_id = str_extract(attributes_2, "gene_id .*?;") %>% str_remove_all("gene_id |;"))
}
```

#### read\_gtf function

``` r
read_gtf <- function(gff_file) {
  read_tsv(file = gff_file,
           comment = "#",
           col_names = c("seqid", "source", "type", "start", "end",
                         "score", "strand", "phase", "attributes")) %>%
    mutate(attributes_2 = str_remove_all(attributes, "\\\\|\"")) %>% 
    mutate(transcript_id = str_extract(attributes_2, "transcript_id .*?; ") %>% str_remove_all("transcript_id |; ")) %>%
    mutate(gene_id = str_extract(attributes_2, "gene_id .*?;") %>% str_remove_all("gene_id |;"))
}
```

#### simplify\_gff3 function

``` r
simplify_gff3 <- function(gff, prefix) {
  gff %>% 
    mutate(seqid_2 = str_remove(seqid, "\\.\\d$")) %>%
    mutate(chrom = str_c(prefix, str_extract(seqid_2,"\\d+$"), sep = "_"), .before = seqid) %>% 
    filter(type == "mRNA") %>%  #gene
    mutate(gene_id = str_c(prefix, name, sep = "__"), .before = name) %>% #id
    select(chrom, start, end, gene_id)
}
#mutate(name_2 == if_else(str_detect(name, "^XM_"), 
#                              paste0(get_cds_name(gff, name)), 
#                              name)) %>% 
```

#### simplify\_gff3\_1.21(gff, prefix) function

``` r
simplify_gff3_1.21 <- function(gff, prefix) {
  gff %>% 
    mutate(seqid_2 = str_remove(seqid, "\\.\\d$")) %>%
    mutate(chrom = str_c(prefix, str_extract(seqid_2,"\\d+$"), sep = "_"), .before = seqid) %>% 
    filter(type %in% c("mRNA", "CDS")) %>% 
    mutate(protein_id = if_else(str_c(id) == lead(parent),
                                lead(name), NA_character_)) %>% 
    mutate(gene_id = str_c(prefix, protein_id, sep = "__"), .before = name) %>%
    filter(type %in% c("mRNA")) %>% 
    select(chrom, start, end, gene_id)
}
```

#### simplify\_gff\_simple function

``` r
simplify_gff_simple <- function(gff, prefix) {
  gff %>% 
    mutate(seqid_2 = str_remove(seqid, "\\.\\d$")) %>%
    mutate(chrom = str_c(prefix, str_extract(seqid_2,"\\d+$"), sep = "_"), .before = seqid) %>% 
    filter(type == "mRNA") %>%  #gene
    mutate(gene_id = str_c(prefix, id, sep = "__"), .before = id) %>% #id
    select(chrom, start, end, gene_id)
}
#mutate(name_2 == if_else(str_detect(name, "^XM_"), 
#                              paste0(get_cds_name(gff, name)), 
#                              name)) %>% 
```

#### simplify\_gff function

``` r
simplify_gff <- function(gff, prefix) {
  gff %>%
    mutate(seqid_2 = str_remove(seqid, "\\.\\d$")) %>%
    mutate(chrom = str_c(prefix, str_extract(seqid_2,"\\d+$"), sep = "_"), .before = seqid) %>% 
    mutate(gene_id = str_c(prefix,attributes, sep = "__"), .before = attributes) %>% 
    filter(type == "transcript") %>%  #double check attributes or transcript_id????
    select(chrom, start, end, gene_id)
}
```

#### simplify\_gtf function

``` r
simplify_gtf <- function(gff, prefix) {
  gff %>%
    mutate(seqid_2 = str_remove(seqid, "\\.\\d$")) %>%
    mutate(chrom = str_c(prefix, str_extract(seqid_2,"\\d+$"), sep = "_"), .before = seqid) %>% 
    mutate(gene_id = str_c(prefix,transcript_id, sep = "__"), .before = transcript_id) %>% 
    filter(type == "transcript") %>% 
    select(chrom, start, end, gene_id)
}
```

#### write\_short\_bed function

``` r
write_short_bed <- function(species_name) {
  print(species_name)
  gff_raw <- str_c("gff/", input$gff_raw[input$species == species_name] %>% str_remove(".gz"))
  prefix <- input$acronym[input$species == species_name]
  gff_version <- read_lines(gff_raw, n_max = 2)
  if (str_sub(gff_version[1], start = 1, end=1) != "#") {
    print("simple_gff_version------------------")
    gff <- read_gff_simple(gff_raw)
    print(gff)
    bed <- simplify_gff_simple(gff, prefix)
    print(bed)
    write_tsv(bed, str_c("short_bed/", prefix, ".bed"), col_names = FALSE) 
  } else if (gff_version[1] == "##gff-version 3") {
    gff <- read_gff3(gff_raw)
    print(gff)
    if (gff_version[2] == "#!gff-spec-version 1.21") {
      print("version 1.21------------------")
      bed <- simplify_gff3_1.21(gff, prefix)
      } else {
      bed <- simplify_gff3(gff, prefix)
      }
    print(bed)
    write_tsv(bed, str_c("short_bed/", prefix, ".bed"), col_names = FALSE)  
  } else if (str_detect(gff_raw, ".gff$")) {
      gff <- read_gff(gff_raw)
      print(gff)
      bed <- simplify_gff(gff, prefix)
      print(bed)
      write_tsv(bed, str_c("short_bed/", prefix, ".bed"), col_names = FALSE)
  } else if (str_detect(gff_raw, ".gtf$")) {
      print(".gtf------------------")
      gff <- read_gtf(gff_raw)
      print(gff)
      bed <- simplify_gtf(gff, prefix)
      print(bed)
      write_tsv(bed, str_c("short_bed/", prefix, ".bed"), col_names = FALSE)  
  } else {
 print(" ERROR: please use .gff3, .gff or .gtf -----------------------------------------------------")  
  }
}
```

### walk input

``` r
dir.create("short_bed")
map(input$species, write_short_bed)
```

## simplify .peps

### rename .peps

``` r
dir.create("short_pep")
for (i in 1:nrow(input)) {
 print(input$pep_raw[[i]])
 print(paste0("pep/",str_remove(input$pep_raw[[i]], ".gz")))
 print(paste0("short_pep/",input$acronym[[i]], ".pep"))
 file.copy(paste0("pep/",str_remove(input$pep_raw[[i]], ".gz")), paste0("short_pep/",input$acronym[[i]], ".pep"))
 }
list.files("short_pep/")
#setwd("short_pep/")
```

### single line & shorten fasta header

``` bash
cd short_pep
for PEP in *.pep
do 
PREFIX=$(echo $PEP | sed 's/.pep//')
echo "$PREFIX"
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $PEP > ${PREFIX}.single
cut -d " " -f1 ${PREFIX}.single > ${PREFIX}.tmp
sed  's/^>\(.*\)/>'"${PREFIX}"'__\1/' ${PREFIX}.tmp > ${PREFIX}.tmp2
sed 's/\.p$//' ${PREFIX}.tmp2 > ${PREFIX}.tmp3
done
```

## compare fasta to bed header

``` bash
for TMP in *.tmp3;
do 
  grep -h -m 3 "^>" $TMP >> pep.txt
done
cd /cluster/work/users/jonathbr/mcscanx/short_bed/
head -n 3  *.bed | grep -v "^=" | sed '/^$/d'|  cut -f 4 | sed 's/^/>/' > bed.txt
comm bed.txt ../short_pep/pep.txt
```

### small fixes

``` bash
#cut -f 10 Aly.bed | cut -d ";" -f 1,2 | grep "ID=cds" | grep "Parent=gene-" | sort -n | uniq | sed 's/ID=cds-/>/' | sed 's/\.[0-9]\;Parent=gene-/\t/' > Aly_pep_bed.tsv
```

### compare gene count

``` bash
cd ../short_pep/
grep -c "^>" *.tmp3 > ../pep_count.txt
cd ../short_bed/
wc -l *.bed > ../bed_count.txt 
```

## rename files

``` bash
cd /cluster/work/users/jonathbr/mcscanx/short_pep
for t in *tmp3;
do
p=$(echo $t | sed 's/.tmp3/.pep/')
echo "
mv $t $p" 
mv $t $p
done 
cd ..
mkdir cleaned
cp short_pep/*.tmp3 cleaned/
cp short_bed/*.bed cleaned/
wget https://github.com/zhaotao1987/SynNet-Pipeline/raw/master/SynetBuild-X.sh
chmod u+x SynetBuild-X.sh
# Change Line70 according to your genome list
# Aaa Ape Ath Aly
```

``` bash
cd /cluster/work/users/jonathbr/mcscanx/cleaned
module purge
module load R/3.6.2-foss-2019b
R
library(tidyverse)
```

### sort\_bed function

``` r
sort_exclude_transcripts_bed <- function(bed) {
  read_tsv(bed, col_names = c("chrom", "start", "end", "id")) %>% 
  select(chrom, id, start, end) %>% 
  mutate(id = str_remove(id, "\\.\\d{1,2}$")) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  arrange(chrom, start, end) %>% 
  write_tsv(paste0("sorted/",bed), col_names = FALSE)
}
#sort_exclude_transcripts_bed("Dst.bed")
files <- list.files(".", pattern = ".bed")
dir.create("sorted")
walk(files, sort_exclude_transcripts_bed)
```

### detect splice variants in .pep

``` r
library(seqinr) 
  
longest_transcript <- function(fasta_file) {
  fasta <- read.fasta(fasta_file, seqtype ="AA")
  fa <- tibble(label = names(fasta), 
             length =getLength(fasta)) %>% 
    mutate(new_label = str_remove(label, "\\.\\d{1,2}$")) %>% 
    arrange(new_label, length) %>% 
    distinct(new_label, .keep_all = T)
  write.fasta(fasta[fa$label],
            names = fa$new_label, 
            file.out = paste0("longest_pep/", fasta_file), 
            nbchar = max(fa$length), as.string = F)
}
files <- list.files(".", pattern = ".pep")
dir.create("longest_pep")
walk(files, longest_transcript)
```

``` bash
rm *.bed
rm *.pep
mv sorted/*.bed .
mv longest_pep/*.pep .
rmdir sorted
rmdir longest_pep

cp ../SynetBuild-X.sh .
sbatch synetbuild.slurm #job 2281144  00:17:56   00:02:45 175 612K
```

# test plots

``` bash
ssh -Y jonathbr@saga.sigma2.no

cd /cluster/work/users/jonathbr/mcscanx/
mkdir test_plots
cd test_plots

PATH=$PATH:/cluster/projects/nn9525k/Programs/MCScanX/

module load Java/11.0.2

cut ../cleaned/Aaa.bed -f1 | sort -n | uniq
cut ../cleaned/Alp.bed -f1 | sort -n | uniq
cut ../cleaned/Aly.bed -f1 | sort -n | uniq
cut ../cleaned/Ath.bed -f1 | sort -n | uniq

wget https://github.com/wyp1125/MCScanX/raw/master/downstream_analyses/dot.ctl #modify chromosomes

java -Xmx2G /cluster/projects/nn9525k/Programs/MCScanX/downstream_analyses/dot_plotter.java -g ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Aaa_Alp.gff -s ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Aaa_Alp.collinearity -c Aaa_Alp.ctr -o Aaa_Alp_dotplot.png

#java dual_synteny_plotter -g gff_file -s collinearity_file -c control_file -o output_PNG_file

java -Xmx2G /cluster/projects/nn9525k/Programs/MCScanX/downstream_analyses/dual_synteny_plotter.java -g ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Aaa_Alp.gff -s ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Aaa_Alp.collinearity -c Aaa_Alp.ctr -o Aaa_Alp_dualsyn.png
```

``` bash
java -Xmx2G /cluster/projects/nn9525k/Programs/MCScanX/downstream_analyses/dual_synteny_plotter.java -g ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Aaa_Ath.gff -s ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Aaa_Ath.collinearity -c Aaa_Ath.ctl -o Aaa_Ath_dualsyn.png
```

``` bash
java -Xmx2G /cluster/projects/nn9525k/Programs/MCScanX/downstream_analyses/dual_synteny_plotter.java -g ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Alp_Ath.gff -s ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Alp_Ath.collinearity -c Alp_Ath.ctl -o Alp_Ath_dualsyn.png
```

``` bash
java -Xmx2G /cluster/projects/nn9525k/Programs/MCScanX/downstream_analyses/dual_synteny_plotter.java -g ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Ath_Aly.gff -s ../cleaned/SynNetBuild20210302_1543-SynNet-k6s5m25/Ath_Aly.collinearity -c Ath_Aly.ctl -o Ath_Aly.dualsyn.png
```
