# Arabidopsis arenosa and A. lyrata genome assemblies
This is a GitHub repository with data associated to the paper:   
   

**Structural evidence for MADS-box type I family expansion seen in new assemblies of <i>Arabidopsis arenosa</i> and <i>A. lyrata</i>**

_The Plant Journal_

​​**Authors:** Jonathan Bramsiepe, Anders K. Krabberød, Katrine N. Bjerkan, Renate M. Alling, Ida M. Nielsen. Karina S. Hornslien, Jason R. Miller, Anne K. Brysting, and Paul E. Grini   
[https://onlinelibrary.wiley.com/doi/10.1111/tpj.16401](https://onlinelibrary.wiley.com/doi/10.1111/tpj.16401)

**Under Construction**
## 1. <i>Arabidopsis arenosa</i> 
The assembly pipeline with scripts and additional information can be found  [here (under construction)](01_arenosa_assembly/).  
The genome can be accessed at [the NCBI Genome page](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_026151155.1/)   
Or you can download it directly from github:
- The full genome: [Arabidopsis_arenosa_genome.fna.gz](01_arenosa_assembly/Arabidopsis_arenosa_genome.fna.gz)
- The genome with repeats softmasked according to RepeatMasker: [Arabidopsis_arenosa_genome.softmasked.fna.gz](01_arenosa_assembly/Arabidopsis_arenosa_genome.softmasked.fna.gz)
- The genome annotation in gtf format with the predicted genes from braker2 [Arabidopsis_arenosa_genome.annotation.gtf.gz](01_arenosa_assembly/Arabidopsis_arenosa_genome.annotation.gtf.gz).

## 2. <i>Arabidopsis lyrata</i> 
The assembly pipeline with scrips and additional information can be found  [here (under construction)](02_lyrata_assembly/)
Link to the [NCBI Genome Page](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_026151145.1/)
- The full genome: [Arabidopsis_lyrata_petraea_genome.fna.gz](02_lyrata_assembly/Arabidopsis_lyrata_petraea_genome.fna.gz)
- The genome with repeats softmasked according to RepeatMasker: [...](02_lyrata_assembly/)
- The genome annotation in gtf format with the predicted genes from braker2 [Arabidopsis_lyrata_petraea_genome.annotation.gtf.gz](02_lyrata_assembly/Arabidopsis_lyrata_petraea_genome.annotation.gtf.gz).

## 3. Analysis
In the paper, we compared the new assemblies with previously published genomes from the Arabidopsis genus, as well as the evolution of the MADS-box gene family in the genus Arabidopsis.  
**Content of the analysis** [here](./03_analysis/)
- [3.1 Assembly and scaffolding stats](./03_analysis/01_assembly_and_scaffolding_stats/)
- [3.2 Orthology analysis](./03_analysis/02_ortholog_prediction/)
- [3.3 Synteny](./03_analysis/03_synteny/)
- [3.5 Phylogeny og MADS-boxes](./03_analysis/05_MADS_phylogeny/)
- [3.6 MEME motif analysis og MADS-box genes](./03_analysis/06_MADS_MEME/)
MEME analysis of  MADS-box genes type I and type II in <i>A. lyrata</i> ssp. <i>petraea</i>, <i>A. arenosa</i> Pusté Pole, <i>A. lyrata </i>ssp. <i>lyrata,</i> <i>A. arenosa</i> Strecno, <i>A. halleri</i>,  <i>Capsella rubella</i> and <i>C. grandiflora</i>.  View the full [MEME results as html](https://htmlpreview.github.io/?https://github.com/krabberod/html_test/blob/main/meme.html) . View/download the results as [a text file](03_analysis/06_MADS_MEME/meme_results.txt). 

*
----
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FPaulGrini%2FArabidopsis_assemblies&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)