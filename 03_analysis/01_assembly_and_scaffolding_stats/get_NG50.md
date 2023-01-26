get_NG50
================
2021-02-01

-   [read google sheet](#read-google-sheet)
-   [add BUSCO](#add-busco)

``` bash
ssh -Y jonathbr@saga.sigma2.no
srun --cpus-per-task=2 --mem-per-cpu=2G --time=04:00:00 --account=nn9525k --x11 --pty bash -i
module purge
module load R/3.6.2-foss-2019b
R
```

``` r
library(tidyverse)
library(fs)
library(googlesheets4)
```

``` r
dir_create("compare_assemblies")
setwd("compare_assemblies")
getwd()
dir_create("assemblies")
```

## read google sheet

``` r
sheets <- c("data", "assemblies", "purged", "scaffolds", "published")
gs4_deauth()
google_list <- map(sheets, read_sheet, ss="https://docs.google.com/spreadsheets/d/1Y-HlibIyD7i_aFtaiN2K53iRynaMvUUzxKhce3amFvc/edit?usp=sharing") %>% set_names(sheets)
```

``` r
#check names
#google_list[-1] %>% map(names) %>% unlist() %>% unique()

google_tb <- google_list[-1] %>% bind_rows(.id = "stage")
```

``` r
#check unique Assembly_name
#google_tb[duplicated(google_tb$Assembly_name),1:3]
```

``` r
#check for existing file under Path 
copy <- google_tb %>% 
  #filter(stage == "assemblies") %>% 
  mutate(f_exists = file_exists(Path),
         f_type = if_else(str_detect(Path, "\\.gz$"), ".fasta.gz", ".fasta"),
         cp_to = paste0("assemblies/",Assembly_name, f_type)) %>% 
  filter(f_exists) %>% 
    select(Path, cp_to)

walk2(copy$Path, copy$cp_to, file_copy)
```

``` r
#gunzip("assemblies/02_falcon_p_h_ctg.fasta.gz")
system("gunzip assemblies/*.gz")
dir_info("assemblies")[c("path","size")]
```

``` r
get_NG50 <- function(ass_path) {
  sysout <- system(paste0("~/Programs/calN50.js -L200000000 ",ass_path), intern = TRUE)
  sysout[15] %>% str_remove("NL\t50\t")
}

NG50_list <- dir_map("assemblies", get_NG50) 

NG50_tb <- tibble(a_name = str_remove_all(dir_ls("assemblies"), ".*assemblies/|\\.fasta"),
                  NL50 = unlist(NG50_list)) %>% 
  separate(NL50, c("NG50", "LG50"), sep = "\t", convert = TRUE)
```

``` r
new_csv <- google_tb %>% 
  select(-c("NG50", "LG50")) %>% 
  left_join(NG50_tb, by = c("Assembly_name" = "a_name"))

write_csv2(new_csv, "new_google_NG50_tmp.csv")

dir_ls()
```

## add BUSCO

``` bash
cd /cluster/projects/nn9525k/
find . -name "short_summary_*.txt" -exec realpath {} \; > busco_summary.fofn
mv busco_summary.fofn /cluster/work/users/jonathbr/compare_assemblies
```

``` r
read_csv("busco_summary.fofn", col_names = "path") %>%
  mutate(file = basename(path),
         dir = dirname(path) %>% str_remove("/cluster/projects/nn9525k/"),
         name_v1 = str_remove_all(file, "short_summary_|\\.txt")) %>% 
write_csv2("buco_summary.csv")
```

``` r
busco <- read_csv2("buco_summary_v2.csv")

busco %>% 
  filter(busco_type == "genomic") %>% 
  count(a_name, sort = TRUE)

busco_sum <- busco %>% 
  filter(busco_type == "genomic") %>% 
  distinct(a_name, .keep_all = TRUE) %>% 
  mutate(cp_to = paste0("busco_sum/", a_name, "_busco.txt")) %>% 
  select(a_name, busco_path = path, cp_to)
```

``` r
dir_create("busco_sum")

walk2(busco_sum$busco_path, busco_sum$cp_to, file_copy)
```

``` r
ng_csv <- read_csv2("new_google_NG50_tmp.csv")

ng_busco <- left_join(ng_csv, busco_sum, by = c("Assembly_name" = "a_name"))
```

``` r
sysout <- system('grep "C:" busco_sum/*.txt', intern = TRUE)

sysout_tb <- tibble(out = sysout) %>% 
  separate(out, c("aname", "busco_string"), sep = "\t") %>% 
  mutate(aname = str_remove(aname, ":"))

ng_busco_sum <- left_join(ng_busco, sysout_tb, by = c("cp_to" = "aname"))

write_csv2(ng_busco_sum, "buco_summary_v3.csv")  
```
