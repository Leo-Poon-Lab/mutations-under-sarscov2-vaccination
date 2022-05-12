# Data selection
## we should include: 1. VOC and B.1 cases unvaccinated cases; 2. all data with known vaccinated history.

library(lubridate)
library(tidyverse)
library(ssh.utils)

data_meta_raw <- read_csv("../../../2021/2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 20000)
data_meta <- data_meta_raw %>% filter(sequenced_by_us) %>% filter(!is.na(`Report date`))

## Unvac VOC and B.1 cases in HK 
data_meta_unvac <- data_meta %>% filter(is.na(`Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)) %>% filter(lineage!="None") %>% filter(lineage!="Unassigned") %>% filter(coverage>=0.85) %>% filter(Ct_value<=25) %>% filter(Ct_value!=0)
sort(table(data_meta_unvac$lineage))
sort(table(data_meta_unvac$variant))
quantile(data_meta_unvac$coverage)
data_meta_unvac <- data_meta_unvac %>% filter(lineage %in% c("B.1.36", "B.1.1.63", "B.1.36.27") | grepl("Alpha", variant) | grepl("Delta", variant) | grepl("Omicron", variant)) # selected

## Vaccinated cases
data_meta_vac <- data_meta %>% filter(!is.na(`Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)) %>% filter(lineage!="None") %>% filter(lineage!="Unassigned") %>% filter(coverage>=0.85) %>% filter(Ct_value<=25) %>% filter(Ct_value!=0)
sort(table(data_meta_vac$lineage))
sort(table(data_meta_vac$variant))
quantile(data_meta_vac$coverage)
data_meta_vac <- data_meta_vac %>% filter(grepl("Delta", variant) | grepl("Omicron", variant) ) # selected
### clean the vaccination section
table(data_meta_vac$`Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)
table(data_meta_vac$`Name of 2nd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)
table(data_meta_vac$`Name of 3rd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)

## Merge metadata
## all selected cases
df_meta <- bind_rows(data_meta_vac, data_meta_unvac)
samples_primer_old <- readLines("../data/samples_primer1.txt")
### primer and Ct value to metadata
df_meta$primer <- ifelse(df_meta$Sample %in% samples_primer_old, "old", "new")

## Upload data to server and re-analyse SNPs/iSNVs
samples_toupload_ori <- unique(df_meta$Sample)
samples_alrd_done <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
samples_alrd_done <- gsub("-trimmed.+", "", samples_alrd_done)
samples_toupload <- samples_toupload_ori[!samples_toupload_ori %in% samples_alrd_done] # we alredy analysis some of the data
sort(samples_toupload)
df_meta %>% filter(Sample %in% samples_toupload) %>% .$lineage %>% table()
# samples_toupload <- samples_toupload[samples_toupload>="WHP4675"]

files_fastq_all <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive/", full.names = T)
samples_all <- sapply(files_fastq_all, function(x){
	tmp <- strsplit(x, "/", fixed=T)[[1]]
	sample_t <- tmp[length(tmp)]
	gsub(".fastq.gz", "", sample_t)
})
files_fastq1 <- files_fastq_all[samples_all %in% paste0(samples_toupload, "_1")]
files_fastq2 <- files_fastq_all[samples_all %in% paste0(samples_toupload, "_2")]

samples_not_available <- samples_toupload[!samples_toupload %in% gsub(".+/", "", gsub("_1.fastq.gz", "", files_fastq1))]

# paste0(samples, "_1")[!paste0(samples, "_1") %in% samples_all] # On 2021-11-03, there is still one sample not available

## copy the files to remote
sapply(c(rbind(files_fastq1, files_fastq2)), function(x){
	cp.remote(path.src = x, remote.src = "", remote.dest = "hggu@147.8.72.238", path.dest = "~/work/2020-09-01_COVID_NGS_pipeline/NGS_data_input/", verbose = T)
})

df_meta <- df_meta %>% filter(!Sample %in% samples_not_available)

write_csv(df_meta, "../data/df_samples.csv")
