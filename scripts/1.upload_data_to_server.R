# Data selection
## we should include: 1. Alpha and B.1 cases in HK before 2021-02-22 as unvaccinated cases; 2. all data with known vaccinated history.

library(lubridate)
library(tidyverse)
library(ssh.utils)

## Alpha and B.1 cases in HK before 2020-12-01
data_hk_seq <- read_csv("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/cov_and_lin.csv")
data_hk_seq <- data_hk_seq %>% filter(!grepl("iseq", Sample))
data_hk_seq_bf <- data_hk_seq %>% filter(dmy(`Report date`)<="2021-02-22") 

system("pangolin ../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome_caseid.fasta --outfile ../results/lineage_hkall.csv")
df_lin_hkall <- read_csv("../results/lineage_hkall.csv")
df_lin_hkall$Sample <- sapply(df_lin_hkall$taxon, function(x){
	strsplit(x, "_")[[1]][1]
})
data_hk_seq_bf <- left_join(data_hk_seq_bf %>% select(Sample:`Analysis date`), df_lin_hkall, "Sample")
sort(table(data_hk_seq_bf$scorpio_call))
data_b1 <- data_hk_seq_bf %>% filter(lineage=="B.1") %>% filter(`HK/Non-HK resident`=="HK resident")
data_alpha <- data_hk_seq_bf %>% filter(grepl("Alpha", scorpio_call)) %>% filter(`HK/Non-HK resident`=="HK resident")

## known vaccinated history
data_selected <- readxl::read_excel("../data/Vaccine Delta SNP list.xlsx", skip = 1)
data_vaccinated_add <- readxl::read_excel("../data/Vaccination case for breakthrough infection .xlsx", skip=1)

data_vaccinated_add$case_id <- data_vaccinated_add$`HK case No.`
data_vaccinated_add <- left_join(data_vaccinated_add, data_hk_seq, "case_id")
unique(data_vaccinated_add$`WHP ID`)
samples_add <- unique(data_vaccinated_add$Sample)

samples_vac <- unique(na.omit(c(data_selected$WHP_ID, samples_add)))
data_selected <- data_selected[!is.na(data_selected$WHP_ID),]
data_selected <- data_selected[!data_selected$WHP_ID=="WHP_ID",]

samples_toupload_ori <- unique(c(samples_vac, data_b1$Sample, data_alpha$Sample))
samples_toupload <- samples_toupload_ori[!samples_toupload_ori %in% data_selected$WHP_ID] # we alredy analysis some of the data

## Merge metadata
data_selected$`Date of 2nd Dose` <- sapply(data_selected$`Date of vaccination`, function(x){
	tmp <- strsplit(x, " and ")[[1]][2]
	as.character(as_date(paste0("2021 ", tmp)))
})
data_selected$`Report Date` <- as_date(as.numeric(data_selected$`Report Date`), origin = "1899-12-30")
data_selected$days_since_last_dose <- as.character(as.numeric(difftime(data_selected$`Date of 2nd Dose`, data_selected$`Report Date`, units = "days")))
df_vac_1 <- data_selected %>% mutate(Sample=WHP_ID, case_id=`Case No.`) %>% select(Sample, case_id, Vaccine, Doses, days_since_last_dose)

data_vaccinated_add$Doses <- ifelse(is.na(data_vaccinated_add$`Date of 2nd Dose`), 1, 2)
df_vac_2 <- data_vaccinated_add %>% mutate(days_since_last_dose=`No. of Days between 2nd dose and confirmed date`) %>% select(Sample, case_id, Vaccine, Doses, days_since_last_dose) %>% filter(!is.na(Sample))

df_unvac <- tibble(Sample=c(data_b1$Sample, data_alpha$Sample), case_id=c(data_b1$case_id, data_alpha$case_id), Vaccine="/")

df_meta <- bind_rows(df_vac_1, df_vac_2 %>% mutate_all(as.character), df_unvac %>% mutate_all(as.character))
df_meta$Vaccine[df_meta$Vaccine=="/"] <- "Non-vaccinated"
samples_primer_old <- readLines("../data/samples_primer1.txt")

### primer and Ct value to metadata
df_meta$primer <- ifelse(df_meta$Sample %in% samples_primer_old, "old", "new")
df_meta <- df_meta %>% filter(Sample!="WHP_ID")
df_meta <- df_meta %>% arrange(case_id) 
df_meta <- df_meta[!duplicated(df_meta$Sample),]

## TODO add Ct values
df_meta <- left_join(df_meta, data_selected  %>% mutate(Sample=WHP_ID) %>% select(Sample, `Ct value`))

write_csv(df_meta, "../data/df_samples.csv")
write_csv(df_meta, "../results/df_samples_clean.csv")

## Upload data to server and re-analyse SNPs/iSNVs
files_fastq_all <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive/", full.names = T)
samples_all <- sapply(files_fastq_all, function(x){
	tmp <- strsplit(x, "/", fixed=T)[[1]]
	sample_t <- tmp[length(tmp)]
	gsub(".fastq.gz", "", sample_t)
})
files_fastq1 <- files_fastq_all[samples_all %in% paste0(samples_toupload, "_1")]
files_fastq2 <- files_fastq_all[samples_all %in% paste0(samples_toupload, "_2")]

# paste0(samples, "_1")[!paste0(samples, "_1") %in% samples_all] # On 2021-11-03, there is still one sample not available

## copy the files to remote
sapply(c(files_fastq1, files_fastq2), function(x){
	cp.remote(path.src = x, remote.src = "", remote.dest = "hggu@147.8.70.135", path.dest = "~/work/2020-09-01_COVID_NGS_pipeline/archive/", verbose = T)
})


