# Data selection
## we should include: 1. Alpha and B.1 cases in HK before 2021-02-22 as unvaccinated cases; 2. all data with known vaccinated history.

library(lubridate)
library(tidyverse)
library(ssh.utils)

## Alpha and B.1 cases in HK before 2021-02-22
data_hk_seq <- read_csv("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/cov_and_lin.csv")
data_hk_seq <- data_hk_seq %>% filter(!grepl("iseq", Sample))
data_hk_seq_bf <- data_hk_seq %>% filter(dmy(`Report date`)<="2021-02-22") 
df_meta_gov <- tryCatch(read_csv("https://www.chp.gov.hk/files/misc/enhanced_sur_covid_19_eng.csv", col_types = cols(.default = "c")), error = function(err){"error"}, warning = function(x){"warning"})
df_meta_gov$case_id <- df_meta_gov$`Case no.`

# system("conda run -n pangolin pangolin ../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome_caseid.fasta --outfile ../results/lineage_hkall.csv")
df_lin_hkall <- read_csv("../results/lineage_hkall.csv")
df_lin_hkall$Sample <- sapply(df_lin_hkall$taxon, function(x){
	strsplit(x, "_")[[1]][1]
})
data_hk_seq_bf <- left_join(data_hk_seq_bf %>% select(Sample:`Analysis date`), df_lin_hkall, "Sample")
sort(table(data_hk_seq_bf$scorpio_call))
sort(table(data_hk_seq_bf$lineage))

## known vaccinated history
data_selected <- readxl::read_excel("../data/Vaccine Delta SNP list.xlsx", skip = 1)
data_selected <- data_selected[!is.na(data_selected$WHP_ID),]
data_selected <- data_selected[!data_selected$WHP_ID=="WHP_ID",]
data_selected$`Date of 2nd Dose` <- sapply(data_selected$`Date of vaccination`, function(x){
	tmp <- strsplit(x, " and ")[[1]][2]
	as.character(as_date(paste0("2021 ", tmp)))
})
data_selected$`Report Date` <- as_date(as.numeric(data_selected$`Report Date`), origin = "1899-12-30")
data_selected$days_since_last_dose <- as.character(as.numeric(difftime(data_selected$`Date of 2nd Dose`, data_selected$`Report Date`, units = "days")))
df_vac_1 <- data_selected %>% mutate(Sample=WHP_ID, case_id=`Case No.`) %>% select(Sample, case_id, Vaccine, Doses, days_since_last_dose)

data_vaccinated_add <- readxl::read_excel("../data/Vaccination case for breakthrough infection .xlsx", skip=1)
data_vaccinated_add$case_id <- data_vaccinated_add$`HK case No.`
data_vaccinated_add <- left_join(data_vaccinated_add, data_hk_seq, "case_id")
data_vaccinated_add$Doses <- ifelse(is.na(data_vaccinated_add$`Date of 2nd Dose`), 1, 2)
df_vac_2 <- data_vaccinated_add %>% mutate(days_since_last_dose=`No. of Days between 2nd dose and confirmed date`) %>% select(Sample, case_id, Vaccine, Doses, days_since_last_dose) %>% filter(!is.na(Sample))

data_annot_raw <- readxl::read_excel("../data/HK case annotation.xlsx")
data_annot <- data_annot_raw %>% filter(grepl("vaccine", tolower(Annotation_DH)) | grepl("vaccine", tolower(Annotation_2)) | grepl("vaccine", tolower(Annotation_1))) %>% mutate(case_id=`HK case`, Sample=`WHP_number`) %>% select(case_id,Sample, Annotation_1:Annotation_DH)
data_master_raw <- readxl::read_excel("../data/Leo WGS master file.xlsx", guess_max = 100000)
data_master <- data_master_raw %>%  filter(grepl("vaccine", tolower(Annotation_DH)) | grepl("vaccine", tolower(Annotation_2)) | grepl("vaccine", tolower(Annotation_1)) | grepl("vaccine", tolower(Remark))) %>% mutate(case_id=`HK Case No.`, Sample=`Lab ID`) %>% select(Sample, case_id, Annotation_1:Annotation_DH, Remark)
df_vac_3 <- bind_rows(data_master, data_annot) %>% filter(!is.na(Sample)) %>% filter(!duplicated(Sample)) %>% filter(!Sample %in% c(data_selected$WHP_ID, data_vaccinated_add$`WHP ID`)) %>% pivot_longer(Annotation_1:Remark) %>% filter(!is.na(value)) %>% mutate(value=tolower(value)) %>% filter(grepl("vaccine", value))
df_vac_3$Vaccine <- sapply(df_vac_3$value, function(x) {
	if(grepl("biontech", x)){return("BioNTech")}
	if(grepl("biontch", x)){return("BioNTech")}
	if(grepl("sinovac", x)){return("Sinovac")}
	tmp <- strsplit(x, " ")[[1]]
	tmp[which(grepl("vaccine", tmp))-1]
})
table(df_vac_3$Vaccine)
df_vac_3$value <- gsub("odses", "doses", df_vac_3$value)
df_vac_3$Doses <- sapply(df_vac_3$value, function(x) {
	tmp <- strsplit(x, " ")[[1]]
	if(!any(grepl("dose", tmp))){
		tmp <- tmp[which(grepl("does", tmp))-1]
	} else {
		tmp <- tmp[which(grepl("dose", tmp))-1]
	}
	if(length(tmp)==0){return(NA)}
	if(any(is.na(as.numeric(tmp)))){
		tmp <- gsub("one","1",tmp)
		tmp <- gsub("first","1",tmp)
		tmp <- gsub("1st","1",tmp)
		tmp <- gsub("two","2",tmp)
		tmp <- gsub("both","2",tmp)
		tmp <- gsub("three","3",tmp)
	}
	max(as.numeric(tmp), na.rm=T)
})
table(df_vac_3$Doses)
df_vac_3$value[is.na(df_vac_3$Doses)]
df_vac_3$Doses[df_vac_3$value=="sample contains n501y variant virus gene_ received both shots for biontech vaccine"] <- 2
df_vac_3$Doses[df_vac_3$value=="patient received at biontech vaccine shot on 15/3"] <- 1

df_vac_3$value[grepl("2020", df_vac_3$value)]
df_vac_3$value[grepl("2021", df_vac_3$value)] ## All cases were Vaccinated on 2021
df_vac_3$value[grepl("2022", df_vac_3$value)]
df_vac_3$Date_vac<- sapply(df_vac_3$value, function(x) {
	tmp <- strsplit(x, " on ")[[1]]
	if(is.na(tmp[2])){
		tmp <- strsplit(x, " in ")[[1]]
	}
	tmp[length(tmp)]
}, USE.NAMES = F)
df_vac_3$value[is.na(df_vac_3$Date_vac)]
df_vac_3$Date_vac_parse <- sapply(df_vac_3$Date_vac, function(x) {
	tmp <- gsub("[[:punct:]]", " ", x)
	tmp <- strsplit(tmp, "\\s")[[1]]
	check_num <- !is.na(as.numeric(tmp))
	if(!any(check_num)){return(NA)}
	idx <- which(check_num)
	dates_x <- sapply(idx, function(i) {
		paste0("2021 ", tmp[i-1], " ", tmp[i])
	})
	max(as.character(ymd(dates_x)), na.rm=T)
}, USE.NAMES = F)
df_vac_3$Date_vac_parse[df_vac_3$Date_vac=="15/3"] <- "2021-03-15"
df_vac_3$Date_vac_parse[df_vac_3$Date_vac=="jan and feb"] <- "2021-02-01"
df_vac_3$Date_vac_parse[df_vac_3$value=="received 2 doses of biontech vaccines on march 19 and april 24; positive on day 3 in nina hotel"] <- "2021-04-24"
df_vac_3$Date_vac_parse[df_vac_3$value=="received 2 doses of biontech vaccines in the uk on august 11 and oct 7; staying in room 1646 and tested positive on 4 nov (day 3)"] <- "2021-10-07"
df_vac_3$Date_vac_parse[df_vac_3$value=="received 2 doses of astrazeneca in the uk on feb 13 and april 27, received one dose of biontech vaccine in the uk on oct 31; positive on day 5 at quarantine hotel"] <- "2021-04-27"
df_vac_3$value[is.na(df_vac_3$Date_vac_parse)]
df_vac_3$Date_vac[is.na(df_vac_3$Date_vac_parse)]

df_vac_3 <- left_join(df_vac_3, df_meta_gov %>% mutate(case_id=as.character(case_id)) %>% filter(!is.na(case_id)) %>% select(case_id, `Report date`) %>% unique(), "case_id")
df_vac_3$`Report date` <- dmy(df_vac_3$`Report date`)
df_vac_3$days_since_last_dose <- as.character(as.numeric(difftime(df_vac_3$`Date_vac_parse`, df_vac_3$`Report date`, units = "days")))
df_vac_3 <- df_vac_3 %>% select(Sample, case_id, Vaccine, Doses, days_since_last_dose)
df_vac_3 <- df_vac_3 %>% filter(Sample %in% data_hk_seq$Sample) # sequenced

df_vac <- bind_rows(df_vac_1%>% mutate_all(as.character), df_vac_2%>% mutate_all(as.character), df_vac_3%>% mutate_all(as.character))
df_vac$lineage <- sapply(df_vac$Sample, function(x){
	tmp <- df_lin_hkall$lineage[df_lin_hkall$Sample==x]
	if(length(tmp)>1){print(tmp)}
	tmp
})

## Non-vaccinated
df_unvac_meta <- data_hk_seq_bf %>% filter(lineage%in%c("B.1.36", "B.1.36.27", "B.1.1.63", "B.1") | grepl("Alpha", scorpio_call)) %>% filter(`HK/Non-HK resident`=="HK resident")
## samples after 2021-07-29, non-vaccinated Delta
data_hk_seq_after729 <- data_hk_seq %>% filter(dmy(`Report date`)>"2021-07-29") 
df_unvac_2 <- data_hk_seq_after729 %>% filter(!Sample %in% df_vac$Sample)

data_annot_raw <- readxl::read_excel("../data/HK case annotation.xlsx")
data_annot_unvac <- data_annot_raw %>%  mutate(case_id=`HK case`, Sample=`WHP_number`) %>% select(case_id,Sample, Annotation_1:Annotation_DH)
data_master_raw <- readxl::read_excel("../data/Leo WGS master file.xlsx", guess_max = 100000)
data_master_unvac <- data_master_raw %>% mutate(case_id=`HK Case No.`, Sample=`Lab ID`) %>% select(Sample, case_id, Annotation_1:Annotation_DH, Remark)
df_annt_info <- bind_rows(data_annot_unvac, data_master_unvac) %>% filter(!is.na(Sample)) %>% filter(!duplicated(Sample))

df_annt_info %>% filter(Sample %in% df_unvac_2$Sample) %>% .$Annotation_1 ## Manual check
df_annt_info %>% filter(Sample %in% df_unvac_2$Sample) %>% .$Annotation_2 ## Manual check
df_annt_info %>% filter(Sample %in% df_unvac_2$Sample) %>% .$Annotation_DH ## Manual check
df_annt_info %>% filter(Sample %in% df_unvac_2$Sample) %>% .$Remark ## Manual check

df_unvac <- tibble(Sample=df_unvac_meta$Sample, case_id=df_unvac_meta$case_id, Vaccine="/", lineage=df_unvac_meta$lineage)
df_unvac_2 <- df_unvac_2 %>% select(Sample, case_id, lineage) %>% mutate(Vaccine="/")

df_unvac_combined <- bind_rows(df_unvac, df_unvac_2)

## Merge metadata
df_meta <- bind_rows(df_vac, df_unvac_combined %>% mutate_all(as.character))
df_meta$Vaccine[df_meta$Vaccine=="/"] <- "Non-vaccinated"
samples_primer_old <- readLines("../data/samples_primer1.txt")

### primer and Ct value to metadata
df_meta$primer <- ifelse(df_meta$Sample %in% samples_primer_old, "old", "new")
df_meta <- df_meta %>% filter(Sample!="WHP_ID")
df_meta <- df_meta %>% arrange(case_id) 
df_meta <- df_meta[!duplicated(df_meta$Sample),]

## add Ct values
df_ct <- readxl::read_excel("../data/WHP Sample List updated as of 20210520.xlsx")
df_ct$`Ct value` <- df_ct$`nCoV-2019 info (CHP/AFCD Ct)`
df_ct$sample_sim <- df_ct$`Lab No.`

df_meta$sample_sim <- gsub("-.+$", "", df_meta$Sample)
df_meta <- left_join(df_meta, df_ct %>% select(sample_sim, `Ct value`), "sample_sim")
df_meta$Vaccine <- stringr::str_to_title(df_meta$Vaccine)
table(df_meta$Vaccine)
df_meta$Vaccine[df_meta$Vaccine=="V"] <- "Sputnik V"
df_meta <- df_meta %>% filter(Vaccine!="Of")

write_csv(df_meta, "../data/df_samples.csv")
write_csv(df_meta, "../results/df_samples_clean.csv")

## Upload data to server and re-analyse SNPs/iSNVs
samples_vac <- unique(na.omit(df_meta$Sample))
samples_vac <- samples_vac[samples_vac!="WHP_ID"]
samples_toupload_ori <- unique(c(samples_vac, df_meta$Sample))
samples_alrd_done <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
samples_alrd_done <- gsub("-t.+", "", samples_alrd_done)
samples_toupload <- samples_toupload_ori[!samples_toupload_ori %in% data_selected$WHP_ID] # we alredy analysis some of the data
samples_toupload <- samples_toupload[!samples_toupload %in% samples_alrd_done] # we alredy analysis some of the data
sort(samples_toupload)
# samples_toupload <- samples_toupload[samples_toupload>="WHP4675"]

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
sapply(c(rbind(files_fastq1, files_fastq2)), function(x){
	cp.remote(path.src = x, remote.src = "", remote.dest = "hggu@147.8.70.135", path.dest = "~/work/2020-09-01_COVID_NGS_pipeline/NGS_data_input/", verbose = T)
})
