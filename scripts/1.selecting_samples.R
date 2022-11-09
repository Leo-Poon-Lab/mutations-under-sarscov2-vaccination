# Data selection
## we should include: 1. Local cases; 2. Dilution samples; 3. Serial samples;
## 1.1. VOC and B.1 cases unvaccinated cases; 1.2. all data with known vaccinated history.

library(lubridate)
library(tidyverse)
library(ssh.utils)
library(ggsci)

## 1. Local cases;
# data_meta_study1 <- read_csv("../data/df_samples_first_version.csv") 
data_meta_study1 <- read_csv("../results/df_samples_clean_first_version.csv") # in the first version of the manuscript, we have 2053 samples passed 
data_meta_raw <- read_csv("../../../2021/2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 2000000)

data_meta <- data_meta_raw %>% filter(sequenced_by_us) %>% filter(!is.na(`Report date`)) %>% filter(!is.na(case_id))
## Unvac VOC and B.1 cases in HK 
data_meta_unvac <- data_meta %>% filter(is.na(`Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)) %>% filter(lineage!="None") %>% filter(lineage!="Unassigned") %>% filter(coverage>=0.9)

# 20I(Alpha, V1) or B.1.1.7
# 20H(Beta, V2) or B.1.351
# 20J(Gamma, V3) or P.1
# 21A(Delta) or B.1.617.2
# 21M(Omicron): 21K(Omicron) or BA.1, 21L(Omicron) or BA.2, 22A(Omicron) or BA.4, 22B(Omicron) or BA.5, 22C(Omicron) or BA.2.12.1, 22D(Omicron) or BA.2.75

data_meta_unvac <- data_meta_unvac %>% filter(!(grepl("Omicron", nextstrain_lineage) & Classification=="Imported")) # filter out Omicron imported cases, to reduce the bias from re-infection. "Reinfection was found in 26 (0.46%) of 5554 Alpha, 209 (1.16%) of 17,941 Delta, and 520 (13.0%) of 3992 Omicron variants." -- https://link.springer.com/article/10.1007/s11845-022-03060-4
data_meta_unvac <- data_meta_unvac %>% filter(lineage %in% c("B.1.36", "B.1.1.63", "B.1.36.27") | nextstrain_lineage=="20I (Alpha,V1)" | nextstrain_lineage=="21J (Delta)" | nextstrain_lineage=="21M (Omicron)" | nextstrain_lineage=="22B (Omicron)") 

source("./helper/get_100reads_coverage.R")
data_meta_unvac$coverage_100reads <- get_100reads_coverage(data_meta_unvac$Sample)
data_meta_unvac <- data_meta_unvac %>% filter(coverage_100reads>=0.9)

sort(table(data_meta_unvac$lineage))
sort(table(data_meta_unvac$nextstrain_lineage))
table(data_meta_unvac$nextstrain_lineage, data_meta_unvac$lineage) # Finally, we have 20I (Alpha, B.1.1.7), 21J (Delta, B.1.617.2.*), Omicron (Local, BA.2.*), and 20B (B.1.1.63), 20A (B.1.36.*) for unvaccinated cases.
table(data_meta_unvac$nextstrain_lineage, data_meta_unvac$Classification=="Imported") 

## Vaccinated cases
data_meta_vac <- data_meta %>% filter(!is.na(`Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)) %>% filter(lineage!="None") %>% filter(lineage!="Unassigned") %>% filter(coverage>=0.9) 
data_meta_vac <- data_meta_vac %>% filter(!(grepl("Omicron", nextstrain_lineage) & Classification=="Imported")) # filter out Omicron imported cases, to reduce the bias from re-infection.

data_meta_vac$coverage_100reads <- get_100reads_coverage(data_meta_vac$Sample)
data_meta_vac <- data_meta_vac %>% filter(coverage_100reads>=0.9)

data_meta_vac <- data_meta_vac %>% filter(nextstrain_lineage=="21J (Delta)" | nextstrain_lineage=="21M (Omicron)" | nextstrain_lineage=="22B (Omicron)") 
# sort(table(data_meta_vac$lineage))
# sort(table(data_meta_vac$nextstrain_lineage))
# table(data_meta_vac$nextstrain_lineage, data_meta_vac$lineage) # Finally, we have 21J (Delta, B.1.617.2.*), Omicron (Local, BA.2.*), and 22B (Local, Omicron, BA.5.*) cases.
table(data_meta_vac$nextstrain_lineage, data_meta_vac$Classification=="Imported") # Make sure that none of the Omicron cases are imported.


# ### clean the vaccination section
# table(data_meta_vac$`Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)
# table(data_meta_vac$`Name of 2nd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)
# table(data_meta_vac$`Name of 3rd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`)

## Merge metadata
## all selected cases
df_meta <- bind_rows(data_meta_vac, data_meta_unvac)
samples_primer_old <- readLines("../data/samples_primer1.txt")
### primer and Ct value to metadata
df_meta$primer <- ifelse(df_meta$Sample %in% samples_primer_old, "old", "new")
## clean metadata
df_meta$Classification_sim <- ifelse(df_meta$Classification=="Imported", "Imported", "Local")
df_meta$sample <- df_meta$Sample
df_meta$lineage_sim <- df_meta$nextstrain_lineage
df_meta$lineage_sim <- gsub(",V1", "", df_meta$lineage_sim)
df_meta$lineage_sim[df_meta$lineage_sim=="20A"] <- "20A (B.1.36.*)"
df_meta$lineage_sim[df_meta$lineage_sim=="20B"] <- "20B (B.1.1.63)"
df_meta$lineage_sim[df_meta$lineage_sim=="21J (Delta)"] <- "21M (Delta, B.1.617.2.*)"
df_meta$lineage_sim[df_meta$lineage_sim=="21M (Omicron)"] <- "21M (Omicron, BA.2.*)"
df_meta$lineage_sim[df_meta$lineage_sim=="22B (Omicron)"] <- "22B (Omicron, BA.5.*)"
table(df_meta$lineage_sim)

df_meta$Doses <- (!is.na(df_meta$`Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`))+
(!is.na(df_meta$`Name of 2nd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`))+
(!is.na(df_meta$`Name of 3rd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`))

table(df_meta$lineage_sim, df_meta$Doses)

df_meta <- df_meta %>% filter(Doses==0 | Doses==2 | (Doses==3 & grepl("Omicron", lineage_sim))) # remove partial vaccination

df_meta$Vaccine <- df_meta %>% 
	select(c("Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)", 
	"Name of 2nd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)",
	"Name of 3rd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)")) %>% 
	apply(1, function(x){
		tmp <- unique(tolower(x[!is.na(x)]))
		if(length(tmp)==0){return("Non-vaccinated")}
		paste(tmp, collapse=" | ")
	})
sort(table(df_meta$Vaccine))
df_meta$Vaccine[df_meta$Vaccine=="biontech/fosun"] <- "biontech"
df_meta <- df_meta %>% filter(Vaccine %in% c("biontech", "sinovac", "Non-vaccinated")) # the number of others cases are too few
vaccine_of_interest <- c("BioNTech", "Sinovac", "Non-vaccinated")
df_meta$Vaccine <- factor(df_meta$Vaccine, levels=c("biontech", "sinovac", "Non-vaccinated"), labels=c("BioNTech", "Sinovac", "Unvaccinated"))

df_date_lag <- df_meta %>% group_by(Sample) %>% select(`Report date`, contains("(date)")) %>% unique() %>% ungroup() 
df_date_lag$days_since_last_dose <- df_date_lag %>% select(-Sample) %>% apply(1, function(x){
	tmp <- lubridate::ymd(x)
	if(all(is.na(tmp[2:length(tmp)]))){
		return(NA)
	} else {
		check <- which(!is.na(tmp))
		abs(as.integer(tmp[1]-tmp[check[length(check)]]))
	}
})
df_meta <- left_join(df_meta, df_date_lag %>% select(Sample, days_since_last_dose), "Sample")

# filtering by throughput
# add stats from fastqc
source("./helper/get_fastqc_stat.R")
df_meta$total_seqs_after_trimming <- get_total_seqs_after_trimming(df_meta$Sample)

sum(as.numeric(df_meta$total_seqs_after_trimming)<50000, na.rm=T)
df_meta <- df_meta %>% filter(Ct_value<=28 | is.na(Ct_value)) # remove all samples with Ct_value > 28

## compare to version 1
data_meta_study1$Vaccine <- factor(data_meta_study1$Vaccine, levels=c("BioNTech", "Sinovac", "Non-vaccinated"))

table(df_meta$lineage_sim, df_meta$Vaccine)
table(df_meta$lineage_sim, df_meta$Vaccine, df_meta$Doses)
sum(table(df_meta$lineage_sim, df_meta$Vaccine, df_meta$Doses))
table(data_meta_study1$lineage_sim, data_meta_study1$Vaccine, data_meta_study1$Doses)

# difference in number of analysed samples
sum(table(df_meta$lineage_sim, df_meta$Vaccine, df_meta$Doses))-sum(table(data_meta_study1$lineage_sim, data_meta_study1$Vaccine))
# number of newly added vaccinated samples
sum(table(df_meta$lineage_sim, df_meta$Vaccine)[,1:2])-sum(table(data_meta_study1$lineage_sim, data_meta_study1$Vaccine)[,1:2])
# number of newly added 3-doses vaccinated samples
sum(df_meta$Doses==0)
sum(df_meta$Doses==2)
sum(df_meta$Doses==3)
# number of newly added vaccinated Omicron samples
sum(table(df_meta$lineage_sim, df_meta$Vaccine)[5:6,1:2])-sum(table(data_meta_study1$lineage_sim, data_meta_study1$Vaccine)[6,1:2])
# number of newly added vaccinated Delta samples
sum(table(df_meta$lineage_sim, df_meta$Vaccine)[4,1:2])-sum(table(data_meta_study1$lineage_sim, data_meta_study1$Vaccine)[5,1:2])
# number of newly added Alpha samples
sum(table(df_meta$lineage_sim, df_meta$Vaccine)[3,])-sum(table(data_meta_study1$lineage_sim, data_meta_study1$Vaccine)[1,])
# difference in number of B.1.36.* samples
sum(table(df_meta$lineage_sim, df_meta$Vaccine)[1,])-sum(table(data_meta_study1$lineage_sim, data_meta_study1$Vaccine)[3:4,])
# difference in number of B.1.1.63 samples
sum(table(df_meta$lineage_sim, df_meta$Vaccine)[2,])-sum(table(data_meta_study1$lineage_sim, data_meta_study1$Vaccine)[2,])

write_csv(df_meta, "../results/df_samples.csv")
