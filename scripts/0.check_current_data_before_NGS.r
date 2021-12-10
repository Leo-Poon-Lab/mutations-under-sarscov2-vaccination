library(tidyverse)
library(readxl)
library(ssh.utils)

data_cov <- read_csv("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/coverage.csv")
df_choose <- read_excel("../data/Vaccine Delta SNP list.xlsx", skip=1)
WHP_id <- df_choose$WHP_ID
WHP_id <- WHP_id[!is.na(WHP_id)]
WHP_id <- WHP_id[WHP_id!="WHP_ID"]

length(WHP_id)
check <- sapply(data_cov$Sample, function(x){
	any(grepl(x, WHP_id)) | any(sapply(WHP_id, function(y){grepl(y, x)}))
})
data_cov_filter <- data_cov[check,] %>% arrange(Sample, genome_coverage_over_5)
write_csv(data_cov_filter, "../results/data_cov_pre-exisit.csv")
