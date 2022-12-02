library(tidyverse)
library(parallel)
library(ggsci)

data_meta_raw_ori <- read_csv("../../../2021/2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 200000)
data_meta_raw <- data_meta_raw_ori %>% filter(sequenced_by_us) %>% filter(!is.na(`Report date`))

# data_meta_study <- read_csv("../data/df_samples.csv", guess_max = 20000)
data_meta_study <- read_csv("../results/df_samples.csv", guess_max=100000)
data_meta_study <- data_meta_study %>% filter(lineage_sim != "22B (Omicron, BA.5.*)")

# logitudinal samples 1
df_voc_logitudinal <- readxl::read_excel("../data/VOC longi samples for Iseq.xlsx")
df_voc_logitudinal$Lab_ID <- gsub("VOC0", "VOC", df_voc_logitudinal$Lab_ID)
df_voc_logitudinal <- df_voc_logitudinal %>% filter(Lab_ID %in% paste0("VOC", c(226, 257, 256, 254, 247, 243)))
df_voc_logitudinal$sample <- paste0(df_voc_logitudinal$Lab_ID, "-H2-iseq")
df_voc_logitudinal <- df_voc_logitudinal %>% group_by(Case) %>% filter(n()!=1) %>% ungroup()

# logitudinal samples 2
samples_controls <- data_meta_raw$samples_all
samples_controls <- unlist(strsplit(samples_controls, ", "))
samples_controls <- samples_controls[grepl("WHP", samples_controls)]
samples_controls <- gsub("_\\d+$", "", samples_controls)

samples_controls_whp_id <- gsub("RE-", "", samples_controls)
samples_controls_whp_id <- gsub("-.+", "", samples_controls_whp_id)
duplicated_samples_whpids <- names(table(samples_controls_whp_id)[table(samples_controls_whp_id)>1])
samples_controls <- samples_controls[samples_controls_whp_id %in% duplicated_samples_whpids]

df_samples_controls <- tibble(sample=samples_controls, whp_id=gsub("[-_].+", "", samples_controls))
df_samples_controls$sequencer <- ifelse(grepl("-S",df_samples_controls$sample), "iSeq", "Novaseq")
df_vm_id <- read_tsv("../../../2021/2021-06-24_merge_metadata/results/data_seqed_vm_id.tsv")
df_samples_controls <- left_join(df_samples_controls, df_vm_id, "whp_id")

df_samples_controls_summary <- df_samples_controls %>% group_by(case_id) %>% summarise(n=n(), n_sequencer=length(unique(sequencer)), diff_collection_date=max(collection_date)-min(collection_date)) 
df_samples_controls_summary %>% filter(n_sequencer>1)
df_samples_controls_summary %>% filter(diff_collection_date>0) %>% left_join(data_meta_study  %>% mutate(case_id = as.character(case_id)) %>% select(case_id, Vaccine, lineage_sim))
df_samples_controls_summary %>% filter(diff_collection_date>0) %>% left_join(data_meta_study  %>% mutate(case_id = as.character(case_id)) %>% select(case_id, Vaccine, lineage_sim)) %>% filter(diff_collection_date<=2)
case_id_longitudinal <- df_samples_controls_summary %>% filter(diff_collection_date>0) %>% left_join(data_meta_study  %>% mutate(case_id = as.character(case_id)) %>% select(case_id, Vaccine, lineage_sim)) %>% filter(!is.na(Vaccine)) %>% .$case_id

df_whp_logitudinal <- df_samples_controls %>% filter(case_id %in% case_id_longitudinal) %>% filter(sequencer=="Novaseq")

samples_longitudinal <- c(df_voc_logitudinal %>% group_by(Case) %>% summarise(sample_list=list(sample)) %>% .$sample_list, df_whp_logitudinal %>% group_by(case_id) %>% summarise(sample_list=list(sample)) %>% .$sample_list)

source("./helper/pysamstats.r")
# read bamreadcounts of the control samples

samples_toupload_ori <- unlist(samples_longitudinal)
samples_alrd_done <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
samples_alrd_done <- gsub("-trimmed.+", "", samples_alrd_done)
samples_toupload <- samples_toupload_ori[!samples_toupload_ori %in% samples_alrd_done] # we alredy analysis some of the data

if(length(samples_toupload)>0){
	files_fastq_all <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive/", full.names = T)
	samples_all <- sapply(files_fastq_all, function(x){
		tmp <- strsplit(x, "/", fixed=T)[[1]]
		sample_t <- tmp[length(tmp)]
		gsub(".fastq.gz", "", sample_t)
	})
	files_fastq1 <- files_fastq_all[samples_all %in% paste0(samples_toupload, "_1")]
	files_fastq2 <- files_fastq_all[samples_all %in% paste0(samples_toupload, "_2")]

	samples_not_available <- samples_toupload[!samples_toupload %in% gsub(".+/", "", gsub("_1.fastq.gz", "", files_fastq1))]
	if(length(samples_not_available)>0){samples_longitudinal <- samples_longitudinal[!sapply(samples_longitudinal, function(x) {any(x%in%samples_not_available)})]}
	## copy the files to remote
	sapply(c(rbind(files_fastq1, files_fastq2)), function(x){
		cp.remote(path.src = x, remote.src = "", remote.dest = "hggu@147.8.70.166", path.dest = "~/work/2020-09-01_COVID_NGS_pipeline/NGS_data_input/", verbose = T)
	})
}

# after running on the server
files_bam <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
files_bam_full <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$", full.names = T)
files_bam_rst_full <- c(list.files("../results/pysamstats/", "tsv$", full.names = T), list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", "tsv$", full.names = T), list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/", "tsv$", full.names = T))
idx <- sapply(paste0("/", unlist(samples_longitudinal), ".tsv"), function(x) {
	tmp <- grep(x, files_bam_rst_full)[1]
})
sum(is.na(idx))
if(sum(is.na(idx))>0){
	samples_no_bamreadcount <- unlist(samples_longitudinal)[is.na(idx)]
	files_bam_int <- files_bam_full[(gsub("-trimmed.+", "", files_bam)) %in% samples_no_bamreadcount]
	process_pysamstats(files_bam_int, n_cores=8)
}
stopifnot(sum(is.na(idx))==0)

df_bam_rst <- mclapply(files_bam_rst_full[idx], read_pysamstats, pp=FALSE, mc.cores=16)
df_bam_rst <- bind_rows(df_bam_rst[!is.na(df_bam_rst)])
unq_samples_bam <- unique(df_bam_rst$sample)
check_completed <- sapply(samples_longitudinal, function(x){ # check whether both samples have enough coverage (default: genome coverage > 90%)
	sum(x%in%unq_samples_bam)>1
})
samples_longitudinal <- samples_longitudinal[check_completed] # we identified 6 cases with longitudinal samples > 1

df_bam_rst <- df_bam_rst %>% filter(sample %in% unlist(samples_longitudinal))
df_bam_rst$muAF <- df_bam_rst$matches/df_bam_rst$reads_all
save(df_bam_rst, file="../results/df_bam_rst_full_notpp_longitudinal.rdata")
# load("../results/df_bam_rst_full_notpp_longitudinal.rdata")

# compare consensus_base
df_con_prop_shared <- lapply(samples_longitudinal, function(whp_id_i) { # calculate proportion of the shared mutations
	# whp_id_i=c("WHP4977", "WHP4988")
	print(whp_id_i)
	df_tmp <- df_bam_rst %>% filter(sample %in% whp_id_i) %>% filter(reads_all>=10)
	if(nrow(df_tmp)==0){return(NA)}
	unq_samples <- unique(df_tmp$sample)	
	df_tmp <- extra_info_from_pysamstats(df_tmp)
	if(length(unq_samples)==1){return(NA)}
	df_tmp <- left_join(df_tmp %>% filter(sample==unq_samples[1]) %>% select(pos, con_base) %>% unique(), df_tmp %>% filter(sample==unq_samples[2]) %>% select(pos, con_base) %>% unique(), "pos")
	
	return(tibble(sample_1=unq_samples[1], sample_2=unq_samples[2], con_base_eq=sum(df_tmp$con_base.x==df_tmp$con_base.y,na.rm=T), con_base_diff=sum(df_tmp$con_base.x!=df_tmp$con_base.y,na.rm=T)))
})
df_con_prop_shared <- bind_rows(df_con_prop_shared)

# MAF_threshold <- 0.05
MAF_threshold <- 0.025
df_bam_rst_filter <- df_bam_rst %>% filter((muAF<(1-MAF_threshold)) & (muAF>MAF_threshold))

names(df_voc_logitudinal)[1:4] <- c("case_id", "whp_id", "days", "Ct_value")
df_voc_logitudinal$case_id <- as.character(df_voc_logitudinal$case_id)
df_voc_logitudinal$sequencer <- "iSeq"
func_days <- function(x){as.numeric(x-min(x))}
df_whp_logitudinal <- df_whp_logitudinal %>% group_by(case_id) %>% mutate(days=func_days(collection_date))
df_all_logitudinal <- bind_rows(df_voc_logitudinal, df_whp_logitudinal)

df_extra_info <- extra_info_from_pysamstats(df_bam_rst_filter)
df_extra_info <- df_extra_info %>% mutate_at(vars(!contains(c("base", "sample"))), as.numeric)
df_extra_info$whp_id <- gsub("-.+$", "", df_extra_info$sample)
stopifnot(all(df_extra_info$whp_id %in% df_all_logitudinal$whp_id))
df_extra_info <- left_join(df_extra_info, df_all_logitudinal %>% select(whp_id, case_id, days, Ct_value), "whp_id")
df_extra_info <- left_join(df_extra_info, data_meta_study %>% select(case_id, lineage_sim, Vaccine), "case_id")

df_extra_info$primer <- "new"

### filtering (QC) by position
#### 1. positions 1:100 and (29903-99):29903 should be removed; 
#### 2. exclude all positions in the PCR primer binding regions
source("./helper/isnv_position_filter.R")
df_extra_info <- filter_by_pos(df_extra_info)

### MAF threshold
df_extra_info <- df_extra_info %>% filter(sec_freq>=MAF_threshold) # MAF>=0.025

### Depth threshold
df_extra_info <- df_extra_info %>% filter(depth>=100) # Depth

### filter structural variants
df_extra_info <- df_extra_info %>% filter(nchar(con_base)==1) # only single-nt mutations
df_extra_info <- df_extra_info %>% filter(nchar(sec_base)==1) # only single-nt mutations
length(unique(df_extra_info$sample))

#### The identified mutations/variants should be supported by at least one forward and one reverse read); 
# df_extra_info <- df_extra_info %>% filter((as.numeric(sec_fwd)+as.numeric(sec_rev))>=5 & as.numeric(sec_fwd)>=1 & as.numeric(sec_rev)>=1)
df_extra_info <- df_extra_info %>% filter(as.numeric(sec_fwd)>=1 & as.numeric(sec_rev)>=1)
#### filter strand bias
df_extra_info <- df_extra_info %>% filter(strand_bias<10) 

#### filter serial adjacent disjoint mutations
fn_check_serial_mut <- function(x){
	if(length(x)<3) return(rep(FALSE, length(x)))
	rst <- rep(FALSE, length(x))
	for (i in seq_len(length(x)-2)) {
		if(x[(i+2)] <= (x[i]+30-1)) {
			rst[i:(i+2)] <- TRUE
		}
	}
	rst
}
# fn_check_serial_mut(c(7102,7123,7142,7150,7219,7234,7319,7384))
df_extra_info <- df_extra_info %>% arrange(pos) %>% group_by(sample) %>% mutate(check_serial_mut=fn_check_serial_mut(pos)) %>% ungroup() %>% filter(!check_serial_mut) # remove serial adjacent separate mutations

df_prop_shared <- lapply(samples_longitudinal, function(whp_id_i) { # calculate proportion of the shared mutations
	# whp_id_i="WHP4845"
	df_tmp <- df_extra_info %>% filter(whp_id %in% gsub("-.+$", "", whp_id_i)) %>% arrange(days)
	if(nrow(df_tmp)==0){return(1)}
	unq_samples <- unique(df_tmp$sample)	
	if(length(unq_samples)==1){return(0)}
	mut_1 <- df_tmp %>% filter(sample==unq_samples[1]) %>% transmute(mutations=paste0(con_base, pos, sec_base)) %>% .$mutations
	day_1 <- df_all_logitudinal %>% filter(sample==unq_samples[1]) %>% .$days %>% unique()
	ct_1 <- df_all_logitudinal %>% filter(sample==unq_samples[1]) %>% .$Ct_value %>% unique()
	mut_2 <- df_tmp %>% filter(sample==unq_samples[2]) %>% transmute(mutations=paste0(con_base, pos, sec_base)) %>% .$mutations
	day_2 <- df_all_logitudinal %>% filter(sample==unq_samples[2]) %>% .$days %>% unique()
	ct_2 <- df_all_logitudinal %>% filter(sample==unq_samples[2]) %>% .$Ct_value %>% unique()

	num_total_muts <- length(unique(c(mut_1, mut_2)))
	num_shared_muts <- sum(table(c(mut_1, mut_2))==2)
	return(tibble(sample_1=unq_samples[1], sample_2=unq_samples[2], day_1=day_1, day_2=day_2, days_diff=day_2-day_1, ct_1=ct_1, ct_2=ct_2, num_shared_muts=num_shared_muts, num_mut_1=length(mut_1), num_mut_2=length(mut_2), num_shared_prop=num_shared_muts/num_total_muts))
})
df_prop_shared <- bind_rows(df_prop_shared)

df_prop_shared # 5 out of 6 cases does not have shared iSNVs
write_csv(df_prop_shared, "../results/df_prop_shared_longitudinal.csv")
