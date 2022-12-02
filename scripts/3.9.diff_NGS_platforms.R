library(tidyverse)
library(parallel)
library(ssh.utils)
source("./helper/pysamstats.r")

data_meta_raw_ori <- read_csv("../../../2021/2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 200000)
data_meta_raw <- data_meta_raw_ori %>% filter(sequenced_by_us) %>% filter(!is.na(`Report date`))

# data_meta_study <- read_csv("../data/df_samples.csv", guess_max = 20000)
data_meta_study <- read_csv("../results/df_samples.csv", guess_max=100000)
data_meta_study <- data_meta_study %>% filter(lineage_sim != "22B (Omicron, BA.5.*)")
table(grepl("iseq", data_meta_study$sample))/nrow(data_meta_study)
table(grepl("iseq", data_meta_study$sample), data_meta_study$lineage_sim)

load("../results/df_plot_n_gene_adj.rdata")
df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, data_meta_study %>% select(sample, Doses), "sample")
df_plot_n_gene_meta_adj$vaccine_doses <- paste0(df_plot_n_gene_meta_adj$Vaccine, "\n(Doses=", df_plot_n_gene_meta_adj$Doses, ")")
df_plot_n_gene_meta_adj$vaccine_doses[df_plot_n_gene_meta_adj$Vaccine=="Unvaccinated"] <- "Unvaccinated"
df_plot_n_gene_meta_adj <- df_plot_n_gene_meta_adj %>% mutate(iseq=grepl("iseq", sample))
df_tmp <- df_plot_n_gene_meta_adj %>% filter(lineage_sim == "21M (Omicron, BA.2.*)") %>% filter(gene=="Full genome")
quantile(df_tmp$n_per_kb_adj[df_tmp$iseq])
quantile(df_tmp$n_per_kb_adj[!df_tmp$iseq])

table(df_tmp$vaccine_doses, df_tmp$iseq)

(unq_vaccine_doses <- sort(unique(df_tmp$vaccine_doses)))

wilcox.test(df_tmp$n_per_kb_adj[df_tmp$iseq & (df_tmp$Vaccine==unq_vaccine_doses[5])], df_tmp$n_per_kb_adj[(!df_tmp$iseq) & df_tmp$Vaccine==unq_vaccine_doses[5]])
wilcox.test(df_tmp$n_per_kb_adj[df_tmp$iseq & (df_tmp$vaccine_doses==unq_vaccine_doses[1])], df_tmp$n_per_kb_adj[(!df_tmp$iseq) & df_tmp$vaccine_doses==unq_vaccine_doses[1]])
wilcox.test(df_tmp$n_per_kb_adj[df_tmp$iseq & (df_tmp$vaccine_doses==unq_vaccine_doses[2])], df_tmp$n_per_kb_adj[(!df_tmp$iseq) & df_tmp$vaccine_doses==unq_vaccine_doses[2]])
wilcox.test(df_tmp$n_per_kb_adj[df_tmp$iseq & (df_tmp$vaccine_doses==unq_vaccine_doses[3])], df_tmp$n_per_kb_adj[(!df_tmp$iseq) & df_tmp$vaccine_doses==unq_vaccine_doses[3]])
wilcox.test(df_tmp$n_per_kb_adj[df_tmp$iseq & (df_tmp$vaccine_doses==unq_vaccine_doses[4])], df_tmp$n_per_kb_adj[(!df_tmp$iseq) & df_tmp$vaccine_doses==unq_vaccine_doses[4]])

ggplot(df_tmp)+geom_boxplot(aes(x=vaccine_doses, y=n_per_kb_adj, color=iseq))
ggsave("../results/tmp.pdf")

# ggplot(df_tmp)+geom_boxplot(aes(x=vaccine_doses, y=n_per_kb, color=iseq))

# select samples sequenced by both platforms
data_meta_raw <- data_meta_raw %>% filter(Sample %in% data_meta_study$Sample)
data_meta_raw <- data_meta_raw %>% filter(grepl(",", samples_all, fixed=T))
data_meta_raw$num_whp <- sapply(data_meta_raw$samples_all, function(x){
	tmp <- strsplit(x, ", ")[[1]]
	sum(grepl("WHP", tmp))
})
data_meta_raw <- data_meta_raw %>% filter(num_whp>1)

# samples different platforms
check_both_platforms <- mapply(function(x, y){
	tmp <- strsplit(x, ", ")[[1]]
	tmp <- gsub("_\\d+$", "", tmp)
	tmp <- tmp[grepl("WHP", tmp)]

	whp_id_ref <-  gsub("_\\d+$", "", y)
	whp_id_ref <- gsub("RE-", "", whp_id_ref)
	whp_id_ref <- gsub("-.+", "", whp_id_ref)

	tmp <- tmp[grepl(whp_id_ref, tmp)]

	check_iseq <- grepl(paste0(whp_id_ref, "-S\\d+-iseq"), tmp)
	check_novaseq <- !grepl("-iseq", tmp)

	if(any(check_iseq) && any(check_novaseq)){return(c(tmp[check_novaseq][1], tmp[check_iseq][1]))}else{return(NA)}
}, data_meta_raw$samples_all, data_meta_raw$Sample, USE.NAMES=F)
samples_both_platforms <- check_both_platforms[!is.na(check_both_platforms)]

# read bamreadcounts
samples_toupload_ori <- unlist(samples_both_platforms)
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
	if(length(samples_not_available)>0){samples_both_platforms <- samples_both_platforms[!sapply(samples_both_platforms, function(x) {any(x%in%samples_not_available)})]}
	## copy the files to remote
	sapply(c(rbind(files_fastq1, files_fastq2)), function(x){
		cp.remote(path.src = x, remote.src = "", remote.dest = "hggu@147.8.70.166", path.dest = "~/work/2020-09-01_COVID_NGS_pipeline/NGS_data_input/", verbose = T)
	})
}

# after running on the server
files_bam <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
files_bam_full <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$", full.names = T)
files_bam_rst_full <- c(list.files("../results/pysamstats/", "tsv$", full.names = T), list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", "tsv$", full.names = T), list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/", "tsv$", full.names = T))
idx <- sapply(paste0("/", unlist(samples_both_platforms), ".tsv"), function(x) {
	tmp <- grep(x, files_bam_rst_full)[1]
})
sum(is.na(idx))
if(sum(is.na(idx))>0){
	samples_no_bamreadcount <- unlist(samples_both_platforms)[is.na(idx)]
	files_bam_int <- files_bam_full[(gsub("-trimmed.+", "", files_bam)) %in% samples_no_bamreadcount]
	process_pysamstats(files_bam_int, n_cores=8)
}
stopifnot(sum(is.na(idx))==0)

df_bam_rst <- mclapply(files_bam_rst_full[idx], read_pysamstats, pp=FALSE, mc.cores=16)
df_bam_rst <- bind_rows(df_bam_rst[!is.na(df_bam_rst)])
unq_samples_bam <- unique(df_bam_rst$sample)
check_completed <- sapply(samples_both_platforms, function(x){ # check whether both samples have enough coverage (default: genome coverage > 90%)
	sum(x%in%unq_samples_bam)==2
})
samples_both_platforms <- samples_both_platforms[check_completed] # we identified 77 samples sequenced by both Novaseq and iSeq platforms having adequate quality for downstream comparison.

df_bam_rst <- df_bam_rst %>% filter(sample %in% unlist(samples_both_platforms))
df_bam_rst$muAF <- df_bam_rst$matches/df_bam_rst$reads_all
save(df_bam_rst, file="../results/df_bam_rst_full_notpp_platform_comp.rdata")
# load("../results/df_bam_rst_full_notpp_platform_comp.rdata")

MAF_threshold <- 0.05
MAF_threshold <- 0.025
df_bam_rst_filter <- df_bam_rst %>% filter((muAF<(1-MAF_threshold)) & (muAF>MAF_threshold))

df_extra_info <- extra_info_from_pysamstats(df_bam_rst_filter)
df_extra_info <- df_extra_info %>% mutate_at(vars(!contains(c("base", "sample"))), as.numeric)
df_extra_info$whp_id <- gsub("-.+$", "", df_extra_info$sample)
df_extra_info <- left_join(df_extra_info, data_meta_study %>% mutate(whp_id = gsub("-.+$", "", sample)) %>% select(whp_id, primer, lineage_sim, Vaccine), "whp_id")

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

prop_shared <- sapply(sort(sapply(samples_both_platforms, function(x)x[1])), function(whp_id_i) { # calculate proportion of the shared mutations
	# whp_id_i="WHP2633"
	df_tmp <- df_extra_info %>% filter(whp_id == whp_id_i)
	if(nrow(df_tmp)==0){return(1)}
	unq_samples <- unique(df_tmp$sample)	
	if(length(unq_samples)==1){return(0)}
	mut_1 <- df_tmp %>% filter(sample==unq_samples[1]) %>% transmute(mutations=paste0(con_base, pos, sec_base)) %>% .$mutations
	mut_2 <- df_tmp %>% filter(sample==unq_samples[2]) %>% transmute(mutations=paste0(con_base, pos, sec_base)) %>% .$mutations

	num_total_muts <- length(unique(c(mut_1, mut_2)))
	num_shared_muts <- sum(table(c(mut_1, mut_2))==2)
	return(num_shared_muts/num_total_muts)
})

quantile(prop_shared)
mean(prop_shared) # 0.60 for MAF0.05; 0.37 for MAF0.025
median(prop_shared) # 0.64 for MAF0.05; 0.40 for MAF0.025


df_extra_info$iseq <- grepl("iseq", df_extra_info$sample)

df_extra_info_n <- df_extra_info %>% group_by(sample, iseq, lineage_sim, Vaccine) %>% summarise(N=n()) 
table(df_extra_info_n$lineage_sim, df_extra_info_n$Vaccine)

ggplot(df_extra_info_n)+
	geom_boxplot(aes(x=lineage_sim,y=N, color=iseq))
