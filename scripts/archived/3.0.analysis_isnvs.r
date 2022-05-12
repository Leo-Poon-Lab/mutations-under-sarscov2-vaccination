# analyse the mutations found in the selected samples
## 1. read the snvs identified by ivar.
## 2. filter out the mutations with low quality/support
## 3. determine the regions with no data / no support for each sample.
## 4. analyse the number of mutations, the frequency of mutations, diversity.

library(tidyverse)
library(Biostrings)
library(ggsci)
library(patchwork)
library(parallel)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")

df_meta <- read_csv("../results/df_samples.csv")

table(df_meta$lineage)
table(df_meta$variant)
df_meta$lineage_sim <- df_meta$lineage
df_meta$lineage_sim[grepl("^AY\\.", df_meta$lineage_sim)] <- "Delta"
df_meta$lineage_sim[df_meta$lineage_sim == "B.1.617.2"] <- "Delta"
df_meta$lineage_sim[grepl("^B\\.1$", df_meta$lineage_sim)] <- "B.1"
df_meta$lineage_sim[grepl("^B\\.1\\.1\\.7$", df_meta$lineage_sim)] <- "Alpha"
df_meta$lineage_sim[grepl("^Q\\.", df_meta$lineage_sim)] <- "Alpha"
df_meta$lineage_sim[grepl("Omicron", df_meta$variant)] <- "Omicron"
df_meta$lineage_sim[grepl("Beta", df_meta$variant)] <- "Beta"
table(df_meta$lineage_sim)

lineage_int <- c("Alpha", "Delta", "Omicron", "B.1.1.63", "B.1.36.27", "B.1.36")
# lineage_int <- c("Alpha", "Delta", "Omicron", "B.1.36")
df_meta <- df_meta %>% filter(lineage_sim %in% lineage_int)

colors_vaccine=pal_jco()(3)
names(colors_vaccine)=c("BioNTech", "Sinovac", "Non-vaccinated")

df_meta$Doses <- (!is.na(df_meta$`Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`))+
(!is.na(df_meta$`Name of 2nd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`))+
(!is.na(df_meta$`Name of 3rd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)`))
df_meta$Vaccine <- df_meta %>% 
	select(c("Name of 1st vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)", 
	"Name of 2nd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)",
	"Name of 3rd vaccine (AstraZeneca, BioNTech, Sinovac, others pls specify)")) %>% 
	apply(1, function(x){
		tmp <- unique(tolower(x[!is.na(x)]))
		if(length(tmp)==0){return("Non-vaccinated")}
		paste(tmp, collapse=" ")
	})
sort(table(df_meta$Vaccine))
df_meta <- df_meta %>% filter(Vaccine %in% c("biontech", "sinovac", "Non-vaccinated")) # the number of others cases are too few
vaccine_of_interest <- c("BioNTech", "Sinovac", "Non-vaccinated")
df_meta$Vaccine <- factor(df_meta$Vaccine, levels=c("biontech", "sinovac", "Non-vaccinated"), labels=c("BioNTech", "Sinovac", "Non-vaccinated"))

df_meta <- df_meta %>% filter(Ct_value<25)
table(df_meta$Vaccine, df_meta$lineage_sim, df_meta$Doses)
table(df_meta$Vaccine)
table(df_meta$lineage_sim, df_meta$Vaccine)
df_date_lag <- df_meta %>% group_by(Sample) %>% select(`Report date`, contains("(date)")) %>% unique() %>% ungroup() 
df_date_lag$days_since_last_dose <- df_date_lag %>% select(-Sample) %>% apply(1, function(x){
	tmp <- lubridate::dmy(x)
	if(all(is.na(tmp[2:length(tmp)]))){
		return(NA)
	} else {
		check <- which(!is.na(tmp))
		abs(as.integer(tmp[1]-tmp[check[length(check)]]))
	}
})
df_meta <- left_join(df_meta, df_date_lag %>% select(Sample, days_since_last_dose), "Sample")
write_csv(df_meta, "../results/df_samples_clean.csv")

## 1. read the snvs identified by ivar.
files_snvs <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/variant_caller/ivar/")
files_snvs_full <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/variant_caller/ivar/", full.names = T)
files_alrd_done <- gsub("\\.tsv\\.vcf$", ".tsv", list.files("../results/vcf/", "\\.tsv\\.vcf$"))
check1 <- files_snvs %in% paste0("ivar_", df_meta$Sample, ".tsv")
check2 <- !files_snvs %in% files_alrd_done
files_snvs_full_int <- files_snvs_full[check1 & check2]
source("./helper/convert_to_vcf.r")

convert_to_vcf(files_snvs_full_int)
files_vcfs <- list.files("../results/vcf/", "vcf$", full.names = T)
files_csv <- list.files("../results/vcf/", "csv$", full.names = T)
files_vcfs_todo <- files_vcfs[!gsub("\\.tsv\\.vcf$", "", files_vcfs) %in% gsub("\\.tsv\\.vcf\\.csv$", "", files_csv)]
annotate_snpeff(files_vcfs_todo)
files_vcfs_csv <- list.files("../results/vcf/", "csv$", full.names = T)
idx <- sapply(paste0(df_meta$Sample, ".tsv"), function(x) {
	tmp <- grep(x, files_vcfs_csv)[1]
})
files_vcfs_csv <- files_vcfs_csv[idx]
df_snvs <- mclapply(files_vcfs_csv, read_snpeff, mc.cores=16)
df_snvs <- bind_rows(df_snvs)
 
## 2. filter out the mutations with low quality/support
## 3. determine the regions with no data / no support for each sample.
### using pysamstats
system("conda run -n base pysamstats --help")
source("./helper/pysamstats.r")
files_bam <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
files_bam_full <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$", full.names = T)
samples_alrd_done <- gsub("\\.tsv$", "", list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/"))
samples_alrd_done_1 <- gsub("\\.tsv$", "", list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/"))
samples_alrd_done_2 <- gsub("\\.tsv$", "", list.files("../results/pysamstats/"))
check1 <- files_bam %in% paste0(df_meta$Sample, "-trimmed.masked.bam")
check2 <- !files_bam %in% paste0(c(samples_alrd_done, samples_alrd_done_1, samples_alrd_done_2), "-trimmed.masked.bam")
files_bam_int <- files_bam_full[check1 & check2]

mclapply(files_bam_int, function(x) {
	process_pysamstats(x)
	return("")
}, mc.cores=16)
files_bam_rst_full <- c(list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/", "tsv$", full.names = T), list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", "tsv$", full.names = T), list.files("../results/pysamstats/", "tsv$", full.names = T))
idx <- sapply(paste0(df_meta$Sample, ".tsv"), function(x) {
	tmp <- grep(x, files_bam_rst_full)[1]
})
files_bam_rst_full <- files_bam_rst_full[idx]
stopifnot(sum(is.na(files_bam_rst_full))==0)
### only nucleotide positions with at least 100 properly paired reads were considered
### only samples with 27000 sites of sequencing depth >= 100 were considered
df_bam_rst <- mclapply(files_bam_rst_full, read_pysamstats, mc.cores=16)
df_bam_rst <- bind_rows(df_bam_rst[!is.na(df_bam_rst)])
save(df_bam_rst, file="../results/df_bam_rst_full.rdata")
# plot(df_bam_rst$n_high_depth)

names(df_snvs)[names(df_snvs)=="sample"] <- "Sample"
df_snvs_meta <- left_join(df_snvs, df_meta, "Sample")
df_snvs_meta <- df_snvs_meta %>% filter(Sample %in% samples_hq)
write_csv(df_snvs_meta, "../results/df_mut_meta.csv")
# df_snvs_meta <- read_csv("../results/df_mut_meta.csv", guess_max=600000)··

### add readcounts from pysamstats results
bases <- c("deletions", "insertions", "A", "C", "T", "G")

df_snvs_meta_append <- mclapply(seq_len(nrow(df_snvs_meta)), function(i){
	print(i)
	df_snvs_meta_i <- df_snvs_meta[i,]
	sample_i <- df_snvs_meta_i$Sample
	pos_i <- df_snvs_meta_i$X2
	ref_i <- df_snvs_meta_i$X4
	alt_i <- df_snvs_meta_i$X5
	df_bam_rst_i <- df_bam_rst[df_bam_rst$sample==sample_i & df_bam_rst$pos==pos_i,]
	
	depth_i <- df_bam_rst_i$reads_all
	depth_pp_i <- df_bam_rst_i$reads_pp # count on paired reads·

	if(nchar(ref_i) > nchar(alt_i)){ # deletion
		alt_i <- "deletions"
	} else if(nchar(ref_i) < nchar(alt_i)){ # insertion
		alt_i <- "insertions"
	} 
	alt_fwd_i <- df_bam_rst_i[[paste0(alt_i, "_pp_fwd")]]
	alt_rev_i <- df_bam_rst_i[[paste0(alt_i, "_pp_rev")]]

	num_bases <- c(df_bam_rst_i$deletions_pp, df_bam_rst_i$insertions_pp, df_bam_rst_i$A_pp, df_bam_rst_i$C_pp, df_bam_rst_i$T_pp, df_bam_rst_i$G_pp)
	if(nchar(alt_i)==1){
		bases <- bases[-(1:2)]
		num_bases <- num_bases[-(1:2)]
	}
	check_consensus_base <- order(num_bases)[length(bases)]
	check_secondary_base <- order(num_bases)[length(bases)-1]
	consensus_base <- bases[check_consensus_base]
	secondary_base <- bases[check_secondary_base]

	con_fwd_i <- df_bam_rst_i[[paste0(consensus_base, "_pp_fwd")]]
	con_rev_i <- df_bam_rst_i[[paste0(consensus_base, "_pp_rev")]]
	sec_fwd_i <- df_bam_rst_i[[paste0(secondary_base, "_pp_fwd")]]
	sec_rev_i <- df_bam_rst_i[[paste0(secondary_base, "_pp_rev")]]

	alle_freq_i <- (alt_fwd_i+alt_rev_i)/depth_pp_i
	con_freq_i <- (con_fwd_i+con_rev_i)/depth_pp_i
	sec_freq_i <- (sec_fwd_i+sec_rev_i)/depth_pp_i

	c(type_alt=alt_i, depth_all=depth_i, depth_pp=depth_pp_i, depth_pp_fwd=df_bam_rst_i$reads_pp_fwd, depth_pp_rev=df_bam_rst_i$reads_pp_rev, alt_fwd=alt_fwd_i, alt_rev=alt_rev_i, alle_freq=alle_freq_i, con_base=consensus_base, con_fwd=con_fwd_i, con_rev=con_rev_i, con_freq=con_freq_i, sec_base=secondary_base, sec_fwd=sec_fwd_i, sec_rev=sec_rev_i, sec_freq=sec_freq_i)
}, mc.cores = 16)
df_snvs_meta_append <- bind_rows(df_snvs_meta_append)
df_mut_meta_add <- bind_cols(df_snvs_meta, df_snvs_meta_append)
write_csv(df_mut_meta_add, "../results/df_mut_meta_add.csv")
df_mut_meta_add <- read_csv("../results/df_mut_meta_add.csv", guess_max=100000)
df_mut_meta_add <- df_mut_meta_add %>% filter(Sample %in% df_meta$Sample)

df_mut_meta_add %>% filter(sec_freq>0.03) %>% group_by(Sample) %>% summarise(n=n()) %>% arrange(desc(n))
df_summary <- df_mut_meta_add %>% filter(sec_freq>0.03) %>% group_by(Sample, lineage_sim, coverage, Ct_value, `Report date`, Vaccine, Doses, `Condition (ever)`) %>% summarise(n=n(), mean_alle_freq = mean(alle_freq, na.rm=T), mean_con_freq=mean(con_freq, na.rm=T), mean_sec_freq=mean(sec_freq, na.rm=T)) %>% arrange(desc(n)) 
df_summary %>% write_csv("../results/num_snvs_meata.csv")

writeLines(df_summary$Sample, "../results/samples_hq.txt")

# filter the outlier samples
#### 0. Remove the samples with extreme number of SNPs (outliers by the 99% percentile cut-off); 
df_tmp <- df_mut_meta_add %>% filter(sec_freq>0.03) %>% group_by(Sample, lineage_sim, Vaccine) %>% summarise(n=n()) %>% arrange(desc(n))
df_tmp %>% arrange(desc(n))
quantile(df_tmp$n,1:100/100) ## remove the extreme values (outliers) by the 99% percentile cut-off.
check <- df_tmp$n>quantile(df_tmp$n,1:100/100)[99]
(samples_extreme <- df_tmp$Sample[check]) # the samples being removed due to extreme number of SNPs
write_csv(df_tmp %>% filter(Sample %in% samples_extreme), "../results/samples_extreme_high_num_of_SNPs.csv") 
df_mut_meta_add <- df_mut_meta_add %>% filter(!Sample %in% samples_extreme)

### filtering (QC)
#### 1. The identified mutations/variants should have at least 100 properly paired reads at the genomic positions (minimum depth of minor allele of 5 required); 
df_snvs_meta_add_qc <- df_mut_meta_add %>% filter(depth_pp>=100) %>% filter((sec_fwd+sec_rev)>=5 & sec_fwd>=1 & sec_rev>=1) %>% filter(nchar(type_alt)==1)
df_indels_meta_add_qc <- df_mut_meta_add %>% filter(depth_pp>=100) %>% filter((sec_fwd+sec_rev)>=5 & sec_fwd>=1 & sec_rev>=1) %>% filter(nchar(type_alt)!=1)
#### 2. positions 1:100 and (29903-99):29903 should be removed; 
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(X2>100 & X2<(29903-99))
#### 3. exclude all positions in the PCR primer binding regions
df_primer_new <- read_tsv("../../2021-10-11_nCoV_primers/results/primers_20211011.bed", col_names=F)
df_primer_old <- read_tsv("../../2021-10-11_nCoV_primers/results/primers_old.bed", col_names=F)

pos_all <- unique(paste(df_snvs_meta_add_qc$X2, df_snvs_meta_add_qc$primer))
check <- sapply(pos_all, function(x) {
	pos_x <- as.numeric(strsplit(x, " ")[[1]][1])
	primer_x <- strsplit(x, " ")[[1]][2]
	if(primer_x=="new"){
		any((df_primer_new$X2 <= pos_x) & (df_primer_new$X3 >= pos_x))	
	} else {
		any((df_primer_old$X2 <= pos_x) & (df_primer_old$X3 >= pos_x))	
	}	
})

df_snvs_meta_add_qc <- df_snvs_meta_add_qc[!paste(df_snvs_meta_add_qc$X2, df_snvs_meta_add_qc$primer) %in% pos_all[check],]

#### test different MAF cutoffs
#### The minor allele frequency (secondary most base frequency) cutoff for determining iSNVs
cut_offs <- 3:10/100
df_plot_maf <- lapply(cut_offs, function(x) {
	df_snvs_meta_add_qc %>% group_by(Sample, lineage_sim, Vaccine, Doses, days_since_last_dose) %>% filter(sec_freq>=x) %>% summarise(n=n(), Ct_value=Ct_value[1], alle_freq=as.character(x))
})
df_plot_maf <- bind_rows(df_plot_maf)

ggplot(df_plot_maf) + # iSNV threshold 
	geom_point(aes(x=as.numeric(Ct_value), y=n, color=Vaccine), alpha=0.6)+
	facet_wrap(vars(alle_freq), ncol=2)+
	xlab("Ct value")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine)
ggsave("../results/MAF-cutoff_Ct.pdf", height = 8, width = 6)

ggplot(df_plot_maf) + 
	geom_point(aes(x=abs(days_since_last_dose), y=n, color=Vaccine), alpha=0.6)+
	facet_wrap(vars(alle_freq), ncol=2)+
	xlab("Days since last dose")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine)
ggsave("../results/MAF-cutoff_doses.pdf", height = 8, width = 6)

thredshold_maf <- 0.05
df_plot_maf_filter <- df_plot_maf %>% filter(alle_freq==thredshold_maf) 
df_plot_maf_filter %>% arrange(desc(n))
cor.test(as.numeric(df_plot_maf_filter$Ct_value), df_plot_maf_filter$n) # correlation between Ct and number of iSNVs
cor.test(abs(df_plot_maf_filter$days_since_last_dose), df_plot_maf_filter$n) # correlation between days post vaccination and number of iSNVs

ggplot(df_plot_maf_filter) + 
	geom_point(aes(x=as.numeric(Ct_value), y=n, color=lineage_sim), alpha=0.6)+
	xlab("Ct value")+
	ylab("No. of identified iSNVs")
ggsave("../results/Ct_vs_lineage.pdf")

p_t_1 <- ggplot(df_plot_maf_filter) + 
	geom_point(aes(x=as.numeric(Ct_value), y=n, color=Vaccine), alpha=0.6)+
	xlab("Ct value")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine)

p_t_2 <- ggplot(df_plot_maf_filter) + 
	geom_point(aes(x=abs(days_since_last_dose), y=n, color=Vaccine), alpha=0.6)+
	xlab("Days since last dose")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine)
p_out <- p_t_1/p_t_2 + plot_layout(guides="collect")
ggsave("../results/supp_fig1.pdf", plot=p_out )

#### 2. The mutations should be supported by at least 5% of the total reads (minor allele frequency >=5%).
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(sec_freq>=thredshold_maf) 
df_snvs_meta_add_qc$sample <- df_snvs_meta_add_qc$Sample
(df_snvs_meta_add_qc %>% select(Sample:lineage_sim,Vaccine) %>% unique() %>% mutate(type=paste(Vaccine, lineage_sim)) %>% .$type %>% table() %>% sort()) # number of cases in each group
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(!X2 %in% c(15494, 15489, 25381, 10194, 22422)) # excluding primer/homoplasy sites, https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28420-7/MediaObjects/41467_2022_28420_MOESM1_ESM.pdf

# remove serial adjacent separate mutations
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
df_snvs_meta_add_qc %>% arrange(X2) %>% group_by(sample) %>% mutate(check_serial_mut=fn_check_serial_mut(X2)) %>% ungroup() %>% filter(check_serial_mut)
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% arrange(X2) %>% group_by(sample) %>% mutate(check_serial_mut=fn_check_serial_mut(X2)) %>% ungroup() %>% filter(!check_serial_mut)

write_csv(df_snvs_meta_add_qc, "../results/df_snvs_meta_add_qc_ivar.csv")

# tmp_list <- tmp$sample[!paste0(tmp$sample, tmp$X2) %in% paste0(df_snvs_meta_add_qc$sample, df_snvs_meta_add_qc$X2)] %>% unique()
# writeLines(tmp_list, "../results/tmp.redolist.txt")
