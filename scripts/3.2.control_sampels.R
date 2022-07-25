library(tidyverse)
library(parallel)
library(ggsci)

data_meta_raw <- read_csv("../../../2021/2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 20000)
data_meta_raw <- data_meta_raw %>% filter(sequenced_by_us) %>% filter(!is.na(`Report date`))

data_meta_study <- read_csv("../results/df_samples_clean.csv")

data_meta_raw <- data_meta_raw %>% filter(Sample %in% data_meta_study$sample)
data_meta_raw <- data_meta_raw %>% filter(grepl(",", samples_all, fixed=T))
data_meta_raw$num_whp <- sapply(data_meta_raw$samples_all, function(x){
	tmp <- strsplit(x, ", ")[[1]]
	sum(grepl("WHP", tmp))
})
data_meta_raw <- data_meta_raw %>% filter(num_whp>1)
samples_controls <- data_meta_raw$samples_all
samples_controls <- unlist(strsplit(samples_controls, ", "))
samples_controls <- samples_controls[grepl("WHP", samples_controls)]
samples_controls <- gsub("_\\d+$", "", samples_controls)

samples_controls_whp_id <- gsub("RE-", "", samples_controls)
samples_controls_whp_id <- gsub("-.+", "", samples_controls_whp_id)
duplicated_samples_whpids <- names(table(samples_controls_whp_id)[table(samples_controls_whp_id)>1])
samples_controls <- samples_controls[samples_controls_whp_id %in% duplicated_samples_whpids]

source("./helper/pysamstats.r")
# read bamreadcounts of the control samples

system("conda run -n base pysamstats --help")
source("./helper/pysamstats.r")
files_bam <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
files_bam_full <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$", full.names = T)
samples_alrd_done <- gsub("\\.tsv$", "", list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/"))
samples_alrd_done_1 <- gsub("\\.tsv$", "", list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/"))
samples_alrd_done_2 <- gsub("\\.tsv$", "", list.files("../results/pysamstats/"))
check1 <- files_bam %in% paste0(samples_controls, "-trimmed.masked.bam")
check2 <- !files_bam %in% paste0(c(samples_alrd_done, samples_alrd_done_1, samples_alrd_done_2), "-trimmed.masked.bam")
files_bam_int <- files_bam_full[check1 & check2]
mclapply(files_bam_int, function(x) {
	process_pysamstats(x)
	return("")
}, mc.cores=16)


files_bam_rst_full <- c(list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/", "tsv$", full.names = T), list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", "tsv$", full.names = T), list.files("../results/pysamstats/", "tsv$", full.names = T))
idx <- sapply(paste0(samples_controls, ".tsv"), function(x) {
	tmp <- grep(x, files_bam_rst_full)[1]
})
files_bam_rst_full <- files_bam_rst_full[idx]
files_bam_rst_full <- files_bam_rst_full[!is.na(files_bam_rst_full)]

### only nucleotide positions with at least 100 properly paired reads were considered
### only samples with 27000 sites of sequencing depth >= 100 were considered
df_bam_rst <- mclapply(files_bam_rst_full, read_pysamstats, pp=FALSE, mc.cores=16)
df_bam_rst <- bind_rows(df_bam_rst[!is.na(df_bam_rst)])

samples_pass_filter <- unique(df_bam_rst$sample)
samples_pass_filter_case_id <- sapply(samples_pass_filter, function(x) {
	tmp <- which(grepl(x, data_meta_raw$samples_all))
	data_meta_raw$case_id[tmp]
})
cases_have_duplicate <- table(samples_pass_filter_case_id)[table(samples_pass_filter_case_id)>1]
length(cases_have_duplicate) # 86 cases have high quality replicates
sort(cases_have_duplicate)
samples_high_quality <- unique(df_bam_rst$sample)[samples_pass_filter_case_id %in% names(cases_have_duplicate)] # 179 samples

df_bam_rst <- df_bam_rst %>% filter(sample %in% samples_high_quality)

samples_pass_filter <- unique(df_bam_rst$sample)
df_tmp <- tibble(sample=samples_pass_filter)
df_tmp$whp_id <- gsub("RE-", "", df_tmp$sample)
df_tmp$whp_id <- gsub("-.+", "", df_tmp$whp_id)

# filter isnvs
muAF <- df_bam_rst$matches/df_bam_rst$reads_all
df_bam_rst_filter <- df_bam_rst[which(muAF<0.95 & muAF>0.05),]

### add readcounts from pysamstats results
bases <- c("deletions", "insertions", "A", "C", "T", "G")

df_bam_rst_filter_append <- mclapply(seq_len(nrow(df_bam_rst_filter)), function(i){
	print(i)
	df_bam_rst_i <- df_bam_rst_filter[i,]
	
	depth_i <- df_bam_rst_i$reads_all # count on paired reads·
	num_bases <- c(df_bam_rst_i$deletions, df_bam_rst_i$insertions, df_bam_rst_i$A, df_bam_rst_i$C, df_bam_rst_i$T, df_bam_rst_i$G)
	ord_t <- order(num_bases, decreasing = T)
	consensus_base <- bases[ord_t[1]]
	secondary_base <- bases[ord_t[2]]

	con_fwd_i <- df_bam_rst_i[[paste0(consensus_base, "_fwd")]]
	con_rev_i <- df_bam_rst_i[[paste0(consensus_base, "_rev")]]
	sec_fwd_i <- df_bam_rst_i[[paste0(secondary_base, "_fwd")]]
	sec_rev_i <- df_bam_rst_i[[paste0(secondary_base, "_rev")]]

	strand_bias_i <- sec_fwd_i/sec_rev_i
	if(is.na(strand_bias_i)){strand_bias_i <- Inf} else	if(strand_bias_i<1){strand_bias_i <- 1/strand_bias_i}

	con_freq_i <- (con_fwd_i+con_rev_i)/depth_i
	sec_freq_i <- (sec_fwd_i+sec_rev_i)/depth_i

	c(depth=depth_i, depth_fwd=df_bam_rst_i$reads_fwd, depth_rev=df_bam_rst_i$reads_rev, con_base=consensus_base, con_fwd=con_fwd_i, con_rev=con_rev_i, con_freq=con_freq_i, sec_base=secondary_base, sec_fwd=sec_fwd_i, sec_rev=sec_rev_i, sec_freq=sec_freq_i, strand_bias=strand_bias_i)
}, mc.cores = 16)

df_bam_rst_filter_append <- bind_rows(df_bam_rst_filter_append)
df_bam_rst_filter_append <- bind_cols(df_bam_rst_filter %>% select(sample, pos), df_bam_rst_filter_append)

df_bam_rst_filter_append$case_id <- sapply(df_bam_rst_filter_append$sample, function(x) {
	tmp <- which(grepl(x, data_meta_raw$samples_all))
	data_meta_raw$case_id[tmp]
})

### filtering (QC)
#### MAF threshold
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(sec_freq>=0.05) # MAF>=0.05
#### filter structural variants
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(nchar(sec_base)==1) # only single-nt mutations
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(nchar(con_base)==1) # only single-nt mutations
length(unique(df_bam_rst_filter_append$sample))

#### filter strand bias
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(strand_bias<10) 
length(unique(df_bam_rst_filter_append$sample))

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
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% arrange(pos) %>% group_by(sample) %>% mutate(check_serial_mut=fn_check_serial_mut(pos)) %>% ungroup() %>% filter(!check_serial_mut) # remove serial adjacent separate mutations
length(unique(df_bam_rst_filter_append$sample))

#### 1. The identified mutations/variants should have at least 100 properly paired reads at the genomic positions (minimum depth of minor allele of 5 required); 
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(as.numeric(depth)>=100) %>% filter((as.numeric(sec_fwd)+as.numeric(sec_rev))>=5 & as.numeric(sec_fwd)>=1 & as.numeric(sec_rev)>=1)
length(unique(df_bam_rst_filter_append$sample))

#### 2. positions 1:100 and (29903-99):29903 should be removed; 
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(pos>100 & pos<(29903-99))

#### 3. exclude all positions in the PCR primer binding regions
df_primer_new <- read_tsv("../../2021-10-11_nCoV_primers/results/primers_20211011.bed", col_names=F)
df_primer_old <- read_tsv("../../2021-10-11_nCoV_primers/results/primers_old.bed", col_names=F)

df_bam_rst_filter_append <- left_join(df_bam_rst_filter_append, data_meta_study %>% select(-sample), "case_id")

pos_all <- unique(paste(df_bam_rst_filter_append$pos, df_bam_rst_filter_append$primer))
check <- sapply(pos_all, function(x) {
	pos_x <- as.numeric(strsplit(x, " ")[[1]][1])
	primer_x <- strsplit(x, " ")[[1]][2]
	if(primer_x=="new"){
		any((df_primer_new$X2 <= pos_x) & (df_primer_new$X3 >= pos_x))	
	} else {
		any((df_primer_old$X2 <= pos_x) & (df_primer_old$X3 >= pos_x))	
	}	
})

df_bam_rst_filter_append <- df_bam_rst_filter_append[!paste(df_bam_rst_filter_append$pos, df_bam_rst_filter_append$primer) %in% pos_all[check],]
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(!pos %in% c(15494, 15489, 25381, 10194, 22422)) # excluding primer/homoplasy sites, https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28420-7/MediaObjects/41467_2022_28420_MOESM1_ESM.pdf

#### remove possible contamination
# df_tmp <- df_bam_rst_filter_append %>% group_by(sample, lineage_sim, Vaccine) %>% summarise(n=n()) %>% arrange(desc(n))
# check <- df_tmp$n>quantile(df_tmp$n,1:100/100)[99]
# (samples_extreme <- df_tmp$sample[check]) # the samples being removed due to extreme number of SNPs
# write_csv(df_tmp %>% filter(Sample %in% samples_extreme), "../results/samples_extreme_high_num_of_SNPs.csv") 
# df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(!Sample %in% samples_extreme)
# length(unique(df_bam_rst_filter_append$sample))

df_tmp <- df_bam_rst_filter_append %>% select(pos, con_base, sec_base)
## annotate gene and mutations
# #CHROM POS     ID        REF    ALT     QUAL FILTER INFO                    
# 20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2
# 20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017   
df_vcf <- tibble(`#CHROM` = rep("NC_045512.2", nrow(df_tmp)))
df_vcf$POS = df_tmp$pos
df_vcf$ID = "."
df_vcf$REF = df_tmp$con_base
df_vcf$ALT = df_tmp$sec_base
df_vcf$QUAL = 100
df_vcf$FILTER = "PASS"
df_vcf$INFO = "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ	0:0:0:0:0:0:0:0"

outfile_vcf <- "../results/vcf_bam_isnvs.vcf"
outfile_snpeff <- "../results/vcf_bam_isnvs.vcf.snpeff"
outfile_csv <- "../results/vcf_bam_isnvs.vcf.snpeff.csv"
write_tsv(df_vcf, outfile_vcf)
system(paste0("java -jar ~/softwares/snpEff/snpEff.jar ", "NC_045512.2", " ", outfile_vcf, " > ", outfile_snpeff))
system(paste0("Rscript ./helper/Parse_SnpEff.r ", outfile_snpeff,  " ", outfile_csv))

source("./helper/convert_to_vcf.r")
df_annt <- read_snpeff(outfile_csv)

df_snvs_meta_add_qc <- bind_cols(df_annt %>% select(-sample), df_bam_rst_filter_append)

## number of recurrent iSNVs 
df_snvs_meta_add_qc$effect_sim <- df_snvs_meta_add_qc$effect
df_snvs_meta_add_qc$effect_sim[grepl("stream", df_snvs_meta_add_qc$effect_sim)] <- "UTR"
df_snvs_meta_add_qc$effect_sim[grepl("stop", df_snvs_meta_add_qc$effect_sim)] <- "Nonsynonymous"
df_snvs_meta_add_qc$effect_sim[grepl("missense", df_snvs_meta_add_qc$effect_sim)] <- "Nonsynonymous"
df_snvs_meta_add_qc$effect_sim[grepl("lost", df_snvs_meta_add_qc$effect_sim)] <- "Nonsynonymous"
df_snvs_meta_add_qc$effect_sim[grepl("^synonymous", df_snvs_meta_add_qc$effect_sim)] <- "Synonymous"
df_snvs_meta_add_qc$effect_sim[grepl("&synonymous", df_snvs_meta_add_qc$effect_sim)] <- "Synonymous"
df_snvs_meta_add_qc$effect_sim <- gsub("_", " ", df_snvs_meta_add_qc$effect_sim)

table(df_snvs_meta_add_qc$effect_sim)
df_tmp <- df_snvs_meta_add_qc %>% group_by(X2, effect_sim) %>% summarise(n=n())
p_n_sharing <- ggplot(df_tmp) +
	geom_histogram(aes(x=n, fill=effect_sim), binwidth = 1, color="black", size=0.3)+
	xlab("Number of samples sharing iSNVs")+
	ylab("Number of iSNVs")+
   scale_fill_manual(name="Variant type", values=pal_jama()(3)[c(3,2,1)])+
	scale_x_continuous(breaks=seq(0,60,5))+
	theme_minimal()+
	theme(legend.position = "top")+
	NULL
ggsave("../results/recurrent_isnvs_control.pdf")

sum(df_tmp$n>1) # 120
sum(df_tmp$n==1) # 125

sum(df_tmp$n==1)/nrow(df_tmp)
sum(df_tmp$n>1)/nrow(df_tmp)


120/243
125/1845
chisq.test(matrix(c(120, 125, 243, 1845), nrow=2)) # p<0.01

unique(df_snvs_meta_add_qc$sample)
df_snvs_meta_add_qc %>% group_by(case_id) %>% summarise(N=n()) %>% filter(N>1)

cases_have_duplicate

df_shared_isnvs <- sapply(unique(df_snvs_meta_add_qc$case_id), function(case_id_i) {
	print(case_id_i)

	df_same_patient_i <- df_snvs_meta_add_qc %>% filter(case_id==case_id_i)
	df_diff_patient_i <- df_snvs_meta_add_qc %>% filter(case_id!=case_id_i)
	samples_same_i <- unique(df_same_patient_i$sample)
	if(length(samples_same_i)<2){return(c(NA,NA))}


	rst_same_i <- apply(combn(samples_same_i,2),2,function(pairs) {
		sites_1 <- df_same_patient_i$X2[df_same_patient_i$sample==pairs[1]]
		sites_2 <- df_same_patient_i$X2[df_same_patient_i$sample==pairs[2]]
		sum(sites_1 %in% sites_2)/length(unique(df_same_patient_i$X2))
	})

	samples_diff_i <- unique(df_diff_patient_i$sample)
	rst_diff_i <- sapply(samples_same_i, function(sample_i) {
		sites_1 <- df_same_patient_i$X2[df_same_patient_i$sample==sample_i]
		rst_tmp <- sapply(samples_diff_i, function(sample_j) {
			sites_2 <- df_diff_patient_i$X2[df_diff_patient_i$sample==sample_j]
			sum(sites_1 %in% sites_2)/length(unique(c(df_diff_patient_i$X2, df_same_patient_i$X2)))
		})
		mean(rst_tmp)		
	})

	return(c(mean(rst_same_i), mean(rst_diff_i)))
})

df_shared_isnvs <- t(df_shared_isnvs)
mean(df_shared_isnvs[,1], na.rm=T)
mean(df_shared_isnvs[,2], na.rm=T)
median(df_shared_isnvs[,1], na.rm=T)
median(df_shared_isnvs[,2], na.rm=T)
wilcox.test(df_shared_isnvs[,1], df_shared_isnvs[,2])
