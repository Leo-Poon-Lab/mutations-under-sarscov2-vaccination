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
df_meta$sample <- df_meta$Sample
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

df_meta <- df_meta %>% filter(Doses==0 | Doses==2) # remove partial vaccination

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

# ## 1. read the snvs identified by ivar.
# files_snvs <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/variant_caller/ivar/")
# files_snvs_full <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/variant_caller/ivar/", full.names = T)
# files_alrd_done <- gsub("\\.tsv\\.vcf$", ".tsv", list.files("../results/vcf/", "\\.tsv\\.vcf$"))
# check1 <- files_snvs %in% paste0("ivar_", df_meta$Sample, ".tsv")
# check2 <- !files_snvs %in% files_alrd_done
# files_snvs_full_int <- files_snvs_full[check1 & check2]
# source("./helper/convert_to_vcf.r")

# convert_to_vcf(files_snvs_full_int)
# files_vcfs <- list.files("../results/vcf/", "vcf$", full.names = T)
# files_csv <- list.files("../results/vcf/", "csv$", full.names = T)
# files_vcfs_todo <- files_vcfs[!gsub("\\.tsv\\.vcf$", "", files_vcfs) %in% gsub("\\.tsv\\.vcf\\.csv$", "", files_csv)]
# annotate_snpeff(files_vcfs_todo)
# files_vcfs_csv <- list.files("../results/vcf/", "csv$", full.names = T)
# idx <- sapply(paste0(df_meta$Sample, ".tsv"), function(x) {
# 	tmp <- grep(x, files_vcfs_csv)[1]
# })
# files_vcfs_csv <- files_vcfs_csv[idx]
# df_snvs <- mclapply(files_vcfs_csv, read_snpeff, mc.cores=16)
# df_snvs <- bind_rows(df_snvs)
 
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
df_bam_rst <- mclapply(files_bam_rst_full, read_pysamstats, pp=FALSE, mc.cores=16)
df_bam_rst <- bind_rows(df_bam_rst[!is.na(df_bam_rst)])
# save(df_bam_rst, file="../results/df_bam_rst_full.rdata")
save(df_bam_rst, file="../results/df_bam_rst_full_notpp.rdata")
# load("../results/df_bam_rst_full_notpp.rdata")

df_meta <- df_meta %>% filter(sample %in% unique(df_bam_rst$sample))
write_csv(df_meta, "../results/df_samples_clean.csv")

muAF <- df_bam_rst$matches/df_bam_rst$reads_all
df_bam_rst_filter <- df_bam_rst[muAF<0.95 & muAF>0.05,]

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

df_bam_rst_filter_append <- left_join(df_bam_rst_filter_append, df_meta %>% mutate(sample=Sample))
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
length(unique(df_bam_rst_filter_append$sample))

#### remove possible contamination
df_tmp <- df_bam_rst_filter_append %>% group_by(Sample, lineage_sim, Vaccine) %>% summarise(n=n()) %>% arrange(desc(n))
check <- df_tmp$n>quantile(df_tmp$n,1:100/100)[99]
(samples_extreme <- df_tmp$Sample[check]) # the samples being removed due to extreme number of SNPs
write_csv(df_tmp %>% filter(Sample %in% samples_extreme), "../results/samples_extreme_high_num_of_SNPs.csv") 
df_bam_rst_filter_append <- df_bam_rst_filter_append %>% filter(!Sample %in% samples_extreme)
length(unique(df_bam_rst_filter_append$sample))

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
write_csv(df_snvs_meta_add_qc, "../results/df_snvs_meta_add_qc_bam.csv")


