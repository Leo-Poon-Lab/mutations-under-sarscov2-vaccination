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

df_meta <- read_csv("../results/df_samples.csv", guess_max=100000)
df_meta <- df_meta %>% filter(lineage_sim != "22B (Omicron, BA.5.*)")

# quantile(df_meta$Ct_value, (1:100)/100, na.rm=T)
# df_meta %>% filter(is.na(Ct_value)) %>% group_by(lineage_sim) %>% summarise(n())
stopifnot(max(df_meta$Ct_value, na.rm=T)<=28)

# table(df_meta$lineage_sim, df_meta$Vaccine, df_meta$Doses)
# table(df_meta$Vaccine)
# table(df_meta$lineage_sim, df_meta$Vaccine)

### using pysamstats
# system("conda run -n base pysamstats --help")
source("./helper/pysamstats.r")

files_bam_rst_full <- c(list.files("../results/pysamstats/", "tsv$", full.names = T), list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", "tsv$", full.names = T), list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/", "tsv$", full.names = T))

idx <- sapply(paste0("/", df_meta$sample, ".tsv"), function(x) {
	tmp <- grep(x, files_bam_rst_full)[1]
})
files_bam_rst_full <- files_bam_rst_full[idx]
stopifnot(sum(is.na(files_bam_rst_full))==0)

df_bam_rst <- mclapply(files_bam_rst_full, read_pysamstats, pp=FALSE, mc.cores=16)
df_bam_rst <- bind_rows(df_bam_rst[!is.na(df_bam_rst)])
# save(df_bam_rst, file="../results/df_bam_rst_full.rdata")

df_bam_rst$muAF <- df_bam_rst$matches/df_bam_rst$reads_all
save(df_bam_rst, file="../results/df_bam_rst_full_notpp.rdata")
# load("../results/df_bam_rst_full_notpp.rdata")

MAF_threshold <- 0.025
df_bam_rst_filter <- df_bam_rst %>% filter((muAF<(1-MAF_threshold)) & (muAF>MAF_threshold))

df_extra_info <- extra_info_from_pysamstats(df_bam_rst_filter)
df_extra_info <- df_extra_info %>% mutate_at(vars(!contains(c("base", "sample"))), as.numeric)
df_extra_info <- left_join(df_extra_info, df_meta %>% select(sample, primer, lineage_sim, Vaccine), "sample")
df_extra_info$primer[is.na(df_extra_info$primer)] <- "new" # RE-WHP2800

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

#### remove possible contamination
df_tmp <- df_extra_info %>% group_by(sample, lineage_sim, Vaccine) %>% summarise(n=n()) %>% arrange(desc(n))
quantile(df_tmp$n,1:100/100)

check <- df_tmp$n>quantile(df_tmp$n,1:100/100)[99]
(samples_extreme <- df_tmp$sample[check]) # the samples being removed due to extreme (>99 percentile) number of SNPs
# df_meta %>% select(sample, lineage_sim, Vaccine) %>% filter(sample %in% samples_extreme) %>% group_by(lineage_sim, Vaccine) %>% summarise(n())
# df_meta %>% select(sample, lineage_sim, Vaccine, Ct_value) %>% filter(sample %in% samples_extreme) %>% print(n=50)
write_csv(df_tmp %>% filter(sample %in% samples_extreme), "../results/samples_extreme_high_num_of_SNPs.csv") 
df_extra_info <- df_extra_info %>% filter(!sample %in% samples_extreme)

df_tmp <- df_extra_info %>% select(pos, con_base, sec_base)
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

df_vcf <- unique(df_vcf)

outfile_vcf <- "../results/vcf_bam_isnvs.vcf"
outfile_snpeff <- "../results/vcf_bam_isnvs.vcf.snpeff"
outfile_csv <- "../results/vcf_bam_isnvs.vcf.snpeff.csv"
write_tsv(df_vcf, outfile_vcf)
system(paste0("java -jar ~/softwares/snpEff/snpEff.jar ", "NC_045512.2", " ", outfile_vcf, " > ", outfile_snpeff))
system(paste0("Rscript ./helper/Parse_SnpEff.r ", outfile_snpeff,  " ", outfile_csv))

source("./helper/convert_to_vcf.r")
df_annt <- read_snpeff(outfile_csv)

df_snvs_meta_add_qc <- left_join(df_extra_info, df_annt %>% mutate(pos=X2, con_base=X4, sec_base=X5) %>% select(-X1, -X2, -X4, -X5, -sample), c( "pos", "con_base", "sec_base"))
write_csv(df_snvs_meta_add_qc, "../results/df_snvs_meta_add_qc_bam.csv")


