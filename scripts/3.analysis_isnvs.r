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

colors_lineage=c("#ff7f00", "#33a02c", "#1f78b4", "#6a3d9a", "#b15928", "#e31a1c") 
names(colors_lineage) <- c("Alpha", "B.1", "Delta", "B.1.1.63", "B.1.36.27", "B.1.36")
colors_vaccine_type=pal_jco()(4)
names(colors_vaccine_type)=c("Adeno", "Inactivated", "Non-vaccinated", "mRNA")

df_meta <- read_csv("../results/df_samples_clean.csv")
# df_meta <- read_csv("../data/df_samples_vaccine_bt_study.csv")
df_meta$`Ct value`[df_meta$`Ct value`==0] <- NA
df_meta <- df_meta %>% filter(!is.na(`Ct value`))

df_vaccine_type <- readxl::read_excel("../data/Vaccination case for breakthrough infection .xlsx", skip=1)
df_vaccine_type <- df_vaccine_type %>% select(Vaccine, `Vaccine type`) %>% unique()
vaccine_of_interest <- c("BioNTech", "Sinovac", "Non-vaccinated")

df_meta$lineage_sim <- df_meta$lineage
df_meta$lineage_sim[grepl("^AY\\.", df_meta$lineage_sim)] <- "Delta"
df_meta$lineage_sim[df_meta$lineage_sim == "B.1.617.2"] <- "Delta"
df_meta$lineage_sim[grepl("^B\\.1$", df_meta$lineage_sim)] <- "B.1"
df_meta$lineage_sim[grepl("^B\\.1\\.1\\.7$", df_meta$lineage_sim)] <- "Alpha"
df_meta$lineage_sim[grepl("^Q\\.", df_meta$lineage_sim)] <- "Alpha"
df_meta <- df_meta %>% filter(lineage_sim %in% names(colors_lineage))
### only include samples with Ct value < 25
df_meta <- df_meta %>% filter(`Ct value`<25)
df_meta <- df_meta %>% filter(Sample!="WHP4499")
table(df_meta$Vaccine, df_meta$lineage_sim)

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
files_snpeff <- list.files("../results/vcf/", "snpeff$", full.names = T)
files_vcfs_todo <- files_vcfs[!gsub("\\.tsv\\.vcf$", "", files_vcfs) %in% gsub("\\.tsv\\.vcf\\.snpeff$", "", files_snpeff)]
annotate_snpeff(files_vcfs_todo)
files_vcfs_csv <- list.files("../results/vcf/", "csv$", full.names = T)
idx <- sapply(paste0(df_meta$Sample, ".tsv"), function(x) {
	tmp <- grep(x, files_vcfs_csv)[1]
})
files_vcfs_csv <- files_vcfs_csv[idx]
df_snvs <- read_snpeff(files_vcfs_csv)

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
df_bam_rst <- mclapply(files_bam_rst_full, read_pysamstats, mc.cores=16)
df_bam_rst <- bind_rows(df_bam_rst[!is.na(df_bam_rst)])
df_bam_rst <- df_bam_rst %>% filter(reads_pp>=100)

### only nucleotide positions with at least 100 properly paired reads were considered
### only samples with 27000 sites of sequencing depth >= 100 were considered
df_bam_rst_depth <- df_bam_rst %>% group_by(sample) %>% summarise(n_high_depth=n()) 
df_bam_rst_depth %>% arrange(desc(n_high_depth))
plot(df_bam_rst_depth$n_high_depth)
samples_hq <- df_bam_rst_depth %>% filter(n_high_depth>=27000) %>% .$sample

writeLines(samples_hq, "../results/samples_hq.txt")
samples_hq <- readLines("../results/samples_hq.txt")

names(df_snvs)[names(df_snvs)=="sample"] <- "Sample"
df_snvs_meta <- left_join(df_snvs, df_meta, "Sample")
df_snvs_meta <- df_snvs_meta %>% filter(Sample %in% samples_hq)
df_snvs_meta$Vaccine[df_snvs_meta$Vaccine=="Biontech"] <- "BioNTech"
df_snvs_meta$Vaccine[df_snvs_meta$Vaccine=="Non-Vaccinated"] <- "Non-vaccinated"
df_snvs_meta <- left_join(df_snvs_meta, df_vaccine_type, "Vaccine")
df_snvs_meta$`Vaccine type`[df_snvs_meta$Vaccine=="Non-vaccinated"] <- "Non-vaccinated"
write_csv(df_snvs_meta, "../results/df_mut_meta.csv")
# df_snvs_meta <- read_csv("../results/df_mut_meta.csv", guess_max=600000)

### add readcounts from pysamstats results
df_bam_rst <- df_bam_rst %>% filter(sample %in% unique(df_snvs_meta$Sample))
bases <- c("deletions", "insertions", "A", "C", "T", "G")

?system.time()
df_snvs_meta_append <- mclapply(seq_len(nrow(df_snvs_meta)), function(i){
	print(i)
	df_snvs_meta_i <- df_snvs_meta[i,]
	sample_i <- df_snvs_meta_i$Sample
	pos_i <- df_snvs_meta_i$X2
	ref_i <- df_snvs_meta_i$X4
	alt_i <- df_snvs_meta_i$X5
	df_bam_rst_i <- df_bam_rst[df_bam_rst$sample==sample_i & df_bam_rst$pos==pos_i,]
	
	depth_i <- df_bam_rst_i$reads_all
	depth_pp_i <- df_bam_rst_i$reads_pp # count on paired reads

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
df_mut_meta_add <- read_csv("../results/df_mut_meta_add.csv", guess_max=60000)
table(df_mut_meta_add$Vaccine)
df_mut_meta_add <- left_join(df_mut_meta_add %>% select(-`Vaccine type`), df_vaccine_type, "Vaccine")
df_mut_meta_add$`Vaccine type`[df_mut_meta_add$Vaccine=="Non-vaccinated"] <- "Non-vaccinated"
table(df_mut_meta_add$`Vaccine type`)
sort(table(df_mut_meta_add$Sample))
quantile(table(df_mut_meta_add$Sample),1:100/100) ## remove the extreme values (outliers) by the 99% percentile cut-off.
check <- table(df_mut_meta_add$Sample)>quantile(table(df_mut_meta_add$Sample),1:100/100)[99]
df_mut_meta_add <- df_mut_meta_add %>% filter(!Sample %in% names(check)[check])

# filter the outlier samples

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
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(Vaccine %in% vaccine_of_interest)
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(!(lineage_sim!="Delta" & Vaccine!="Non-vaccinated"))

df_snvs_meta_add_qc %>% select(Sample:lineage_sim) %>% unique()
df_snvs_meta_add_qc %>% select(Sample:lineage_sim) %>% unique() %>% mutate(type=paste(Vaccine, lineage_sim)) %>% .$type %>% table() %>% sort()

#### test different MAF cutoffs
#### The minor allele frequency (secondary most base frequency) cutoff for determining iSNVs
cut_offs <- 3:10/100
df_plot_maf <- lapply(cut_offs, function(x) {
	df_snvs_meta_add_qc %>% group_by(Sample, Vaccine, Doses, days_since_last_dose, `Vaccine type`) %>% filter(sec_freq>=x) %>% summarise(n=n(), Ct_value=`Ct value`[1], alle_freq=as.character(x))
})
df_plot_maf <- bind_rows(df_plot_maf)

ggplot(df_plot_maf) + # iSNV threshold
	geom_point(aes(x=as.numeric(Ct_value), y=n, color=`Vaccine type`), alpha=0.6)+
	facet_wrap(vars(alle_freq), ncol=2)+
	xlab("Ct value")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine_type)
ggsave("../results/MAF-cutoff_Ct.pdf", height = 8, width = 6)

ggplot(df_plot_maf) + 
	geom_point(aes(x=abs(days_since_last_dose), y=n, color=`Vaccine type`), alpha=0.6)+
	facet_wrap(vars(alle_freq), ncol=2)+
	xlab("Days since last dose")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine_type)
ggsave("../results/MAF-cutoff_doses.pdf", height = 8, width = 6)

df_plot_maf_003 <- df_plot_maf %>% filter(alle_freq==0.03)
df_plot_maf_003 %>% arrange(desc(n))

p_t_1 <- ggplot(df_plot_maf_003) + 
	geom_point(aes(x=as.numeric(Ct_value), y=n, color=`Vaccine type`), alpha=0.6)+
	xlab("Ct value")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine_type)
cor.test(as.numeric(df_plot_maf_003$Ct_value), df_plot_maf_003$n) # correlation between Ct and number of iSNVs

p_t_2 <- ggplot(df_plot_maf_003) + 
	geom_point(aes(x=abs(days_since_last_dose), y=n, color=`Vaccine type`), alpha=0.6)+
	xlab("Days since last dose")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine_type)
p_out <- p_t_1/p_t_2
ggsave("../results/supp_fig1.pdf", plot=p_out)
cor.test(abs(df_plot_maf_003$days_since_last_dose), df_plot_maf_003$n) # correlation between days post vaccination and number of iSNVs

#### 2. The mutations should be supported by at least 3% of the total reads (minor allele frequency >=3%).
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(sec_freq>=0.03) 
df_snvs_meta_add_qc$sample <- df_snvs_meta_add_qc$Sample
write_csv(df_snvs_meta_add_qc, "../results/df_snvs_meta_add_qc_ivar.csv")
df_snvs_meta_add_qc <- read_csv("../results/df_snvs_meta_add_qc_ivar.csv", guess_max=600000)

# df_snvs_meta_add_qc_bam <- read_csv("../results/df_snvs_meta_add_qc_bam.csv")
# df_tmp <- df_snvs_meta_add_qc_bam[df_snvs_meta_add_qc_bam$con_base!=df_snvs_meta_add_qc_bam$ref,]
# write_csv(df_tmp, "../results/tmp.csv")

# (check <- paste(df_snvs_meta_add_qc$Sample, df_snvs_meta_add_qc$X2) %in% paste(df_snvs_meta_add_qc_bam$Sample, df_snvs_meta_add_qc_bam$pos))
# as.character(df_snvs_meta_add_qc[!check,])

# (check <- paste(df_snvs_meta_add_qc_bam$Sample, df_snvs_meta_add_qc_bam$pos) %in% paste(df_snvs_meta_add_qc$Sample, df_snvs_meta_add_qc$X2))
# df_tmp <- df_snvs_meta_add_qc_bam[!check,]
# write_csv(df_tmp, "../results/tmp.csv")

## 4. analyse the number of mutations, the frequency of mutations, diversity.
### Number of iSNV mutations 
# df_plot <- df_snvs_meta_add_qc_bam %>% filter(sec_freq>0.03)
df_plot <- df_snvs_meta_add_qc
df_plot$sample <- df_plot$Sample
df_tmp <- df_plot %>% group_by(sample) %>% summarise(n=n())
df_tmp$sample_label <- paste0(df_tmp$sample, " (N=", df_tmp$n, ")")
df_plot <- left_join(df_plot, df_tmp, "sample")
df_plot <- left_join(df_plot, df_bam_rst_depth, "sample")
df_plot$n_per_pos <- df_plot$n/df_plot$n_high_depth
df_plot$n_per_kb <- df_plot$n_per_pos*1000
df_plot_n <- df_plot %>% select(sample, n_per_kb, Vaccine, `Vaccine type`, lineage_sim) %>% unique()
ggplot(df_plot_n, aes(x=`Vaccine type`, y=n_per_kb, color=lineage_sim))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	ylab("Numer of iSNVs per Kb")+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	theme_classic()+
	NULL
ggsave("../results/Num_isnvs_by_vaccine_type.pdf", width=6, height=4)
save_pptx("../results/Num_isnvs_by_vaccine_type.pptx", width=6, height=4)

ggplot(df_plot_n, aes(x=`Vaccine`, y=n_per_kb, color=lineage_sim))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1), size=0.8)+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	ylab("Numer of iSNVs per Kb")+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	theme_classic()+
	# theme_minimal()+
	NULL
ggsave("../results/Num_isnvs_by_vaccine.pdf", width=6, height=4)
save_pptx("../results/Num_isnvs_by_vaccine.pptx", width=6, height=4)

df_orf_sim <- read_csv("../data/ORF_SCoV2_sim.csv")
df_plot_n_gene <- df_plot %>% group_by(sample, Vaccine, lineage_sim, gene) %>% summarise(n=n())
df_bam_rst$sample_fct <- factor(df_bam_rst$sample_fct)

system.time(df_bam_rst$sample_fct == df_plot_n_gene[1,1])
system.time(df_bam_rst$sample == df_plot_n_gene[1,1])

length_gene <- apply(df_plot_n_gene, 1, function(x) {
	print(x)
	df_bam_rst
	df_bam_rst %>% filter(sample==x[1] & pos>=df_orf_sim$start[df_orf_sim$sequence==x[4]] & pos<=df_orf_sim$stop[df_orf_sim$sequence==x[4]]) %>% nrow()
})
df_plot_n_gene$length_gene <- length_gene
df_plot_n_gene$n_per_kb <- df_plot_n_gene$n/df_plot_n_gene$length_gene*1000

df_plot_n$gene <- "Full genome"
df_plot_n_gene <- bind_rows(df_plot_n, df_plot_n_gene)

ggplot(df_plot_n_gene %>% filter(gene %in% c("Full genome", "ORF1ab", "S")), aes(x=`Vaccine`, y=n_per_kb, color=lineage_sim))+ # mutation rate per gene
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_grid(rows=vars(gene))+
	ylab("Numer of iSNVs per Kb")+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	theme_classic()+
	# theme_minimal()+
	NULL	
	
ggsave("../results/Num_isnvs_by_vaccine_gene.pdf", width=8, height=6)
save_pptx("../results/Num_isnvs_by_vaccine_gene.pptx", width=8, height=6)
mean(df_plot_n_gene$n_per_kb[df_plot_n_gene$gene=="S"])
mean(df_plot_n_gene$n_per_kb[df_plot_n_gene$gene=="ORF1ab"])
wilcox.test(df_plot_n_gene$n_per_kb[df_plot_n_gene$gene=="S"], df_plot_n_gene$n_per_kb[df_plot_n_gene$gene=="ORF1ab"])


df_groups <- df_plot_n_gene %>% filter(gene %in% c("Full genome", "ORF1ab", "S")) %>% select(gene) %>% unique() 
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	gene_i <- df_groups$gene[i]
	df_tmp <- df_plot_n_gene %>% filter(gene==gene_i)
	df_tmp_grps <- df_tmp %>% select(Vaccine, lineage_sim) %>% unique()
	df_tmp_grps$x_grps <- df_tmp_grps %>% apply(1,paste0, collapse="_")
	mat_pairs <- combn(df_tmp_grps$x_grps, 2)
	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		value1 <- df_tmp$n_per_kb[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var1] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var1]]
		value2 <- df_tmp$n_per_kb[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var2] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var2]]
		rst <- wilcox.test(value1, value2)
		tibble(var1=var1, var2=var2, var1_mean=mean(value1,na.rm=T), var2_mean=mean(value2,na.rm=T), median_var1=median(value1,na.rm=T), median_var2=median(value2,na.rm=T), p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	df_rst <- df_rst %>% mutate(gene=gene_i)
	df_rst %>% arrange(p_value)
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test %>% filter(p_value<0.05) # four pairs of comparision are significant (P<0.05)
write_csv(df_wilc_test , "../results/df_test_num_isnvs_by_vaccine_gene.csv")

### Allele frequency of iSNV mutations 
df_plot_full <- bind_rows(df_plot, df_plot %>% mutate(gene="Full genome"))
ggplot(df_plot, aes(x=Vaccine, y=sec_freq, color=lineage_sim))+ # freq of isnvs
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	ylab("Minor allele frequency (MAF)")+
	xlab("Group")+
	theme_classic()+
	NULL
ggsave("../results/Freq_isnvs_by_vaccine.pdf", width=6, height=4)
ggplot(df_plot_full %>% filter(gene %in% c("Full genome", "ORF1ab", "S")), aes(x=Vaccine, y=sec_freq, color=lineage_sim))+ # freq of isnvs
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_grid(rows=vars(gene))+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	ylab("Minor allele frequency (MAF)")+
	xlab("Group")+
	theme_classic()+
	NULL
ggsave("../results/Freq_isnvs_by_vaccine_gene.pdf", width=8, height=6)

df_groups <- df_plot_full %>% filter(gene %in% c("Full genome", "ORF1ab", "S")) %>% select(gene) %>% unique() 
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	gene_i <- df_groups$gene[i]
	df_tmp <- df_plot_full %>% filter(gene==gene_i)
	df_tmp_grps <- df_tmp %>% select(Vaccine, lineage_sim) %>% unique()
	df_tmp_grps$x_grps <- df_tmp_grps %>% apply(1,paste0, collapse="_")
	mat_pairs <- combn(df_tmp_grps$x_grps, 2)
	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		value1 <- df_tmp$sec_freq[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var1] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var1]]
		value2 <- df_tmp$sec_freq[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var2] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var2]]
		rst <- wilcox.test(value1, value2)
		tibble(var1=var1, var2=var2, var1_mean=mean(value1,na.rm=T), var2_mean=mean(value2,na.rm=T), median_var1=median(value1,na.rm=T), median_var2=median(value2,na.rm=T), p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	df_rst <- df_rst %>% mutate(gene=gene_i)
	df_rst %>% arrange(p_value)
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test %>% filter(p_value<0.05) # four pairs of comparision are significant (P<0.05)
write_csv(df_wilc_test , "../results/df_test_freq_isnvs_by_vaccine_gene.csv")


### Diversity of samples
#### we use nucleotide diversity pi here (Ref: https://academic.oup.com/ve/article/5/1/vey041/5304643, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4316684/)
source("./helper/cal_nu_diveristy_pi.r")
# all(unique(df_plot$gene) %in% df_orf_sim$sequence)

list_samples <- df_plot$sample %>% unique()
mclapply(list_samples, function(sample_t) {
	print(sample_t)
	# sample_t <- list_samples[1]
	df_mut_t <- df_plot %>% filter(sample==sample_t)
	Generate_cvf_files(df_mut_t)
	run_snpgenie(sample_t)
}, mc.cores=8)

df_diveristy_pop <-  read_snpgenie_rst(list.files("../results/snpgenie/", pattern="population_summary", recursive=T, full.names=T))
df_bam_t <- df_bam_rst %>% filter(reads_pp>=100) %>% group_by(sample) %>% summarise(n=n())
df_diveristy_pop$sample <- gsub("\\.vcf", "", df_diveristy_pop$file)
df_diveristy_pop <- left_join(df_diveristy_pop, df_bam_t, "sample")
df_diveristy_pop$scale_fct <- 29903/df_diveristy_pop$n # scale
df_diveristy_pop <- left_join(df_diveristy_pop, df_plot %>% select(lineage_sim, sample, Vaccine) %>% unique(), "sample")
df_diveristy_pop <- bind_cols(df_diveristy_pop %>% select(!contains("pi")), df_diveristy_pop %>% select(contains("pi")) %>% mutate_all(as.numeric))
df_diveristy_pop$pi_scale <- df_diveristy_pop$pi*df_diveristy_pop$scale_fct

df_gene_div <-  read_snpgenie_rst(list.files("../results/snpgenie/", pattern="product", recursive=T, full.names=T))
df_scale_gene <- lapply(seq_len(nrow(df_orf_sim)), function(i) {
	# print(i)
	gene_i <- df_orf_sim$sequence[i]
	start_i <- df_orf_sim$start[i]
	stop_i <- df_orf_sim$stop[i]
	df_bam_rst %>% filter(reads_pp>=100) %>% filter(pos>=start_i & pos<=stop_i) %>% group_by(sample) %>% summarise(n=n()) %>% mutate(gene=gene_i, scale_fct=(stop_i-start_i+1)/n)
})
df_scale_gene <- bind_rows(df_scale_gene)
df_gene_div$gene <- gsub("gene-", "", df_gene_div$product)
df_gene_div$gene <- gsub("orf", "ORF", df_gene_div$gene)
df_gene_div$sample <- gsub("\\.vcf", "", df_gene_div$file)
df_gene_div <- left_join(df_gene_div, df_scale_gene, by=c("gene", "sample"))
df_gene_div <- left_join(df_gene_div, df_plot %>% select(lineage_sim, sample, Vaccine) %>% unique(), "sample")
df_gene_div <- bind_cols(df_gene_div %>% select(!contains("pi")), df_gene_div %>% select(contains("pi")) %>% mutate_all(as.numeric))
df_gene_div$pi <- (as.numeric(df_gene_div$N_diffs) + as.numeric(df_gene_div$S_diffs))/(as.numeric(df_gene_div$N_sites) + as.numeric(df_gene_div$S_sites))
df_gene_div$pi_scale <- df_gene_div$pi*df_gene_div$scale_fct

# plot for pi
df_diveristy_pop$gene <- "Full genome"
df_plot_pi <- bind_rows(df_gene_div, df_diveristy_pop)

ggplot(df_plot_pi %>% filter(gene %in% c("Full genome", "ORF1ab", "S")), aes(x=`Vaccine`, y=pi_scale, color=lineage_sim))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_grid(rows=vars(gene))+
	# ylab(expression(pi[N]~"/"~pi["S"]))+
	ylab(expression("Nucleotide diversity"~pi))+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	theme_classic()+
	# theme_minimal()+
	NULL
ggsave("../results/diveristy_pi_gene.pdf", width=8, height=6)

df_groups <- df_plot_pi %>% filter(gene %in% c("Full genome", "ORF1ab", "S")) %>% select(gene) %>% unique() 
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	gene_i <- df_groups$gene[i]
	df_tmp <- df_plot_pi %>% filter(gene==df_groups$gene[i])
	df_tmp_grps <- df_tmp %>% select(Vaccine, lineage_sim) %>% unique()
	df_tmp_grps$x_grps <- df_tmp_grps %>% apply(1,paste0, collapse="_")
	mat_pairs <- combn(df_tmp_grps$x_grps, 2)
	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		value1 <- df_tmp$pi_scale[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var1] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var1]]
		value2 <- df_tmp$pi_scale[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var2] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var2]]
		rst <- wilcox.test(value1, value2)
		tibble(var1=var1, var2=var2, var1_mean=mean(value1,na.rm=T), var2_mean=mean(value2,na.rm=T), median_var1=median(value1,na.rm=T), median_var2=median(value2,na.rm=T), p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	df_rst <- df_rst %>% mutate(gene=gene_i)
	df_rst %>% arrange(p_value)
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test %>% filter(p_value<0.05) # four pairs of comparision are significant (P<0.05)
write_csv(df_wilc_test , "../results/df_test_diversity_pi_by_vaccine_gene.csv")

# plot for piN_piS
df_gene_div$piN_piS <- df_gene_div$piN/df_gene_div$piS
df_diveristy_pop$piN_piS <- df_diveristy_pop$piN/df_diveristy_pop$piS
df_plot_pi <- bind_rows(df_gene_div, df_diveristy_pop)
ggplot(df_plot_pi %>% filter(gene %in% c("Full genome", "ORF1ab", "S")), aes(x=`Vaccine`, y=piN_piS, color=lineage_sim))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_grid(rows=vars(gene))+
	ylab(expression(pi[N]~"/"~pi["S"]))+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	theme_classic()+
	# theme_minimal()+
	NULL
ggsave("../results/diveristy_piN_piS_gene.pdf", width=8, height=6)
save_pptx("../results/diveristy_piN_piS_gene.pptx", width=8, height=6)

df_groups <- df_plot_pi %>% filter(gene %in% c("Full genome", "ORF1ab", "S")) %>% select(gene) %>% unique() 
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	gene_i <- df_groups$gene[i]
	df_tmp <- df_plot_pi %>% filter(gene==df_groups$gene[i])
	df_tmp_grps <- df_tmp %>% select(Vaccine, lineage_sim) %>% unique()
	df_tmp_grps$x_grps <- df_tmp_grps %>% apply(1,paste0, collapse="_")
	mat_pairs <- combn(df_tmp_grps$x_grps, 2)
	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		value1 <- df_tmp$piN_piS[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var1] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var1]]
		value2 <- df_tmp$piN_piS[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var2] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var2]]
		rst <- wilcox.test(value1, value2)
		tibble(var1=var1, var2=var2, var1_mean=mean(value1,na.rm=T), var2_mean=mean(value2,na.rm=T), median_var1=median(value1,na.rm=T), median_var2=median(value2,na.rm=T), p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	df_rst <- df_rst %>% mutate(gene=gene_i)	
	df_rst %>% arrange(p_value)
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test %>% filter(p_value<0.05) # four pairs of comparision are significant (P<0.05)
write_csv(df_wilc_test, "../results/df_test_diversity_piN_piS_by_vaccine_gene.csv")

#### mutations by gene and epitopes
df_epitope <- readxl::read_excel("../data/df_t_cell_epitope.xlsx")
df_epitope$found_in <- sapply(strsplit(df_epitope$`Reference(s)*`, ", "), length)
source("./helper/annotate_gene.r")

df_plot$effect_sim <- df_plot$effect
df_plot$effect_sim[grepl("stream", df_plot$effect_sim)] <- "UTR"
df_plot$effect_sim[grepl("stop", df_plot$effect_sim)] <- "missense_variant"
df_plot$effect_sim <- gsub("_", " ", df_plot$effect_sim)
table(df_plot$effect_sim)

df_plot <- bind_cols(df_plot, get_orf(df_plot$X2))

##### annotate epitope https://www.sciencedirect.com/science/article/pii/S1931312821002389?via%3Dihub#app2

table(df_plot$gene, df_plot$gene_nsp)
df_plot$epitope_type <- NA
df_plot$epitope_found_in <- NA

sapply(seq_len(nrow(df_plot)), function(i) {
	# print(i)
	gene_i <- toupper(df_plot$gene_nsp[i])
	pos_nsp_i <- df_plot$pos_orf[i]
	df_epitope_i <- df_epitope %>% filter(toupper(Antigen)==gene_i & Start<=pos_nsp_i & End >=pos_nsp_i)
	if(nrow(df_epitope_i)==0){
		return("")
	} else {
		df_plot$epitope_type[i] <<- paste0(sort(unique(df_epitope_i$`Restriction`)), collapse=" and ")
		df_plot$epitope_found_in[i] <<- max(df_epitope_i$found_in)
		return("")
	}	
})

color_epitope <- pal_aaas()(3)
names(color_epitope) <- sort(unique(df_plot$epitope_type))

# TODO
#### The mutations on spike
df_plot$pos_aa <- as.numeric(sapply(strsplit(df_plot$pos_aa, "\\/"), function(x) {x[1]}))
assign_spike_aa_pos <- function(pos_spike) {
	df_s_r <- read_tsv("./helper/spike_region.tsv")
	out <- sapply(pos_spike, function(x) {
		tmp <- df_s_r$Region[df_s_r$start<=x & df_s_r$stop>=x]
		if(length(tmp)!=1){return(NA)}else{return(tmp)}
	})	
	factor(out, levels=df_s_r$Region)
	# out
}

df_plot2 <- df_plot %>% filter(gene=="S" & effect_sim!="UTR") %>% mutate(spike_region=assign_spike_aa_pos(pos_aa)) %>% group_by(X2, mut_aa, Vaccine, lineage_sim, spike_region, effect_sim, epitope_type, epitope_found_in) %>% mutate(mut_aa=gsub("p\\.", "", mut_aa)) %>% mutate(mut_aa=factor(mut_aa, levels=unique(mut_aa))) %>% summarise(n=n(), n_grp=sum(df_plot$Vaccine==Vaccine[1] & df_plot$lineage_sim==lineage_sim[1]), pect=n/n_grp) %>% ungroup()
df_plot2$epitope_type[df_plot2$epitope_found_in<2] <- NA
mut_common <- names(table(df_plot2$mut_aa)[table(df_plot2$mut_aa)>1]) # commonly seen in more than one group
df_plot2 <- df_plot2 %>% filter(mut_aa %in% mut_common)

tmp <- as.character(unique(df_plot2$spike_region))
names(tmp) <- tmp
tmp[grepl("-", tmp)] <- ""

facet_labeller <- function(variable, value) {
	tmp
}

p_spike <- ggplot(df_plot2)+
 	geom_tile(aes(x=mut_aa, y=paste(Vaccine, lineage_sim), fill=epitope_type))+
	facet_grid(effect_sim ~ spike_region, scales="free_x",space="free_x", labeller=labeller(spike_region=as_labeller(facet_labeller)))+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
	xlab("Position")+
	ylab("Group")+
	# scale_fill_viridis_c(name="Proportion. of\nsamples")+
	scale_fill_manual(name="T cell epitope", values=color_epitope, na.value="grey")+
	ggtitle("Spike")+
	NULL
ggsave("../results/mut_isnvs_by_vaccine_spike.pdf", width=12, height=4)

# df_plot %>% filter(gene=="S" & effect_sim!="UTR") %>% mutate(spike_region=assign_spike_aa_pos(pos_aa)) %>% arrange(X2) %>% mutate(mut_nt=factor(paste0(X4, X2, X5), levels=unique(paste0(X4, X2, X5)))) %>% group_by(mut_nt, mut_aa, Vaccine, lineage_sim, spike_region, effect_sim) %>% summarise(n=n(), n_grp=sum(df_plot$Vaccine==Vaccine[1] & df_plot$lineage_sim==lineage_sim[1]), pect=n/n_grp) %>% 
# 	ggplot()+
#  	geom_tile(aes(x=mut_nt, y=paste(Vaccine, lineage_sim), fill=round(pect,2)))+
# 	facet_grid(effect_sim ~ spike_region, scales="free_x",space="free_x", labeller=labeller(spike_region=as_labeller(facet_labeller)))+
# 	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
# 	xlab("Position")+
# 	ylab("Group")+
# 	scale_fill_viridis_c(name="Proportion. of\nsamples")+
# 	NULL
# ggsave("../results/mut_isnvs_by_vaccine_spike_nt.pdf", width=12)

### The mutations on Other genes
df_plot3 <- df_plot %>% filter(gene=="ORF1ab" & effect_sim!="UTR") %>% group_by(X2, mut_aa, Vaccine, lineage_sim, effect_sim, epitope_type, epitope_found_in) %>% mutate(mut_aa=gsub("p\\.", "", mut_aa)) %>% mutate(mut_aa=factor(mut_aa, levels=unique(mut_aa))) %>% summarise(n=n(), n_grp=sum(df_plot$Vaccine==Vaccine[1] & df_plot$lineage_sim==lineage_sim[1]), pect=n/n_grp) %>% ungroup()
df_plot3$epitope_type[df_plot3$epitope_found_in<2] <- NA
df_plot3 <- df_plot3 %>% filter(!grepl("Ter", mut_aa))
mut_common <- names(table(df_plot3$mut_aa)[table(df_plot3$mut_aa)>1]) # commonly seen in more than one group
df_plot3 <- df_plot3 %>% filter(mut_aa %in% mut_common)

p_orf1ab <- ggplot(df_plot3)+
 	geom_tile(aes(x=mut_aa, y=paste(Vaccine, lineage_sim), fill=epitope_type))+
	facet_grid(effect_sim ~ ., scales="free_x",space="free_x")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
	xlab("Position")+
	ylab("Group")+
	# scale_fill_viridis_c(name="Proportion. of\nsamples")+
	scale_fill_manual(name="T cell epitope", values=color_epitope, na.value="grey")+
	ggtitle("ORF1ab")+
  	# theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	NULL
ggsave("../results/mut_isnvs_by_vaccine_orf1ab.pdf", width=12, height=4)

p_orf1ab/p_spike + plot_layout(guides = "collect")
ggsave("../results/mut_isnvs_by_vaccine_combine.pdf", width=10, height=8)

df_plot %>% filter(grepl("D6304G", mut_aa)) %>% .$gene_nsp
df_plot %>% filter(grepl("D6304G", mut_aa)) %>% .$pos_orf
