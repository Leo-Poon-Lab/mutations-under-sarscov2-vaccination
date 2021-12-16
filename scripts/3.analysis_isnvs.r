# analyse the mutations found in the selected samples
## 1. read the snvs identified by ivar.
## 2. filter out the mutations with low quality/support
## 3. determine the regions with no data / no support for each sample.
## 4. analyse the number of mutations, the frequency of mutations, diversity.

library(tidyverse)
library(Biostrings)
library(ggsci)
library(parallel)

df_epitope <- readxl::read_excel("../data/df_t_cell_epitope.xlsx")
df_epitope$found_in <- sapply(strsplit(df_epitope$`Reference(s)*`, ", "), length)

colors_lineage=c("#ff7f00", "#33a02c", "#1f78b4")
names(colors_lineage) <- c("Alpha", "B.1", "Delta")
colors_vaccine_type=pal_jco()(4)
names(colors_vaccine_type)=c("Adeno", "Inactivated", "Non-vaccinated", "mRNA")

df_meta <- read_csv("../data/df_samples_vaccine_bt_study.csv")
df_meta$`Ct value`[df_meta$`Ct value`==0] <- NA

df_vaccine_type <- readxl::read_excel("../data/Vaccination case for breakthrough infection .xlsx", skip=1)
df_vaccine_type <- df_vaccine_type %>% select(Vaccine, `Vaccine type`) %>% unique()
vaccine_of_interest <- c("BioNTech", "Sinovac", "Non-vaccinated")

df_lin <- read_csv("../results/lineage.csv")
df_lin$Sample <- sapply(df_lin$taxon, function(x) {
	strsplit(x, "_")[[1]][2]
})
df_lin <- df_lin %>% select(Sample, lineage, scorpio_call)
df_lin$lineage_sim <- df_lin$lineage
df_lin$lineage_sim[grepl("^AY\\.", df_lin$lineage_sim)] <- "Delta"
df_lin$lineage_sim[df_lin$lineage_sim == "B.1.617.2"] <- "Delta"
df_lin$lineage_sim[grepl("^B\\.1$", df_lin$lineage_sim)] <- "B.1"
df_lin$lineage_sim[grepl("^B\\.1\\.1\\.7$", df_lin$lineage_sim)] <- "Alpha"
df_lin$lineage_sim[grepl("^Q\\.", df_lin$lineage_sim)] <- "Alpha"

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
df_snvs <- read_snpeff(files_vcfs_csv)

## 2. filter out the mutations with low quality/support
## 3. determine the regions with no data / no support for each sample.
### using pysamstats
system("conda run -n base pysamstats --help")
source("./helper/pysamstats.r")
files_bam <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
files_bam_full <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$", full.names = T)
samples_alrd_done <- gsub("\\.tsv$", "", list.files("../results/pysamstats/"))
check1 <- files_bam %in% paste0(df_meta$Sample, "-trimmed.masked.bam")
check2 <- !files_bam %in% paste0(samples_alrd_done, "-trimmed.masked.bam")
files_bam_int <- files_bam_full[check1 & check2]

process_pysamstats(files_bam_int)
files_bam_rst_full <- list.files("../results/pysamstats/", "tsv$", full.names = T)
df_bam_rst <- mclapply(files_bam_rst_full, read_pysamstats, mc.cores=8)
df_bam_rst <- bind_rows(df_bam_rst)

### only nucleotide positions with at least 100 properly paired reads were considered
df_bam_rst_depth <- df_bam_rst %>% group_by(sample) %>% filter(reads_pp>=100) %>% summarise(n_high_depth=n()) 
plot(df_bam_rst_depth$n_high_depth)
samples_hq <- df_bam_rst_depth %>% filter(n_high_depth>=27000) %>% .$sample

### only include samples with Ct value < 25
samples_hq <- samples_hq[samples_hq %in% df_meta$Sample[round(as.numeric(df_meta$`Ct value`)) < 25]]
writeLines(samples_hq, "../results/samples_hq.txt")
samples_hq <- readLines("../results/samples_hq.txt")

### only samples with 27000 sites of sequencing depth >= 100 were considered
names(df_snvs)[names(df_snvs)=="sample"] <- "Sample"
df_snvs_meta <- left_join(df_snvs, df_meta, "Sample")
df_snvs_meta <- df_snvs_meta %>% filter(Sample %in% samples_hq)
df_snvs_meta <- left_join(df_snvs_meta, df_vaccine_type)
df_snvs_meta$`Vaccine type`[df_snvs_meta$Vaccine=="Non-vaccinated"] <- "Non-vaccinated"
df_snvs_meta <- left_join(df_snvs_meta, df_lin, "Sample")
write_csv(df_snvs_meta, "../results/df_mut_meta.csv")

### add readcounts from pysamstats results
df_snvs_meta_append <- mclapply(seq_len(nrow(df_snvs_meta)), function(i){
	print(i)
	df_snvs_meta_i <- df_snvs_meta[i,]
	sample_i <- df_snvs_meta_i$Sample
	pos_i <- df_snvs_meta_i$X2
	ref_i <- df_snvs_meta_i$X4
	alt_i <- df_snvs_meta_i$X5
	df_bam_rst_i <- df_bam_rst %>% filter(sample==sample_i, pos==pos_i)
	depth_i <- df_bam_rst_i$reads_all
	depth_pp_i <- df_bam_rst_i$reads_pp # count on paired reads

	if(nchar(ref_i) > nchar(alt_i)){ # deletion
		alt_i <- "deletions"
	} else if(nchar(ref_i) < nchar(alt_i)){ # insertion
		alt_i <- "insertions"
	} 
	alt_fwd_i <- df_bam_rst_i[[paste0(alt_i, "_pp_fwd")]]
	alt_rev_i <- df_bam_rst_i[[paste0(alt_i, "_pp_rev")]]

	bases <- c("deletions", "insertions", "A", "C", "T", "G")
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

	tibble(type_alt=alt_i, depth_all=depth_i, depth_pp=depth_pp_i, depth_pp_fwd=df_bam_rst_i$reads_pp_fwd, depth_pp_rev=df_bam_rst_i$reads_pp_rev, alt_fwd=alt_fwd_i, alt_rev=alt_rev_i, alle_freq=alle_freq_i, con_base=consensus_base, con_fwd=con_fwd_i, con_rev=con_rev_i, con_freq=con_freq_i, sec_base=secondary_base, sec_fwd=sec_fwd_i, sec_rev=sec_rev_i, sec_freq=sec_freq_i)
}, mc.cores = 8)
df_snvs_meta_append <- bind_rows(df_snvs_meta_append)
df_mut_meta_add <- bind_cols(df_snvs_meta, df_snvs_meta_append)
write_csv(df_mut_meta_add, "../results/df_mut_meta_add.csv")
# df_mut_meta_add <- read_csv("../results/df_mut_meta_add.csv")

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

#### test different MAF cutoffs
#### The minor allele frequency (secondary most base frequency) cutoff for determining iSNVs
cut_offs <- 3:10/100
df_plot <- lapply(cut_offs, function(x) {
	df_snvs_meta_add_qc %>% group_by(Sample, Vaccine, `Vaccine type`) %>% filter(sec_freq>=x) %>% summarise(n=n(), Ct_value=`Ct value`[1], alle_freq=as.character(x))
})
df_plot <- bind_rows(df_plot)

ggplot(df_plot) + # iSNV threshold
	geom_point(aes(x=as.numeric(Ct_value), y=n, color=`Vaccine type`), alpha=0.6)+
	facet_wrap(vars(alle_freq), ncol=2)+
	xlab("Ct value")+
	ylab("No. of identified iSNVs")+
	scale_color_manual(values=colors_vaccine_type)
ggsave("../results/MAF-cutoff.pdf", height = 8, width = 6)

#### 2. The mutations should be supported by at least 3% of the total reads (minor allele frequency >=3%).
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(sec_freq>=0.03) 
df_snvs_meta_add_qc$sample <- df_snvs_meta_add_qc$Sample
write_csv(df_snvs_meta_add_qc, "../results/df_snvs_meta_add_qc_ivar.csv")

df_snvs_meta_add_qc_bam <- read_csv("../results/df_snvs_meta_add_qc_bam.csv")

# df_tmp <- df_snvs_meta_add_qc_bam[df_snvs_meta_add_qc_bam$con_base!=df_snvs_meta_add_qc_bam$ref,]
# write_csv(df_tmp, "../results/tmp.csv")

# (check <- paste(df_snvs_meta_add_qc$Sample, df_snvs_meta_add_qc$X2) %in% paste(df_snvs_meta_add_qc_bam$Sample, df_snvs_meta_add_qc_bam$pos))
# as.character(df_snvs_meta_add_qc[!check,])

# (check <- paste(df_snvs_meta_add_qc_bam$Sample, df_snvs_meta_add_qc_bam$pos) %in% paste(df_snvs_meta_add_qc$Sample, df_snvs_meta_add_qc$X2))
# df_tmp <- df_snvs_meta_add_qc_bam[!check,]
# write_csv(df_tmp, "../results/tmp.csv")

## 4. analyse the number of mutations, the frequency of mutations, diversity.
### Number of iSNV mutations 
df_plot <- df_snvs_meta_add_qc_bam %>% filter(sec_freq>0.05)
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
ggsave("../results/Num_isnvs_by_vaccine_type.pdf")

ggplot(df_plot_n, aes(x=`Vaccine`, y=n_per_kb, color=lineage_sim))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	ylab("Numer of iSNVs per Kb")+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	theme_classic()+
	# theme_minimal()+
	NULL
ggsave("../results/Num_isnvs_by_vaccine.pdf")

##### Comparison within unvaccinated
wilcox.test(df_plot_n$n_per_kb[df_plot_n$lineage_sim=="Alpha" & df_plot_n$Vaccine=="Non-vaccinated"], df_plot_n$n_per_kb[df_plot_n$lineage_sim=="Delta" & df_plot_n$Vaccine=="Non-vaccinated"])
wilcox.test(df_plot_n$n_per_kb[df_plot_n$lineage_sim=="Alpha" & df_plot_n$Vaccine=="Non-vaccinated"], df_plot_n$n_per_kb[df_plot_n$lineage_sim=="B.1" & df_plot_n$Vaccine=="Non-vaccinated"]) 
wilcox.test(df_plot_n$n_per_kb[df_plot_n$lineage_sim=="Delta" & df_plot_n$Vaccine=="Non-vaccinated"], df_plot_n$n_per_kb[df_plot_n$lineage_sim=="B.1" & df_plot_n$Vaccine=="Non-vaccinated"]) 

##### Comparison between unvaccinated and vaccinated
wilcox.test(df_plot_n$n_per_kb[df_plot_n$lineage_sim=="Delta" & df_plot_n$Vaccine=="Non-vaccinated"], df_plot_n$n_per_kb[df_plot_n$lineage_sim=="Delta" & df_plot_n$Vaccine!="Non-vaccinated"]) # Vac Delta and Unvac Delta
wilcox.test(df_plot_n$n_per_kb[df_plot_n$lineage_sim=="B.1" & df_plot_n$Vaccine=="Non-vaccinated"], df_plot_n$n_per_kb[df_plot_n$lineage_sim=="Delta" & df_plot_n$Vaccine!="Non-vaccinated"]) # Vac Delta and Unvac B.1
wilcox.test(df_plot_n$n_per_kb[df_plot_n$Vaccine=="Non-vaccinated" & df_plot_n$lineage_sim=="Delta"], df_plot_n$n_per_kb[df_plot_n$Vaccine=="BioNTech" & df_plot_n$lineage_sim=="Delta"]) # Vac BioNTech Delta and Unvac Delta
wilcox.test(df_plot_n$n_per_kb[df_plot_n$Vaccine=="Non-vaccinated" & df_plot_n$lineage_sim=="Delta"], df_plot_n$n_per_kb[df_plot_n$Vaccine=="Sinovac" & df_plot_n$lineage_sim=="Delta"])

wilcox.test(df_plot_n$n_per_kb[df_plot_n$`Vaccine type`=="Non-vaccinated"], df_plot_n$n_per_kb[df_plot_n$`Vaccine type`=="Inactivated"])
wilcox.test(df_plot_n$n_per_kb[df_plot_n$`Vaccine type`=="Non-vaccinated" & df_plot_n$lineage_sim=="Delta"], df_plot_n$n_per_kb[df_plot_n$`Vaccine type`=="Inactivated" & df_plot_n$lineage_sim=="Delta"])
wilcox.test(df_plot_n$n_per_kb[df_plot_n$`Vaccine type`=="Non-vaccinated"], df_plot_n$n_per_kb[df_plot_n$`Vaccine type`=="mRNA"]) # Unvac all and mRNA Delta/All

##### Comparison within vaccinated
wilcox.test(df_plot_n$n_per_kb[df_plot_n$Vaccine=="Sinovac"], df_plot_n$n_per_kb[df_plot_n$Vaccine=="BioNTech"])

### Diversity of samples
#### we use nucleotide diversity pi here (Ref: https://academic.oup.com/ve/article/5/1/vey041/5304643, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4316684/)
source("./helper/cal_nu_diveristy_pi.r")
df_orf_sim <- read_csv("../data/ORF_SCoV2_sim.csv")
# all(unique(df_plot$gene) %in% df_orf_sim$sequence)

list_samples <- df_plot$sample %>% unique()
df_diveristy <- mclapply(list_samples, function(sample_t) {
	print(sample_t)
	# sample_t <- list_samples[1]
	df_mut_t <- df_plot %>% filter(sample==sample_t)
	df_bam_t <- df_bam_rst %>% filter(sample==sample_t)	
	df_out <- tibble(gene="Full", pi_all=NA, pi_s=NA, pi_n=NA)

	## pi over the full genome
	dls <- sapply(df_mut_t$X2, function(pos_i) {
		df_bam_t_i <- df_bam_t %>% filter(pos==pos_i)
		Cal_dl(c(df_bam_t_i$A_pp, df_bam_t_i$C_pp, df_bam_t_i$T_pp, df_bam_t_i$G_pp ))
	})
	length_full <- df_bam_t %>% filter(reads_pp>=100) %>% summarise(n_high_depth=n())  %>% .$n_high_depth
	df_out$pi_all <- Cal_pi(dls, length_full)
	## pi s
	df_out$pi_s <- Cal_pi(dls[df_mut_t$effect=="synonymous_variant"], length_full)
	## pi n
	df_out$pi_n <- Cal_pi(dls[df_mut_t$effect=="missense_variant"], length_full)

	## pi in different genes (pi n, pi s)
	df_out_gene <- lapply(seq_len(nrow(df_orf_sim)), function(i) {
		# print(i)
		gene_i <- df_orf_sim$sequence[i]
		df_mut_t_i <- df_mut_t %>% filter(gene==gene_i)
		df_out_i <- tibble(gene=gene_i, pi_all=NA, pi_s=NA, pi_n=NA)
		if(nrow(df_mut_t_i)==0){return(df_out_i)}
		# pi all
		dls <- sapply(df_mut_t_i$X2, function(pos_i) {
			df_bam_t_i <- df_bam_t %>% filter(pos==pos_i)
			Cal_dl(c(df_bam_t_i$A_pp, df_bam_t_i$C_pp, df_bam_t_i$T_pp, df_bam_t_i$G_pp ))
		})
		df_bam_t_i <- df_bam_t %>% filter(pos>=df_orf_sim$start[i] & pos<=df_orf_sim$stop[i])
		length_full <- df_bam_t_i %>% filter(reads_pp>=100) %>% summarise(n_high_depth=n())  %>% .$n_high_depth
		df_out_i$pi_all <- Cal_pi(dls, length_full)
		# pi s
		df_out_i$pi_s <- Cal_pi(dls[df_mut_t_i$effect=="synonymous_variant"], length_full)
		# pi n
		df_out_i$pi_n <- Cal_pi(dls[df_mut_t_i$effect=="missense_variant"], length_full)
		return(df_out_i)
	})
	df_out <- bind_rows(df_out, df_out_gene)
	df_out$sample=sample_t
	return(df_out)
}, mc.cores=8)
df_diveristy <- bind_rows(df_diveristy)
df_diveristy$pi_n_pi_s <- df_diveristy$pi_n/df_diveristy$pi_s
df_diveristy <- left_join(df_diveristy, df_plot %>% select(lineage_sim, sample, Vaccine) %>% unique(), "sample")

ggplot(df_diveristy %>% filter(gene %in% c("Full", "ORF1ab", "S")), aes(x=`Vaccine`, y=pi_n_pi_s, color=lineage_sim))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_wrap(vars(gene), scales="free_x", ncol=1)+
	ylab(expression(pi[N]~"/"~pi["S"]))+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	theme_classic()+
	# theme_minimal()+
	NULL
ggsave("../results/diveristy_pin_pis.pdf")

ggplot(df_diveristy %>% filter(gene %in% c("Full", "ORF1ab", "S")), aes(x=`Vaccine`, y=pi_all, color=lineage_sim))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_wrap(vars(gene), scales="free_x", ncol=1)+
	ylab("Nucleotide diversity")+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	theme_classic()+
	# theme_minimal()+
	NULL
ggsave("../results/tmp.pdf")

df_groups <- df_diveristy %>% filter(gene %in% c("Full", "ORF1ab", "S")) %>% select(gene) %>% unique() 
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	df_tmp <- df_diveristy %>% filter(gene==df_groups$gene[i])
	df_tmp_grps <- df_tmp %>% select(Vaccine, lineage_sim) %>% unique()
	df_tmp_grps$x_grps <- df_tmp_grps %>% apply(1,paste0, collapse="_")
	mat_pairs <- combn(df_tmp_grps$x_grps, 2)
	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		rst <- wilcox.test(df_tmp$pi_n_pi_s[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var1] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var1]], df_tmp$pi_n_pi_s[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var2] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var2]])
		tibble(var1=var1, var2=var2, p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	bind_cols(df_tmp %>% select(gene) %>% unique(), df_rst)	
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test %>% filter(p_value<0.05) # four pairs of comparision are significant (P<0.05)
write_csv(df_wilc_test %>% filter(p_value<0.05), "../results/df_test_num_isnvs_by_vaccine_gene.csv")



### Allele frequency of iSNV mutations 
ggplot(df_plot, aes(x=Vaccine, y=sec_freq, color=lineage_sim))+ # freq of isnvs
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	ylab("Minor allele frequency (MAF)")+
	xlab("Group")+
	theme_classic()+
	NULL

wilcox.test(df_plot$sec_freq[df_plot$Vaccine=="Non-vaccinated" & df_plot$lineage_sim=="Delta"], df_plot$sec_freq[df_plot$Vaccine=="BioNTech"])
wilcox.test(df_plot$sec_freq[df_plot$Vaccine=="Non-vaccinated" & df_plot$lineage_sim=="Delta"], df_plot$sec_freq[df_plot$Vaccine=="Sinovac"]) # Vac Sinovac Delta and Unvac Delta
wilcox.test(df_plot$sec_freq[df_plot$Vaccine=="Non-vaccinated" & df_plot$lineage_sim=="Delta"], df_plot$sec_freq[df_plot$Vaccine!="Non-vaccinated" & df_plot$lineage_sim=="Delta"]) 

#### Which gene
source("./helper/annotate_gene.r")

df_plot$effect_sim <- df_plot$effect
df_plot$effect_sim[grepl("stream", df_plot$effect_sim)] <- "UTR"
df_plot$effect_sim[grepl("stop", df_plot$effect_sim)] <- "missense_variant"
df_plot$effect_sim <- gsub("_", " ", df_plot$effect_sim)
table(df_plot$effect_sim)

df_plot <- bind_cols(df_plot, get_orf(df_plot$X2))

##### annotate epitope
table(df_plot$gene, df_plot$gene_nsp)

lapply(seq_len(nrow(df_plot)), function(i) {
	gene_i <- toupper(df_plot$gene_nsp[i])
	pos_nsp_i <- df_plot$pos_orf[i]
	df_epitope_i <- df_epitope %>% filter(toupper(Antigen)==gene_i & Start<=pos_nsp_i & End >=pos_nsp_i)
	df_epitope_i$Dominant
})

df_epitope$Antigen



df_plot2 <- df_plot %>% group_by(n_high_depth, Vaccine, lineage_sim, Sample, gene, effect_sim) %>% summarise(n=n()) %>% filter(effect_sim!="UTR") %>% filter(gene %in% c("ORF1ab", "S")) %>% ungroup()
df_plot2$n_per_kb <- df_plot2$n/df_plot2$n_high_depth*1000

df_plot2 %>% ggplot(aes(x=Vaccine, y=n_per_kb, color=lineage_sim))+
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_grid(rows=vars(effect_sim), cols=vars(gene), scales="free_x",space="free_x")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	ylab("Numer of iSNVs per Kb")+
	xlab("Group")+
	theme_classic()+
	# theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=0.5))+
	NULL
ggsave("../results/Num_isnvs_by_vaccine_gene.pdf")

df_groups <- df_plot2 %>% select(gene, effect_sim) %>% unique() 

df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	df_tmp <- df_plot2 %>% filter(gene==df_groups$gene[i], effect_sim==df_groups$effect_sim[i])
	df_tmp_grps <- df_tmp %>% select(Vaccine, lineage_sim) %>% unique()
	df_tmp_grps$x_grps <- df_tmp_grps %>% apply(1,paste0, collapse="_")
	mat_pairs <- combn(df_tmp_grps$x_grps, 2)
	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		rst <- wilcox.test(df_tmp$n_per_kb[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var1] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var1]], df_tmp$n_per_kb[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var2] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var2]])
		tibble(var1=var1, var2=var2, p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	bind_cols(df_tmp %>% select(gene, effect_sim) %>% unique(), df_rst)	
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test %>% filter(p_value<0.05) # four pairs of comparision are significant (P<0.05)
write_csv(df_wilc_test %>% filter(p_value<0.05), "../results/df_test_num_isnvs_by_vaccine_gene.csv")

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

df_plot2 <- df_plot %>% filter(gene=="S" & effect_sim!="UTR") %>% mutate(spike_region=assign_spike_aa_pos(pos_aa)) %>% group_by(X2, mut_aa, Vaccine, lineage_sim, spike_region, effect_sim) %>% mutate(mut_aa=gsub("p\\.", "", mut_aa)) %>% mutate(mut_aa=factor(mut_aa, levels=unique(mut_aa))) %>% summarise(n=n(), n_grp=sum(df_plot$Vaccine==Vaccine[1] & df_plot$lineage_sim==lineage_sim[1]), pect=n/n_grp) %>% ungroup()

tmp <- as.character(unique(df_plot2$spike_region))
names(tmp) <- tmp
tmp[grepl("-", tmp)] <- ""

facet_labeller <- function(variable, value) {
	tmp
}

ggplot(df_plot2)+
 	geom_tile(aes(x=mut_aa, y=paste(Vaccine, lineage_sim), fill=round(pect,2)))+
	facet_grid(effect_sim ~ spike_region, scales="free_x",space="free_x", labeller=labeller(spike_region=as_labeller(facet_labeller)))+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
	xlab("Position")+
	ylab("Group")+
	scale_fill_viridis_c(name="Proportion. of\nsamples")+
	NULL
ggsave("../results/Freq_isnvs_by_vaccine_spike.pdf", width=12)

df_plot %>% filter(gene=="S" & effect_sim!="UTR") %>% mutate(spike_region=assign_spike_aa_pos(pos_aa)) %>% arrange(X2) %>% mutate(mut_nt=factor(paste0(X4, X2, X5), levels=unique(paste0(X4, X2, X5)))) %>% group_by(mut_nt, mut_aa, Vaccine, lineage_sim, spike_region, effect_sim) %>% summarise(n=n(), n_grp=sum(df_plot$Vaccine==Vaccine[1] & df_plot$lineage_sim==lineage_sim[1]), pect=n/n_grp) %>% 
	ggplot()+
 	geom_tile(aes(x=mut_nt, y=paste(Vaccine, lineage_sim), fill=round(pect,2)))+
	facet_grid(effect_sim ~ spike_region, scales="free_x",space="free_x", labeller=labeller(spike_region=as_labeller(facet_labeller)))+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
	xlab("Position")+
	ylab("Group")+
	scale_fill_viridis_c(name="Proportion. of\nsamples")+
	NULL
ggsave("../results/Freq_isnvs_by_vaccine_spike_nt.pdf", width=12)

### The mutations on Other genes



### silent/non-silent mutations, dS/dN


### Diversity
