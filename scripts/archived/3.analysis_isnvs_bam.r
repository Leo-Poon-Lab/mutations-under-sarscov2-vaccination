# analyse the mutations found in the selected samples
## 1. read mpileup summary from pysamstats
## 2. filter out the mutations with low quality/support
## 3. determine the regions with no data / no support for each sample.
## 4. analyse the number of mutations, the frequency of mutations, diversity.

library(tidyverse)
library(Biostrings)
library(ggsci)
library(parallel)

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

## 1. read mpileup summary from pysamstats
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
df_snvs <- df_bam_rst

## 2. filter out the mutations with low quality/support
## 3. determine the regions with no data / no support for each sample.
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

### isnvs were called with filtering (QC)
#### 1. The identified isnvs should have at least 100 properly paired reads at the genomic positions 
df_snvs_meta <- df_snvs_meta %>% filter(reads_pp>=100)
#### 2. positions 1:100 and (29903-99):29903 should be removed; 
df_snvs_meta <- df_snvs_meta %>% filter(pos>100 & pos<(29903-99))
#### 3. exclude all positions in the PCR primer binding regions
df_primer_new <- read_tsv("../../2021-10-11_nCoV_primers/results/primers_20211011.bed", col_names=F)
df_primer_old <- read_tsv("../../2021-10-11_nCoV_primers/results/primers_old.bed", col_names=F)

pos_all <- unique(paste(df_snvs_meta$pos, df_snvs_meta$primer))
check <- sapply(pos_all, function(x) {
	pos_x <- as.numeric(strsplit(x, " ")[[1]][1])
	primer_x <- strsplit(x, " ")[[1]][2]
	if(primer_x=="new"){
		any((df_primer_new$X2 <= pos_x) & (df_primer_new$X3 >= pos_x))	
	} else {
		any((df_primer_old$X2 <= pos_x) & (df_primer_old$X3 >= pos_x))	
	}	
})
df_snvs_meta <- df_snvs_meta[!paste(df_snvs_meta$pos, df_snvs_meta$primer) %in% pos_all[check],]
#### 4. MAF >=0.03 (minimum depth of minor allele of 5 required, at least 1 on forward strand and 1 on reverse strand); 
df_snvs_meta <- df_snvs_meta %>% filter(matches_pp/reads_pp<0.97)
bases <- c("deletions", "insertions", "A", "C", "T", "G")
df_snvs_meta_append <- mclapply(seq_len(nrow(df_snvs_meta)), function(i) {
	num_bases <- c(df_snvs_meta$deletions[i], df_snvs_meta$insertions[i], df_snvs_meta$A_pp[i], df_snvs_meta$C_pp[i], df_snvs_meta$T_pp[i], df_snvs_meta$G_pp[i])
	
	check_consensus_base <- order(num_bases)[length(bases)]
	check_secondary_base <- order(num_bases)[length(bases)-1]
	consensus_base <- bases[check_consensus_base]
	secondary_base <- bases[check_secondary_base]

	con_fwd_i <- df_snvs_meta[[paste0(consensus_base, "_pp_fwd")]][i]
	con_rev_i <- df_snvs_meta[[paste0(consensus_base, "_pp_rev")]][i]
	sec_fwd_i <- df_snvs_meta[[paste0(secondary_base, "_pp_fwd")]][i]
	sec_rev_i <- df_snvs_meta[[paste0(secondary_base, "_pp_rev")]][i]

	con_freq_i <- (con_fwd_i+con_rev_i)/df_snvs_meta$reads_pp[i]
	sec_freq_i <- (sec_fwd_i+sec_rev_i)/df_snvs_meta$reads_pp[i]

	tibble(i=i, con_base=consensus_base, con_fwd=con_fwd_i, con_rev=con_rev_i, con_freq=con_freq_i, sec_base=secondary_base, sec_fwd=sec_fwd_i, sec_rev=sec_rev_i, sec_freq=sec_freq_i)
}, mc.cores=8)
df_snvs_meta_append <- bind_rows(df_snvs_meta_append)
df_mut_meta_add <- bind_cols(df_snvs_meta, df_snvs_meta_append)

df_snvs_meta_add_qc_bam <- df_mut_meta_add %>% filter(reads_pp>=100) %>% filter(sec_freq>=0.03) %>% filter((sec_fwd+sec_rev)>=5 & sec_fwd>=1 & sec_rev>=1) %>% filter(nchar(sec_base)==1 & nchar(con_base)==1)
df_snvs_meta_add_qc_bam <- df_snvs_meta_add_qc_bam %>% filter(Vaccine %in% vaccine_of_interest)
df_snvs_meta_add_qc_bam <- df_snvs_meta_add_qc_bam %>% filter(!(lineage_sim!="Delta" & Vaccine!="Non-vaccinated"))

write_csv(df_snvs_meta_add_qc_bam, "../results/df_snvs_meta_add_qc_bam.csv")

df_snvs_meta_add_qc <- df_snvs_meta_add_qc_bam
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

## 4. analyse the number of mutations, the frequency of mutations, diversity.
### Number of iSNV mutations 
df_plot <- df_snvs_meta_add_qc
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
	ylab("Numer of mutations per position")+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	NULL
ggsave("../results/Num_isnvs_by_vaccine_type.pdf")

ggplot(df_plot_n, aes(x=`Vaccine`, y=n_per_kb, color=lineage_sim))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1))+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	ylab("Numer of mutations per position")+
	xlab("Group")+
	scale_color_manual(name="Lineage", values=colors_lineage)+
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

### Number of iSNV mutations 
df_plot %>% filter(sec_freq>0.5)

ggplot(df_plot, aes(x=Vaccine, y=sec_freq, color=lineage_sim))+ # freq of isnvs
	geom_boxplot()+
	# facet_wrap(vars(Vaccine), scales="free_x", ncol=1)+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	NULL
ggsave("../results/tmp.pdf")

wilcox.test(df_plot$sec_freq[df_plot$Vaccine=="Non-vaccinated"], df_plot$sec_freq[df_plot$Vaccine=="BioNTech"])
wilcox.test(df_plot$sec_freq[df_plot$Vaccine=="Non-vaccinated"], df_plot$sec_freq[df_plot$Vaccine=="Sinovac"])
wilcox.test(df_plot$sec_freq[df_plot$Vaccine=="Non-vaccinated"], df_plot$sec_freq[df_plot$Vaccine!="Non-vaccinated"])


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
	df_bam_t <- df_bam_rst %>% filter(sample==sample_t) %>% filter(reads_pp>=100)
	df_out <- tibble(gene="Full", pi_all=NA, pi_s=NA, pi_n=NA)

	## pi over the full genome
	df_bam_t_i <- df_bam_t %>% select(A_pp, C_pp, T_pp, G_pp)
	dls <- apply(df_bam_t_i, 1, function(x) {Cal_dl(x)})
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

ggplot(df_diveristy, aes(x=`Vaccine`, y=pi_all, color=lineage_sim))+ # mutation rate
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

df_groups <- df_diveristy %>% select(gene) %>% unique() 
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	df_tmp <- df_diveristy %>% filter(gene==df_groups$gene[i])
	df_tmp_grps <- df_tmp %>% select(Vaccine, lineage_sim) %>% unique()
	df_tmp_grps$x_grps <- df_tmp_grps %>% apply(1,paste0, collapse="_")
	mat_pairs <- combn(df_tmp_grps$x_grps, 2)
	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		rst <- wilcox.test(df_tmp$pi_all[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var1] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var1]], df_tmp$pi_all[df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var2] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var2]])
		tibble(var1=var1, var2=var2, p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	bind_cols(df_tmp %>% select(gene) %>% unique(), df_rst)	
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test %>% filter(p_value<0.05) # four pairs of comparision are significant (P<0.05)
write_csv(df_wilc_test %>% filter(p_value<0.05), "../results/df_test_num_isnvs_by_vaccine_gene.csv")




### Which gene
### silent/non-silent mutations, dS/dN
df_plot


### Diversity
