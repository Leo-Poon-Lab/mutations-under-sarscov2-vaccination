library(tidyverse)
library(Biostrings)
library(ggsci)
library(patchwork)
library(parallel)

df_orf_sim <- read_csv("../data/ORF_SCoV2_sim.csv")
df_snvs_meta_add_qc <- read_csv("../results/df_snvs_meta_add_qc_bam.csv", guess_max=60000)
df_meta <- read_csv("../results/df_samples.csv", guess_max=100000)
df_meta <- df_meta %>% filter(lineage_sim != "22B (Omicron, BA.5.*)")
df_orf_sim$length_full <- df_orf_sim$stop-df_orf_sim$start+1

load("../results/df_bam_rst_full_notpp.rdata")

# table(df_meta$lineage_sim)
# dates_omicron_1 <- df_meta %>% filter(lineage_sim=="21M (Omicron, BA.2.*)") %>% .$collection_date
# range(dates_omicron_1, na.rm=T)
# dates_omicron_2 <- df_meta %>% filter(lineage_sim=="22B (Omicron, BA.5.*)") %>% .$collection_date
# range(dates_omicron_2, na.rm=T)

## 1. analyse the number of mutations, the frequency of mutations, diversity.
### Number of iSNV mutations 
df_tmp <- df_snvs_meta_add_qc %>% group_by(sample) %>% summarise(n=n())
df_plot <- left_join(df_snvs_meta_add_qc, df_tmp, "sample")
df_bam_rst_depth <- df_bam_rst %>% group_by(sample) %>% filter(reads_all>=100) %>% summarise(length_gene=n()) 

df_plot_n <- df_plot %>% select(sample, n, Vaccine, lineage_sim) %>% unique()
df_plot_n <- left_join(df_meta %>% select(sample, Vaccine, lineage_sim), df_plot_n)
df_plot_n$n[is.na(df_plot_n$n)] <- 0

df_plot_n <- left_join(df_plot_n, df_bam_rst_depth)
df_plot_n$n_per_pos <- df_plot_n$n/df_plot_n$length_gene
df_plot_n$n_per_kb <- df_plot_n$n_per_pos*1000
mean(df_plot_n$n) # mean number of iSNVs per sample: 9.8
median(df_plot_n$n) # median number of iSNVs per sample: 5
df_plot_n$gene <- "Full genome"
df_plot_n$length_full <- 29903

df_samples_gene <- full_join(df_meta %>% select(sample, Vaccine, lineage_sim), df_orf_sim %>% mutate(gene=sequence) %>% select(gene), by = character())
df_plot_n_gene <- df_snvs_meta_add_qc %>% group_by(sample, Vaccine, lineage_sim, gene) %>% summarise(n=n()) %>% ungroup()

df_plot_n_gene <- left_join(df_samples_gene, df_plot_n_gene)
df_plot_n_gene$n[is.na(df_plot_n_gene$n)] <- 0

pos_lab <- sapply(1:29903, function(x){
	tmp <- df_orf_sim$sequence[which(df_orf_sim$start<=x & df_orf_sim$stop>=x)]
	if(length(tmp)==0){return(NA)}else{tmp}
})
df_tmp <- tibble(pos=1:29903, gene=pos_lab)
df_tmp <- df_bam_rst %>% left_join(df_tmp, "pos") %>% filter(reads_all>=100) %>% group_by(sample,gene) %>% summarise(length_gene=n()) %>% ungroup()

df_plot_n_gene <- left_join(df_plot_n_gene, df_tmp, c("sample", "gene"))
df_plot_n_gene$n_per_kb <- df_plot_n_gene$n/df_plot_n_gene$length_gene*1000

df_plot_n_gene <- left_join(df_plot_n_gene, df_orf_sim %>% mutate(gene=sequence) %>% select(gene, length_full))

df_plot_n_gene <- bind_rows(df_plot_n_gene, df_plot_n %>% select(-n_per_pos))
df_plot_n_gene$n[is.na(df_plot_n_gene$n)] <- 0
df_plot_n_gene$n_per_kb[is.na(df_plot_n_gene$n_per_kb)] <- 0
save(df_plot_n_gene, file="../results/df_plot_n_gene.rdata")

### Diversity of samples
#### we use nucleotide diversity pi here (Ref: https://academic.oup.com/ve/article/5/1/vey041/5304643, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4316684/)
source("./helper/cal_nu_diveristy_pi.r")

## run diversity analysis by group
df_tmp <- df_plot %>% group_by(Vaccine,lineage_sim) %>% summarise(n=n())

dir.create("~/Downloads/snpgenie/")
commands <- mclapply(seq_len(nrow(df_tmp)), function(i){
	# i=1
	vaccine_i <- df_tmp$Vaccine[i]
	lineage_i <<- df_tmp$lineage_sim[i]
	group_i <- paste0(vaccine_i, "-", lineage_i)
	wd_group_i <- paste0("/Users/haogao/Downloads/snpgenie/", group_i)
	wd_group_complete_i <- wd_group_i
	dir.create(wd_group_i)
	
	df_plot_i <- df_plot %>% filter(Vaccine==vaccine_i) %>% filter(lineage_sim==lineage_i)
	commands_i <- sapply(unique(df_plot_i$sample), function(sample_t){
		# sample_t=df_plot_i$sample[1]
		df_vcf_i <- Generate_cvf_files(df_plot_i %>% filter(sample==sample_t))
		file_vcf_i <- paste0(wd_group_i, "/", sample_t, ".vcf")
		write_tsv(df_vcf_i, file_vcf_i)	
		paste0("~/softwares/SNPGenie/snpgenie.pl ", "--fastafile=reference.fasta --gtffile=SCoV2_genes.gtf --snpreport=", sample_t, '.vcf --vcfformat=2 --workdir="', wd_group_complete_i, '/" --outdir=', sample_t, "  2>> parallele.out.err > /dev/null")
	})
	
	file.copy("../data/reference.fasta", paste0(wd_group_i, "/reference.fasta"), overwrite=T)
	file.copy("../data/SCoV2_genes.gtf", paste0(wd_group_i, "/SCoV2_genes.gtf"), overwrite=T)	

	cur_wd <- getwd()
	setwd(wd_group_i)
	# system("rm -r SNPGenie_Results")
	# system("rm temp*")
	# system("~/softwares/SNPGenie/snpgenie.pl --vcfformat=2")
	setwd(cur_wd)
	commands_i
}, mc.cores=10)

### parallel running
# snpgenie.pl --fastafile=chr21.fa --gtffile=chr21_genes.gtf --snpreport=chr21_GATK.vcf -vcfformat=4 --minfreq=0.001 --workdir=/home/ohta/human/data/ --outdir=SNPGenie_Results
commands <- unlist(commands)

# system(commands[grep("WHP3907", commands)])

groups <- cut(seq_len(length(commands)), 10)
batch_cmds <- sapply(seq_along(levels(groups)), function(i) {
	group <- levels(groups)[i]
	idx <- which(groups==group)
	out_name <- paste0("/Users/haogao/Downloads/parallel_diversity_analysis_",i, ".sh")
	writeLines(c("#!/bin/bash", paste0(commands[idx], " &"), "wait"), out_name)
	system(paste0("chmod 755 ", out_name))
	out_name
})
batch_cmds <- unlist(batch_cmds)

writeLines(c("#!/bin/bash", batch_cmds), "./parallel_diversity_analysis.sh")
system("./parallel_diversity_analysis.sh")
# system("cp -r /Users/haogao/Downloads/snpgenie ../results/")

