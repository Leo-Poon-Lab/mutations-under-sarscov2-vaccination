# QC of the NGS data
# Analysis the consensus sequences
## conda activate pangolin
library(tidyverse)
library(Biostrings)
library(ggrepel)
library(slider)
library(ggsci)

df_samples <- read_csv("../data/df_samples.csv")

# QC
files_rc <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/read_count/", full.names = T)
check <- sapply(files_rc, function(x){
	any(sapply(paste0(df_samples$Sample, "_"), function(y){grepl(y, x)}))
})
files_rc_filter <- files_rc[check]
files_rc_filter <- files_rc_filter[!grepl("iseq", files_rc_filter)]
files_rc_filter <- sort(files_rc_filter, decreasing = T)
samples_rc <- sapply(files_rc_filter, function(x){
	tmp <- strsplit(x, "_")[[1]]
	tmp[grepl("WHP", tmp)]
})
files_rc_filter <- files_rc_filter[!duplicated(samples_rc)]

df_rc <- lapply(files_rc_filter, function(x){
	tmp <- read_csv(x)
})
df_rc <- bind_rows(df_rc)
df_rc <- df_rc %>% group_by(Sample) %>% mutate(log10.depth.slide=log10(slide_mean(depth, before=200))) %>% ungroup()

df_samples$Sample <- df_samples$Sample
df_rc_bind <- left_join(df_rc, df_samples)
df_rc_bind %>% filter(log10.depth.slide<2, pos%%100==0) %>% .$pos %>% table() %>% sort()

df_rc_bind %>% filter(is.na(Vaccine))

ggplot(df_rc_bind)+
	geom_line(aes(x=pos, y=log10.depth.slide, group = Sample, color=primer), alpha = 0.8, data=. %>% filter(pos%%50==0))+
	geom_hline(yintercept=2, linetype="dashed")+
	facet_wrap(vars(Vaccine), ncol=1)+
	# geom_text_repel(aes(x=pos, y=log10.depth.slide, label = Sample),  data=. %>% filter(pos%in%c(1300, 8000, 13500, 18000, 19400), log10.depth.slide<log10(100)), size=2, min.segment.length=0, bg.color = "white", segment.color="red")+
	scale_color_jco(name="Primer set")+
	scale_x_continuous(breaks = seq(0, 30000, 2500))+
	ylab("Depth (log10, slide window of 200bp)")+
	xlab("Position")+
	NULL
ggsave("../results/coverage.pdf", width = 10, height = 8)

# Consensus sequence
file_consensus_full <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus/", full.names=T)
check <- sapply(file_consensus_full, function(x){
	any(sapply(df_samples$Sample, function(y){grepl(paste0(y, "_c"), x)}))
})
files_consensus_filter <- file_consensus_full[check]
files_consensus_filter <- files_consensus_filter[!grepl("iseq", files_consensus_filter)]
files_consensus_filter <- sort(files_consensus_filter, decreasing = T)

seqs <- lapply(files_consensus_filter, function(x){
	readDNAStringSet(x)
})
seqs <- do.call(c, seqs)

## align to ref
file_seq <- "../results/seqs.fasta"
file_seq_aln <- "../results/seqs_aln.fasta"
writeXStringSet(seqs, file_seq)
system(paste0("mafft --auto --maxiterate 1000 --thread 8 --keeplength --addfragments ", file_seq, " ../results/ref_seq.fasta > ", file_seq_aln))
system(paste0("pangolin ", file_seq_aln, " --outfile ../results/lineage.csv")) 
seq_aln <- readDNAStringSet("../results/seqs_aln.fasta")

# # filter high-coverage
# df_lin <- read_csv("../results/lineage.csv")
# cov_check_raw <- apply(alphabetFrequency(seq_aln), 1, function(x){
# 	sum(x[1:4])
# })
# cov_check <- round(cov_check_raw/29903*100, 2)
# names(cov_check) <- names(seq_aln)
# head(sort(cov_check))
# seq_aln_hq <- seq_aln[cov_check_raw>27000]
# seq_aln_hq <- c(seq_aln_hq)
# samples_hq <- sapply(names(seq_aln_hq), function(x){
# 	tmp <- strsplit(x, "_")[[1]]
# 	tmp[grepl("WHP", tmp)]
# })
# samples_hq <- unlist(samples_hq)
# writeXStringSet(seq_aln_hq, "../results/seqs_aln_hq.fasta")
# writeLines(samples_hq, "../results/samples_hq.txt")

