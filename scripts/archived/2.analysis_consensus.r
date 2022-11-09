# QC of the NGS data
# Analysis the consensus sequences
## conda activate pangolin
library(tidyverse)
library(Biostrings)
library(ggrepel)
library(slider)
library(ggsci)

df_samples <- read_csv("../results/df_samples.csv")

# QC and plotting the sequencing depth
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

# # Consensus sequence
# file_consensus_full <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus/", full.names=T)
# tmp <- paste0("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus//", df_samples$Sample, "_consensus.fa" )
# check <- sapply(file_consensus_full, function(x){
# 	x %in% tmp
# })
# files_consensus_filter <- file_consensus_full[check]
# files_consensus_filter <- sort(files_consensus_filter, decreasing = T)

# seqs <- lapply(files_consensus_filter, function(x){
# 	readDNAStringSet(x)
# })
# seqs <- do.call(c, seqs)

# names_new <- gsub("^Consensus_", "", names(seqs))
# names_new <- gsub("_consensus.+", "", names_new)
# names(seqs) <- names_new

# ## align to ref
# file_seq <- "../results/seqs.fasta"
# writeXStringSet(seqs, file_seq)

# # file_seq_aln <- "../results/seqs_aln.fasta"
# # system(paste0("mafft --auto --maxiterate 1000 --thread 8 --keeplength --addfragments ", file_seq, " ../results/ref_seq.fasta > ", file_seq_aln))
# # system(paste0("pangolin ", file_seq_aln, " --outfile ../results/lineage.csv")) 
# # seq_aln <- readDNAStringSet("../results/seqs_aln.fasta")

# system("chmod 755 ./helper/run_nextclade.sh ")
# system("./helper/run_nextclade.sh ../results/seqs.fasta ../results/nextclade/")

# df_samples$group <- paste(df_samples$Vaccine, df_samples$lineage_sim)
# df_samples <- df_samples %>% mutate(group=gsub("Non-", "Un", group))

# # sapply(unique(df_samples$group), function(group_i){
# # 	samples_i <- df_samples$sample[df_samples$group==group_i]
# # 	seqs_i <- seqs[names(seqs) %in% samples_i]
# # 	stopifnot(length(seqs_i)==length(samples_i))
# # 	group_i_no_space <- gsub(" ", "_", group_i)
# # 	file_seq_i <- paste0("../results/seqs_", group_i_no_space, ".fasta")
# # 	writeXStringSet(seqs_i, file_seq_i)
# # 	system(paste0("./helper/run_nextclade.sh ", file_seq_i, " ../results/nextclade_", group_i_no_space, "/"))
# # })

# seqs_new_name <- seqs
# names(seqs_new_name) <- sapply(names(seqs_new_name), function(sample_i){
# 	tmp <- paste(df_samples$group[df_samples$sample==sample_i], sample_i)
# 	gsub(" ", "_", tmp)
# })
# file_seq_rename <- "../results/seqs_rename.fasta"
# writeXStringSet(seqs_new_name[order(names(seqs_new_name))], file_seq_rename)
# system(paste0("./helper/run_nextclade.sh ", file_seq_rename, " ../results/nextclade_rename/"))
