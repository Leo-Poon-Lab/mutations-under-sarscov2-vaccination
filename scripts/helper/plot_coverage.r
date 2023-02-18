# QC of the NGS data
# Analysis the consensus sequences
## conda activate pangolin
library(tidyverse)
library(Biostrings)
library(ggrepel)
library(slider)
library(ggsci)

# QC
plot_depth <- function(samples_cur_batch, file_prefix, sample_names=NA){
	files_rc_filter <- list.files("../results/pysamstats/", full.names = T)
	samples <- sapply(files_rc_filter, function(x){
		tmp <- strsplit(x, "\\/")[[1]]
		tmp <- gsub(".tsv", "", tmp[length(tmp)], fixed=T)
	}, USE.NAMES = F)
	files_rc_filter <- files_rc_filter[samples %in% samples_cur_batch]
	# files_rc_filter <- files_rc_filter[grepl("iseq", files_rc_filter)] 
	files_rc_filter <- sort(files_rc_filter, decreasing = T)

	samples_selected <- sapply(files_rc_filter, function(x){
		tmp <- strsplit(x, "\\/")[[1]]
		tmp <- gsub(".tsv", "", tmp[length(tmp)], fixed=T)
	}, USE.NAMES = F)

	files_rc_filter <- files_rc_filter[!duplicated(files_rc_filter)]

	df_rc_slide <- mclapply(files_rc_filter, function(x){
		tmp <- read_tsv(x)
		slide_mean(tmp$reads_all, before=200)
	}, mc.cores = 16)

	df_rc_slide_m <- as.data.frame(df_rc_slide)
	df_rc_slide_mean <- apply(df_rc_slide_m, 1, mean)

	df_rc_slide_all <- tibble(pos=rep(1:29903, length(df_rc_slide)), depth=unlist(df_rc_slide), sample=rep(samples_selected, each=29903))
	df_rc_bind <- tibble(pos=seq_along(df_rc_slide_mean), depth=df_rc_slide_mean)

	if(is.na(sample_names)){df_rc_slide_all$sample_names=df_rc_slide_all$sample}else{df_rc_slide_all$sample_names=factor(df_rc_slide_all$sample, levels=samples_cur_batch, labels=sample_names)}

	p_out <- ggplot(df_rc_slide_all)+
		geom_line(aes(x=pos, y=log10(depth), group=sample, color=sample), data=. %>% filter(pos%%50==0), alpha=0.28, show.legend = F)+
		# geom_line(aes(x=pos, y=log10(depth)), alpha = 0.8, data=df_rc_bind %>% filter(pos%%50==0), color="dark red")+
		# geom_line(aes(x=pos, y=log10(depth), color=sample), alpha = 0.8, data=df_rc_slide_new_primer %>% filter(pos%%50==0))+
		geom_hline(yintercept=2, linetype="dashed")+
		# facet_wrap(vars(Vaccine), ncol=1)+
		geom_text_repel(aes(x=pos, y=log10(depth), label = sample_names, color=sample),  data=. %>% filter(pos%in%c(12500)), size=2, min.segment.length=0, bg.color = "white", segment.color="red", max.overlaps = 100, show.legend = F)+
		scale_color_viridis_d()+
		scale_x_continuous(breaks = seq(0, 30000, 2000))+
		ylab("Depth (log10, slide window of 200bp)")+
		# ylab("Depth (number of reads, slide window of 200bp)")+
		xlab("Position")+
		NULL
	ggsave(paste0(file_prefix, ".pdf"), width = 12, height = 8, plot=p_out)
	write_csv(df_rc_slide_all %>% group_by(sample) %>% summarise(mean_slide_depth=mean(depth)), paste0(file_prefix, ".csv"))
	# writeLines(samples_selected, "../results/samples_selected.txt")
	p_out
}

