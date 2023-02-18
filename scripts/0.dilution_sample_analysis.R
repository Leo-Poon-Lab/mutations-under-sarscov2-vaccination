# Analyze the within-host mutations in the serial dilution samples, # exact Ct values for samples in "data/dilution/VOC0013 PCR Novaseq dilution.xlsx"
# VOC0013 E-1, VOC0013 E-2, VOC0013 E-3, VOC0013 E-4, VOC0013 E-5, VOC0013 E-6 1, VOC0013 E-6 2, VOC0013 E-7 1, VOC0013 E-7 2, VOC0013 E-8 1, VOC0013 E-8 2, VOC0013 E-9 1, VOC0013 E-9 2, VOC0013 E-10 1, VOC0013 E-10 2, VOC0013 E-11 1, VOC0013 E-11 2

library(tidyverse)
library(parallel)
library(naturalsort)
# install.packages("plotROC", repos='http://cran.us.r-project.org')
library(plotROC)
library(pROC)
library(scico)
library(ggrepel)
library(patchwork)
library(ggbreak)
source("./helper/pysamstats.r")

df_ct_info <- readxl::read_excel("../data/dilution/VOC0013 PCR Novaseq dilution.xlsx")
df_ct_info <- df_ct_info %>% filter(Sample!="VOC0013-e-4") # e-4 sample failed
df_ct_info$Ct <- round(as.numeric(df_ct_info$Ct),2)
df_ct_info$info <- df_ct_info$Sample
df_ct_info$info[6:18] <- paste0(gsub("VOC0013-", "", df_ct_info$Sample[6:18]), " (Ct:", df_ct_info$Ct[6:18], ")")
df_ct_info$color <- NA
df_ct_info$color[2:10] <- scico(9, palette = 'batlow')
df_ct_info$color[is.na(df_ct_info$color)] <- "#000000"

files_bam <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
files_bam_full <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$", full.names = T)
files_bam_int <- naturalsort(files_bam_full[grepl("VOC0013-e", files_bam_full)])

# process_pysamstats(files_bam_int, n_cores=8)
files_bam_rst_full <- c(list.files("../results/pysamstats/", "tsv$", full.names = T))
files_bam_rst_full <- naturalsort(files_bam_rst_full[grepl("VOC0013-e", files_bam_rst_full)])
stopifnot(sum(is.na(files_bam_rst_full))==0)

read_pysamstats_dilution <- function(x, pp=TRUE) {
	print(x)
	if(pp){
		df_bam_rst <- read_tsv(x) %>% select(pos, contains("_pp"))
	} else {
		df_bam_rst <- read_tsv(x) %>% select(pos, !contains("_pp"))
	}
	sample_x <- strsplit(x, "/", fixed=T)[[1]]
	sample_x <- sample_x[length(sample_x)]
	sample_x <- gsub(".tsv", "", sample_x, fixed=T)
	df_bam_rst$sample <- sample_x
	# if(sum(df_bam_rst$reads_all>=100)>=27000){return(df_bam_rst)} else{return(NA)}
	return(df_bam_rst)
}

df_bam_rst <- mclapply(files_bam_rst_full, read_pysamstats_dilution, pp=FALSE, mc.cores=16)
stopifnot(all(sapply(df_bam_rst, nrow)==29903))
df_bam_rst <- bind_rows(df_bam_rst[!is.na(df_bam_rst)])
samples_all <- unique(df_bam_rst$sample)
samples_ref <- "VOC0013-e-1"

# df_bam_rst$sample <- factor(df_bam_rst$sample, samples_all)
# levels(df_bam_rst$sample)

# list difference between original samples and dilution samples
df_extra_info <- extra_info_from_pysamstats(df_bam_rst)
df_extra_info <- df_extra_info %>% mutate_at(vars(!contains(c("base", "sample"))), as.numeric)
df_extra_info <- full_join(df_bam_rst %>% mutate(depth=reads_all) %>% select(pos,sample,depth), df_extra_info %>% select(-depth), c("pos", "sample"))

# df_extra_info <- left_join(df_extra_info, df_bam_rst, c("pos", "sample"))

# df_bam_rst %>% filter(pos==18181) %>% select(A, C, G, T, N, insertions, deletions, sample) # checking
# df_extra_info %>% filter(pos==18181) %>% select(con_base, con_freq, sec_base, sec_freq, sample) %>% print(n=20) # checking
# df_bam_rst %>% filter(pos==24091) %>% select(A, C, G, T, N, insertions, deletions, sample) # checking
# df_extra_info %>% filter(pos==24091) %>% select(con_base, con_freq, sec_base, sec_freq, sample) %>% print(n=20) # checking
# df_bam_rst %>% filter(pos==25381) %>% select(A, C, G, T, N, insertions, deletions, sample) # checking
# df_extra_info %>% filter(pos==25381) %>% select(con_base, con_freq, sec_base, sec_freq, sample) %>% print(n=20) # checking

# sequencing depth
source("./helper/plot_coverage.r")
# plot_depth(samples_all, "../results/dilution_depth", sample_names=NA)
p_depth <- plot_depth(df_ct_info$Sample[2:17], "../results/dilution_depth_Ct", sample_names=df_ct_info$info[2:17])

## mismatch of consensus_base
df_extra_info_ref <- df_extra_info %>% filter(sample == samples_ref)
df_extra_info_others <- df_extra_info %>% filter(sample != samples_ref)
samples_others <- unique(df_extra_info_others$sample)
df_con_base_compare <- lapply(samples_others, function(sample_i){
	print(sample_i)
	df_extra_info_others_i <- df_extra_info_others %>% filter(sample == sample_i)
	df_compare_i <- left_join(df_extra_info_ref %>% select(pos, con_base), df_extra_info_others_i %>% select(pos, con_base), "pos")
	check_i <- which(df_compare_i$con_base.x!=df_compare_i$con_base.y)
	diff_con_base_i <- length(check_i)
	if(diff_con_base_i>0){
		return(tibble(num_diff=diff_con_base_i, pos=paste(df_compare_i$pos[check_i], collapse = ","), sample=sample_i))
	} else {
		return(tibble(num_diff=0, pos=NA, sample=sample_i))
	}
})
df_con_base_compare <- bind_rows(df_con_base_compare)
df_con_base_diff_details <- lapply(seq_len(nrow(df_con_base_compare)), function(i){
	pos_i <- df_con_base_compare$pos[i]
	if(is.na(pos_i)){return(NA)}
	pos_i <- as.numeric(unlist(strsplit(pos_i, ",")))
	sample_i <- df_con_base_compare$sample[i]
	df_extra_info %>% filter(sample==sample_i, pos%in%pos_i)
})
df_con_base_diff_details <- bind_rows(df_con_base_diff_details[!is.na(df_con_base_diff_details)])
quantile(df_con_base_diff_details$depth, seq(1,100)/100) # All the consensus_base errors occur in positions with depth <= 4

## mismatch of minor allele
## Our aim is to determine the thresholds for excluding fasle positive iSNVs.
### First to construct a set of true-positive minor variants. (defined as minor variants that are found in at least two of the three highest-viral-load samples)
### Then plot the RCO curve, to study effects of CT values, and to identify optimal thresholds
# basic filter
df_indels_true_positive_major <- df_extra_info %>% filter(depth>=10) %>% filter(sample %in% c("VOC0013-e-1", "VOC0013-e-2", "VOC0013-e-3")) %>% group_by(pos, con_base) %>% summarise(n=n()) %>% filter(n>1) %>% filter(con_base %in% c("deletions", "insertions")) %>% print(n=50)
df_indels_true_positive_minor <- df_extra_info %>% filter(sec_freq>=0.01) %>% filter(depth>=10) %>% filter(sample %in% c("VOC0013-e-1", "VOC0013-e-2", "VOC0013-e-3")) %>% group_by(pos, sec_base) %>% summarise(n=n()) %>% filter(n>1) %>% filter(sec_base %in% c("deletions", "insertions")) %>% print(n=50)

source("./helper/isnv_position_filter.R")
df_extra_info$primer <- "new"
pos_indels <- c(df_indels_true_positive_major$pos, df_indels_true_positive_minor$pos)
if(length(pos_indels)>0){df_extra_info <- filter_by_pos(df_extra_info, pos_indels=pos_indels)}

# df_extra_info %>% filter(sec_freq>=0.01) %>% filter(depth>=10) %>% group_by(sample) %>% summarise(n=n())
# df_extra_info %>% filter(sec_freq>=0.02) %>% filter(depth>=10) %>% group_by(sample) %>% summarise(n=n())
# df_extra_info %>% filter(sec_freq>=0.03) %>% filter(depth>=10) %>% group_by(sample) %>% summarise(n=n())
# df_extra_info %>% filter(sec_freq>=0.01) %>% filter(depth>=100) %>% group_by(sample) %>% summarise(n=n())
# df_extra_info %>% filter(sec_freq>=0.05) %>% filter(depth>=100) %>% group_by(sample) %>% summarise(n=n())
# df_extra_info %>% filter(sec_freq>=0.01) %>% filter(depth>=10) %>% filter(sample %in% c("VOC0013-e-1", "VOC0013-e-2", "VOC0013-e-3")) %>% group_by(pos, sec_base) %>% summarise(n=n()) %>% filter(n>1) %>% print(n=50)


## Combination test of depth cutoff and allele-freq cutoff
## plot sensitivity test for thresholds
df_isnvs_true_positive <- df_extra_info %>% filter(sec_freq>=0.01) %>% filter(depth>=10) %>% filter(sample %in% c("VOC0013-e-1", "VOC0013-e-2","VOC0013-e-3")) %>% group_by(pos, sec_base) %>% summarise(n=n()) %>% filter(n>1)  %>% filter (!sec_base %in% c("deletions", "insertions")) %>% print(n=50)

df_roc <- full_join(df_extra_info %>% filter(sec_freq>=0.01) %>% select("sample", "pos", "sec_base", "sec_freq", "depth"), df_isnvs_true_positive, c("pos", "sec_base"))
df_roc$tp <- !is.na(df_roc$n)
df_roc <- left_join(df_roc %>% mutate(Sample=sample), df_ct_info %>% mutate(info=gsub("VOC0013-", "", info)), "Sample")

p1 <- ggplot(df_roc %>% filter(sample %in% samples_all[1:8]), aes(d=tp, m=depth, color=info, fill=info))+
	geom_roc(n.cuts = 5, labelsize = 2, alpha=0.8) + 
	geom_rocci(ci.at = c(100), labels=F, linetype=1, size=0.5, show.legend = FALSE)+
	# geom_rocci(ci.at = c(10), labels=F, linetype=2, size=0.5, show.legend = FALSE)+
	style_roc()+
	scale_color_manual(name="Sample", values=df_ct_info$color[df_ct_info$Sample %in% samples_all[1:9]])+
	scale_fill_manual(name="Sample", values=df_ct_info$color[df_ct_info$Sample %in% samples_all[1:9]])+
	facet_wrap(vars(info))+
	# scale_x_continuous(guide = guide_axis(n.dodge=1))+
	NULL
ggsave("../results/Dilution_threshold_depth.pdf", width=10, height=10)

df_roc_depth_100 <- full_join(df_extra_info %>% filter(depth>=100, sec_freq>=0.01) %>% select("sample", "pos", "sec_base", "sec_freq", "depth"), full_join(tibble(sample=samples_all), df_isnvs_true_positive, by = character()), c("pos", "sec_base", "sample"))
df_roc_depth_100$sec_freq[is.na(df_roc_depth_100$sec_freq)] <- 0
df_roc_depth_100$depth[is.na(df_roc_depth_100$depth)] <- 0
df_roc_depth_100$tp <- !is.na(df_roc_depth_100$n)
df_roc_depth_100 <- left_join(df_roc_depth_100 %>% mutate(Sample=sample), df_ct_info %>% mutate(info=gsub("VOC0013-", "", info)), "Sample")
df_roc_depth_100 %>% filter(sample=="VOC0013-e-4") %>% filter(tp)

p2 <- ggplot(df_roc_depth_100 %>% filter(sample %in% samples_all[1:8]), aes(d=tp, m=sec_freq))+
	geom_roc(cutoffs.at = 1:5/100, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_rocci(ci.at = c(0.025), labels=F, linetype=1, size=0.5, show.legend = FALSE, fill="dark red")+
	# geom_rocci(ci.at = c(0.02), labels=F, linetype=2, size=0.5, show.legend = FALSE, fill="dark blue")+
	style_roc()+
	# scale_color_manual(name="Sample", values=df_ct_info$color[df_ct_info$Sample %in% samples_all[1:9]])+
	# scale_fill_manual(name="Sample", values=df_ct_info$color[df_ct_info$Sample %in% samples_all[1:9]])+
	facet_wrap(vars(info))+
	# scale_x_continuous(guide = guide_axis(n.dodge=1))+
	NULL
ggsave("../results/Dilution_threshold_maf.pdf", width=10, height=10)

df_auc <- calc_auc(p2) %>% print()
df_thred_optimal <- lapply(samples_all[1:8], function(sample_i){
	df_tmp <- df_roc_depth_100 %>% filter(sample %in% sample_i)
	my_curve <- roc(df_tmp$tp, df_tmp$sec_freq, direction="<")
	df_thred <- coords(my_curve, x=1:20/200)
	df_thred$sample=sample_i
	df_thred
})
df_thred_optimal <- bind_rows(df_thred_optimal)

df_thred_optimal <- df_thred_optimal %>% mutate(diff=sensitivity-(1-specificity)) %>% arrange(desc(diff)) %>% group_by(sample) %>% filter(row_number()==1) %>% arrange(sample)
df_thred_optimal <- left_join(df_thred_optimal, df_ct_info %>% mutate(sample=Sample) %>% select(-Sample, -color))
df_thred_optimal$info <- gsub("VOC0013-", "", df_thred_optimal$info)
df_thred_optimal <- df_thred_optimal %>% select(-sample, -Ct)
df_thred_optimal <- left_join(df_thred_optimal, df_auc %>% select(-PANEL, -group))
writexl::write_xlsx(df_thred_optimal, "../results/df_thred_optimal.xlsx")

depth_cuts <- c(10, 1:5*20)
df_roc_thresholds <- lapply(depth_cuts, function(depth_cut_i){
	df_tmp <- full_join(df_extra_info %>% filter(depth>=depth_cut_i, sec_freq>=0.01) %>% select("sample", "pos", "sec_base", "sec_freq", "depth"), full_join(tibble(sample=samples_all), df_isnvs_true_positive, by = character()), c("pos", "sec_base", "sample"))
	df_tmp$depth_cut <- depth_cut_i
	df_tmp$sec_freq[is.na(df_tmp$sec_freq)] <- 0
	df_tmp$depth[is.na(df_tmp$depth)] <- 0
	df_tmp$tp <- !is.na(df_tmp$n)
	df_tmp <- left_join(df_tmp %>% mutate(Sample=sample), df_ct_info %>% mutate(info=gsub("VOC0013-", "", info)), "Sample")
	df_tmp
})
df_roc_thresholds <- bind_rows(df_roc_thresholds)
df_roc_thresholds$depth_cut <- factor(df_roc_thresholds$depth_cut, depth_cuts)

p3_001 <- ggplot(df_roc_thresholds %>% filter(sample %in% samples_all[1:8]), aes(d=tp, m=sec_freq, color=depth_cut, fill=depth_cut))+
	# geom_roc(cutoffs.at = 1:5/100, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_roc(n.cuts = 0, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_rocci(ci.at = c(0.01), labels=F, linetype=2, size=0.5, show.legend = FALSE)+
	# geom_rocci(ci.at = c(0.05), labels=F, linetype=1, size=0.5, show.legend = FALSE)+
	style_roc()+
	scale_color_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	scale_fill_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	facet_wrap(vars(info))+
	ggtitle(paste0("A (MAF cutoff:", 0.01, ")"))+
	# scale_x_continuous(guide = guide_axis(n.dodge=1))+
	NULL
p3_002 <- ggplot(df_roc_thresholds %>% filter(sample %in% samples_all[1:8]), aes(d=tp, m=sec_freq, color=depth_cut, fill=depth_cut))+
	# geom_roc(cutoffs.at = 1:5/100, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_roc(n.cuts = 0, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_rocci(ci.at = c(0.02), labels=F, linetype=2, size=0.5, show.legend = FALSE)+
	# geom_rocci(ci.at = c(0.05), labels=F, linetype=1, size=0.5, show.legend = FALSE)+
	style_roc()+
	scale_color_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	scale_fill_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	facet_wrap(vars(info))+
	ggtitle(paste0("B (MAF cutoff:", 0.02, ")"))+
	# scale_x_continuous(guide = guide_axis(n.dodge=1))+
	NULL
p3_003 <- ggplot(df_roc_thresholds %>% filter(sample %in% samples_all[1:8]), aes(d=tp, m=sec_freq, color=depth_cut, fill=depth_cut))+
	# geom_roc(cutoffs.at = 1:5/100, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_roc(n.cuts = 0, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_rocci(ci.at = c(0.03), labels=F, linetype=2, size=0.5, show.legend = FALSE)+
	# geom_rocci(ci.at = c(0.05), labels=F, linetype=1, size=0.5, show.legend = FALSE)+
	style_roc()+
	scale_color_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	scale_fill_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	facet_wrap(vars(info))+
	ggtitle(paste0("C (MAF cutoff:", 0.03, ")"))+
	# scale_x_continuous(guide = guide_axis(n.dodge=1))+
	NULL
p3_004 <- ggplot(df_roc_thresholds %>% filter(sample %in% samples_all[1:8]), aes(d=tp, m=sec_freq, color=depth_cut, fill=depth_cut))+
	# geom_roc(cutoffs.at = 1:5/100, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_roc(n.cuts = 0, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_rocci(ci.at = c(0.04), labels=F, linetype=2, size=0.5, show.legend = FALSE)+
	# geom_rocci(ci.at = c(0.05), labels=F, linetype=1, size=0.5, show.legend = FALSE)+
	style_roc()+
	scale_color_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	scale_fill_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	facet_wrap(vars(info))+
	ggtitle(paste0("D (MAF cutoff:", 0.04, ")"))+
	# scale_x_continuous(guide = guide_axis(n.dodge=1))+
	NULL
p3_005 <- ggplot(df_roc_thresholds %>% filter(sample %in% samples_all[1:8]), aes(d=tp, m=sec_freq, color=depth_cut, fill=depth_cut))+
	# geom_roc(cutoffs.at = 1:5/100, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_roc(n.cuts = 0, labelsize = 2, labelround = 2, alpha=0.88) + 
	geom_rocci(ci.at = c(0.05), labels=F, linetype=2, size=0.5, show.legend = FALSE)+
	# geom_rocci(ci.at = c(0.05), labels=F, linetype=1, size=0.5, show.legend = FALSE)+
	style_roc()+
	scale_color_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	scale_fill_manual(name="Depth threshold", values=rev(scico(length(depth_cuts), palette = 'batlow', end=0.8)))+
	facet_wrap(vars(info))+
	# scale_x_continuous(guide = guide_axis(n.dodge=1))+
	ggtitle(paste0("D (MAF cutoff:", 0.05, ")"))+
	NULL
p3_out <- p3_001+p3_002+p3_003+p3_004+plot_layout(ncol=2, guides='collect') & theme(legend.position='bottom') & guides(colour = guide_legend(nrow = 1, override.aes = list(size=10)))
ggsave("../results/Dilution_threshold_maf_2.pdf", width=18, height=18)

# genome coverage
df_coverage <- df_bam_rst %>% group_by(sample) %>% summarise(genome_cov=sum(reads_all>=10), average_depth=mean(reads_all)) %>% ungroup()
df_coverage <- left_join(df_coverage %>% mutate(Sample=sample), df_ct_info %>% mutate(info=gsub("VOC0013-", "", info)), "Sample")
df_coverage$info <- naturalfactor(df_coverage$info)
p4_1 <- ggplot(df_coverage, aes(x=log10(average_depth), y=genome_cov, label=info, fill=color)) + 
	geom_point(size=3, shape=21, color="black")+
	xlab("Average depth (Log10)")+
	ylab("Genome coverage")+
	geom_text_repel(bg.color = "white", bg.r=0.1, alpha=0.8, size=2.5, min.segment.length=0.1, max.overlaps=1000)+
	scale_fill_identity()+
	scale_y_continuous(breaks=0:6*5000)+
	theme_minimal()+
	NULL

# Plot MAF fluctuation between samples
df_num_isnvs <- df_roc_thresholds %>% filter(sample %in% samples_all[1:9]) %>% group_by(sample)%>% mutate(n_before_filter=sum(depth_cut==10 & sec_freq>=0.01)) %>% filter(depth_cut==100, sec_freq>=0.025) %>% summarise(n_before_filter=n_before_filter[1], n=n(), n_tp=sum(tp))
df_num_isnvs <- df_num_isnvs %>% pivot_longer(n_before_filter:n_tp)
df_num_isnvs$name[df_num_isnvs$name=="n_before_filter"] <- "Total (Depth>=10, MAF>=0.01)"
df_num_isnvs$name[df_num_isnvs$name=="n"] <- "Total (Depth>=100, MAF>=0.025)"
df_num_isnvs$name[df_num_isnvs$name=="n_tp"] <- "True positives"
df_num_isnvs <- left_join(df_num_isnvs, df_ct_info %>% mutate(sample=Sample) %>% select(-Sample, -color))
df_num_isnvs$info <- gsub("VOC0013-", "", df_num_isnvs$info)

p4_2 <- ggplot(df_num_isnvs)+
	geom_col(aes(x=info, y= value, fill=name), position="dodge", color="black")+
	scale_fill_manual(name="iSNVs", values=c("Dark red", "Dark blue", "Dark green"))+
	scale_y_break(breaks = c(20, 21), scales=0.5, expand=F)+
	scale_y_break(breaks = c(80, 100), scales=0.5, expand=F)+
	# scale_x_discrete(guide = guide_axis(n.dodge=2))+
	xlab("Sample")+
	ylab("Number")+
	coord_flip()+
	theme_minimal()+
	theme(legend.position='bottom')+
	NULL
p4_btm <- ((p4_1+ggtitle("B"))/(p4_2+ggtitle("C"))+plot_layout(heights = c(1, 1.2)))
p4 <- (p_depth+theme_minimal()+ggtitle("A"))/p4_btm+plot_layout(heights = c(0.6, 1))
ggsave("../results/Dilution_genome_coverage.pdf", height=10, width=8)



