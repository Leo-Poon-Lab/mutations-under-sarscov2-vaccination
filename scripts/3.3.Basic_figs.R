library(tidyverse)
library(ggsci)
library(lubridate)
library(patchwork)
library(ggrepel)
library(writexl)
library(slider)
library(gggenes)
library(ggpubr)
library(boot)
library(parallel)
library(ggbreak)

source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")

colors_lineage=c("#e41a1c", "#33a02c", "#1f78b4", "#ff7f00", "#f781bf", "#666666") 
names(colors_lineage) <- c("Alpha", "Delta", "Omicron", "B.1.36", "B.1.36.27", "B.1.1.63")
colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine)=c("BioNTech", "Sinovac", "Non-vaccinated")

load("../results/df_plot_n_gene.rdata")
load("../results/df_bam_rst_full_notpp.rdata")

df_meta <- read_csv("../results/df_samples_clean.csv", guess_max = 60000)
df_meta$lineage_sim <- factor(df_meta$lineage_sim, levels = names(colors_lineage))
df_meta$detection_lag <- as.numeric(dmy(df_meta$`Report date`) - dmy(df_meta$`Onset date`))

df_plot_n_gene$lineage_sim <- factor(df_plot_n_gene$lineage_sim, levels = names(colors_lineage))

df_snvs_meta_add_qc <- read_csv("../results/df_snvs_meta_add_qc_bam.csv", guess_max=600000)
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(sample %in% df_meta$Sample)
df_snvs_meta_add_qc$lineage_sim <- factor(df_snvs_meta_add_qc$lineage_sim, levels = names(colors_lineage))

df_snvs_meta_add_qc$effect_sim <- df_snvs_meta_add_qc$effect
df_snvs_meta_add_qc$effect_sim[grepl("stream", df_snvs_meta_add_qc$effect_sim)] <- "UTR"
df_snvs_meta_add_qc$effect_sim[grepl("stop", df_snvs_meta_add_qc$effect_sim)] <- "Nonsynonymous"
df_snvs_meta_add_qc$effect_sim[grepl("missense", df_snvs_meta_add_qc$effect_sim)] <- "Nonsynonymous"
df_snvs_meta_add_qc$effect_sim[grepl("lost", df_snvs_meta_add_qc$effect_sim)] <- "Nonsynonymous"
df_snvs_meta_add_qc$effect_sim[grepl("^synonymous", df_snvs_meta_add_qc$effect_sim)] <- "Synonymous"
df_snvs_meta_add_qc$effect_sim[grepl("&synonymous", df_snvs_meta_add_qc$effect_sim)] <- "Synonymous"
df_snvs_meta_add_qc$effect_sim <- gsub("_", " ", df_snvs_meta_add_qc$effect_sim)

# df_snvs_meta_add_qc %>% group_by(effect_sim) %>% summarise(n=n(), prop=n/nrow(df_snvs_meta_add_qc))

# Basic stats on iSNVs
## number of iSNVs per sample
df_n_isnvs <- df_snvs_meta_add_qc %>% group_by(sample) %>% summarise(n=n()) %>% arrange(desc(n))
df_tmp <- df_meta %>% select(sample,Ct_value) %>% left_join(df_n_isnvs)
df_tmp$n[is.na(df_tmp$n)] <- 0
sum(df_tmp$n==0)
(n_isnvs_mean <- mean(df_tmp$n))
df_tmp$Ct_value_break <- cut(df_tmp$Ct_value, 3)

p_n_isnvs <- ggplot(df_tmp) +
	geom_histogram(aes(x=n, fill=Ct_value_break), binwidth = 1, color="black", size=0.3)+
   geom_vline(aes(xintercept=n_isnvs_mean), linetype="dashed")+
	xlab("Number of iSNV sites in a sample")+
	ylab("Number of samples")+
	scale_x_continuous(breaks=seq(0,30,5))+
	scale_fill_manual(name="Ct value range", values=pal_jama()(6)[4:6])+
	theme_minimal()+
	theme(legend.position = "top")+
	NULL

# ## mutation rate of nt pairs
# df_tmp <- left_join(df_snvs_meta_add_qc %>% select(sample, pos, con_base, sec_base), df_bam_rst[1:29903,] %>% select(pos, ref))

# df_tmp <- df_tmp %>% mutate(mutation=paste0(ref,">",sec_base)) %>% select(mutation, sample) %>%  group_by(sample, mutation) %>% summarise(n_sample=n()) %>% ungroup() %>% group_by(mutation) %>% summarise(mean_n=mean(n_sample), sd=sd(n_sample))
# ggplot(df_tmp) +
# 	geom_boxplot(aes(x=mutation, y=mean_n))

## number of recurrent iSNVs 
table(df_snvs_meta_add_qc$effect_sim)
df_tmp <- df_snvs_meta_add_qc %>% group_by(X2, effect_sim) %>% summarise(n=n())
sum(df_tmp$n>1) # 243
sum(df_tmp$n==1) # 1845
sum(df_tmp$n==1)/nrow(df_tmp)

p_n_sharing <- ggplot(df_tmp %>% ungroup() %>% group_by(effect_sim,n) %>% summarise(n_group=n(), n_group_log10=log10(n_group))) +
	# geom_histogram(aes(x=n, fill=effect_sim), binwidth = 1, color="black", size=0.3)+
	geom_col(aes(x=n, y=n_group, fill=effect_sim), color="black", size=0.3)+
	xlab("Number of samples sharing iSNVs")+
	ylab("Number of iSNVs")+
   scale_fill_manual(name="Variant type", values=pal_jama()(3)[c(3,2,1)])+
	theme_minimal()+
	theme(legend.position = "top")+
	scale_x_continuous(breaks=seq(0,60,5))+
   scale_y_break(breaks=c(30, 65), scales=0.3, expand=FALSE, ticklabels=c(65, 70, 75))+
   scale_y_break(breaks=c(75, 95), scales=0.8, expand=FALSE, ticklabels=c(100, 500, 1000, 1500))+
	NULL

## sequencing depth
unique(df_bam_rst$sample)
tmp <- df_bam_rst %>% group_by(sample) %>% summarise(meadian_depth=median(reads_all)) %>% .$meadian_depth
range(tmp)

df_rc_slide <- df_bam_rst %>% group_by(sample) %>% summarise(slide_depth = slide_mean(reads_all, before=200), pos=pos) %>% filter(pos%%50==0)
df_rc_bind_mean <- df_rc_slide %>% group_by(pos) %>% summarise(mean_depth=mean(slide_depth))

p_depth <- ggplot(df_rc_slide)+
	geom_line(aes(x=pos, y=log10(slide_depth), group=sample), data=. %>% filter(pos%%50==0), alpha=0.01)+
	geom_line(aes(x=pos, y=log10(mean_depth)), alpha = 0.8, data=df_rc_bind_mean %>% filter(pos%%50==0), color="dark red")+
	geom_hline(yintercept=2, linetype="dashed")+
	# facet_wrap(vars(Vaccine), ncol=1)+
	scale_x_continuous(breaks = seq(0, 30000, 2000))+
	ylab("Depth (log10, slide window of 200bp)")+
	ylim(0,5)+
	theme_minimal()+
	# ylab("Depth (number of reads, slide window of 200bp)")+
   xlab("Genomic positions")+
	NULL
ggsave("../results/depth.pdf", width=10, height=6, plot=p_depth)

# iSNVs incidence across genome
### Average number of iSNVs per Kb
df_plot <- df_snvs_meta_add_qc
length(unique(df_snvs_meta_add_qc$X2))
df_tmp <- df_snvs_meta_add_qc %>% select(sample, Vaccine, lineage_sim, X2) %>% unique() %>% group_by(Vaccine, lineage_sim) %>% summarise(N_sites=n())
df_meta %>% group_by(Vaccine, lineage_sim) %>% summarise(N_samples=n()) %>% left_join(df_tmp)

df_tmp <- df_snvs_meta_add_qc %>% mutate(site=X2, consensus_base=con_base, secondary_base=sec_base) %>% select(sample, Vaccine, lineage_sim, site, gene, consensus_base, secondary_base, mut_aa) %>% group_by(Vaccine, lineage_sim, site, gene, consensus_base, secondary_base, mut_aa) %>% summarise(N_samples_with_mutation_in_group=n()) 
df_out <- df_meta %>% mutate(group=paste(Vaccine, lineage_sim)) %>% group_by(Vaccine, lineage_sim, group) %>% summarise(N_samples_in_group=n()) %>% left_join(df_tmp %>% select(N_samples_with_mutation_in_group, everything())) %>% ungroup() %>% select(-Vaccine, -lineage_sim)
df_out$mut_aa <- gsub("^p\\.", "", df_out$mut_aa)
write_tsv(df_out, "../results/identified_minor_variants_in_groups.tsv")

df_samples_high_n <- df_snvs_meta_add_qc %>% group_by(sample) %>% summarise(n=n()) %>% filter(n>10)
df_samples_high_n_mut <- df_snvs_meta_add_qc %>% mutate(group=paste(Vaccine, lineage_sim)) %>% filter(sample %in% df_samples_high_n$sample) %>% mutate(site=X2, consensus_base=con_base, secondary_base=sec_base) %>% select(sample, group, site, gene, consensus_base, secondary_base, mut_aa, sec_freq) %>% arrange(sample, site)
unique(df_samples_high_n_mut$sample)
range(df_samples_high_n_mut$sec_freq)
df_samples_high_n_mut$mut_aa <- gsub("^p\\.", "", df_samples_high_n_mut$mut_aa)
write_tsv(df_samples_high_n_mut, "../results/df_samples_high_n_mut.tsv")

df_plot$sample <- df_plot$Sample
breaks_geneome <- seq(0,30000, 1000)
n_breaks <- length(breaks_geneome)
df_plot$pos_cut = cut(df_plot$X2, breaks_geneome, dig.lab = 10)
df_bam_rst$pos_cut = cut(df_bam_rst$pos, breaks_geneome, dig.lab = 10)

df_tmp <- df_plot %>% group_by(sample, Vaccine, lineage_sim, effect_sim, pos_cut) %>% summarise(n=n())
df_bam_rst_depth <- df_bam_rst %>% filter(sample %in% unique(df_tmp$sample)) %>% group_by(sample, pos_cut) %>% filter(reads_all>=100) %>% summarise(n_high_depth=n()) 

df_full <- full_join(df_meta %>% select(sample, Vaccine, lineage_sim) %>% arrange(Vaccine, lineage_sim, sample), df_plot %>% select(pos_cut) %>% unique(), by = character())  ## add samples without iSNVs
df_full$group <- paste(df_full$Vaccine, df_full$lineage_sim)
df_full <- df_full %>% mutate(group=gsub("Non-", "Un", group))
df_full$group <- factor(df_full$group, c("BioNTech Delta","BioNTech Omicron","Sinovac Delta","Sinovac Omicron","Unvaccinated Alpha","Unvaccinated Delta","Unvaccinated Omicron","Unvaccinated B.1.36","Unvaccinated B.1.36.27","Unvaccinated B.1.1.63"))

df_full <- full_join(df_full, tibble(effect_sim=c("Nonsynonymous", "Synonymous")), by = character())  ## add samples without iSNVs
df_plot_join <- left_join(df_full, df_tmp)
df_plot_join$n[is.na(df_plot_join$n)] <- 0

df_plot_join <- left_join(df_plot_join, df_bam_rst_depth, c("sample", "pos_cut"))
df_plot_join$n_per_pos <- df_plot_join$n/df_plot_join$n_high_depth
df_plot_join$n_per_pos[is.na(df_plot_join$n_per_pos)] <- 0
df_plot_join$n_per_kb <- df_plot_join$n_per_pos*1000
write_xlsx(df_plot_join, "../results/iSNVs_incidence_across_genome_more.xlsx")

df_plot_join_n <- df_plot_join %>% select(sample, pos_cut, n_per_kb, Vaccine, lineage_sim, effect_sim) %>% group_by(pos_cut, effect_sim) %>% summarise(n=n(), n_per_kb_mean=mean(n_per_kb))

df_n_s_sites <- read_tsv("../results/snpgenie/BioNTech-Delta/WHP4592/codon_results.txt")
df_n_s_sites <- df_n_s_sites %>% arrange(site)
df_plot_join_n$sites <- sapply(seq_len(nrow(df_plot_join_n)), function(i) {
   pos_stop <- gsub(".+,", "", df_plot_join_n$pos_cut[i])
   pos_stop <- as.numeric(gsub("]", "", pos_stop))
   pos_start <- pos_stop-999
   if(df_plot_join_n$effect_sim[i]=="Nonsynonymous"){
      return(sum(df_n_s_sites$N_sites[(df_n_s_sites$site>=pos_start) & (df_n_s_sites$site<=pos_stop)]))
   } else {
      return(sum(df_n_s_sites$S_sites[(df_n_s_sites$site>=pos_start) & (df_n_s_sites$site<=pos_stop)]))
   }
})
df_plot_join_n$n_per_sites <- df_plot_join_n$n_per_kb_mean*1000 / df_plot_join_n$sites

p1 <- ggplot(df_plot_join_n) +
	geom_col(aes(x=pos_cut, y=n_per_sites, fill=effect_sim), position=position_dodge(preserve="total"), color="black", size=0.3)+
	# facet_wrap(vars(lineage_sim), ncol=1)+
	scale_x_discrete(breaks = levels(df_plot$pos_cut)[seq(1,30,2)], expand = c(0.05,0.05), guide = guide_axis(n.dodge = 2))+
   scale_fill_manual(name="Variant type", values=pal_jama()(3)[c(3,2,1)])+
   theme_minimal()+
	# theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1)) +
	theme(legend.position = "top")+
	xlab("Genomic positions")+
	ylab("Average number of iSNVs per\nsynonymous/non-synonymous site")+
	NULL
p1

df_plot_n_more <- df_plot_join %>% select(sample, pos_cut, n_per_kb, Vaccine, lineage_sim, effect_sim, group) %>% group_by(pos_cut, Vaccine, lineage_sim, effect_sim, group) %>% summarise(n=n(), n_per_kb_mean=mean(n_per_kb))
write_xlsx(df_plot_n_more, "../results/iSNVs_incidence_across_genome_more.xlsx")

p1_more <- ggplot(df_plot_n_more) +
	geom_col(aes(x=pos_cut, y=n_per_kb_mean, fill=effect_sim), position=position_dodge(preserve="total"), color="black", size=0.3)+
	facet_wrap(vars(group), ncol=1)+
	scale_x_discrete(breaks = levels(df_plot$pos_cut)[seq(1,30,2)], expand = c(0.05,0.05), guide = guide_axis(n.dodge = 2))+
   scale_fill_manual(name="Variant type", values=pal_jama()(3)[c(3,2,1)])+
   theme_minimal()+
	# theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1)) +
	theme(legend.position = "top")+
	xlab("Genomic positions")+
	ylab("Average number of iSNVs per Kb")+
	NULL
p1_more
ggsave("../results/iSNVs_incidence_across_genome_more.pdf", height=12, width=8)

### mutation type frequency: e.g. C->U, A->C
N_sites_total <- sum(df_n_s_sites$N_sites)
S_sites_total <- sum(df_n_s_sites$S_sites)

df_full <- full_join(tibble(REF=c("A", "C", "G", "T")), tibble(ALT=c("A", "C", "G", "T")), by=character())
df_full <- df_full[df_full$REF!=df_full$ALT,]
df_full <- full_join(df_meta %>% select(sample), df_full, by=character())
df_full <- full_join(df_full, tibble(effect_sim = c("Nonsynonymous", "Synonymous")), by=character())
df_tmp <- df_snvs_meta_add_qc %>% mutate(REF=X4, ALT=X5) %>% group_by(sample, REF, ALT, effect_sim) %>% summarise(n=n())
df_tmp <- left_join(df_full, df_tmp)
df_tmp$n[is.na(df_tmp$n)] <- 0
df_tmp$mutation <- paste0(df_tmp$REF, ">", df_tmp$ALT)
(df_tmp_summary <- df_tmp %>% group_by(effect_sim) %>% summarise(mean=mean(n), max=max(n), n_samples=length(unique(sample))) %>% arrange(mean))
df_tmp_summary$mean_adj <- sapply(seq_len(nrow(df_tmp_summary)), function(i) {
   num_adjust <- ifelse(df_tmp_summary$effect_sim[i]=="Nonsynonymous", N_sites_total, S_sites_total)
   df_tmp_summary$mean[i]/num_adjust
})

df_tmp <- df_tmp %>% select(-REF, -ALT)

df_boot <- mclapply(unique(df_tmp_summary$mutation), function(mutation_i) {
   # mutation_i <- df_tmp_summary$mutation[1]
   print(mutation_i)
   rsts <- c()
   for (effect_i in c("Nonsynonymous", "Synonymous")) {
      num_adjust <- ifelse(effect_i=="Nonsynonymous", N_sites_total, S_sites_total)
      mut_function <- function(D, indices, num_adjust) {
         df_tmp <- D[indices,]
         mean(df_tmp$n[df_tmp$mutation==mutation_i & df_tmp$effect_sim==effect_i], na.rm = TRUE)/num_adjust
      }
      boot_mean <- boot(data = df_tmp, R = 10000, statistic = mut_function, parallel = 'multicore', ncpus = 2, num_adjust=num_adjust)
      rst <- c(mutation_i, effect_i, boot_mean$t0, sd(boot_mean$t))
      print(effect_i)
      print(rst)
      rsts <- c(rsts, rst)
   }
   return(rsts)
}, mc.cores=8)

df_boot_m <- matrix(unlist(df_boot), ncol=4, byrow=T)
colnames(df_boot_m) <- c("Mutation", "Type", "Mean", "SD")
df_boot_m <- as_tibble(df_boot_m)
df_boot_m <- df_boot_m %>% mutate_at(vars(`Mean`:`SD`), as.numeric)
df_boot_m$Mutation <- gsub("T", "U", df_boot_m$Mutation)
p_mut_type_boot <- ggplot(df_boot_m, aes(x=Mutation)) +
   geom_errorbar(mapping =  aes(ymin = Mean-SD, ymax = Mean+SD, color=Type), position = position_dodge(width = 0.5), width = 0, size = 2) +
    geom_point(mapping =  aes(y = Mean-SD, color=Type ), 
               position = position_dodge(width = 0.5), size = 1.) + # size = 1
    geom_point(mapping =  aes(y = Mean+SD, color=Type), 
               position = position_dodge(width = 0.5), size = 1.2) +
   geom_point(aes(y=Mean, fill=Type), position = position_dodge(width = 0.5), size=3, shape=21, alpha=0.9) +
   geom_hline(aes(yintercept=mean_adj, color=effect_sim), data=df_tmp_summary, linetype="dashed")+
    scale_color_manual(name="Variant type", values=pal_jama()(3)[c(3,2)])+   
    scale_fill_manual(name="Variant type", values=pal_jama()(3)[c(3,2)])+   
   theme_minimal()+
	theme(legend.position = "top")+
   ylab("iSNVs frequency per sample per\nsynonymous/non-synonymous site")+
    NULL
p_mut_type_boot


### Proportion of samples share iSNVs
df_plot <- df_snvs_meta_add_qc
df_plot$mutation <- paste0(df_plot$X4, df_plot$X2, df_plot$X5)
df_plot_n_mut <- df_plot %>% select(effect_sim, sample, gene, X2, mutation, mut_aa,  Vaccine, lineage_sim) %>% group_by(X2, mutation, gene, mut_aa, effect_sim) %>% summarise(n=n()) %>% ungroup()
df_tmp <- df_meta %>% summarise(n_group=n()) %>% ungroup()
df_plot_n_mut$prop <- df_plot_n_mut$n/df_tmp$n_group
df_plot_n_mut$mutation <- paste0(df_plot_n_mut$mutation, " (", gsub("p.", "", df_plot_n_mut$mut_aa, fixed=T), ", N=", df_plot_n_mut$n, ")")

p2 <- ggplot(df_plot_n_mut) +
	geom_point(aes(x=X2, y=prop, color=effect_sim), position="identity")+
	geom_point(aes(x=X2, y=prop), color="black", position="identity", shape=21, stroke=0.2)+
	# facet_wrap(vars(lineage_sim), ncol=1)+
	scale_x_continuous(breaks = breaks_geneome[seq(1,n_breaks,2)])+
	theme_minimal()+
	xlab("Genomic positions")+
   scale_color_manual(name="Variant type", values=pal_jama()(3)[c(3,2,1)])+   
	ylab("Proportion of samples sharing iSNVs")+
	theme(legend.position = "top")+
	geom_text_repel(aes(x=X2, y=prop, label=mutation), data = . %>% filter(prop>0.01))+
	# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	NULL
p2

df_tmp <- df_plot %>% filter(mut_aa %in% (df_plot_n_mut %>% filter(prop>0.01) %>% .$mut_aa)) %>% group_by(mut_aa) %>% mutate(n_mut = n()) %>% group_by(lineage_sim,Vaccine, gene, mut_aa) %>% summarise(prop_mut=n()/n_mut[1]) %>% arrange(mut_aa, desc(prop_mut)) %>% ungroup() %>% mutate(mut_aa=gsub("p.", "", mut_aa, fixed=T))
write_csv(df_tmp, "../results/shared_isnvs_in_different_groups.csv")

### genes, mutation counts and silent mutations.
sars2_ALL_genes <- tibble(
  molecule = "SARS-CoV2",
  gene = factor(c("ORF1a", "ORF1b","S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"),
    levels=c("ORF1a", "ORF1b","S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")),
  start = c(266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558),
  end = c(13483, 21555,25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674),
  strand = "forward",
  direction = 1
)
p_genes <- ggplot(
    sars2_ALL_genes,
    aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)
  ) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "left") +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(name = "Coding region", palette = "Set3") +
  scale_x_continuous(breaks = breaks_geneome[seq(1,n_breaks,2)])+
  ylab("")+
  theme_genes()+
  theme(legend.position = "top")+
  guides(fill=guide_legend(nrow=1))


### combine figures
p_out_top <- (p_n_isnvs+ggtitle("A"))|(print(p_n_sharing+ggtitle("B"))) 
p_out_bottom <- (p_mut_type_boot+ggtitle("C"))/(p1+ggtitle("D"))/(p2+ggtitle("E"))/(p_genes+ggtitle("F")) + plot_layout(heights = c(1, 1, 1, 0.3)) 
# p_out <- ((p_n_isnvs+ggtitle("A"))|(p_n_sharing+ggtitle("B")))/((p_depth+ggtitle("C"))/(p1+ggtitle("D"))/(p2+ggtitle("E"))) + plot_layout(heights = c(1, 3))
p_out <- p_out_top/p_out_bottom + plot_layout(heights = c(1, 3.5)) 
ggsave("../results/Figure 1.pdf",height=15,width=15/sqrt(2))
# ggsave("../results/tmp.pdf",height=15,width=15/sqrt(2))


# Data and quality control
unique(df_meta$Sample)
range(lubridate::dmy(df_meta$`Report date`))
table(df_meta$lineage_sim)
table(df_meta$Vaccine)
table(df_meta$lineage_sim, df_meta$Vaccine)

df_meta_filterd <- df_meta %>% filter(sample %in% unique(df_snvs_meta_add_qc$sample))
table(df_meta_filterd$lineage_sim, df_meta_filterd$Vaccine)

## Correlation to Ct values and detection lag
df_plot_n_gene_meta <- left_join(df_plot_n_gene, df_meta %>% select(Sample, detection_lag, Ct_value, days_since_last_dose) %>% mutate(sample=Sample), "sample") 

add_residuals <- function(df,x_var,y_var,out_var){
	genes = unique(df$gene)
	df_out <- lapply(seq_along(genes), function(i) {
		gene_i <- genes[i]
		df_tmp <- df %>% filter(gene==gene_i)
		lm_i <- lm(paste0(y_var,"~",x_var), data=df_tmp)
		df_tmp[[out_var]] <- lm_i$residuals
		df_tmp
	})
	bind_rows(df_out)
}

df_plot_n_gene_meta_adj <- add_residuals(df_plot_n_gene_meta, x_var = "Ct_value", y_var = "n_per_kb", out_var = "n_per_kb_adj")
df_snvs_meta_add_qc_adj <- add_residuals(bind_rows(df_snvs_meta_add_qc %>% mutate(gene="Full genome"), df_snvs_meta_add_qc), x_var = "Ct_value", y_var = "sec_freq", out_var = "sec_freq_adj")

save(df_plot_n_gene_meta_adj, file="../results/df_plot_n_gene_adj.rdata")

### Ct value
load(file="../results/df_plot_n_gene_adj.rdata")
p_pre <- ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x = "Ct_value", y = "n_per_kb",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   alpha=0.8,
   size=0.8,
   conf.int = TRUE # Add confidence interval
   )+
   xlab("Ct value")+
   ylab("Number of iSNVs per Kb")+
   stat_cor(method = "pearson", label.x = 10, label.y = 1.3)+
   NULL
p_aft <- ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x = "Ct_value", y = "n_per_kb_adj",
   add = "reg.line",  # Add regressin line
   alpha=0.8,
   size=0.8,
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+
   xlab("Ct value")+
   ylab("Number of iSNVs per Kb (adjusted)")+
   stat_cor(method = "pearson", label.x = 10, label.y = 1.3)+
   NULL
p_pre2 <- ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x = "detection_lag", y = "n_per_kb",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   alpha=0.8,
   size=0.8,
   conf.int = TRUE # Add confidence interval
   )+
   xlab("Detection lag (days)")+
   ylab("Number of iSNVs per Kb")+
   stat_cor(method = "pearson", label.x = 0, label.y = 1.3)+
   NULL
ggsave("../results/detection_lag_vs_isnvs.pdf", plot=p_pre2)
p_aft2 <- ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x = "detection_lag", y = "n_per_kb_adj",
   add = "reg.line",  # Add regressin line
   alpha=0.8,
   size=0.8,
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+
   xlab("Detection lag (days)")+
   ylab("Number of iSNVs per Kb (adjusted)")+
   stat_cor(method = "pearson", label.x = 0, label.y = 1.3)+
   NULL


p_out_lag <- ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x = "Ct_value", y = "detection_lag",
   add = "reg.line",  # Add regressin line
   size=0.8,
   alpha=0.8,
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+
   xlab("Ct value")+
   ylab("Detection lag (days)")+
   stat_cor(method = "pearson", label.x = 10, label.y = 20)+
   NULL

#### MAF
p_maf <- ggscatter(df_snvs_meta_add_qc_adj %>% filter(gene=="Full genome"), x = "Ct_value", y = "sec_freq",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   size = 0.8
   )+
   xlab("Ct value")+
   ylab("MAF")+
   stat_cor(method = "pearson", label.x = 0, label.y = 0.6)+
   NULL

p_out <- ((p_pre+ggtitle("A"))/(p_maf+ggtitle("C"))/(p_aft+ggtitle("E")))|((p_pre2+ggtitle("B"))/(p_out_lag+ggtitle("D"))/(p_aft2+ggtitle("F")))
ggsave("../results/n_per_kb_adjusted.pdf", width=10, height=10)

# ggscatter(df_snvs_meta_add_qc_adj %>% filter(gene=="Full genome"), x = "Ct_value", y = "sec_freq_adj",
#    add = "reg.line",  # Add regressin line
#    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#    conf.int = TRUE, # Add confidence interval
#    size = 0.8
#    )+
#    xlab("Ct value")+
#    ylab("MAF (adjusted)")+
#    stat_cor(method = "pearson", label.x = 10, label.y = 0.6)+
#    NULL

# ggscatter(df_plot_pi_meta %>% filter(gene=="Full genome"), x = "Ct_value", y = "piN_piS",
#    add = "reg.line",  # Add regressin line
#    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#    conf.int = TRUE, # Add confidence interval
#    size = 0.8
#    )+
#    xlab("Ct value")+
#    ylab(expression(pi[N]~"/"~pi["S"]))+
#    stat_cor(method = "pearson", label.x = 10, label.y = 3)+
#    NULL

# p_tmp <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x_var="Vaccine", y_var="Ct_value", color_var="lineage_sim", y_lab="Ct value")
# p_tmp <- p_tmp + scale_color_manual(name="Lineage", values=colors_lineage)
# ggsave("../results/Ct_value_more.pdf", width=8, height=10)

# df_wilc_test <- cal_wilc_test(df_meta %>% mutate(gene="Full genome"), "Ct_value", genes="Full genome")
# df_wilc_test <- highlight_diff(df_wilc_test)
# write_xlsx(df_wilc_test , "../results/df_test_Ct_value_more.xlsx")


# ### detection lag
# p1_1_detection_lag <- plot_box(df_plot=df_meta, x_var="Vaccine", y_var="detection_lag", color_var="lineage_sim", y_lab="Detection lag (days)")
# p1_1_detection_lag <- p1_1_detection_lag + scale_color_manual(name="Lineage", values=colors_lineage)
# ggsave("../results/days_detection_lag_more.pdf", width=8, height=6)
# df_wilc_test <- cal_wilc_test(df_meta %>% mutate(gene="Full genome"), "detection_lag", genes="Full genome")
# df_wilc_test <- highlight_diff(df_wilc_test)
# write_xlsx(df_wilc_test , "../results/df_test_detection_lag_more.xlsx")

# ### days since last dose
# p_tmp <- plot_box(df_plot=df_meta %>% filter(Vaccine!="Non-vaccinated"), x_var="Vaccine", y_var="days_since_last_dose", color_var="lineage_sim", y_lab="Days since last dose")
# p_tmp <- p_tmp + scale_color_manual(name="Lineage", values=colors_lineage)
# ggsave("../results/days_since_last_dose_more.pdf", width=8, height=6)
# df_wilc_test <- cal_wilc_test(df_meta %>% filter(Vaccine!="Non-vaccinated") %>% mutate(gene="Full genome"), "days_since_last_dose", genes="Full genome")
# df_wilc_test <- highlight_diff(df_wilc_test)
# write_xlsx(df_wilc_test , "../results/df_test_days_since_last_dose_more.xlsx")

# ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x = "days_since_last_dose", y = "Ct_value",
#    add = "reg.line",  # Add regressin line
#    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#    conf.int = TRUE, # Add confidence interval
#    size = 0.8
#    )+
#    xlab("Days since last dose")+
#    ylab("Ct value")+
#    stat_cor(method = "pearson", label.x = 3, label.y = 1)+
#    NULL

# p_out <- (p_pre+ggtitle("A"))|(p_out_lag+ggtitle("B"))
# ggsave("../results/ct_value_detection_lag.pdf", width=8, height=4)

# ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome") %>% filter(Vaccine!="Non-vaccinated"), x = "days_since_last_dose", y = "n_per_kb_adj",
#    color = "Vaccine",
#    add = "reg.line",  # Add regressin line
# #    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#    conf.int = TRUE, # Add confidence interval
#    size = 0.8
#    )+
#    xlab("Days since last dose")+
#    ylab("Number of iSNVs per Kb (adjusted)")+
#    facet_wrap(vars(as.character(lineage_sim)),ncol=1)+
#    scale_color_manual(name="Vaccine", values=colors_vaccine[1:2])+
#    scale_fill_manual(name="Vaccine", values=colors_vaccine[1:2])+
#    stat_cor(aes(color = Vaccine), label.x = 3)+
#    NULL
# ggsave("../results/cor_last_dose_isnvs.pdf", width=8, height=6)

# df_tmp_0 <- df_snvs_meta_add_qc_adj %>% filter(gene=="Full genome") %>% group_by(sample, gene) %>% summarise(median_sec_freq_adj=median(sec_freq_adj)) %>% ungroup %>% select(sample,median_sec_freq_adj) 
# df_tmp <- left_join(df_tmp_0, df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), "sample")
# ggscatter(df_tmp %>% filter(gene=="Full genome"), x = "median_sec_freq_adj", y = "n_per_kb_adj",
#    color = "Vaccine",
#    add = "reg.line",  # Add regressin line
# #    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#    conf.int = TRUE, # Add confidence interval
#    size = 0.5
#    )+
#    xlab("Median adjusted MAF")+
#    ylab("Number of iSNVs per Kb (adjusted)")+
#    facet_wrap(vars(as.character(lineage_sim)),ncol=1)+
#    scale_color_manual(name="Vaccine", values=colors_vaccine)+
#    scale_fill_manual(name="Vaccine", values=colors_vaccine)+
#    stat_cor(aes(color = Vaccine), label.x = 0.1)+
#    NULL
# ggsave("../results/cor_maf_isnvs.pdf", width=8, height=12)