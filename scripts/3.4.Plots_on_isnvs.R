library(tidyverse)
library(Biostrings)
library(ggsci)
library(patchwork)
library(parallel)
library(lubridate)
library(writexl)
library(ggpubr)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")
source('./helper/plot_box.R')
source("./helper/df_test.R")

colors_lineage=c("#e41a1c", "#33a02c", "#1f78b4", "#ff7f00", "#f781bf", "#666666") 
names(colors_lineage) <- c("Alpha", "Delta", "Omicron", "B.1.36", "B.1.36.27", "B.1.1.63")
colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine)=c("BioNTech", "Sinovac", "Non-vaccinated")
colors_vaccine_new <- colors_vaccine
names(colors_vaccine_new)=c("BioNTech", "Sinovac", "Unvaccinated")
df_orf_sim <- read_csv("../data/ORF_SCoV2_sim.csv")

# load("../results/df_bam_rst.rdata")
df_meta <- read_csv("../results/df_samples_clean.csv", guess_max = 60000)
df_meta$lineage_sim <- factor(df_meta$lineage_sim, levels = names(colors_lineage))
df_meta$detection_lag <- as.numeric(dmy(df_meta$`Report date`) - dmy(df_meta$`Onset date`))

# for iSNVs incidence
load("../results/df_plot_n_gene_adj.rdata")
df_plot_n_gene_meta_adj$lineage_sim <- factor(df_plot_n_gene_meta_adj$lineage_sim, levels = names(colors_lineage))

# for MAF
df_snvs_meta_add_qc <- read_csv("../results/df_snvs_meta_add_qc_bam.csv", guess_max=600000)
# df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% filter(sample %in% df_meta$Sample)
df_snvs_meta_add_qc$lineage_sim <- factor(df_snvs_meta_add_qc$lineage_sim, levels = names(colors_lineage))
df_snvs_meta_add_qc <- bind_rows(df_snvs_meta_add_qc %>% mutate(gene = "Full genome"), df_snvs_meta_add_qc)
write_csv(df_snvs_meta_add_qc, "../results/df_snvs_meta_add_qc_ivar_clean.csv")


# General plots
df_plot_n_gene_meta_adj$Vaccine <- factor(df_plot_n_gene_meta_adj$Vaccine, levels=names(colors_vaccine), labels=names(colors_vaccine_new))
df_snvs_meta_add_qc$Vaccine <- factor(df_snvs_meta_add_qc$Vaccine, levels=names(colors_vaccine), labels=names(colors_vaccine_new))

## plot for incidence
p1_1_ori <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x_var="Vaccine", y_var="n_per_kb", color_var="lineage_sim", y_lab="Number of iSNVs per Kb")
p1_1_ori <- p1_1_ori + scale_color_manual(name="Lineage", values=colors_lineage)
p1_1 <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x_var="Vaccine", y_var="n_per_kb_adj", color_var="lineage_sim", y_lab="Number of iSNVs per Kb (adjusted)")
p1_1 <- p1_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p_out <- p1_1_ori/p1_1 + plot_layout(guides='collect')
ggsave("../results/Num_isnvs_by_vaccine_more.pdf", width=6, height=6)
save_pptx("../results/Num_isnvs_by_vaccine_more.pptx", width=6, height=6)

p1_2 <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x_var="lineage_sim", y_var="n_per_kb_adj", color_var="Vaccine", y_lab="Number of iSNVs per Kb (adjusted)")
p1_2 <- p1_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_new)
ggsave("../results/Num_isnvs_by_vaccine_2_more.pdf", width=6, height=4)
save_pptx("../results/Num_isnvs_by_vaccine_2_more.pptx", width=6, height=4)

p1_3 <- plot_box(df_plot_n_gene_meta_adj %>% filter(gene %in% c("ORF1ab", "S")), x_var="Vaccine", y_var="n_per_kb_adj", color_var="lineage_sim", y_lab="Number of iSNVs per Kb (adjusted)")
p1_3 <- p1_3 + scale_color_manual(name="Lineage", values=colors_lineage)+ facet_grid(rows=vars(gene))
ggsave("../results/Num_isnvs_by_vaccine_gene_more.pdf", width=8, height=6)
save_pptx("../results/Num_isnvs_by_vaccine_gene_more.pptx", width=8, height=6)

p1_4 <- plot_box(df_plot_n_gene_meta_adj %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="n_per_kb_adj", color_var="Vaccine", y_lab="Number of iSNVs per Kb (adjusted)")
p1_4 <- p1_4 + scale_color_manual(name="Vaccine", values=colors_vaccine_new)+facet_grid(rows=vars(gene))
ggsave("../results/Num_isnvs_by_vaccine_gene_2_more.pdf", width=8, height=6)
save_pptx("../results/Num_isnvs_by_vaccine_gene_2_more.pptx", width=8, height=6)

df_wilc_test <- cal_wilc_test(df_plot_n_gene_meta_adj, "n_per_kb_adj", genes = c("Full genome", "ORF1ab", "S"))
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test , "../results/df_test_num_isnvs_by_vaccine_gene_more.xlsx")

ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome") %>% filter(Vaccine!="Non-vaccinated"), x = "days_since_last_dose", y = "n_per_kb_adj",
   color = "Vaccine",
   add = "reg.line",  # Add regressin line
#    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   size = 0.8
   )+
   xlab("Days since last dose")+
   ylab("Number of iSNVs per Kb (adjusted)")+
   facet_wrap(vars(as.character(lineage_sim)),ncol=1)+
   scale_color_manual(name="Vaccine", values=colors_vaccine_new[1:2])+
   scale_fill_manual(name="Vaccine", values=colors_vaccine_new[1:2])+
   stat_cor(aes(color = Vaccine), label.x = 3)+
   NULL
ggsave("../results/cor_last_dose_isnvs.pdf", width=8, height=6)

### non-vaccinated density
p_den <- ggplot(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome")) +
	geom_density(aes(x=n_per_kb_adj, color=lineage_sim, fill=lineage_sim), alpha=0.8)+
	facet_grid(rows=vars(lineage_sim), cols=vars(Vaccine))+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	scale_fill_manual(name="Lineage", values=colors_lineage)+
	xlab("Nubmer of iSNVs per Kb (adjusted)")+
	theme(axis.text.x = element_text(angle = 45, hjust=1))+
	NULL
ggsave("../results/Density_isnvs_by_vaccine_more.pdf", width=8, height=6)
save_pptx("../results/Density_isnvs_by_vaccine_more.pptx", width=8, height=6)

## plot for MAF
p2_1 <- plot_box(df_snvs_meta_add_qc %>% filter(gene=="Full genome"), x_var="Vaccine", y_var="sec_freq", color_var="lineage_sim", y_lab="MAF (adjusted)")
p2_1 <- p2_1 + scale_color_manual(name="Lineage", values=colors_lineage)
ggsave("../results/MAF_more.pdf", width=8, height=6)
save_pptx("../results/MAF_more.pptx", width=8, height=6)

p2_2 <- plot_box(df_snvs_meta_add_qc %>% filter(gene=="Full genome"), x_var="lineage_sim", y_var="sec_freq", color_var="Vaccine", y_lab="MAF (adjusted)")
p2_2 <- p2_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_new)
ggsave("../results/MAF_2_more.pdf", width=8, height=6)
save_pptx("../results/MAF_2_more.pptx", width=8, height=6)

p2_3 <- plot_box(df_snvs_meta_add_qc %>% filter(gene %in% c("ORF1ab", "S")), x_var="Vaccine", y_var="sec_freq", color_var="lineage_sim", y_lab="MAF (adjusted)")
p2_3 <- p2_3 + scale_color_manual(name="Lineage", values=colors_lineage)+facet_grid(rows=vars(gene))
ggsave("../results/MAF_gene_more.pdf", width=8, height=6)
save_pptx("../results/MAF_gene_more.pptx", width=8, height=6)

p2_4 <- plot_box(df_snvs_meta_add_qc %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="sec_freq", color_var="Vaccine", y_lab="MAF (adjusted)")
p2_4 <- p2_4 + scale_color_manual(name="Vaccine", values=colors_vaccine_new)+
	facet_grid(rows=vars(gene))
ggsave("../results/MAF_gene_2_more.pdf", width=8, height=6)
save_pptx("../results/MAF_gene_2_more.pptx", width=8, height=6)

df_wilc_test <- cal_wilc_test(df_snvs_meta_add_qc, "sec_freq")
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test , "../results/df_MAF_by_vaccine_gene_more.xlsx")

# plot for pi
df_plot_pi <- read_csv("../results/codon_results_NOL_bySample.csv")
df_plot_pi$lineage_sim <- factor(df_plot_pi$lineage_sim, levels = names(colors_lineage))
df_plot_pi$Vaccine <- factor(df_plot_pi$Vaccine, levels = names(colors_vaccine), labels=names(colors_vaccine_new))

mean(df_plot_pi$piN, na.rm=T)
mean(df_plot_pi$piS, na.rm=T)
# boxplot(df_plot_pi$piN-df_plot_pi$piS)
wilcox.test(df_plot_pi$piN-df_plot_pi$piS)

p3_1 <- plot_box(df_plot_pi, x_var="Vaccine", y_var="pi", color_var="lineage_sim", y_lab=expression("Nucleotide diversity ("~pi~")"))
p3_1 <- p3_1 + scale_color_manual(name="Lineage", values=colors_lineage)
ggsave("../results/diveristy_pi_more.pdf", width=8, height=6)
save_pptx("../results/diveristy_pi_more.pptx", width=8, height=6)

p3_2 <- plot_box(df_plot_pi, x_var="lineage_sim", y_var="pi", color_var="Vaccine", y_lab=expression("Nucleotide diversity ("~pi~")"))
p3_2 <- p3_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_new)
ggsave("../results/diveristy_pi_2_more.pdf", width=8, height=6)
save_pptx("../results/diveristy_pi_2_more.pptx", width=8, height=6)

# p3_3 <- plot_box(df_plot_pi %>% filter(gene %in% c("ORF1ab", "S")), x_var="Vaccine", y_var="pi", color_var="lineage_sim", y_lab=expression("Nucleotide diversity ("~pi~")"))
# p3_3 <- p3_3 + scale_color_manual(name="Lineage", values=colors_lineage)+
# 	facet_grid(rows=vars(gene)) 
# ggsave("../results/diveristy_pi_gene_more.pdf", width=8, height=6)
# save_pptx("../results/diveristy_pi_gene_more.pptx", width=8, height=6)

# p3_4 <- plot_box(df_plot_pi %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="pi", color_var="Vaccine", y_lab=expression("Nucleotide diversity ("~pi~")"))
# p3_4 <- p3_4 + scale_color_manual(name="Vaccine", values=colors_vaccine_new)+
# 	facet_grid(rows=vars(gene)) 
# ggsave("../results/diveristy_pi_gene_2_more.pdf", width=8, height=6)
# save_pptx("../results/diveristy_pi_gene_2_more.pptx", width=8, height=6)

df_wilc_test <- cal_wilc_test(df_plot_pi %>% mutate(gene="Full genome"), "pi", genes=c("Full genome"))
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test, "../results/df_test_diversity_pi_by_vaccine_gene_more.xlsx")

### combine figures
p_out_1 <- (p1_1+ggtitle("A")) + (p2_1+ggtitle("B")) + (p3_1+ggtitle("C")) + plot_layout(ncol=1, guides="collect") & theme(legend.position='bottom')
p_out_2 <- (p1_2+ggtitle("D")) + (p2_2+ggtitle("E")) + (p3_2+ggtitle("F")) + plot_layout(ncol=1, guides="collect") & theme(legend.position='bottom')
p_out <- p_out_1 | p_out_2
ggsave("../results/plot_comb_by_variant_full_more.pdf", width=10, height=8)
save_pptx("../results/plot_comb_by_variant_full_more.pptx", width=8, height=8)



# difference between genes/positions
(df_tmp <- df_plot_n_gene_meta_adj %>% group_by(gene) %>% summarise(mean_n_per_kb=mean(n_per_kb,na.rm=T), mean_n_per_kb_adj=mean(n_per_kb_adj,na.rm=T)) %>% mutate(gene=factor(gene, levels=c("Full genome", df_orf_sim$sequence))) %>% arrange(gene))
write_csv(df_tmp, "../results/mean_isnvs_per_kb.csv")
df_tmp %>% arrange(mean_n_per_kb_adj)

df_tmp <- df_plot_n_gene_meta_adj %>% filter(gene%in%c("ORF1ab", "S")) %>% select(gene, n_per_kb_adj, Vaccine, lineage_sim, sample)
median(df_tmp$n_per_kb_adj[df_tmp$gene=="ORF1ab"])
median(df_tmp$n_per_kb_adj[df_tmp$gene=="S"])
wilcox.test(df_tmp$n_per_kb_adj[df_tmp$gene=="ORF1ab"], df_tmp$n_per_kb_adj[df_tmp$gene=="S"])

df_tmp2 <- df_tmp %>% pivot_wider(names_from="gene", values_from="n_per_kb_adj")
p_out <- ggscatter(df_tmp2, x = "ORF1ab", y = "S", color = "Vaccine",
   add = "reg.line",  # Add regressin line
   conf.int = TRUE, # Add confidence interval
   size = 0.8
   )+
   scale_color_manual(name="Vaccine", values=colors_vaccine_new)+
   scale_fill_manual(name="Vaccine", values=colors_vaccine_new)+
   stat_cor(aes(color = Vaccine), label.x = 0.1)+
   facet_wrap(vars(as.character(lineage_sim)),ncol=2)+
   NULL
p_out <- p_out+geom_abline(intercept = 0, slope = 1, linetype="dashed")
ggsave("../results/cor_isnvs_orf1ab_s.pdf", width=8, height=12)

p_length_isnv <- ggscatter(df_plot_n_gene_meta_adj %>% filter(gene!="Full genome"), x = "length_full", y = "n_per_kb_adj", color="gene",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   alpha=0.8,
   size=0.8,
   conf.int = TRUE # Add confidence interval
   )+
   # xlab("Ct value")+
   # ylab("Number of iSNVs per Kb")+
   stat_cor(method = "pearson", label.x = 10, label.y = 1.3)+
   NULL
p_length_isnv


df_tmp <- df_snvs_meta_add_qc %>% filter(gene%in%c("ORF1ab", "S")) %>% select(gene, sec_freq, Vaccine, lineage_sim, sample) %>% group_by(gene,sample,Vaccine, lineage_sim) %>% summarise(sec_freq=mean(sec_freq))
median(df_tmp$sec_freq[df_tmp$gene=="ORF1ab"])
median(df_tmp$sec_freq[df_tmp$gene=="S"])
wilcox.test(df_tmp$sec_freq[df_tmp$gene=="ORF1ab"], df_tmp$sec_freq[df_tmp$gene=="S"])

df_tmp2 <- df_tmp %>% pivot_wider(names_from="gene", values_from="sec_freq")
p_out <- ggscatter(df_tmp2, x = "ORF1ab", y = "S", color = "Vaccine",
   add = "reg.line",  # Add regressin line
   conf.int = TRUE, # Add confidence interval
   size = 0.8
   )+
   scale_color_manual(name="Vaccine", values=colors_vaccine_new)+
   scale_fill_manual(name="Vaccine", values=colors_vaccine_new)+
   stat_cor(aes(color = Vaccine), label.x = 0.1)+
   facet_wrap(vars(as.character(lineage_sim)),ncol=2)+
   NULL
p_out2 <- p_out+geom_abline(intercept = 0, slope = 1, linetype="dashed")
ggsave("../results/cor_MAF_orf1ab_s.pdf", width=8, height=12)

