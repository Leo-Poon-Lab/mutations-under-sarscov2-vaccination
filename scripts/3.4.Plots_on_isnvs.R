library(tidyverse)
library(Biostrings)
library(ggsci)
library(patchwork)
library(parallel)
library(lubridate)
library(writexl)
library(ggpubr)
library(ggsignif)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")
source('./helper/plot_box.R')
source("./helper/df_test.R")

df_meta <- read_csv("../results/df_samples.csv", guess_max=100000)
df_meta <- df_meta %>% filter(lineage_sim != "22B (Omicron, BA.5.*)")
lineages_all <- sort(unique(df_meta$lineage_sim))
colors_lineage=rev(c("#4daf4a", "#984ea3", "#ff7f00", "#f781bf", "#666666"))
names(colors_lineage) <- lineages_all
vaccine_all <- sort(unique(df_meta$Vaccine))
colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine)=vaccine_all
colors_vaccine_doses = c("#a65628", "#633318", "#7570b3", "#46436b", "#999999")
names(colors_vaccine_doses)=c("Comirnaty\n(Doses=2)", "Comirnaty\n(Doses=3)", "CoronaVac\n(Doses=2)", "CoronaVac\n(Doses=3)", "Unvaccinated")

df_meta$lineage_sim <- factor(df_meta$lineage_sim, levels = names(colors_lineage))
# collection lag
df_meta$detection_lag <- as.numeric(df_meta$`Report date` - ymd(df_meta$`Onset date`))
df_meta$collection_lag <- as.numeric(df_meta$collection_date - ymd(df_meta$`Onset date`))
df_meta$collection_lag[df_meta$collection_lag>100] <- NA # one ourlier
df_meta$collection_lag[df_meta$collection_lag<0] <- NA # one ourlier

colors_vaccine_new <- colors_vaccine
names(colors_vaccine_new)=c("Comirnaty", "CoronaVac", "Unvaccinated")
df_orf_sim <- read_csv("../data/ORF_SCoV2_sim.csv")

# for iSNVs incidence
load("../results/df_plot_n_gene_adj.rdata")
# unique(df_plot_n_gene_meta_adj$lineage_sim)

df_plot_n_gene_meta_adj$lineage_sim <- factor(df_plot_n_gene_meta_adj$lineage_sim, levels = names(colors_lineage))
df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, df_meta %>% select(sample, Doses), "sample")
df_plot_n_gene_meta_adj$Vaccine <- factor(df_plot_n_gene_meta_adj$Vaccine, levels=names(colors_vaccine), labels=names(colors_vaccine_new))
df_plot_n_gene_meta_adj$vaccine_doses <- paste0(df_plot_n_gene_meta_adj$Vaccine, "\n(Doses=", df_plot_n_gene_meta_adj$Doses, ")")
df_plot_n_gene_meta_adj$vaccine_doses[df_plot_n_gene_meta_adj$Vaccine=="Unvaccinated"] <- "Unvaccinated"
df_plot_n_gene_meta_adj$vaccine_doses <- factor(df_plot_n_gene_meta_adj$vaccine_doses, levels = names(colors_vaccine_doses))

# for MAF
df_snvs_meta_add_qc <- read_csv("../results/df_snvs_meta_add_qc_bam_adj.csv", guess_max=600000)
df_snvs_meta_add_qc$lineage_sim <- factor(df_snvs_meta_add_qc$lineage_sim, levels = names(colors_lineage))
df_snvs_meta_add_qc <- bind_rows(df_snvs_meta_add_qc %>% mutate(gene = "Full genome"), df_snvs_meta_add_qc)
write_csv(df_snvs_meta_add_qc, "../results/df_snvs_meta_add_qc_ivar_clean.csv")
df_snvs_meta_add_qc$Vaccine <- factor(df_snvs_meta_add_qc$Vaccine, levels=names(colors_vaccine), labels=names(colors_vaccine_new))
df_snvs_meta_add_qc$vaccine_doses <- paste0(df_snvs_meta_add_qc$Vaccine, "\n(Doses=", df_snvs_meta_add_qc$Doses, ")")
df_snvs_meta_add_qc$vaccine_doses[df_snvs_meta_add_qc$Vaccine=="Unvaccinated"] <- "Unvaccinated"
df_snvs_meta_add_qc$vaccine_doses <- factor(df_snvs_meta_add_qc$vaccine_doses, levels = names(colors_vaccine_doses))

# for pi
df_plot_pi <- read_csv("../results/codon_results_NOL_bySample.csv")
df_plot_pi$lineage_sim <- factor(df_plot_pi$lineage_sim, levels = names(colors_lineage))
df_plot_pi$Vaccine <- factor(df_plot_pi$Vaccine, levels = names(colors_vaccine), labels=names(colors_vaccine_new))

df_plot_pi$vaccine_doses <- paste0(df_plot_pi$Vaccine, "\n(Doses=", df_plot_pi$Doses, ")")
df_plot_pi$vaccine_doses[df_plot_pi$Vaccine=="Unvaccinated"] <- "Unvaccinated"
df_plot_pi$vaccine_doses <- factor(df_plot_pi$vaccine_doses, levels = names(colors_vaccine_doses))

# Plot adjusted values vs original values
source("./helper/plot_reg.R")
p_pre <- plot_reg(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), "Ct_value", "n_per_kb", x_lab="Ct value", y_lab="Number of iSNVs per Kb")
p_aft <- plot_reg(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), "Ct_value", "n_per_kb_adj", x_lab="Ct value", y_lab="Number of iSNVs per Kb (adjusted)")
p_pre2 <- plot_reg(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), "collection_lag", "n_per_kb", x_lab="Collection lag (days)", y_lab="Number of iSNVs per Kb")
p_aft2 <- plot_reg(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), "collection_lag", "n_per_kb_adj", x_lab="Collection lag (days)", y_lab="Number of iSNVs per Kb (adjusted)")
ggsave("../results/collection_lag_vs_isnvs_adj.pdf", plot=p_aft2)

p_out_lag <- plot_reg(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), "Ct_value", "collection_lag", x_lab="Ct value", y_lab="Collection lag (days)")
#### MAF
p_maf <- plot_reg(df_snvs_meta_add_qc %>% filter(gene=="Full genome"), "Ct_value", "sec_freq", x_lab="Ct value", y_lab="MAF") # R value very small: -0.054, therefore no need to use the adjusted MAF
p_maf_adj <- plot_reg(df_snvs_meta_add_qc %>% filter(gene=="Full genome"), "Ct_value", "sec_freq_adj", x_lab="Ct value", y_lab="MAF (adjusted)")
#### Pi
p_pi <- plot_reg(df_plot_pi, "Ct_value", "pi_ori", x_lab="Ct value", y_lab=expression("Nucleotide diversity ("~pi~")"))
p_pi_adj <- plot_reg(df_plot_pi, "Ct_value", "pi", x_lab="Ct value", y_lab=expression("Nucleotide diversity ("~pi~") (adjusted)")) 

# p_out <- ((p_pre+ggtitle("A"))|(p_aft+ggtitle("B")))/(p_pre2+ggtitle("C")|(p_aft2+ggtitle("D")))/((p_maf+ggtitle("E"))|p_maf_adj+ggtitle("F"))/((p_pi+ggtitle("G"))|p_pi_adj+ggtitle("H"))
p_out <- ((p_pre+ggtitle("A"))|(p_aft+ggtitle("B")))/(p_pre2+ggtitle("C")|(p_aft2+ggtitle("D")))/((p_pi+ggtitle("E"))|p_pi_adj+ggtitle("F"))
ggsave("../results/n_per_kb_adjusted.pdf", width=10, height=10)


# test significance
df_wilc_test <- cal_wilc_test(df_plot_n_gene_meta_adj, "n_per_kb_adj", genes = c("Full genome", "ORF1ab", "S"))
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test , "../results/df_test_num_isnvs_by_vaccine_gene_more.xlsx")

df_wilc_test <- cal_wilc_test(df_snvs_meta_add_qc, "sec_freq")
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test , "../results/df_MAF_by_vaccine_gene_more.xlsx")

df_wilc_test <- cal_wilc_test(df_plot_pi %>% mutate(gene="Full genome"), "pi", genes=c("Full genome"))
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test, "../results/df_test_diversity_pi_by_vaccine_gene_more.xlsx")

# General plots
## plot for incidence
df_wilc_test <- cal_wilc_test(df_plot_n_gene_meta_adj, "n_per_kb", genes = c("Full genome", "ORF1ab", "S"))
df_wilc_test <- highlight_diff(df_wilc_test)
df_wilc_test %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()
p1_1_ori <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x_var="vaccine_doses", y_var="n_per_kb", color_var="lineage_sim", y_lab="Number of iSNVs per Kb", x_lab="Vaccine")
p1_1_ori <- p1_1_ori + scale_color_manual(name="Lineage", values=colors_lineage)
p1_1_ori_p <- p1_1_ori + # comirnaty 2-doses
   # ggtitle("A")+
   geom_signif(y_position = 2.5+seq(0,0)/5, xmin = c(0.8), xmax = c(1.2), annotation = c("*"), color="black", vjust=0.65, tip_length = 0.01)
p1_1_ori_p <- p1_1_ori_p+ # unvaccinated
   geom_signif(y_position = 2.5+seq(0,6)/5, xmin = c(4.7,4.7,4.7,4.85,4.85,4.85,4.85), xmax = c(5,5.15,5.3,4.7,5,5.15,5.3), annotation = c("**","**","**","**","**","**","**"), color="black", vjust=0.65, tip_length = 0.01)

df_test_isnvs <- readxl::read_excel("../results/df_test_num_isnvs_by_vaccine_gene_more.xlsx")
df_test_isnvs %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()

p1_1 <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x_var="vaccine_doses", y_var="n_per_kb_adj", color_var="lineage_sim", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Vaccine")
p1_1 <- p1_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p1_1_p <- p1_1+ # comirnaty 2-doses
   # ggtitle("A")+
   geom_signif(y_position = 2.5+seq(0,0)/5, xmin = c(0.8), xmax = c(1.2), annotation = c("**"), color="black", vjust=0.65, tip_length = 0.01)
p1_1_p <- p1_1_p+ # unvaccinated
   geom_signif(y_position = 2.5+seq(0,4)/5, xmin = c(4.7,4.7,4.85,4.85,4.85), xmax = c(5,5.15,5,5.15,5.3), annotation = c("*","**","**","**", "*"), color="black", vjust=0.65, tip_length = 0.01)

p_out <- p1_1_ori_p/p1_1_p + plot_layout(guides='collect') & theme(legend.position='bottom') & guides(col = guide_legend(nrow = 2))
ggsave("../results/Num_isnvs_by_vaccine_more.pdf", width=8, height=8)
save_pptx("../results/Num_isnvs_by_vaccine_more.pptx", width=8, height=8)

df_test_isnvs <- readxl::read_excel("../results/df_test_num_isnvs_by_vaccine_gene_more.xlsx")
df_test_isnvs %>% filter(gene=="Full genome", check_p_adj, same_lineage) %>% select(var1, var2, notation_adj) %>% t()
p1_2 <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), x_var="lineage_sim", y_var="n_per_kb_adj", color_var="vaccine_doses", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Lineage")
p1_2 <- p1_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses) + scale_x_discrete(guide = guide_axis(n.dodge = 2))

p1_2_p <- p1_2 + # Delta
   geom_signif(y_position = 2.5+seq(0,0)/5, xmin = c(3.8), xmax = c(4.2), annotation = c("*"), color="black", vjust=0.65, tip_length = 0.01)
p1_2_p <- p1_2_p + # Omicron
   geom_signif(y_position = 2.5+seq(0,1)/5, xmin = c(4.7, 5), xmax = c(4.85, 4.85), annotation = c("**","**"), color="black", vjust=0.65, tip_length = 0.01)
p1_2_p <- p1_2_p + theme(legend.position='bottom') + guides(col = guide_legend(nrow = 2))
ggsave("../results/Num_isnvs_by_vaccine_2_more.pdf", width=8, height=6)
save_pptx("../results/Num_isnvs_by_vaccine_2_more.pptx", width=8, height=6)

p1_3 <- plot_box(df_plot_n_gene_meta_adj %>% filter(gene %in% c("ORF1ab", "S")), x_var="vaccine_doses", y_var="n_per_kb_adj", color_var="lineage_sim", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Vaccine")
p1_3 <- p1_3 + scale_color_manual(name="Lineage", values=colors_lineage)+ facet_grid(rows=vars(gene))
p1_3 <- p1_3 + theme(legend.position='bottom') + guides(col = guide_legend(nrow = 2))
# ggsave("../results/Num_isnvs_by_vaccine_gene_more.pdf", width=8, height=8)
# save_pptx("../results/Num_isnvs_by_vaccine_gene_more.pptx", width=8, height=6)

p1_4 <- plot_box(df_plot_n_gene_meta_adj %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="n_per_kb_adj", color_var="vaccine_doses", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Lineage")
p1_4 <- p1_4 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses)+facet_grid(rows=vars(gene))+ scale_x_discrete(guide = guide_axis(n.dodge = 2))
# ggsave("../results/Num_isnvs_by_vaccine_gene_2_more.pdf", width=8, height=6)
# save_pptx("../results/Num_isnvs_by_vaccine_gene_2_more.pptx", width=8, height=6)

ggscatter(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome") %>% filter(Vaccine!="Unvaccinated"), x = "days_since_last_dose", y = "n_per_kb_adj",
   color = "vaccine_doses",
   add = "reg.line",  # Add regressin line
#    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   size = 0.8
   )+
   xlab("Days since last dose")+
   ylab("Number of iSNVs per Kb (adjusted)")+
   facet_wrap(vars(as.character(lineage_sim)),ncol=1)+
   scale_color_manual(name="Vaccine", values=colors_vaccine_doses[1:4])+
   scale_fill_manual(name="Vaccine", values=colors_vaccine_doses[1:4])+
   stat_cor(aes(color = vaccine_doses), label.x = 3)+
   NULL
ggsave("../results/cor_last_dose_isnvs.pdf", width=8, height=12)

### non-vaccinated density
p_den <- ggplot(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome")) +
	geom_density(aes(x=n_per_kb_adj, color=lineage_sim, fill=lineage_sim), alpha=0.8)+
	facet_grid(rows=vars(lineage_sim), cols=vars(vaccine_doses))+
	scale_color_manual(name="Lineage", values=colors_lineage)+
	scale_fill_manual(name="Lineage", values=colors_lineage)+
	xlab("Nubmer of iSNVs per Kb (adjusted)")+
	theme(axis.text.x = element_text(angle = 45, hjust=1))+
	NULL
ggsave("../results/Density_isnvs_by_vaccine_more.pdf", width=8, height=6)
save_pptx("../results/Density_isnvs_by_vaccine_more.pptx", width=8, height=6)

## plot for MAF
df_test_maf <- readxl::read_excel("../results/df_MAF_by_vaccine_gene_more.xlsx")
df_test_maf %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()

p2_1 <- plot_box(df_snvs_meta_add_qc %>% filter(gene=="Full genome"), x_var="vaccine_doses", y_var="sec_freq", color_var="lineage_sim", y_lab="MAF", x_lab="Vaccine")
p2_1 <- p2_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p2_1_p <- p2_1 + # CoronaVac 2-doses
   geom_signif(y_position = 0.5+seq(0,0)/5, xmin = c(2.8), xmax = c(3.2), annotation = c("*"), color="black", vjust=0.65, tip_length = 0.01)
ggsave("../results/MAF_more.pdf", width=8, height=6)
save_pptx("../results/MAF_more.pptx", width=8, height=6)

df_test_maf %>% filter(gene=="Full genome", check_p_adj, same_lineage) %>% select(var1, var2, notation_adj) %>% t()
p2_2 <- plot_box(df_snvs_meta_add_qc %>% filter(gene=="Full genome"), x_var="lineage_sim", y_var="sec_freq", color_var="vaccine_doses", y_lab="MAF", x_lab="Lineage")
p2_2 <- p2_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave("../results/MAF_2_more.pdf", width=8, height=6)
save_pptx("../results/MAF_2_more.pptx", width=8, height=6) 

df_test_maf %>% filter(gene!="Full genome", check_p_adj, same_vaccine) %>% select(gene, var1, var2, notation_adj) %>% t()
p2_3 <- plot_box(df_snvs_meta_add_qc %>% filter(gene %in% c("ORF1ab", "S")), x_var="vaccine_doses", y_var="sec_freq", color_var="lineage_sim", y_lab="MAF", x_lab="Vaccine")
p2_3 <- p2_3 + scale_color_manual(name="Lineage", values=colors_lineage)+facet_grid(rows=vars(gene))
# ggsave("../results/MAF_gene_more.pdf", width=8, height=6)
# save_pptx("../results/MAF_gene_more.pptx", width=8, height=6)

p2_4 <- plot_box(df_snvs_meta_add_qc %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="sec_freq", color_var="vaccine_doses", y_lab="MAF", x_lab="Lineage")
p2_4 <- p2_4 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses)+
	facet_grid(rows=vars(gene)) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
# ggsave("../results/MAF_gene_2_more.pdf", width=8, height=6)
# save_pptx("../results/MAF_gene_2_more.pptx", width=8, height=6)

# plot for pi
df_test_pi <- readxl::read_excel("../results/df_test_diversity_pi_by_vaccine_gene_more.xlsx")
df_test_pi %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()

mean(df_plot_pi$piN, na.rm=T)
mean(df_plot_pi$piS, na.rm=T)
# boxplot(df_plot_pi$piN-df_plot_pi$piS)
wilcox.test(df_plot_pi$piN-df_plot_pi$piS)

p3_1 <- plot_box(df_plot_pi, x_var="vaccine_doses", y_var="pi", color_var="lineage_sim", y_lab="Nucleotide diversity (adjusted)", x_lab="Vaccine")
p3_1 <- p3_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p3_1_p <- p3_1+ # comirnaty 2-doses
   # ggtitle("A")+
   geom_signif(y_position = 0.0008+seq(0,0)/5000, xmin = c(0.8), xmax = c(1.2), annotation = c("**"), color="black", vjust=0.65, tip_length = 0.01)
p3_1_p <- p3_1_p+ # unvaccinated
   geom_signif(y_position = 0.0008+seq(0,3)/18000, xmin = c(4.7,4.7,4.85,4.85), xmax = c(5,5.15,5,5.15), annotation = c("*","**","**","**"), color="black", vjust=0.65, tip_length = 0.01)

ggsave("../results/diveristy_pi_more.pdf", width=8, height=6)
save_pptx("../results/diveristy_pi_more.pptx", width=8, height=6)

df_test_pi <- readxl::read_excel("../results/df_test_diversity_pi_by_vaccine_gene_more.xlsx")
df_test_pi %>% filter(gene=="Full genome", check_p_adj, same_lineage) %>% select(var1, var2, notation_adj) %>% t()

p3_2 <- plot_box(df_plot_pi, x_var="lineage_sim", y_var="pi", color_var="vaccine_doses", y_lab="Nucleotide diversity (adjusted)", x_lab="Lineage")
p3_2 <- p3_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
p3_2_p <- p3_2 + # Delta
   geom_signif(y_position = 0.0008+seq(0,0)/18000, xmin = c(3.8), xmax = c(4.2), annotation = c("*"), color="black", vjust=0.65, tip_length = 0.01)
p3_2_p <- p3_2_p + # Omicron
   geom_signif(y_position = 0.0008+seq(0,3)/18000, xmin = c(4.85, 4.85, 4.85, 4.85), xmax = c(4.7, 5, 5.15, 5.3), annotation = c("**","**","*","*"), color="black", vjust=0.65, tip_length = 0.01)

ggsave("../results/diveristy_pi_2_more.pdf", width=8, height=6)
save_pptx("../results/diveristy_pi_2_more.pptx", width=8, height=6)
 
# p3_3 <- plot_box(df_plot_pi %>% filter(gene %in% c("ORF1ab", "S")), x_var="Vaccine", y_var="pi", color_var="lineage_sim", y_lab="Nucleotide diversity (adjusted)")
# p3_3 <- p3_3 + scale_color_manual(name="Lineage", values=colors_lineage)+
# 	facet_grid(rows=vars(gene)) 
# ggsave("../results/diveristy_pi_gene_more.pdf", width=8, height=6)
# save_pptx("../results/diveristy_pi_gene_more.pptx", width=8, height=6)

# p3_4 <- plot_box(df_plot_pi %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="pi", color_var="Vaccine", y_lab="Nucleotide diversity (adjusted)")
# p3_4 <- p3_4 + scale_color_manual(name="Vaccine", values=colors_vaccine_new)+
# 	facet_grid(rows=vars(gene)) 
# ggsave("../results/diveristy_pi_gene_2_more.pdf", width=8, height=6)
# save_pptx("../results/diveristy_pi_gene_2_more.pptx", width=8, height=6)


# update on 2022-07-20
## Figure 2 upper
df_test_isnvs %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()
p1_1 <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome", Vaccine=="Unvaccinated"), x_var="vaccine_doses", y_var="n_per_kb_adj", color_var="lineage_sim", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Vaccine")+theme(axis.title.x = element_blank())
p1_1 <- p1_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p1_1_p <- p1_1+ # comirnaty 2-doses
   ggtitle("A")+
   geom_signif(y_position = 2.5+seq(0,4)/5, xmin = c(4.7,4.7,4.85,4.85,4.85)-4, xmax = c(5,5.15,5,5.15,5.3)-4, annotation = c("*","**","**","**", "*"), color="black", vjust=0.65, tip_length = 0.01)

df_test_maf %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()
p2_1 <- plot_box(df_snvs_meta_add_qc %>% filter(gene=="Full genome", Vaccine=="Unvaccinated"), x_var="vaccine_doses", y_var="sec_freq", color_var="lineage_sim", y_lab="MAF", x_lab="Vaccine")+theme(axis.title.x = element_blank())
p2_1 <- p2_1 + scale_color_manual(name="Lineage", values=colors_lineage)

p2_1_p <- p2_1 # no significant difference found

df_test_pi %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()
p3_1 <- plot_box(df_plot_pi %>% filter(Vaccine=="Unvaccinated"), x_var="vaccine_doses", y_var="pi", color_var="lineage_sim", y_lab="Nucleotide diversity (adjusted)", x_lab="Vaccine")+theme(axis.title.x = element_blank())
p3_1 <- p3_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p3_1_p <- p3_1+
   ggtitle("C")+ # unvaccinated
   geom_signif(y_position = 0.0009+seq(0,3)/18000, xmin = c(4.7,4.7,4.85,4.85)-4, xmax = c(5,5.15,5,5.15)-4, annotation = c("*","**","**","**"), color="black", vjust=0.65, tip_length = 0.01)

p_out_upper_1 <- (p1_1_p) + (p2_1_p+ggtitle("B")) + (p3_1_p) + plot_layout(nrow=1, guides="collect") & theme(legend.position='bottom') & guides(col = guide_legend(nrow = 1))
ggsave("../results/Figure 2_upper.pdf", width=8, height=4, plot= p_out_upper_1)

## Figure 3 upper
df_test_isnvs %>% filter(gene=="Full genome", check_p_adj, same_lineage) %>% select(var1, var2, notation_adj) %>% t()
p1_2 <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome", lineage_sim %in% lineages_all[4:5]), x_var="lineage_sim", y_var="n_per_kb_adj", color_var="vaccine_doses", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Lineage")+theme(axis.title.x = element_blank())
p1_2 <- p1_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses, limits = force)
p1_2_p <- p1_2 + # Delta
   ggtitle("A")+
   geom_signif(y_position = 2.5+seq(0,0)/5, xmin = c(3.8)-3, xmax = c(4.2)-3, annotation = c("*"), color="black", vjust=0.65, tip_length = 0.01)
p1_2_p <- p1_2_p + # Omicron
   geom_signif(y_position = 2.5+seq(0,1)/5, xmin = c(4.7, 5)-3, xmax = c(4.85, 4.85)-3, annotation = c("**", "**"), color="black", vjust=0.65, tip_length = 0.01)
p1_2_p <- p1_2_p + theme(legend.position='bottom') + guides(col = guide_legend(nrow = 1))+ scale_x_discrete(guide = guide_axis(n.dodge = 2))

df_test_maf %>% filter(gene=="Full genome", check_p_adj, same_lineage) %>% select(var1, var2, notation_adj) %>% t()
p2_2 <- plot_box(df_snvs_meta_add_qc %>% filter(gene=="Full genome", lineage_sim %in% lineages_all[4:5]) , x_var="lineage_sim", y_var="sec_freq", color_var="vaccine_doses", y_lab="MAF", x_lab="Vaccine")+theme(axis.title.x = element_blank())
p2_2 <- p2_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses, limits = force)+ scale_x_discrete(guide = guide_axis(n.dodge = 2))

df_test_pi %>% filter(gene=="Full genome", check_p_adj, same_lineage) %>% select(var1, var2, notation_adj) %>% t()

p3_2 <- plot_box(df_plot_pi %>% filter(lineage_sim %in% lineages_all[4:5]) , x_var="lineage_sim", y_var="pi", color_var="vaccine_doses", y_lab="Nucleotide diversity (adjusted)", x_lab="Vaccine")+theme(axis.title.x = element_blank())
p3_2 <- p3_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses, limits = force)+ scale_x_discrete(guide = guide_axis(n.dodge = 2))
p3_2_p <- p3_2 + # Delta
   ggtitle("C")+
   geom_signif(y_position = 0.0009+seq(0,0)/18000, xmin = c(3.8)-3, xmax = c(4.2)-3, annotation = c("*"), color="black", vjust=0.65, tip_length = 0.01)
p3_2_p <- p3_2_p + # Omicron
    geom_signif(y_position = 0.0009+seq(0,3)/18000, xmin = c(4.85, 4.85, 4.85, 4.85)-3, xmax = c(4.7, 5, 5.15, 5.3)-3, annotation = c("**","**","*","*"), color="black", vjust=0.65, tip_length = 0.01)

p_out_upper_2 <- (p1_2_p) + (p2_2+ggtitle("B")) + (p3_2_p) + plot_layout(nrow=1, guides="collect") & theme(legend.position='bottom') & guides(col = guide_legend(nrow = 1))
ggsave("../results/Figure 3_upper.pdf", width=8, height=4, plot= p_out_upper_2)

## Figure S3
df_test_isnvs %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()
p1_1 <- plot_box(df_plot=df_plot_n_gene_meta_adj %>% filter(gene=="Full genome", lineage_sim %in% lineages_all[4:5]), x_var="vaccine_doses", y_var="n_per_kb_adj", color_var="lineage_sim", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Vaccine")+theme(axis.title.x = element_blank())
p1_1 <- p1_1 + scale_color_manual(name="Lineage", values=colors_lineage, limits = force)
p1_1_p <- p1_1+
   ggtitle("A")+ # Comirnaty_Doses2
   geom_signif(y_position = 2.5, xmin = c(0.8), xmax = c(1.2), annotation = c("**"), color="black", vjust=0.65, tip_length = 0.01)

df_test_maf %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()
p2_1 <- plot_box(df_snvs_meta_add_qc %>% filter(gene=="Full genome", lineage_sim %in% lineages_all[4:5]), x_var="vaccine_doses", y_var="sec_freq", color_var="lineage_sim", y_lab="MAF", x_lab="Vaccine")+theme(axis.title.x = element_blank())
p2_1 <- p2_1 + scale_color_manual(name="Lineage", values=colors_lineage, limits = force)
p2_1_p <- p2_1+
   ggtitle("B")+
   geom_signif(y_position = 0.5, xmin = c(0.8)+2, xmax = c(1.2)+2, annotation = c("*"), color="black", vjust=0.65, tip_length = 0.01)

df_test_pi %>% filter(gene=="Full genome", check_p_adj, same_vaccine) %>% select(var1, var2, notation_adj) %>% t()
p3_1 <- plot_box(df_plot_pi %>% filter(lineage_sim %in% lineages_all[4:5]), x_var="vaccine_doses", y_var="pi", color_var="lineage_sim", y_lab="Nucleotide diversity (adjusted)", x_lab="Vaccine")+theme(axis.title.x = element_blank())
p3_1 <- p3_1 + scale_color_manual(name="Lineage", values=colors_lineage, limits = force)
p3_1_p <- p3_1+
   ggtitle("C")+ # Comirnaty_Doses2
   geom_signif(y_position = 0.0009+seq(0,0)/18000, xmin = c(0.8), xmax = c(1.2), annotation = "**", color="black", vjust=0.65, tip_length = 0.01)

p_out_2 <- (p1_1_p) + (p2_1_p) + (p3_1_p) + plot_layout(ncol=1, guides="collect") & theme(legend.position='bottom') & guides(col = guide_legend(nrow = 1))

ggsave("../results/Figure 2_supp.pdf", width=8/1.618, height=8, plot= p_out_2)

# difference between genes/positions
(df_tmp <- df_plot_n_gene_meta_adj %>% group_by(gene) %>% summarise(mean_n_per_kb=mean(n_per_kb,na.rm=T), mean_n_per_kb_adj=mean(n_per_kb_adj,na.rm=T)) %>% mutate(gene=factor(gene, levels=c("Full genome", df_orf_sim$sequence))) %>% arrange(gene))
write_csv(df_tmp %>% arrange(mean_n_per_kb_adj), "../results/mean_isnvs_per_kb.csv")

df_tmp <- df_plot_n_gene_meta_adj %>% filter(gene%in%c("ORF1ab", "S")) %>% select(gene, n_per_kb_adj, vaccine_doses, lineage_sim, sample)
median(df_tmp$n_per_kb_adj[df_tmp$gene=="ORF1ab"])
median(df_tmp$n_per_kb_adj[df_tmp$gene=="S"])
wilcox.test(df_tmp$n_per_kb_adj[df_tmp$gene=="ORF1ab"], df_tmp$n_per_kb_adj[df_tmp$gene=="S"]) # P<0.001

df_tmp2 <- df_tmp %>% pivot_wider(names_from="gene", values_from="n_per_kb_adj")
df_tmp2 <- left_join(df_tmp2, (df_tmp2 %>% group_by(lineage_sim) %>% summarise(n=n()) %>% mutate(lineage_sim_n=paste0(lineage_sim, " (N=", n, ")" )) %>% select(-n)))

p_out <- ggscatter(df_tmp2, x = "ORF1ab", y = "S", color = "vaccine_doses",
   add = "reg.line",  # Add regressin line
   conf.int = TRUE, # Add confidence interval
   size = 0.8
   )+
   scale_color_manual(name="Vaccine", values=colors_vaccine_doses)+
   scale_fill_manual(name="Vaccine", values=colors_vaccine_doses)+
   stat_cor(aes(color = vaccine_doses), label.x = 0.1)+
   facet_wrap(vars(as.character(lineage_sim_n)),ncol=2)+
   NULL
p_out <- p_out+geom_abline(intercept = 0, slope = 1, linetype="dashed")
ggsave("../results/cor_isnvs_orf1ab_s.pdf", width=8, height=12) # interesting

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

df_tmp <- df_snvs_meta_add_qc %>% filter(gene%in%c("ORF1ab", "S")) %>% select(gene, sec_freq, vaccine_doses, lineage_sim, sample) %>% group_by(gene,sample,vaccine_doses, lineage_sim) %>% summarise(sec_freq=mean(sec_freq))
median(df_tmp$sec_freq[df_tmp$gene=="ORF1ab"])
median(df_tmp$sec_freq[df_tmp$gene=="S"])
wilcox.test(df_tmp$sec_freq[df_tmp$gene=="ORF1ab"], df_tmp$sec_freq[df_tmp$gene=="S"])

df_tmp2 <- df_tmp %>% pivot_wider(names_from="gene", values_from="sec_freq")
p_out <- ggscatter(df_tmp2, x = "ORF1ab", y = "S", color = "vaccine_doses",
   add = "reg.line",  # Add regressin line
   conf.int = TRUE, # Add confidence interval
   size = 0.8
   )+
   scale_color_manual(name="Vaccine", values=colors_vaccine_doses)+
   scale_fill_manual(name="Vaccine", values=colors_vaccine_doses)+
   stat_cor(aes(color = vaccine_doses), label.x = 0.1)+
   facet_wrap(vars(as.character(lineage_sim)),ncol=2)+
   NULL
p_out2 <- p_out+geom_abline(intercept = 0, slope = 1, linetype="dashed")
ggsave("../results/cor_MAF_orf1ab_s.pdf", width=8, height=12)

