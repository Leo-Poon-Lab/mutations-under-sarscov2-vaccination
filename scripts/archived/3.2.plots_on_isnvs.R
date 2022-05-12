library(tidyverse)
library(Biostrings)
library(ggsci)
library(patchwork)
library(parallel)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")
source('./helper/plot_box.R')
source("./helper/df_test.R")

colors_lineage=c("#e41a1c", "#33a02c", "#1f78b4", "#ff7f00") 
names(colors_lineage) <- c("Alpha", "Delta", "Omicron", "B.1.36")
colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine)=c("BioNTech", "Sinovac", "Non-vaccinated")

load("../results/df_bam_rst.rdata")
load("../results/df_plot_n_gene.rdata")
load("../results/df_plot_pi.rdata")
df_meta <- read_csv("../results/df_samples_clean.csv", guess_max = 60000)
df_meta <- df_meta %>% filter(Sample %in% unique(df_plot_n_gene$sample))

df_meta <- df_meta %>% filter(Doses==0 | Doses==2) # remove partial vaccination
df_plot_n_gene <- df_plot_n_gene %>% filter(sample %in% df_meta$Sample)
df_plot_pi <- df_plot_pi %>% filter(sample %in% df_meta$Sample)

# Data and quality control
unique(df_meta$Sample)
range(lubridate::dmy(df_meta$`Report date`))
table(df_meta$lineage_sim)
table(df_meta$Vaccine)
table(df_meta$lineage_sim, df_meta$Vaccine)
table(df_meta$lineage_sim, df_meta$Doses, df_meta$Vaccine)

df_tmp <- df_plot_n_gene %>% filter(gene=="Full genome") %>% select(sample, Vaccine, n_per_kb)
df_tmp$Ct_value <- df_meta$Ct_value[match(df_tmp$sample, df_meta$Sample)]
df_tmp$age <- df_meta$`Age (yr)`[match(df_tmp$sample, df_meta$Sample)]
df_tmp$days_since_last_dose <- df_meta$days_since_last_dose[match(df_tmp$sample, df_meta$Sample)]

cor.test(as.numeric(df_tmp$Ct_value), df_tmp$n_per_kb) # correlation between Ct and number of iSNVs
plot(as.numeric(df_tmp$Ct_value), df_tmp$n_per_kb)
cor.test(df_tmp$age, df_tmp$n_per_kb) 
cor.test(df_tmp$days_since_last_dose, df_tmp$n_per_kb) 

# General plots
## plot for incidence
p1_1 <- plot_box(df_plot=df_plot_n_gene %>% filter(gene=="Full genome"), x_var="Vaccine", y_var="n_per_kb", color_var="lineage_sim", y_lab="Numer of iSNVs per Kb")
p1_1 <- p1_1 + scale_color_manual(name="Lineage", values=colors_lineage)
ggsave("../results/Num_isnvs_by_vaccine.pdf", width=6, height=4)
save_pptx("../results/Num_isnvs_by_vaccine.pptx", width=6, height=4)

p1_2 <- plot_box(df_plot=df_plot_n_gene %>% filter(gene=="Full genome"), x_var="lineage_sim", y_var="n_per_kb", color_var="Vaccine", y_lab="Numer of iSNVs per Kb")
p1_2 <- p1_2 + scale_color_manual(name="Vaccine", values=colors_vaccine)
ggsave("../results/Num_isnvs_by_vaccine_2.pdf", width=6, height=4)
save_pptx("../results/Num_isnvs_by_vaccine_2.pptx", width=6, height=4)

p1_3 <- plot_box(df_plot_n_gene %>% filter(gene %in% c("ORF1ab", "S")), x_var="Vaccine", y_var="n_per_kb", color_var="lineage_sim", y_lab="Numer of iSNVs per Kb")
p1_3 <- p1_3 + scale_color_manual(name="Lineage", values=colors_lineage)+
	facet_grid(rows=vars(gene))
ggsave("../results/Num_isnvs_by_vaccine_gene.pdf", width=8, height=6)
save_pptx("../results/Num_isnvs_by_vaccine_gene.pptx", width=8, height=6)

p1_4 <- plot_box(df_plot_n_gene %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="n_per_kb", color_var="Vaccine", y_lab="Numer of iSNVs per Kb")
p1_4 <- p1_4 + scale_color_manual(name="Vaccine", values=colors_vaccine)+
	facet_grid(rows=vars(gene))
ggsave("../results/Num_isnvs_by_vaccine_gene_2.pdf", width=8, height=6)
save_pptx("../results/Num_isnvs_by_vaccine_gene_2.pptx", width=8, height=6)

df_wilc_test <- cal_wilc_test(df_plot_n_gene, "n_per_kb")
write_csv(df_wilc_test , "../results/df_test_num_isnvs_by_vaccine_gene.csv")

## plot for pi
p2_1 <- plot_box(df_plot_pi %>% filter(gene %in% c("Full genome")), x_var="Vaccine", y_var="pi_scale", color_var="lineage_sim", y_lab=expression("Nucleotide diversity"~pi))
p2_1 <- p2_1 + scale_color_manual(name="Lineage", values=colors_lineage)
ggsave("../results/diveristy_pi.pdf", width=8, height=6)
save_pptx("../results/diveristy_pi.pptx", width=8, height=6)

p2_2 <- plot_box(df_plot_pi %>% filter(gene %in% c("Full genome")), x_var="lineage_sim", y_var="pi_scale", color_var="Vaccine", y_lab=expression("Nucleotide diversity"~pi))
p2_2 <- p2_2 + scale_color_manual(name="Vaccine", values=colors_vaccine)
ggsave("../results/diveristy_pi_2.pdf", width=8, height=6)
save_pptx("../results/diveristy_pi_2.pptx", width=8, height=6)

p2_3 <- plot_box(df_plot_pi %>% filter(gene %in% c("ORF1ab", "S")), x_var="Vaccine", y_var="pi_scale", color_var="lineage_sim", y_lab=expression("Nucleotide diversity"~pi))
p2_3 <- p2_3 + scale_color_manual(name="Lineage", values=colors_lineage)+
	facet_grid(rows=vars(gene))
ggsave("../results/diveristy_pi_gene.pdf", width=8, height=6)
save_pptx("../results/diveristy_pi_gene.pptx", width=8, height=6)

p2_4 <- plot_box(df_plot_pi %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="pi_scale", color_var="Vaccine", y_lab=expression("Nucleotide diversity"~pi))
p2_4 <- p2_4 + scale_color_manual(name="Vaccine", values=colors_vaccine)+
	facet_grid(rows=vars(gene))
ggsave("../results/diveristy_pi_gene_2.pdf", width=8, height=6)
save_pptx("../results/diveristy_pi_gene_2.pptx", width=8, height=6)

df_wilc_test <- cal_wilc_test(df_plot_pi, "pi_scale")
write_csv(df_wilc_test , "../results/df_test_diversity_pi_by_vaccine_gene.csv")

# plot for piN_piS
p3_1 <- plot_box(df_plot_pi %>% filter(gene %in% c("Full genome")), x_var="Vaccine", y_var="piN_piS", color_var="lineage_sim", y_lab=expression(pi[N]~"/"~pi["S"]))
p3_1 <- p3_1 + scale_color_manual(name="Lineage", values=colors_lineage)
ggsave("../results/diveristy_piN_piS.pdf", width=8, height=6)
save_pptx("../results/diveristy_piN_piS.pptx", width=8, height=6)

p3_2 <- plot_box(df_plot_pi %>% filter(gene %in% c("Full genome")), x_var="lineage_sim", y_var="piN_piS", color_var="Vaccine", y_lab=expression(pi[N]~"/"~pi["S"]))
p3_2 <- p3_2 + scale_color_manual(name="Vaccine", values=colors_vaccine)
ggsave("../results/diveristy_piN_piS_2.pdf", width=8, height=6)
save_pptx("../results/diveristy_piN_piS_2.pptx", width=8, height=6)

p3_3 <- plot_box(df_plot_pi %>% filter(gene %in% c("ORF1ab", "S")), x_var="Vaccine", y_var="piN_piS", color_var="lineage_sim", y_lab=expression(pi[N]~"/"~pi["S"]))
p3_3 <- p3_3 + scale_color_manual(name="Lineage", values=colors_lineage)+
	facet_grid(rows=vars(gene))
ggsave("../results/diveristy_piN_piS_gene.pdf", width=8, height=6)
save_pptx("../results/diveristy_piN_piS_gene.pptx", width=8, height=6)

p3_4 <- plot_box(df_plot_pi %>% filter(gene %in% c("ORF1ab", "S")), x_var="lineage_sim", y_var="piN_piS", color_var="Vaccine", y_lab=expression(pi[N]~"/"~pi["S"]))
p3_4 <- p3_4 + scale_color_manual(name="Vaccine", values=colors_vaccine)+
	facet_grid(rows=vars(gene))
ggsave("../results/diveristy_piN_piS_gene_2.pdf", width=8, height=6)
save_pptx("../results/diveristy_piN_piS_gene_2.pptx", width=8, height=6)

df_wilc_test <- cal_wilc_test(df_plot_pi, "piN_piS", genes=c("Full genome", "ORF1ab", "S"))
write_csv(df_wilc_test, "../results/df_test_diversity_piN_piS_by_vaccine_gene.csv")

### combine figures
p_out_1 <- (p1_1+ggtitle("A")) + (p2_1+ggtitle("B")) + (p3_1+ggtitle("C")) + plot_layout(ncol=1, guides="collect") & theme(legend.position='bottom')
p_out_2 <- (p1_2+ggtitle("D")) + (p2_2+ggtitle("E")) + (p3_2+ggtitle("F")) + plot_layout(ncol=1, guides="collect") & theme(legend.position='bottom')
p_out <- p_out_1 | p_out_2
ggsave("../results/plot_comb_by_variant_full.pdf", width=10, height=8)
save_pptx("../results/plot_comb_by_variant_full.pptx", width=8, height=8)


# difference between variants
df_tmp <- df_plot_n_gene %>% filter(Vaccine=="Non-vaccinated", gene=="Full genome")
kruskal.test(df_tmp$n_per_kb, df_tmp$lineage_sim)

df_tmp <- df_plot_pi %>% filter(Vaccine=="Non-vaccinated", gene=="Full genome")
kruskal.test(df_tmp$pi_scale, df_tmp$lineage_sim)

# difference between vaccination

# difference between genes/positions, 

# age, patient status