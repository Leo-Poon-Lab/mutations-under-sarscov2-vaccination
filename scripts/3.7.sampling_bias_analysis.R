library(tidyverse)
library(lubridate)
library(patchwork)
library(ggsignif)
source("./helper/df_test.R")
source('./helper/plot_box.R')

df_meta <- read_csv("../results/df_samples.csv", guess_max=100000)
df_meta <- df_meta %>% filter(lineage_sim != "22B (Omicron, BA.5.*)")
df_meta$Vaccine <- gsub("BioNTech", "Comirnaty", df_meta$Vaccine)
df_meta$Vaccine <- gsub("Sinovac", "CoronaVac", df_meta$Vaccine)
sort(table(df_meta$specimen_type))
df_meta$specimen_type[df_meta$specimen_type=="specimen"] <- NA
df_meta$specimen_type[df_meta$specimen_type=="Swab"] <- NA

lineages_all <- sort(unique(df_meta$lineage_sim))
colors_lineage=rev(c("#4daf4a", "#984ea3", "#ff7f00", "#f781bf", "#666666"))
names(colors_lineage) <- lineages_all
vaccine_all <- sort(unique(df_meta$Vaccine))
colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine)=vaccine_all
colors_vaccine_new <- colors_vaccine
names(colors_vaccine_new)=c("Comirnaty", "CoronaVac", "Unvaccinated")
colors_vaccine_doses = c("#a65628", "#633318", "#7570b3", "#46436b", "#999999")
names(colors_vaccine_doses)=c("Comirnaty\n(Doses=2)", "Comirnaty\n(Doses=3)", "CoronaVac\n(Doses=2)", "CoronaVac\n(Doses=3)", "Unvaccinated")

df_pcr_tech <- readxl::read_excel("../data/PCR_tech_data.xlsx")
names(df_pcr_tech)[1] <- "vm_id"
df_meta <- left_join(df_meta, df_pcr_tech, "vm_id")

load("../results/df_plot_n_gene_adj.rdata")
df_plot_n_gene_meta_adj$lineage_sim <- factor(df_plot_n_gene_meta_adj$lineage_sim, levels = names(colors_lineage))
df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, df_meta %>% select(sample, Doses), "sample")
df_plot_n_gene_meta_adj$Vaccine <- gsub("BioNTech", "Comirnaty", df_plot_n_gene_meta_adj$Vaccine)
df_plot_n_gene_meta_adj$Vaccine <- gsub("Sinovac", "CoronaVac", df_plot_n_gene_meta_adj$Vaccine)
df_plot_n_gene_meta_adj$Vaccine <- factor(df_plot_n_gene_meta_adj$Vaccine, levels=names(colors_vaccine), labels=names(colors_vaccine_new))
df_plot_n_gene_meta_adj$vaccine_doses <- paste0(df_plot_n_gene_meta_adj$Vaccine, "\n(Doses=", df_plot_n_gene_meta_adj$Doses, ")")
df_plot_n_gene_meta_adj$vaccine_doses[df_plot_n_gene_meta_adj$Vaccine=="Unvaccinated"] <- "Unvaccinated"
df_plot_n_gene_meta_adj$vaccine_doses <- factor(df_plot_n_gene_meta_adj$vaccine_doses, levels = names(colors_vaccine_doses))

df_meta$lineage_sim <- factor(df_meta$lineage_sim, levels = names(colors_lineage))
# collection lag
df_meta$detection_lag <- as.numeric(df_meta$`Report date` - ymd(df_meta$`Onset date`))
df_meta$collection_lag <- as.numeric(df_meta$collection_date - ymd(df_meta$`Onset date`))
df_meta$collection_lag[df_meta$collection_lag>100] <- NA # one ourlier
df_meta$collection_lag[df_meta$collection_lag<0] <- NA # one ourlier

df_meta$vaccine_doses <- paste0(df_meta$Vaccine, "\n(Doses=", df_meta$Doses, ")")
df_meta$vaccine_doses[df_meta$Vaccine=="Unvaccinated"] <- "Unvaccinated"

df_meta$vaccine_doses <- factor(df_meta$vaccine_doses, levels = names(colors_vaccine_doses))
# names(df_meta)

# Viral loads between groups
which((is.na(df_meta$Ct_value)) & (!is.na(df_meta$Ct)))
source("./helper/df_test.R")
df_test <- cal_wilc_test(df_meta, "Ct_value", genes=NA)
highlight_diff(df_test) %>% filter(check_p_adj) %>% select(var1, var2, p_value_adj, notation_adj) %>% t()

# wilcox.test(df_meta %>% filter(Vaccine=="Unvaccinated") %>% filter(lineage_sim %in% names(colors_lineage)[1:2]) %>% .$Ct_value, df_meta %>% filter(Vaccine=="Unvaccinated") %>% filter(!lineage_sim %in% names(colors_lineage)[1:2]) %>% .$Ct_value)

p1_1 <- plot_box(df_plot=df_meta %>% filter(Vaccine=="Unvaccinated"), x_var="vaccine_doses", y_var="Ct_value", color_var="lineage_sim", y_lab="Ct value", x_lab="Vaccine")
p1_1 <- p1_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p1_1_p <- p1_1+ # unvaccinated
   # ggtitle("A")+
   geom_signif(y_position = 28.5+seq(0,3)/1.5, xmin = c(0.85, 0.85, 0.85, 0.85), xmax = c(0.7, 1, 1.15, 1.3), annotation = c("**", "**", "**", "**"), color="black", vjust=0.65, tip_length = 0.01)

p1_2 <- plot_box(df_plot=df_meta %>% filter(lineage_sim%in%names(colors_lineage)[4:5]), x_var="lineage_sim", y_var="Ct_value", color_var="vaccine_doses", y_lab="Ct value", x_lab="Lineage")
p1_2 <- p1_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses) + 
	# scale_x_discrete(guide = guide_axis(n.dodge = 2))+
	NULL
p1_2_p <- p1_2

p_1_out <- (p1_1_p+ggtitle("A"))/(p1_2_p+ggtitle("B"))+
	# plot_layout(guides='collect') & theme(legend.position='bottom') & guides(col = guide_legend(nrow = 3))+
	NULL
ggsave("../results/Ct_Value_by_group.pdf", width=8, height=8)

df_meta$report_year_month <- paste0(year(df_meta$`Report date`), "-", month(df_meta$`Report date`))
df_meta$report_year_month <- naturalsort::naturalfactor(df_meta$report_year_month)

p1_3 <- ggplot(df_meta %>% filter(!is.na(Ct_value)) %>% filter(Vaccine=="Unvaccinated") %>% filter(lineage_sim %in% names(colors_lineage)[1:2]), aes(x=report_year_month, y=Ct_value, color=lineage_sim)) +
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1), size=0.8)+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_wrap(vars(lineage_sim), ncol=1)+
	scale_color_manual(name="Lineage", values=colors_lineage[1:2])+
	xlab("Report dates")+
	ylab("Ct value")
ggsave("../results/Ct_value_by_time.pdf", width=8, height=6)

# Collection timepoints during acute infection
source("./helper/df_test.R")
df_test <- cal_wilc_test(df_meta, "collection_lag", genes=NA)
highlight_diff(df_test) %>% filter(check_p_adj) %>% select(var1, var2, median_var1, median_var2, notation_adj) %>% t()

p2_1 <- plot_box(df_plot=df_meta %>% filter(Vaccine=="Unvaccinated"), x_var="vaccine_doses", y_var="collection_lag", color_var="lineage_sim", y_lab="Collection lag (days)", x_lab="Vaccine")
p2_1 <- p2_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p2_1_p <- p2_1+ # unvaccinated
   # ggtitle("A")+
   geom_signif(y_position = 20+seq(0,0)/1.5, xmin = c(0.85), xmax = c(1.15), annotation = c("*"), color="black", vjust=0.65, tip_length = 0.01)

p2_2 <- plot_box(df_plot=df_meta %>% filter(lineage_sim%in%names(colors_lineage)[4:5]), x_var="lineage_sim", y_var="collection_lag", color_var="vaccine_doses", y_lab="Collection lag (days)", x_lab="Lineage")
p2_2 <- p2_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses) 
p2_2_p <- p2_2

p_2_out <- (p2_1_p+ggtitle("A"))/(p2_2_p+ggtitle("B"))
ggsave("../results/Collection_lag_by_group.pdf", width=8, height=8)

#names(df_meta)

p2_3 <- ggplot(df_meta %>% filter(!is.na(collection_lag)) %>% filter(Vaccine=="Unvaccinated") %>% filter(lineage_sim %in% names(colors_lineage)[1:2]), aes(x=report_year_month, y=collection_lag, color=lineage_sim)) +
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1), size=0.8)+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	# facet_wrap(vars(Classification_sim))+
	scale_color_manual(name="Lineage", values=colors_lineage[1:2])+
	xlab("Report dates")+
	ylab("Collection lag (days)")+
	NULL
ggsave("../results/Collection_lag_by_time.pdf", width=8, height=6)


p2_4 <- ggplot(df_meta %>% filter(!is.na(specimen_type)) %>% filter(Vaccine=="Unvaccinated") %>% filter(lineage_sim %in% names(colors_lineage)[1:2]), aes(x=report_year_month, y=collection_lag, color=specimen_type)) +
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1), size=0.8)+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_grid(cols=vars(lineage_sim), rows=vars(specimen_type))+
	# scale_color_manual(name="Lineage", values=colors_lineage[1:2])+
	xlab("Report dates")+
	ylab("Collection lag (days)")+
	scale_x_discrete(guide = guide_axis(n.dodge = 2))+
	NULL
ggsave("../results/Collection_lag_by_specimen_type.pdf", width=10, height=12, plot=p2_4)


p2_5 <- ggplot(df_meta %>% filter(!is.na(Ct_value)) %>% filter(Vaccine=="Unvaccinated") %>% filter(lineage_sim %in% names(colors_lineage)[1:2]), aes(x=report_year_month, y=Ct_value, color=PCR)) +
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1), size=0.8)+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	facet_wrap(vars(lineage_sim), nrow=2)+
	xlab("Report dates")+
	ylab("Ct value")+
	NULL

# p_wave3_4 <- (p1_3+ggtitle("A"))/(p2_4+ggtitle("B")) & theme(legend.position='bottom')
ggsave("../results/Ct_Value_wave3_4_by_time.pdf", width=8, height=6, plot=p2_5)

# quantile(df_meta$Ct_value[(df_meta$PCR=="NCOV-10") & (df_meta$lineage_sim==names(colors_lineage)[2])], na.rm=T)
# sum((df_meta$PCR=="NCOV-10") & (df_meta$lineage_sim==names(colors_lineage)[2]), na.rm=T)
# quantile(df_meta$Ct_value[(df_meta$PCR!="NCOV-10") & (df_meta$lineage_sim==names(colors_lineage)[2])], na.rm=T)
# quantile(df_meta$Ct_value[(df_meta$lineage_sim==names(colors_lineage)[2])], na.rm=T)

# specimen_type
df_meta$specimen_type_fct <- factor(df_meta$specimen_type, rev(names(sort(table(df_meta$specimen_type)))))
table(df_meta$specimen_type_fct, df_meta$vaccine_doses)
table(df_meta$specimen_type_fct, df_meta$vaccine_doses, df_meta$lineage_sim)

p3_1 <- ggplot(df_meta)+
	geom_histogram(aes(x=specimen_type_fct, fill=vaccine_doses), stat="count")+
	scale_fill_manual(name="Vaccine", values=colors_vaccine_doses) +
	# scale_x_discrete(guide = guide_axis(n.dodge = 4))+
	theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
	xlab("Specimen type")+
	ylab("Count")+
	NULL
p3_2 <- ggplot(df_meta)+
	geom_histogram(aes(x=specimen_type_fct, fill=lineage_sim), stat="count")+
	scale_fill_manual(name="Lineage", values=colors_lineage)+
	# scale_x_discrete(guide = guide_axis(n.dodge = 4))+
	theme_classic()+
	theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
	xlab("Specimen type")+
	ylab("Count")+
	facet_wrap(vars(vaccine_doses), ncol=1, scale="free_y")+
	NULL
p_3_out <- (p3_2+theme(legend.position='bottom'))
ggsave("../results/Specimen_type_by_group.pdf", width=8, height=10)
## Apparently, there are some difference in lineage/vaccine composition among specimen types, this is inevitable, because as the pandemic progress, the sample collection methods are changing. However, we can test whether different specimen type is associated with different number of iSNVs. This can serve as a reference to estimated the bias introduced by sample collection methods.
df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, df_meta %>% select(sample, specimen_type), "sample")
df_plot_n_gene_meta_adj_filter <- df_plot_n_gene_meta_adj %>% filter(gene=="Full genome") %>% filter(!is.na(specimen_type)) %>% filter(specimen_type!="Stool")
df_meta$specimen_type[df_meta$specimen_type=="Stool"] <- NA
p_3_3 <- plot_box(df_plot=df_plot_n_gene_meta_adj_filter, x_var="specimen_type", y_var="n_per_kb_adj", color_var="lineage_sim", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Specimen type")+facet_grid(rows=vars(vaccine_doses), cols=vars(lineage_sim), scale="free_x", space="free")
p_3_3 <- p_3_3 + scale_color_manual(name="Lineage", values=colors_lineage)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='bottom')

df_input <- df_plot_n_gene_meta_adj_filter
df_groups <- df_input %>% select(lineage_sim) %>% unique()
df_input$specimen_type <- gsub(" ", ".", df_input$specimen_type)
var <- "n_per_kb_adj"
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	# i=2
	lineage_i <- df_groups$lineage_sim[i]
	df_tmp <- df_input %>% filter(lineage_sim==lineage_i)
	df_tmp <- df_tmp %>% mutate(x_grps=paste(vaccine_doses, specimen_type, lineage_sim, sep="_"))
	mat_pairs <- combn(unique(df_tmp$x_grps), 2)

	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		value1 <- df_tmp[[var]][df_tmp$x_grps==var1]
		value2 <- df_tmp[[var]][df_tmp$x_grps==var2]
		
		rst <- wilcox.test(value1, value2)
		tibble(var1=var1, var2=var2, var1_mean=mean(value1,na.rm=T), var2_mean=mean(value2,na.rm=T), median_var1=median(value1,na.rm=T), median_var2=median(value2,na.rm=T), p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	df_rst <- df_rst %>% mutate(gene=lineage_i)
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test$p_value_adj <- p.adjust(df_wilc_test$p_value, "BH")

highlight_diff <- function(df){
	var1_split <- strsplit(df$var1, "_", fixed=T)
	var2_split <- strsplit(df$var2, "_", fixed=T)
	check1 <- sapply(var1_split, function(x)x[1]) == sapply(var2_split, function(x)x[1]) # Vaccine
	check2 <- sapply(var1_split, function(x)x[2]) == sapply(var2_split, function(x)x[2]) # specimen type
	check3 <- sapply(var1_split, function(x)x[3]) == sapply(var2_split, function(x)x[3]) # lineage
	
	df$same_vaccine <- check1
	df$same_lineage <- check3
	df$within_group <- df$same_vaccine & df$same_lineage
	medians_low <- sapply(seq_len(nrow(df)), function(i){
		tmp <- abs(c(df$median_var1[i], df$median_var2[i]))
		tmp[which.min(tmp)]
	})
	df$median_diff <- abs(df$median_var1-df$median_var2)/medians_low
	df$median_diff_over_10percent <- df$median_diff>0.1
	df$check_three <- df$median_diff_over_10percent & df$within_group & (df$p_value<0.05)
	
	df$check_p_adj <- df$median_diff_over_10percent & df$within_group & (df$p_value_adj<0.05)

	df$notation <- NA
	df$notation[df$p_value<0.1] <- "^"
	df$notation[df$p_value<0.05] <- "*"
	df$notation[df$p_value<0.01] <- "**"

	df$notation_adj <- NA
	df$notation_adj[df$p_value_adj<0.1] <- "^"
	df$notation_adj[df$p_value_adj<0.05] <- "*"
	df$notation_adj[df$p_value_adj<0.01] <- "**"

	df %>% arrange(gene, check_p_adj, check_three, same_vaccine, same_lineage, var1, var2)
}
df_wilc_test <- highlight_diff(df_wilc_test)
df_wilc_test %>% filter(check_p_adj) %>% select(var1, var2, median_var1, median_var2, median_diff, p_value_adj, notation_adj) %>% t()

annotation_df <- data.frame(
  vaccine_doses = c("Unvaccinated"),
  lineage_sim = c("20A (B.1.36.*)"),
  start = c("Throat saliva"),
  end = c("Throat and nasal swab"),
  y = c(3.8),
  label = c("*")
)

p_3_3 +
  geom_signif(
    data = annotation_df,
    aes(xmin = start, xmax = end, annotations = label, y_position = y), color="black",
    # textsize = 3, vjust = -0.2,
    manual = TRUE
  )
ggsave("../results/num_isnvs_by_Specimen_type_by_group.pdf", width=12, height=12)

p_3_4 <- plot_box(df_plot=df_meta, x_var="specimen_type", y_var="Ct_value", color_var="lineage_sim", y_lab="Ct value", x_lab="Specimen type")+facet_grid(rows=vars(vaccine_doses), cols=vars(lineage_sim), scale="free_x", space="free")
p_3_4 <- p_3_4 + scale_color_manual(name="Lineage", values=colors_lineage)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='bottom')
ggsave("../results/Ct_value_by_Specimen_type_by_group.pdf", width=12, height=8)

p_3_5 <- plot_box(df_plot=df_meta %>% filter(lineage_sim%in%names(colors_lineage)[1:2]), x_var="specimen_type", y_var="Ct_value", color_var="lineage_sim", y_lab="Ct value", x_lab="Specimen type")
p_3_5 <- p_3_5 + scale_color_manual(name="Lineage", values=colors_lineage)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='bottom')
ggsave("../results/Ct_value_by_Specimen_type_by_group_wave3_4.pdf", width=12, height=8)

p_3_6 <- plot_box(df_plot=df_meta %>% filter(lineage_sim%in%names(colors_lineage)[1:2]), x_var="specimen_type", y_var="Ct_value", color_var="lineage_sim", y_lab="Ct value", x_lab="Specimen type")+facet_grid(rows=vars(report_year_month), scale="free", space="free")
p_3_6 <- p_3_6 + scale_color_manual(name="Lineage", values=colors_lineage)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='bottom')
ggsave("../results/Ct_value_by_Specimen_type_by_group_wave3_4_by_time.pdf", width=12, height=20)

# days_since_last_dose
names(df_meta)
source("./helper/df_test.R")
df_test <- cal_wilc_test(df_meta %>% filter(Vaccine!="Unvaccinated"), "days_since_last_dose", genes=NA)
highlight_diff(df_test) %>% filter(p_value_adj<0.05) %>% filter(same_lineage) %>% arrange(var1, var2) %>% select(var1, var2, p_value_adj, notation_adj) %>% t()

p_4_1 <- plot_box(df_plot=df_meta %>% filter(Vaccine!="Unvaccinated"), x_var="lineage_sim", y_var="days_since_last_dose", color_var="vaccine_doses", y_lab="Days since last dose", x_lab="Lineage")
p_4_1 <- p_4_1 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses[1:4])

p_4_1_p <- p_4_1 +
	geom_signif(y_position = 350+seq(0,4)*10, xmin = c(0.7, 0.7, 0.7, 1.1, 1.1)+1, xmax = c(0.9, 1.1, 1.3, 0.9, 1.3)+1, annotation = c("**", "**", "**", "**", "**"), color="black", vjust=0.65, tip_length = 0.01)
ggsave("../results/days_since_last_dose_by_lineage.pdf", width=8, height=6)

# Age, chronic disease, condition
df_meta$Age <- as.numeric(df_meta$"Age (yr)")
df_test <- cal_wilc_test(df_meta, "Age", genes=NA)
highlight_diff(df_test) %>% filter(check_p_adj) %>% arrange(var1, var2) %>% select(var1, var2, p_value_adj, notation_adj) %>% t()

source("./helper/plot_reg.R")
df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, df_meta %>% select(sample, Age), "sample")
p5_01 <- plot_reg(df_meta, "Age", "Ct_value", y_lab="Ct value")
p5_02 <- plot_reg(df_plot_n_gene_meta_adj %>% filter(gene=="Full genome"), "Age", "n_per_kb_adj", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Age")
p5_out <- p5_01/p5_02
ggsave("../results/cor_age_vs_n_per_kd_adj.pdf", width=8, height=8)


p5_1 <- plot_box(df_plot=df_meta %>% filter(Vaccine=="Unvaccinated"), x_var="vaccine_doses", y_var="Age", color_var="lineage_sim", y_lab="Age", x_lab="Vaccine")+xlab("")
p5_1 <- p5_1 + scale_color_manual(name="Lineage", values=colors_lineage)

p5_2 <- plot_box(df_plot=df_meta %>% filter(lineage_sim%in%names(colors_lineage)[4:5]), x_var="lineage_sim", y_var="Age", color_var="vaccine_doses", y_lab="Age", x_lab="Lineage")+xlab("")
p5_2 <- p5_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses) + 
	# scale_x_discrete(guide = guide_axis(n.dodge = 2))+
	NULL

p_5_out <- (p5_1+ggtitle("A"))/(p5_2+ggtitle("B"))+
	# plot_layout(guides='collect') & theme(legend.position='bottom') & guides(col = guide_legend(nrow = 3))+
	NULL
ggsave("../results/Age_by_group_no_signf.pdf", width=8, height=8)

names(df_meta)
sort(table(tolower(df_meta$`chronic disease`)))
df_meta$chronic_disease <- NA
sum(tolower(df_meta$`chronic disease`) %in% c("nil\r\a", "nil"))
sum(is.na(df_meta$`chronic disease`))
df_meta$chronic_disease[which(tolower(df_meta$`chronic disease`) %in% c("nil\r\a", "nil"))] <- "No"
df_meta$chronic_disease[which(!tolower(df_meta$`chronic disease`) %in% c("nil\r\a", "nil"))] <- "Yes"
df_meta$chronic_disease[is.na(df_meta$`chronic disease`)] <- NA

df_meta$`chronic disease`[which(df_meta$chronic_disease=="No")]
df_meta$`chronic disease`[which(df_meta$chronic_disease=="Yes")]

# df_meta$chronic_disease[grep("HT", df_meta$`chronic disease`)] <- "HT"
# df_meta$chronic_disease[grep("DM", df_meta$`chronic disease`)] <- "DM"
# df_meta$chronic_disease[grep("hyperlipidaemia", tolower(df_meta$`chronic disease`))] <- "hyperlipidaemia"

apply(table(df_meta$chronic_disease[df_meta$Vaccine=="Unvaccinated"], df_meta$lineage_sim[df_meta$Vaccine=="Unvaccinated"]), 2, function(x){x[2]/sum(x)})
chisq.test(table(df_meta$chronic_disease[df_meta$Vaccine=="Unvaccinated"], df_meta$lineage_sim[df_meta$Vaccine=="Unvaccinated"]))

df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, df_meta %>% select(sample, chronic_disease), "sample")
df_plot_n_gene_meta_adj_filter <- df_plot_n_gene_meta_adj %>% filter(gene=="Full genome") %>% filter(Vaccine=="Unvaccinated") %>% filter(!is.na(chronic_disease))

p_5_3 <- plot_box(df_plot=df_plot_n_gene_meta_adj_filter, x_var="chronic_disease", y_var="n_per_kb_adj", color_var="chronic_disease", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="chronic disease")+facet_grid(rows=vars(vaccine_doses), cols=vars(lineage_sim), scale="free_x", space="free")
p_5_4 <- plot_box(df_plot=df_meta[df_meta$Vaccine=="Unvaccinated",], x_var="chronic_disease", y_var="Ct_value", color_var="chronic_disease", y_lab="Ct value", x_lab="chronic disease")+facet_grid(rows=vars(vaccine_doses), cols=vars(lineage_sim), scale="free_x", space="free")

df_input <- df_plot_n_gene_meta_adj_filter
df_groups <- df_input %>% select(lineage_sim) %>% unique()
var <- "n_per_kb_adj"
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	# i=2
	lineage_i <- df_groups$lineage_sim[i]
	df_tmp <- df_input %>% filter(lineage_sim==lineage_i)
	df_tmp <- df_tmp %>% mutate(x_grps=paste(vaccine_doses, chronic_disease, lineage_sim, sep="_"))
	mat_pairs <- combn(unique(df_tmp$x_grps), 2)

	df_rst <- apply(mat_pairs,2,function(y) {
		var1 <- y[1]
		var2 <- y[2]
		value1 <- df_tmp[[var]][df_tmp$x_grps==var1]
		value2 <- df_tmp[[var]][df_tmp$x_grps==var2]
		
		rst <- wilcox.test(value1, value2)
		tibble(var1=var1, var2=var2, var1_mean=mean(value1,na.rm=T), var2_mean=mean(value2,na.rm=T), median_var1=median(value1,na.rm=T), median_var2=median(value2,na.rm=T), p_value=rst$p.value)
	})
	df_rst <- bind_rows(df_rst)
	df_rst <- df_rst %>% mutate(gene=lineage_i)
})
df_wilc_test <- bind_rows(df_wilc_test)
df_wilc_test$p_value_adj <- p.adjust(df_wilc_test$p_value, "BH")
df_wilc_test$var1

highlight_diff <- function(df){
	var1_split <- strsplit(df$var1, "_", fixed=T)
	var2_split <- strsplit(df$var2, "_", fixed=T)
	check1 <- sapply(var1_split, function(x)x[1]) == sapply(var2_split, function(x)x[1]) # Vaccine
	check2 <- sapply(var1_split, function(x)x[2]) == sapply(var2_split, function(x)x[2]) # chronic_disease
	check3 <- sapply(var1_split, function(x)x[3]) == sapply(var2_split, function(x)x[3]) # lineage
	
	df$same_vaccine <- check1
	df$same_lineage <- check3
	df$within_group <- df$same_vaccine & df$same_lineage
	medians_low <- sapply(seq_len(nrow(df)), function(i){
		tmp <- abs(c(df$median_var1[i], df$median_var2[i]))
		tmp[which.min(tmp)]
	})
	df$median_diff <- abs(df$median_var1-df$median_var2)/medians_low
	df$median_diff_over_10percent <- df$median_diff>0.1
	df$check_three <- df$median_diff_over_10percent & df$within_group & (df$p_value<0.05)
	
	df$check_p_adj <- df$median_diff_over_10percent & df$within_group & (df$p_value_adj<0.05)

	df$notation <- NA
	df$notation[df$p_value<0.1] <- "^"
	df$notation[df$p_value<0.05] <- "*"
	df$notation[df$p_value<0.01] <- "**"

	df$notation_adj <- NA
	df$notation_adj[df$p_value_adj<0.1] <- "^"
	df$notation_adj[df$p_value_adj<0.05] <- "*"
	df$notation_adj[df$p_value_adj<0.01] <- "**"

	df %>% arrange(gene, check_p_adj, check_three, same_vaccine, same_lineage, var1, var2)
}
df_wilc_test <- highlight_diff(df_wilc_test)
df_wilc_test %>% filter(check_p_adj) %>% select(var1, var2, median_var1, median_var2, median_diff, p_value_adj, notation_adj) %>% t()

df_meta_sub <- df_meta %>% filter(lineage_sim %in% names(colors_lineage)[1:2])
df_plot_n_gene_meta_adj_filter_sub <- df_plot_n_gene_meta_adj_filter %>% filter(lineage_sim %in% names(colors_lineage)[1:2])
wilcox.test(df_meta_sub$Ct_value[df_meta_sub$chronic_disease=="Yes"], df_meta_sub$Ct_value[df_meta_sub$chronic_disease=="No"])
wilcox.test(df_plot_n_gene_meta_adj_filter$n_per_kb[df_plot_n_gene_meta_adj_filter$chronic_disease=="Yes"], df_plot_n_gene_meta_adj_filter$n_per_kb[df_plot_n_gene_meta_adj_filter$chronic_disease=="No"])
wilcox.test(df_plot_n_gene_meta_adj_filter$n_per_kb_adj[df_plot_n_gene_meta_adj_filter$chronic_disease=="Yes"], df_plot_n_gene_meta_adj_filter$n_per_kb_adj[df_plot_n_gene_meta_adj_filter$chronic_disease=="No"])

mean(df_plot_n_gene_meta_adj_filter$n_per_kb_adj[df_plot_n_gene_meta_adj_filter$chronic_disease=="Yes"])
mean(df_plot_n_gene_meta_adj_filter$n_per_kb_adj[df_plot_n_gene_meta_adj_filter$chronic_disease=="No"])
# cases with chronic disease infected by B.1.36.* seem to have higher adj num of iSNVs, also, the proportion of chronic-disease cases in non-VOC groups are higher. In this case, the adj num of iSNVs in non-VOC groups should be higher, however we observed the opposite (Figure 2 of the paper), this further suggests that VOC cases have higher adj num of iSNVs than non-VOC cases.

df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, df_meta %>% select(sample, Sex), "sample")
df_plot_n_gene_meta_adj_filter <- df_plot_n_gene_meta_adj %>% filter(gene=="Full genome") %>% filter(Vaccine=="Unvaccinated")

wilcox.test(df_meta$Ct_value[df_meta$Sex=="Male"], df_meta$Ct_value[df_meta$Sex=="Female"])
wilcox.test(df_plot_n_gene_meta_adj_filter$n_per_kb_adj[df_plot_n_gene_meta_adj_filter$Sex=="Male"], df_plot_n_gene_meta_adj_filter$n_per_kb_adj[df_plot_n_gene_meta_adj_filter$Sex=="Female"])
# no diff in Ct/num_isnvs_adj by gender