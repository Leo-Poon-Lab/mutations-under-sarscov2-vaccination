library(tidyverse)
library(lubridate)
library(patchwork)
library(ggsignif)
source("./helper/df_test.R")
source('./helper/plot_box.R')

df_meta <- read_csv("../results/df_samples.csv", guess_max=100000)
df_meta <- df_meta %>% filter(lineage_sim != "22B (Omicron, BA.5.*)")
lineages_all <- sort(unique(df_meta$lineage_sim))
colors_lineage=rev(c("#4daf4a", "#984ea3", "#ff7f00", "#f781bf", "#666666"))
names(colors_lineage) <- lineages_all
vaccine_all <- sort(unique(df_meta$Vaccine))
colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine)=vaccine_all
colors_vaccine_new <- colors_vaccine
names(colors_vaccine_new)=c("Comirnaty", "CoronaVac", "Unvaccinated")

load("../results/df_plot_n_gene_adj.rdata")
df_plot_n_gene_meta_adj$lineage_sim <- factor(df_plot_n_gene_meta_adj$lineage_sim, levels = names(colors_lineage))
df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, df_meta %>% select(sample, Doses), "sample")
df_plot_n_gene_meta_adj$Vaccine <- factor(df_plot_n_gene_meta_adj$Vaccine, levels=names(colors_vaccine), labels=names(colors_vaccine_new))
df_plot_n_gene_meta_adj$vaccine_doses <- paste0(df_plot_n_gene_meta_adj$Vaccine, "\n(Doses=", df_plot_n_gene_meta_adj$Doses, ")")
df_plot_n_gene_meta_adj$vaccine_doses[df_plot_n_gene_meta_adj$Vaccine=="Unvaccinated"] <- "Unvaccinated"

df_meta$lineage_sim <- factor(df_meta$lineage_sim, levels = names(colors_lineage))
# collection lag
df_meta$detection_lag <- as.numeric(df_meta$`Report date` - ymd(df_meta$`Onset date`))
df_meta$collection_lag <- as.numeric(df_meta$collection_date - ymd(df_meta$`Onset date`))
df_meta$collection_lag[df_meta$collection_lag>100] <- NA # one ourlier
df_meta$collection_lag[df_meta$collection_lag<0] <- NA # one ourlier

df_meta$vaccine_doses <- paste0(df_meta$Vaccine, "\n(Doses=", df_meta$Doses, ")")
df_meta$vaccine_doses[df_meta$Vaccine=="Unvaccinated"] <- "Unvaccinated"

colors_vaccine_doses = c("#a65628", "#844420", "#7570b3", "#5d598f", "#999999")
names(colors_vaccine_doses)=sort(unique(df_meta$vaccine_doses))
# names(df_meta)

# Viral loads between groups
df_test <- cal_wilc_test(df_meta, "Ct_value", genes=NA)
highlight_diff(df_test) %>% filter(check_three) %>% select(var1, var2, notation) %>% t()

p1_1 <- plot_box(df_plot=df_meta, x_var="vaccine_doses", y_var="Ct_value", color_var="lineage_sim", y_lab="Ct value", x_lab="Vaccine")
p1_1 <- p1_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p1_1_p <- p1_1+ # unvaccinated
   # ggtitle("A")+
   geom_signif(y_position = 28.5+seq(0,3)/1.5, xmin = c(0.7, 0.7, 0.7, 0.7)+4, xmax = c(0.85, 1, 1.15, 1.3)+4, annotation = c("**", "**", "**", "**"), color="black", vjust=0.65, tip_length = 0.01)

p1_2 <- plot_box(df_plot=df_meta, x_var="lineage_sim", y_var="Ct_value", color_var="vaccine_doses", y_lab="Ct value", x_lab="Lineage")
p1_2 <- p1_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
p1_2_p <- p1_2

p_1_out <- (p1_1_p+theme(legend.position='bottom'))/(p1_2_p+theme(legend.position='bottom'))+
	# plot_layout(guides='collect') & theme(legend.position='bottom') & guides(col = guide_legend(nrow = 3))+
	NULL
ggsave("../results/Ct_Value_by_group.pdf", width=8, height=8)

# specimen_type
sort(table(df_meta$specimen_type))
df_meta$specimen_type[df_meta$specimen_type=="specimen"] <- NA
df_meta$specimen_type[df_meta$specimen_type=="Swab"] <- NA

df_meta$specimen_type_fct <- factor(df_meta$specimen_type, rev(names(sort(table(df_meta$specimen_type)))))
table(df_meta$specimen_type_fct, df_meta$vaccine_doses)
table(df_meta$specimen_type_fct, df_meta$vaccine_doses, df_meta$lineage_sim)

p2_1 <- ggplot(df_meta)+
	geom_histogram(aes(x=specimen_type_fct, fill=vaccine_doses), stat="count")+
	scale_fill_manual(name="Vaccine", values=colors_vaccine_doses) +
	# scale_x_discrete(guide = guide_axis(n.dodge = 4))+
	theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
	xlab("Specimen type")+
	ylab("Count")+
	NULL
p2_2 <- ggplot(df_meta)+
	geom_histogram(aes(x=specimen_type_fct, fill=lineage_sim), stat="count")+
	scale_fill_manual(name="Lineage", values=colors_lineage)+
	# scale_x_discrete(guide = guide_axis(n.dodge = 4))+
	theme_classic()+
	theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
	xlab("Specimen type")+
	ylab("Count")+
	facet_wrap(vars(vaccine_doses), ncol=1, scale="free_y")+
	NULL
# p_2_out <- (p2_1+theme(legend.position='bottom'))/(p2_2+theme(legend.position='bottom'))
p_2_out <- (p2_2+theme(legend.position='bottom'))
ggsave("../results/Specimen_type_by_group.pdf", width=8, height=10)
## Apparently, there are some difference in lineage/vaccine composition among specimen types, this is inevitable, because as the pandemic progress, the sample collection methods are changing. However, we can test whether different specimen type is associated with different number of iSNVs. This can serve as a reference to estimated the bias introduced by sample collection methods.
df_plot_n_gene_meta_adj <- left_join(df_plot_n_gene_meta_adj, df_meta %>% select(sample, specimen_type), "sample")
	scale_fill_manual(name="Lineage", values=colors_lineage)
df_plot_n_gene_meta_adj_filter <- df_plot_n_gene_meta_adj %>% filter(gene=="Full genome") %>% filter(Vaccine=="Unvaccinated") %>% filter(!is.na(specimen_type))
p_2_3 <- plot_box(df_plot=df_plot_n_gene_meta_adj_filter, x_var="specimen_type", y_var="n_per_kb_adj", color_var="lineage_sim", y_lab="Number of iSNVs per Kb (adjusted)", x_lab="Specimen type")+facet_grid(rows=vars(vaccine_doses), cols=vars(lineage_sim), scale="free", space="free")
p_2_3 <- p_2_3 + scale_color_manual(name="Lineage", values=colors_lineage)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='bottom')
ggsave("../results/num_isnvs_by_Specimen_type_by_group.pdf", width=12, height=8)

df_groups <- df_plot_n_gene_meta_adj_filter %>% select(lineage_sim) %>% unique()
df_input <- df_plot_n_gene_meta_adj_filter
df_input$specimen_type <- gsub(" ", ".", df_input$specimen_type)
var <- "n_per_kb_adj"
df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
	# i=1
	lineage_i <- df_groups$lineage_sim[i]
	df_tmp <- df_input %>% filter(lineage_sim==lineage_i)
	df_tmp <- df_tmp %>% mutate(x_grps=paste(Vaccine, specimen_type, lineage_sim, sep="_"))
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
df_wilc_test %>% filter(check_p_adj) %>% select(var1, var2, notation_adj) %>% t()

sort(df_wilc_test$p_value_adj)

df_tmp <- df_plot_n_gene_meta_adj_filter %>% filter(lineage_sim==names(colors_lineage)[1])
kruskal.test(x=df_tmp$n_per_kb_adj, g=df_tmp$specimen_type)


# Collection timepoints during acute infection
df_test <- cal_wilc_test(df_meta, "Ct_value", genes=NA)
highlight_diff(df_test) %>% filter(check_three) %>% select(var1, var2, notation) %>% t()

p1_1 <- plot_box(df_plot=df_meta, x_var="vaccine_doses", y_var="Ct_value", color_var="lineage_sim", y_lab="Ct value", x_lab="Vaccine")
p1_1 <- p1_1 + scale_color_manual(name="Lineage", values=colors_lineage)
p1_1_p <- p1_1+ # unvaccinated
   # ggtitle("A")+
   geom_signif(y_position = 28.5+seq(0,3)/1.5, xmin = c(0.7, 0.7, 0.7, 0.7)+4, xmax = c(0.85, 1, 1.15, 1.3)+4, annotation = c("**", "**", "**", "**"), color="black", vjust=0.65, tip_length = 0.01)

p1_2 <- plot_box(df_plot=df_meta, x_var="lineage_sim", y_var="Ct_value", color_var="vaccine_doses", y_lab="Ct value", x_lab="Lineage")
p1_2 <- p1_2 + scale_color_manual(name="Vaccine", values=colors_vaccine_doses) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
p1_2_p <- p1_2


#
df_meta %>% filter(is.na(`onset to report (days)`))

# Use of hospital resources

# Age, chronic disease, condition