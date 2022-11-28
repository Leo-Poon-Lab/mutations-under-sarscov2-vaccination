library(tidyverse)
library(ggsci)
library(patchwork)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")
colors_lineage=c("#e41a1c", "#33a02c", "#1f78b4", "#ff7f00", "#f781bf", "#666666") 
names(colors_lineage) <- c("Alpha", "Delta", "Omicron", "B.1.36", "B.1.36.27", "B.1.1.63")
colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine)=c("BioNTech", "Sinovac", "Non-vaccinated")

df_plot <- read_csv("../results/df_snvs_meta_add_qc_ivar_clean.csv", guess_max=600000)
df_plot <- df_plot %>% filter(!X2 %in% c(15494, 15489, 25381, 10194, 22422)) # excluding primer/homoplasy sites, https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28420-7/MediaObjects/41467_2022_28420_MOESM1_ESM.pdf
df_plot$sample <- df_plot$Sample

#### mutations by gene and epitopes
df_epitope <- readxl::read_excel("../data/df_t_cell_epitope.xlsx")
df_epitope$found_in <- sapply(strsplit(df_epitope$`Reference(s)*`, ", "), length)
source("./helper/annotate_gene.r")

df_orf <- read_csv("../data/ORF_SCoV2.csv")
df_orf$sequence

df_baseline <- lapply(seq_len(nrow(df_orf)), function(i) {
	x <- df_orf[i,]
	df_tmp <- df_epitope %>% filter(toupper(Antigen)==toupper(x[[1]])) %>% filter(found_in>=2) 
	pos_tmp <- 1:((x[[3]]-x[[2]]+1)/3)
	check_epitopes <- sapply(pos_tmp, function(y){
		any((df_tmp$Start<=y) & (df_tmp$End>=y))
	})
	prop=sum(check_epitopes)/length(pos_tmp)
	tibble(gene=x[[1]], prop=prop, n=sum(check_epitopes), length=length(pos_tmp))
})
df_baseline <- bind_rows(df_baseline)
df_baseline_orf1ab <- df_baseline %>% filter(grepl("nsp", gene)) %>% summarise(gene="ORF1ab", n = sum(n), length=sum(length)) %>% mutate(prop=n/length)
df_baseline <- bind_rows(df_baseline_orf1ab, df_baseline %>% filter(!grepl("nsp", gene)))
df_baseline$gene <- factor(df_baseline$gene, levels=df_baseline$gene)


df_plot$effect_sim <- df_plot$effect
df_plot$effect_sim[grepl("stream", df_plot$effect_sim)] <- "UTR"
df_plot$effect_sim[grepl("stop", df_plot$effect_sim)] <- "missense_variant"
df_plot$effect_sim[grepl("missense", df_plot$effect_sim)] <- "missense_variant"
df_plot$effect_sim[grepl("synonymous", df_plot$effect_sim)] <- "synonymous variant"
df_plot$effect_sim <- gsub("_", " ", df_plot$effect_sim)
table(df_plot$effect_sim)

df_plot <- bind_cols(df_plot, get_orf(df_plot$X2))

##### annotate epitope https://www.sciencedirect.com/science/article/pii/S1931312821002389?via%3Dihub#app2

table(df_plot$gene, df_plot$gene_nsp)
df_plot$epitope_type <- NA
df_plot$epitope_found_in <- NA

sapply(seq_len(nrow(df_plot)), function(i) {
	# print(i)
	gene_i <- toupper(df_plot$gene_nsp[i])
	pos_nsp_i <- df_plot$pos_orf[i]
	df_epitope_i <- df_epitope %>% filter(toupper(Antigen)==gene_i & Start<=pos_nsp_i & End >=pos_nsp_i)
	if(nrow(df_epitope_i)==0){
		return("")
	} else {
		df_plot$epitope_type[i] <<- paste0(sort(unique(df_epitope_i$`Restriction`)), collapse=" and ")
		df_plot$epitope_found_in[i] <<- max(df_epitope_i$found_in)
		return("")
	}	
})

color_epitope <- pal_aaas()(3)
names(color_epitope) <- sort(unique(df_plot$epitope_type))

# TODO
#### The mutations on spike
df_plot$pos_aa <- as.numeric(sapply(strsplit(df_plot$pos_aa, "\\/"), function(x) {x[1]}))
assign_spike_aa_pos <- function(pos_spike) {
	df_s_r <- read_tsv("./helper/spike_region.tsv")
	out <- sapply(pos_spike, function(x) {
		tmp <- df_s_r$Region[df_s_r$start<=x & df_s_r$stop>=x]
		if(length(tmp)!=1){return(NA)}else{return(tmp)}
	})	
	factor(out, levels=df_s_r$Region)
	# out
}

df_plot2 <- df_plot %>% filter(gene=="S" & effect_sim!="UTR") %>% mutate(spike_region=assign_spike_aa_pos(pos_aa)) %>% mutate(mut_aa=gsub("p\\.", "", mut_aa)) %>% arrange(X2) %>% mutate(mut_aa=factor(mut_aa, levels=unique(mut_aa)))%>% group_by(X2, mut_aa, Vaccine, lineage_sim, spike_region, effect_sim, epitope_type, epitope_found_in) %>% summarise(n=n(), n_grp=sum(df_plot$Vaccine==Vaccine[1] & df_plot$lineage_sim==lineage_sim[1]), pect=n/n_grp) %>% ungroup()
df_plot2$epitope_type[df_plot2$epitope_found_in<2] <- NA
mut_common <- names(table(df_plot2$mut_aa)[table(df_plot2$mut_aa)>1]) # commonly seen in more than two group
df_plot2 <- df_plot2 %>% filter(mut_aa %in% mut_common)

tmp <- as.character(unique(df_plot2$spike_region))
names(tmp) <- tmp
tmp[grepl("-", tmp)] <- ""

facet_labeller <- function(variable, value) {
	tmp
}
df_plot2$epitope_type <- factor(df_plot2$epitope_type, levels=names(color_epitope))

df_plot2 %>% filter(mut_aa=="F133V")
df_plot %>% filter(mut_aa=="p.F133V")
df_plot2 %>% filter(mut_aa=="L241L")
df_plot %>% filter(mut_aa=="p.L241L")

p_spike <- ggplot(df_plot2)+
 	geom_tile(aes(x=mut_aa, y=paste(Vaccine, lineage_sim), fill=epitope_type))+
	facet_grid(cols=vars(spike_region), scales="free_x",space="free_x", labeller=labeller(spike_region=as_labeller(facet_labeller)))+
	theme_classic()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
	xlab("Position")+
	ylab("Group")+
	# scale_fill_viridis_c(name="Proportion. of\nsamples")+
	scale_fill_manual(name="T cell epitope", values=color_epitope, na.value="grey", drop=F)+
	# scale_color_manual(name="Mutation", values=c("dark blue", "yellow"), na.value="grey")+
	ggtitle("Spike")+
	NULL
ggsave("../results/mut_isnvs_by_vaccine_spike.pdf", width=12, height=4)
 
sum(check_epitopes_spike)/length(pos_spike)
df_tmp1 <- df_plot2 %>% filter(!duplicated(df_plot2$X2))
m1 <- matrix(c(df_baseline$n[df_baseline$gene=="S"], df_baseline$length[df_baseline$gene=="S"], sum(!is.na(df_tmp1$epitope_type)), length(df_tmp1$epitope_type)), ncol=2)
chisq.test(m1)

# df_plot %>% filter(gene=="S" & effect_sim!="UTR") %>% mutate(spike_region=assign_spike_aa_pos(pos_aa)) %>% arrange(X2) %>% mutate(mut_nt=factor(paste0(X4, X2, X5), levels=unique(paste0(X4, X2, X5)))) %>% group_by(mut_nt, mut_aa, Vaccine, lineage_sim, spike_region, effect_sim) %>% summarise(n=n(), n_grp=sum(df_plot$Vaccine==Vaccine[1] & df_plot$lineage_sim==lineage_sim[1]), pect=n/n_grp) %>% 
# 	ggplot()+
#  	geom_tile(aes(x=mut_nt, y=paste(Vaccine, lineage_sim), fill=round(pect,2)))+
# 	facet_grid(effect_sim ~ spike_region, scales="free_x",space="free_x", labeller=labeller(spike_region=as_labeller(facet_labeller)))+
# 	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
# 	xlab("Position")+
# 	ylab("Group")+
# 	scale_fill_viridis_c(name="Proportion. of\nsamples")+
# 	NULL
# ggsave("../results/mut_isnvs_by_vaccine_spike_nt.pdf", width=12)

### The mutations on Other genes
df_plot3 <- df_plot %>% filter(gene=="ORF1ab" & effect_sim!="UTR") %>% mutate(mut_aa=gsub("p\\.", "", mut_aa)) %>% arrange(X2) %>% mutate(mut_aa=factor(mut_aa, levels=unique(mut_aa)))%>% group_by(X2, mut_aa, Vaccine, lineage_sim, effect_sim, epitope_type, epitope_found_in) %>% summarise(n=n(), n_grp=sum(df_plot$Vaccine==Vaccine[1] & df_plot$lineage_sim==lineage_sim[1]), pect=n/n_grp) %>% ungroup()
df_plot3$epitope_type[df_plot3$epitope_found_in<2] <- NA
df_plot3 <- df_plot3 %>% filter(!grepl("Ter", mut_aa))
mut_common <- names(table(df_plot3$mut_aa)[table(df_plot3$mut_aa)>1]) # commonly seen in more than two group
df_plot3 <- df_plot3 %>% filter(mut_aa %in% mut_common)

p_orf1ab <- ggplot(df_plot3)+
 	geom_tile(aes(x=mut_aa, y=paste(Vaccine, lineage_sim), fill=epitope_type))+
	# facet_grid(effect_sim ~ ., scales="free_x",space="free_x")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
	xlab("Position")+
	ylab("Group")+
	# scale_fill_viridis_c(name="Proportion. of\nsamples")+
	scale_fill_manual(name="T cell epitope", values=color_epitope, na.value="grey")+
	ggtitle("ORF1ab")+
  	# theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	coord_flip()+
	NULL
ggsave("../results/mut_isnvs_by_vaccine_orf1ab.pdf", width=4, height=30)

df_tmp2 <- df_plot3 %>% filter(!duplicated(df_plot3$X2))
df_baseline$prop[df_baseline$gene=="ORF1ab"]
sum(!is.na(df_tmp2$epitope_type))/length(df_tmp2$epitope_type)
m2 <- matrix(c(df_baseline$n[df_baseline$gene=="ORF1ab"], df_baseline$length[df_baseline$gene=="ORF1ab"], sum(!is.na(df_tmp2$epitope_type)), length(df_tmp2$epitope_type)), ncol=2)
chisq.test(m2)

p_out <- p_orf1ab/p_spike + plot_layout(guides = "collect")
ggsave("../results/mut_isnvs_by_vaccine_combine.pdf", width=10, height=8)
save_pptx("../results/mut_isnvs_by_vaccine_combine.pptx", width=10, height=8)

# df_plot %>% filter(grepl("D6304G", mut_aa)) %>% .$gene_nsp
# df_plot %>% filter(grepl("D6304G", mut_aa)) %>% .$pos_orf
## across genes
df_plot_n <- df_plot 
breaks_geneome <- seq(1,30001, 1000)
n_breaks <- length(breaks_geneome)

df_plot_n$epitope_type[df_plot_n$epitope_found_in<2] <- NA
df_plot_n$check_epitope <- !is.na(df_plot_n$epitope_type)
df_tmp <- df_plot_n %>% arrange(X2) %>% mutate(gene=factor(gene, levels=unique(gene))) %>% select(sample, gene, check_epitope, Vaccine, lineage_sim) %>% group_by(sample, gene, check_epitope, Vaccine, lineage_sim) %>% summarise(n=n()) %>% ungroup() %>% group_by(gene, check_epitope, Vaccine, lineage_sim) %>% summarise(n_sum=sum(n))

df_tmp$group <- paste0(df_tmp$Vaccine, " ", df_tmp$lineage_sim)

ggplot(df_tmp) +
	geom_col(aes(x=gene, y=n_sum, fill=check_epitope), position="fill", alpha=0.8)+
	geom_segment(aes(x=as.numeric(gene)-0.4, xend=as.numeric(gene)+0.4, y=prop, yend=prop), color="black", data=df_baseline, alpha=0.1)+
	geom_segment(aes(x=as.numeric(gene)-0.4, xend=as.numeric(gene)+0.4, y=prop, yend=prop), color="black", data=df_baseline, linetype="dotted",size=0.3)+
	facet_wrap(vars(group), ncol=2)+
	coord_flip()+
	scale_fill_manual(name = "T-cell epitope", values=c("#ef8a62", "#67a9cf"))+
	theme_classic()+
	ylab("Relative frequency")+
	xlab("Gene")+
	NULL
ggsave("../results/tmp.pdf", height=10, width=8)


## T cell epitopes frequency
source('./helper/plot_box.R')
df_tmp <- df_plot_n %>% select(sample, gene, check_epitope, Vaccine, lineage_sim) %>% group_by(sample, gene, Vaccine, lineage_sim) %>% summarise(n_epitope=sum(check_epitope), n_non_epitope=sum(!check_epitope), n_total=n()) %>% ungroup() %>% mutate(freq_epitope = n_epitope/n_total)
df_tmp <- bind_rows(df_tmp %>% mutate(gene="Full genome"), df_tmp)
p1_2 <- plot_box(df_plot=df_tmp %>% filter(gene=="Full genome"), x_var="Vaccine", y_var="freq_epitope", color_var="lineage_sim", y_lab="freq_epitope")
p1_2 <- p1_2 + scale_color_manual(name="Vaccine", values=colors_lineage) 
ggsave("../results/tmp.pdf", height=6, width=8)

# ggsave("../results/Num_isnvs_by_vaccine_2_more.pdf", width=6, height=4)
# save_pptx("../results/Num_isnvs_by_vaccine_2_more.pptx", width=6, height=4)


# age, patient status
