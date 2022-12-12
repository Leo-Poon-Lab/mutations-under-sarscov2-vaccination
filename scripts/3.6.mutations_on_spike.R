library(tidyverse)
library(trackViewer)

# read data
df_spike_region <- read_tsv("../data/spike_region.tsv")

df_snvs_meta_add_qc <- read_csv("../results/df_snvs_meta_add_qc_bam_adj.csv", guess_max=600000)
df_snvs_meta_add_qc$Vaccine <- gsub("BioNTech", "Comirnaty", df_snvs_meta_add_qc$Vaccine)
df_snvs_meta_add_qc$Vaccine <- gsub("Sinovac", "CoronaVac", df_snvs_meta_add_qc$Vaccine)
df_snvs_meta_add_qc$vaccine_doses <- paste0(df_snvs_meta_add_qc$Vaccine, "(Doses=", df_snvs_meta_add_qc$Doses, ")\n")
df_snvs_meta_add_qc$vaccine_doses[df_snvs_meta_add_qc$Vaccine=="Unvaccinated"] <- "Unvaccinated\n"
df_snvs_meta_add_qc <- df_snvs_meta_add_qc %>% mutate(Group=paste(vaccine_doses, lineage_sim))

df_meta <- read_csv("../results/df_samples.csv", guess_max=100000)
df_meta <- df_meta %>% filter(lineage_sim != "22B (Omicron, BA.5.*)")
df_meta$Vaccine <- gsub("BioNTech", "Comirnaty", df_meta$Vaccine)
df_meta$Vaccine <- gsub("Sinovac", "CoronaVac", df_meta$Vaccine)
df_meta$vaccine_doses <- paste0(df_meta$Vaccine, "(Doses=", df_meta$Doses, ")\n")
df_meta$vaccine_doses[df_meta$Vaccine=="Unvaccinated"] <- "Unvaccinated\n"

df_tmp <- df_meta %>% mutate(Group=paste(vaccine_doses, lineage_sim)) %>% group_by(Group) %>% summarise(N=n(), Group_more=paste0(Group, "\n(N=", N, ")")) %>% unique()

df_snvs_meta_add_qc$Group_more <- factor(df_snvs_meta_add_qc$Group, levels=df_tmp$Group, labels=df_tmp$Group_more)

# spike non-synonymous mutation
unique(df_snvs_meta_add_qc$lineage_sim)
df_spike_nonsyn_muts <- df_snvs_meta_add_qc %>% filter(gene=="S") %>% filter((grepl("Omicron", lineage_sim))|(grepl("Delta", lineage_sim))) %>% filter(effect=="missense_variant") %>% group_by(Group_more, mut_aa, pos_aa) %>% summarise(N=n()) %>% mutate(pos_aa=as.numeric(gsub("/1273", "", pos_aa, fixed=T)),) %>% arrange(pos_aa, mut_aa, N, Group_more)

df_spike_nonsyn_muts$spike_regions <- sapply(df_spike_nonsyn_muts$pos_aa, function(x){
	df_spike_region$Region[df_spike_region$start<=x & df_spike_region$stop>=x]
})


# Non-synonymous mutations on Spike region
df_spike_region_nogap <- df_spike_region %>% filter(!grepl("-", Region))

features <- GRanges("chr1", IRanges(df_spike_region_nogap$start, 
                                    width=(df_spike_region_nogap$stop-df_spike_region_nogap$start+1),
                                    names=df_spike_region_nogap$Region))
colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")
features$fill <- colors
df_spike_nonsyn_muts$color <- factor(df_spike_nonsyn_muts$spike_regions, levels=df_spike_region_nogap$Region, labels=colors)

groups_unique <- sort(as.character(unique(df_spike_nonsyn_muts$Group_more)))

muts_list <- list()
features_list <- list()

for (group_i in groups_unique) {
	df_tmp <- df_spike_nonsyn_muts %>% filter(Group_more==group_i, N>1) # only shows occurrence more than once
	if(nrow(df_tmp)>1){
		SNP <- df_tmp$pos_aa
		Mutations <- GRanges("chr1", IRanges(SNP, width=1, names=gsub("p.", "", df_tmp$mut_aa, fixed=T)))
		Mutations$score <- df_tmp$N
		Mutations$label <- as.character(df_tmp$N)
		Mutations$label.col <- "black"
		# Mutations$color <- as.character(df_tmp$color)
		muts_list[[group_i]] <- Mutations
		features_list[[group_i]] <- features
	}	
}

source("./helper/lolliplot_mod.R")
idx1 <- rev(grep("Delta", names(muts_list)))
pdf("../results/spike_muts_Delta.pdf", width=10, height=10)
lolliplot_mod(muts_list[idx1], features_list[idx1], yaxis=FALSE, cex=0.8)
dev.off()

idx2 <- rev(grep("Omicron", names(muts_list)))
pdf("../results/spike_muts_Omicron.pdf", width=10, height=20)
lolliplot_mod(muts_list[idx2], features_list[idx2], yaxis=FALSE)
dev.off()

idx_all <- c(idx2, idx1)
pdf("../results/spike_muts_all.pdf", width=10, height=20)
lolliplot_mod(muts_list[idx_all], features_list[idx_all], yaxis=TRUE)
dev.off()

df_spike_nonsyn_muts_summary <- df_spike_nonsyn_muts %>% group_by(spike_regions, Group_more) %>% summarise(N=sum(N)) %>% arrange(spike_regions, Group_more, N)

mut_aa_int <- df_spike_nonsyn_muts %>% filter(spike_regions=="RBD") %>% .$mut_aa
df_snvs_meta_add_qc %>% filter(mut_aa %in% mut_aa_int) %>% select(Group, sample, mut_aa) %>% arrange(Group, sample)

df_spike_nonsyn_muts %>% filter(N>1) %>% arrange(Group_more, pos_aa, desc(N)) %>% write_csv("../results/df_spike_nonsyn_muts.csv")
