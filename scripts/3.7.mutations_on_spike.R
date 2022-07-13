library(tidyverse)
library(trackViewer)

# read data
df_spike_region <- read_tsv("../data/spike_region.tsv")

df_snvs_meta_add_qc <- read_csv("../results/df_snvs_meta_add_qc_bam.csv", guess_max=600000)
df_snvs_meta_add_qc$Vaccine[df_snvs_meta_add_qc$Vaccine=="Non-vaccinated"] <- "Unvaccinated"
df_snvs_meta_add_qc$Group <- paste(df_snvs_meta_add_qc$Vaccine, df_snvs_meta_add_qc$lineage_sim)

df_meta <- read_csv("../results/df_samples_clean.csv", guess_max = 60000)
df_meta$Vaccine[df_meta$Vaccine=="Non-vaccinated"] <- "Unvaccinated"
df_tmp <- df_meta %>% mutate(Group=paste(Vaccine, lineage_sim)) %>% group_by(Group) %>% summarise(N=n(), Group_more=paste0(Group, "\n(N=", N, ")")) %>% unique()

df_snvs_meta_add_qc$Group_more <- factor(df_snvs_meta_add_qc$Group, levels=df_tmp$Group, labels=df_tmp$Group_more)

# spike non-synonymous mutation
df_spike_nonsyn_muts <- df_snvs_meta_add_qc %>% filter(gene=="S") %>% filter(lineage_sim %in% c("Alpha", "Delta", "Omicron")) %>% filter(effect=="missense_variant") %>% group_by(Group_more, mut_aa, pos_aa) %>% summarise(N=n()) %>% mutate(pos_aa=as.numeric(gsub("/1273", "", pos_aa, fixed=T)),) %>% arrange(pos_aa, mut_aa, N, Group_more)

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
	df_tmp <- df_spike_nonsyn_muts %>% filter(Group_more==group_i)
	SNP <- df_tmp$pos_aa
	Mutations <- GRanges("chr1", IRanges(SNP, width=1, names=gsub("p.", "", df_tmp$mut_aa, fixed=T)))
	Mutations$score <- df_tmp$N
	# Mutations$color <- as.character(df_tmp$color)
	muts_list[[group_i]] <- Mutations
	features_list[[group_i]] <- features
}

idx <- rev(grep("Delta", names(muts_list)))
pdf("../results/spike_muts_Delta.pdf", width=10, height=10)
lolliplot(muts_list[idx], features_list[idx], yaxis=FALSE)
dev.off()

idx <- rev(grep("Omicron", names(muts_list)))
pdf("../results/spike_muts_Omicron.pdf", width=10, height=10)
lolliplot(muts_list[idx], features_list[idx], yaxis=FALSE)
dev.off()

df_spike_nonsyn_muts_summary <- df_spike_nonsyn_muts %>% group_by(spike_regions, Group_more) %>% summarise(N=sum(N)) %>% arrange(spike_regions, Group_more, N)

mut_aa_int <- df_spike_nonsyn_muts %>% filter(spike_regions=="RBD") %>% .$mut_aa
df_snvs_meta_add_qc %>% filter(mut_aa %in% mut_aa_int) %>% select(Group, sample, mut_aa) %>% arrange(Group, sample)
