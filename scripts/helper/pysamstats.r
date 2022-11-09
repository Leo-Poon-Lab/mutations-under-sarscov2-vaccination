# call the pysamstats results for a list of bam files
library(parallel)
if(!"./helper" %in% list.dirs()){dir.create("./helper")}
if(!"pysamstats" %in% list.dirs(path="../results/", full.names = FALSE, recursive = FALSE)){dir.create("../results/pysamstats")}
system("wget https://raw.githubusercontent.com/Leo-Poon-Lab/VOC-SNPs-checking/main/data/reference.fasta -P ./helper/ -nc")
file_ref_seq <- "./helper/reference.fasta"

process_pysamstats <- function(file_input, n_cores){
	mclapply(file_input, function(x){
		print(x)
		sample_t <- strsplit(x, "/", fixed=T)[[1]]
		sample_t <- sample_t[length(sample_t)]
		sample_sim_t <- strsplit(sample_t, "-trimmed", fixed=T)[[1]][1]
			
		file_out_x <- paste("../results/")
		system(paste0("/usr/local/Caskroom/miniconda/base/bin/pysamstats --type variation_strand -d -D 100000 -f ", file_ref_seq, " ", x, " > ../results/pysamstats/", sample_sim_t, ".tsv"))
	}, mc.cores=n_cores)
}

read_pysamstats <- function(x, pp=TRUE, genome_coverage=0.9) {
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
	if(sum(df_bam_rst$reads_all>=100)>=29903*genome_coverage){return(df_bam_rst)} else{return(NA)}
}

extra_info_from_pysamstats <- function(df_bam_rst_filter){
	bases <- c("deletions", "insertions", "A", "C", "T", "G")
	# bases <- c("A", "C", "T", "G") # ignore INDELs
	tmp <- mclapply(seq_len(nrow(df_bam_rst_filter)), function(i){
		# print(i)
		df_bam_rst_i <- df_bam_rst_filter[i,]
		depth_i <- df_bam_rst_i$reads_all # count on paired reads·
		meta_i <- c(pos=df_bam_rst_i$pos, sample=df_bam_rst_i$sample)

		if(depth_i==0){ # depth is 0
			tmp <- list(c(meta_i, depth="0", depth_fwd="0", depth_rev="0", con_base=NA, con_fwd="0", con_rev="0", con_freq="0", sec_base=NA, sec_fwd="0", sec_rev="0", sec_freq="0", strand_bias=NA))
			return(tmp)
		} else {
			num_bases <- c(df_bam_rst_i$deletions, df_bam_rst_i$insertions, df_bam_rst_i$A, df_bam_rst_i$C, df_bam_rst_i$T, df_bam_rst_i$G)
			freq_bases <- num_bases/depth_i
			ord_t <- order(num_bases, decreasing = T)
			ord_t <- ord_t[freq_bases[ord_t]>0.001] # only keep variants with freq > 1%

			consensus_base <- bases[ord_t[1]]
			con_fwd_i <- df_bam_rst_i[[paste0(consensus_base, "_fwd")]]
			con_rev_i <- df_bam_rst_i[[paste0(consensus_base, "_rev")]]
			con_freq_i <- freq_bases[ord_t[1]]
			
			if(length(ord_t)==1){
				return(list(c(meta_i, depth=depth_i, depth_fwd=df_bam_rst_i$reads_fwd, depth_rev=df_bam_rst_i$reads_rev, con_base=consensus_base, con_fwd=con_fwd_i, con_rev=con_rev_i, con_freq=con_freq_i, sec_base=NA, sec_fwd="0", sec_rev="0", sec_freq="0", strand_bias=NA)))
			} else {
				tmp <- lapply(ord_t[-1], function(ord_t_i){
					secondary_base <- bases[ord_t_i]
					sec_fwd_i <- df_bam_rst_i[[paste0(secondary_base, "_fwd")]]
					sec_rev_i <- df_bam_rst_i[[paste0(secondary_base, "_rev")]]
					sec_freq_i <- freq_bases[ord_t_i]
					
					strand_bias_i <- sort(c(sec_fwd_i, sec_rev_i))
					strand_bias_i <- strand_bias_i[2]/strand_bias_i[1]

					return(c(meta_i, depth=depth_i, depth_fwd=df_bam_rst_i$reads_fwd, depth_rev=df_bam_rst_i$reads_rev, con_base=consensus_base, con_fwd=con_fwd_i, con_rev=con_rev_i, con_freq=con_freq_i, sec_base=secondary_base, sec_fwd=sec_fwd_i, sec_rev=sec_rev_i, sec_freq=sec_freq_i, strand_bias=strand_bias_i))
				})
				return(tmp)
			}
		}
	}, mc.cores = 16)
	bind_rows(tmp)
}
