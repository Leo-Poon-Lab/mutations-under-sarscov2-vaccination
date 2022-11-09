library(parallel)
library(tidyverse)
# dir.create("../results/sam_flagstat/", showWarnings=F)

get_total_seqs_after_trimming <- function(samples_todo){
	# samples_todo <- df_meta$Sample
	files_bam <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
	files_bam_full <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$", full.names = T)
	
	check_bam_todo_no_tsv <- which(files_bam %in% paste0(samples_todo, "-trimmed.masked.bam"))

	stopifnot(all(files_bam[check_bam_todo_no_tsv] == paste0(sort(samples_todo), "-trimmed.masked.bam")))
	
	check_bam_todo_no_tsv_match <- check_bam_todo_no_tsv[match(paste0(samples_todo, "-trimmed.masked.bam"), files_bam[check_bam_todo_no_tsv])]

	stopifnot(all(files_bam[check_bam_todo_no_tsv_match] == paste0(samples_todo, "-trimmed.masked.bam")))
	files_bam_todo_no_tsv <- files_bam_full[check_bam_todo_no_tsv_match]
	stopifnot(length(files_bam_todo_no_tsv)==length(samples_todo))

	seqs_info <- mclapply(seq_along(files_bam_todo_no_tsv), function(i) {
		print(i)
		sample_i <- samples_todo[i]
		file_bam_i <- files_bam_todo_no_tsv[i]
		# command_i <- paste0("/usr/local/Caskroom/miniconda/base/envs/artic-ncov2019/bin/samtools flagstats -@ 4 -O tsv ", file_bam_i, " > ", "../results/sam_flagstat/", sample_i, ".txt")
		file_size_mb <- file.info(file_bam_i, extra_cols = F)$size/1000/1000
		if(file_size_mb > 50){
			return("file_size_mb > 50")
		} else {
			print(paste0("Runing samtools stats for ", sample_i))
			command_i <- paste0("/usr/local/Caskroom/miniconda/base/envs/artic-ncov2019/bin/samtools view -@ 4 -c ", file_bam_i)
			seqs_i <- system(command_i, intern=T)
			return(seqs_i)
		}
	}, mc.cores=4)

	seqs_info <- unlist(seqs_info)
	stopifnot(length(samples_todo) == length(seqs_info))

	return(seqs_info)

}