library(tidyverse)
source("./helper/pysamstats.r")

get_100reads_coverage <- function(samples_todo){
	# samples_todo <- data_meta_vac$Sample
	files_bam <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
	files_bam_full <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$", full.names = T)

	files_bam_rst_full <- c(list.files("../results/pysamstats/", "tsv$", full.names = T), list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/", "tsv$", full.names = T), list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", "tsv$", full.names = T))

	idx <- sapply(paste0("/", samples_todo, ".tsv"), function(x) {
		tmp <- grep(x, files_bam_rst_full)[1]
	})
	sum(is.na(idx))

	samples_todo_no_tsv <- samples_todo[is.na(idx)]
	check_bam_todo_no_tsv <- files_bam %in% paste0(samples_todo_no_tsv, "-trimmed.masked.bam")
	files_bam_todo_no_tsv <- files_bam_full[check_bam_todo_no_tsv]

	## Upload data to server and re-analyse SNPs/iSNVs
	samples_toupload_ori <- samples_todo_no_tsv
	samples_alrd_done <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/", "bam$")
	samples_alrd_done <- gsub("-trimmed.+", "", samples_alrd_done)
	samples_toupload <- samples_toupload_ori[!samples_toupload_ori %in% samples_alrd_done] # we already analysis some of the data
	# sort(samples_toupload)

	
	files_fastq_all <- list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive/", full.names = T)
	samples_all <- sapply(files_fastq_all, function(x){
		tmp <- strsplit(x, "/", fixed=T)[[1]]
		sample_t <- tmp[length(tmp)]
		gsub(".fastq.gz", "", sample_t)
	})
	files_fastq1 <- files_fastq_all[samples_all %in% paste0(samples_toupload, "_1")]
	files_fastq2 <- files_fastq_all[samples_all %in% paste0(samples_toupload, "_2")]

	samples_not_available <- samples_toupload[!samples_toupload %in% gsub(".+/", "", gsub("_1.fastq.gz", "", files_fastq1))]

	if(length(files_fastq1)>0){	
		## copy the files to remote
		print(paste0("Uploading ", length(files_fastq1), " samples!"))
		sapply(c(rbind(files_fastq1, files_fastq2)), function(x){
			cp.remote(path.src = x, remote.src = "", remote.dest = "hggu@147.8.70.166", path.dest = "~/work/2020-09-01_COVID_NGS_pipeline/NGS_data_input/", verbose = T)
		})

		if(!all(samples_todo_no_tsv%in%samples_not_available)){ # make sure all the bam files are ready
			print("Bam files of some samples are still waiting, please sync '~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/archive_bam/' and try again!")
			return("")
		}
	} else {
		process_pysamstats(files_bam_todo_no_tsv, n_cores=8)
	}


	files_bam_rst_full <- c(list.files("../results/pysamstats/", "tsv$", full.names = T), list.files("~/../../Volumes/Backup/Haogao/work/2020/2020-09-01_COVID_NGS_pipeline/results/pysamstats/", "tsv$", full.names = T), list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", "tsv$", full.names = T))
	
	sapply(paste0(samples_todo, ".tsv"), function(x) {
		# x="WHP1335.tsv"
		print(x)
		tmp <- grep(x, files_bam_rst_full)[1]
		if(is.na(tmp)){return(NA)}
		df_tmp <- read_pysamstats(files_bam_rst_full[tmp], pp=FALSE, genome_coverage=0.9)
		if(is.na(df_tmp)){return(NA)}else{return(sum(df_tmp$reads_all>=100)/29903)}
	})
}
