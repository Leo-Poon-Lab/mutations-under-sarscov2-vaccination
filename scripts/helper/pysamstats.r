# call the pysamstats results for a list of bam files

dir.create("./helper")
dir.create("../results/pysamstats")
system("wget https://raw.githubusercontent.com/Leo-Poon-Lab/VOC-SNPs-checking/main/data/reference.fasta -P ./helper/ -nc")
file_ref_seq <- "./helper/reference.fasta"

process_pysamstats <- function(file_input){
	sapply(file_input, function(x){
		print(x)
		sample_t <- strsplit(x, "/", fixed=T)[[1]]
		sample_t <- sample_t[length(sample_t)]
		sample_sim_t <- strsplit(sample_t, "-trimmed", fixed=T)[[1]][1]
			
		file_out_x <- paste("../results/")
		system(paste0("conda run -n base pysamstats --type variation_strand -d -D 100000 -f ", file_ref_seq, " ", x, " > ../results/pysamstats/", sample_sim_t, ".tsv"))
	})
}

read_pysamstats <- function(files_bam_rst_full) {
	df_bam_rst <- lapply(files_bam_rst_full, function(x){
		tmp <- read_tsv(x)
		sample_x <- strsplit(x, "/", fixed=T)[[1]]
		sample_x <- sample_x[length(sample_x)]
		sample_x <- gsub(".tsv", "", sample_x, fixed=T)
		tmp$sample <- sample_x
		return(tmp)
	})
	return(bind_rows(df_bam_rst))
}

