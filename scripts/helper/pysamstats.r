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

read_pysamstats <- function(x, pp=TRUE) {
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
	if(sum(df_bam_rst$reads_all>=100)>=27000){return(df_bam_rst)} else{return(NA)}
}
