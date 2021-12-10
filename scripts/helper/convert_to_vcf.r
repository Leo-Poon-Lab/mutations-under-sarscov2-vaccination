# convert ivar tsv to vcf

dir.create("./helper")
dir.create("../results/vcf")
system("wget https://raw.githubusercontent.com/Leo-Poon-Lab/SNVs-comparison-table/c6b04f29bf4cf209f1755376f90978cc79863b98/scripts/helper/ivar_variants_to_vcf.py -P ./helper/ -nc")
system("wget https://raw.githubusercontent.com/Leo-Poon-Lab/SNVs-comparison-table/ebea548fcc185a129f86b54c6f0538f1576f2c79/scripts/helper/Parse_SnpEff.r -P ./helper/ -nc")
system("chmod 755 ./helper/ivar_variants_to_vcf.py")

convert_to_vcf <- function(files_tsv){
	df_ann_raw <- lapply(seq_along(files_tsv), function(i){
		file_tsv_i <- files_tsv[i]
		sample_t <- strsplit(file_tsv_i, "/", fixed=T)[[1]]
		sample_t <- sample_t[length(sample_t)]
		outfile <- paste0("../results/vcf/", sample_t, ".vcf")
		
		system(paste0("./helper/ivar_variants_to_vcf.py ", file_tsv_i, " ", outfile))
		
	})
}

annotate_snpeff <- function(files_vcf){
	genome <- "NC_045512.2"
	df_out <- lapply(files_vcf, function(x){
		tmp <- readLines(x)
		writeLines(gsub("MN908947_3", genome, tmp), x)
		
		outfile_snpeff <- paste0(x, ".snpeff")
		outfile_snpeff_csq <- paste0(x, ".snpeff.csq")
		outfile_csv <- paste0(x, ".csv")
		system(paste0("java -jar ~/softwares/snpEff/snpEff.jar ", genome, " ", x, " > ", outfile_snpeff))
		# system(paste0("bcftools csq --force --phase a -f ../data/refernece.fasta -g ../data/Sars_cov_2.ASM985889v3.101.primary_assembly.MN908947.3.gff3 ", outfile_snpeff, " > ", outfile_snpeff_csq))
		
		system(paste0("Rscript ./helper/Parse_SnpEff.r ", outfile_snpeff,  " ", outfile_csv))
		
	})
}

read_snpeff <- function(files_snpeff_csv) {
	df_tmp <- lapply(files_snpeff_csv, function(x) {
		df_tmp <- read_csv(x)
		sample_t <- strsplit(gsub("\\.tsv$", "", x), "/", fixed=T)[[1]]
		sample_t <- sample_t[length(sample_t)]
		sample_sim_t <- strsplit(sample_t, ".", fixed=T)[[1]][1]
		df_tmp$sample <- gsub("ivar_", "", sample_sim_t)
		return(df_tmp)
	})
	return(bind_rows(df_tmp))
}
