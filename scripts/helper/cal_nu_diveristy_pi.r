# calculate the nucleotide diversity pi over a genomic region
## The function should be ideally works on pysamstats results, where it can use the number of properly paired reads.
## we use nucleotide diversity pi here (Ref: https://academic.oup.com/ve/article/5/1/vey041/5304643, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4316684/)

## Cal Dl
### At the locus l, where ni copies of the allele i are observed, the pro- portion of pairwise differences between alleles may be calcu- lated as
### Input should be allele frequencies
Cal_dl <- function(freqs){
	N <- sum(freqs)
	sum_tmp <- sum(sapply(freqs, function(ni) {ni*(ni-1)}))
	return(1-(sum_tmp/(N*(N-1))))
}
# # unit test
# Cal_dl(1:4)
# Cal_dl(1:2)
# Cal_dl(1:20)
# Cal_dl(c(16,1,1,2))
# Cal_dl(c(160,1,1,2))
# Cal_dl(rep(4,4))
# Cal_dl(rep(400,4))

## Cal pi
### The statistic p may then be calculated for a genome 
Cal_pi <- function(dls, length_genome) {
	if(length(dls)==0){return(NA)}
	sum(dls)/length_genome
}

# Broadly speaking, excess nonsynonymous polymorphism (πN/πS > 1) points toward diversifying or positive selection while excess synonymous polymorphism (πN/πS < 1) indicates purifying selection. When πN / πS is approximately 1, genetic drift, i.e., stochastic changes in the frequency of viral genotypes over time, can be an important force shaping genetic diversity. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7946358/)

# The below scripts use SNPGenie for calculation of pi
dir.create("../results/snpgenie/")
Generate_cvf_files <- function(df_mut_t, pos_name="pos", ref_name="con_base", alt_name="sec_base") {
	df_vcf <- tibble('#CHROM' = "MN908947.3", POS=df_mut_t[[pos_name]])
	df_vcf$ID <- "."
	df_vcf$REF <- df_mut_t[[ref_name]]
	df_vcf$ALT <- df_mut_t[[alt_name]]
	df_vcf$QUAL <- 30
	df_vcf$FILTER <- "PASS"
	df_vcf$INFO <- paste0("DP=", df_mut_t$con_fwd+df_mut_t$con_rev+df_mut_t$sec_fwd+df_mut_t$sec_rev, ";", "AF=", df_mut_t$sec_freq)
	# df_vcf$FORMAT <- 
	sample_name_t <- df_mut_t$sample[1]
	df_vcf$`<SAMPLE>` <- sample_name_t
	df_vcf
}

prepare_snpgenie <- function(sample_name_t) {
	file.copy("../data/reference.fasta", paste0("../results/snpgenie/", sample_name_t, "/reference.fasta"), overwrite=T)
	file.copy("../data/SCoV2_genes.gtf", paste0("../results/snpgenie/", sample_name_t, "/SCoV2_genes.gtf"), overwrite=T)	
	cur_wd <- getwd()
	setwd(paste0("../results/snpgenie/", sample_name_t))
	system("rm -r SNPGenie_Results")
	# system("~/softwares/SNPGenie/snpgenie.pl --vcfformat=2")
	setwd(cur_wd)
}

read_snpgenie_rst <- function(files) {
	header_tmp <- read_tsv(files[1], col_types=cols(.default = "c")) %>% names()
	read_tsv_self <- function(x){
		tmp <- read_tsv(x, col_types=cols(.default = "c"), col_names=F)
		tmp %>% filter(X1!="file") %>% filter(!grepl("^temp", X1))
	}
	tmp_all <- lapply(files, read_tsv_self)
	tmp_all <- bind_rows(tmp_all) 
	names(tmp_all) <- header_tmp
	tmp_all
}


# perl snpgenie.pl --vcfformat=4 --snpreport='/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data/Cats/SARSCoV2_transmission_in_domestic_cats/data_derived/1_2A.vcf.recode.vcf' --fastafile='/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data/Cats/SARSCoV2_transmission_in_domestic_cats/MW219695.1/MW219695.1.fasta' --gtffile='/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data/Cats/SARSCoV2_transmission_in_domestic_cats/MW219695.1/genes.gtf' --o='/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data/Cats/SARSCoV2_transmission_in_domestic_cats/data_derived/diversity/1_2A'
