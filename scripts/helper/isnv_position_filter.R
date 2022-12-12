# This function is for filtering the low-quality iSNVs by specified positions
library(tidyverse)

pos_head_tail <- c(1:100, (29903-99):29903) # positions 1:100 and (29903-99):29903 should be removed; 

df_primer_new <- read_tsv("../../2021-10-11_nCoV_primers/results/primers_20211011.bed", col_names=F) # primers
df_primer_old <- read_tsv("../../2021-10-11_nCoV_primers/results/primers_old.bed", col_names=F)

filter_by_pos <- function(df_input, pos_indels=NA){
	stopifnot(all(c("pos", "con_base", "sec_base", "sample", "primer") %in% names(df_input)))

	if(!is.numeric(df_input$pos[1])){df_input$pos <- as.numeric(df_input$pos)}

	if(!is.na(pos_indels)){
		df_input <- df_input %>% filter(!pos %in% pos_indels) # 1. remove INDELs
	}
	
	df_input <- df_input %>% filter(!pos %in% pos_head_tail) # 2. remove head and tail 100 bases

	# 3. exclude all positions in the PCR primer binding regions
	df_input <- df_input %>% mutate(pos_combn=paste(pos, primer))
	pos_all <- unique(df_input$pos_combn)
	check <- sapply(pos_all, function(x) {
		pos_x <- as.numeric(strsplit(x, " ")[[1]][1])
		primer_x <- strsplit(x, " ")[[1]][2]
		if(primer_x=="new"){
			any((df_primer_new$X2 <= pos_x) & (df_primer_new$X3 >= pos_x))	
		} else {
			any((df_primer_old$X2 <= pos_x) & (df_primer_old$X3 >= pos_x))	
		}	
	})
	df_input <- df_input %>% filter(!pos_combn %in% pos_all[check])

	df_input <- df_input %>% filter(!pos %in% c(15494, 15489, 25381, 10194, 22422)) # 4. excluding primer/homoplasy sites, https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28420-7/MediaObjects/41467_2022_28420_MOESM1_ESM.pdf
	
	return(df_input)
}
