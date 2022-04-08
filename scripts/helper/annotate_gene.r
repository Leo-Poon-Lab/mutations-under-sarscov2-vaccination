# annotate gene
data_nsp_region <- readr::read_csv("../data/ORF_SCoV2.csv")

get_orf <- function(position) {
	Gene <- sapply(as.numeric(position), function(x){
		tmp <- data_nsp_region$sequence[data_nsp_region$start<=x & data_nsp_region$stop>=x]
		if(length(tmp)==0){return(NA)}else{return(tmp)}
	})
	AmA_pos <- sapply(as.numeric(position), function(x){
		check <- data_nsp_region$start<=x & data_nsp_region$stop>=x
		tmp <- floor((x-data_nsp_region$start[check])/3)+1
		if(length(tmp)==0){return(NA)}else{return(tmp)}
	})
	pos_adjust <- data_nsp_region$stop[data_nsp_region$sequence=="nsp12_1"] - data_nsp_region$start[data_nsp_region$sequence=="nsp12_1"]+1
	AmA_pos[which(Gene=="nsp12_2")] <- AmA_pos[which(Gene=="nsp12_2")]+pos_adjust/3
	Gene[which(Gene=='nsp12_2')] <- "nsp12"
	return(dplyr::tibble(position=position, gene_nsp=Gene, pos_orf=AmA_pos))
}

get_orf_gene <- function(position) {
	Gene <- sapply(as.numeric(position), function(x){
		tmp <- data_nsp_region$sequence[data_nsp_region$start<=x & data_nsp_region$stop>=x]
		if(length(tmp)==0){return(NA)}else{return(tmp)}
	})
}
