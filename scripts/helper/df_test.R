cal_wilc_test <- function(df_input, var, genes=c("Full genome", "ORF1ab", "S", "N")) {
	df_groups <- df_input %>% filter(gene %in% genes) %>% select(gene) %>% unique() 
	df_wilc_test <- lapply(seq_len(nrow(df_groups)), function(i) {
		gene_i <- df_groups$gene[i]
		df_tmp <- df_input %>% filter(gene==gene_i)
		df_tmp$Doses <- paste0("Doses",df_tmp$Doses)
		df_tmp_grps <- df_tmp %>% select(Vaccine, Doses, lineage_sim) %>% unique()
		df_tmp_grps$x_grps <- df_tmp_grps %>% apply(1,paste0, collapse="_")
		mat_pairs <- combn(df_tmp_grps$x_grps, 2)
		df_rst <- apply(mat_pairs,2,function(y) {
			var1 <- y[1]
			var2 <- y[2]
			value1 <- df_tmp[[var]][df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var1] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var1]]
			value2 <- df_tmp[[var]][df_tmp$Vaccine==df_tmp_grps$Vaccine[df_tmp_grps$x_grps==var2] & df_tmp$lineage_sim==df_tmp_grps$lineage_sim[df_tmp_grps$x_grps==var2]]
			rst <- wilcox.test(value1, value2)
			tibble(var1=var1, var2=var2, var1_mean=mean(value1,na.rm=T), var2_mean=mean(value2,na.rm=T), median_var1=median(value1,na.rm=T), median_var2=median(value2,na.rm=T), p_value=rst$p.value)
		})
		df_rst <- bind_rows(df_rst)
		df_rst <- df_rst %>% mutate(gene=gene_i)
		df_rst %>% arrange(p_value)
	})
	bind_rows(df_wilc_test)
}

highlight_diff <- function(df){
	var1_split <- strsplit(df$var1, "_", fixed=T)
	var2_split <- strsplit(df$var2, "_", fixed=T)
	check1 <- sapply(var1_split, function(x)x[1]) == sapply(var2_split, function(x)x[1]) # Vaccine
	check2 <- sapply(var1_split, function(x)x[2]) == sapply(var2_split, function(x)x[2]) # Doses
	check3 <- sapply(var1_split, function(x)x[3]) == sapply(var2_split, function(x)x[3]) # lineage
	
	df$same_vaccine <- (check1 & check2)
	df$same_lineage <- check3
	df$within_group <- df$same_vaccine | df$same_lineage
	medians_low <- sapply(seq_len(nrow(df)), function(i){
		tmp <- abs(c(df$median_var1[i], df$median_var2[i]))
		tmp[which.min(tmp)]
	})
	df$median_diff <- abs(df$median_var1-df$median_var2)/medians_low
	df$median_diff_over_10percent <- df$median_diff>0.1
	df$check_three <- df$median_diff_over_10percent & df$within_group & (df$p_value<0.05)

	df$notation <- NA
	df$notation[df$p_value<0.1] <- "^"
	df$notation[df$p_value<0.05] <- "*"
	df$notation[df$p_value<0.01] <- "**"

	df %>% arrange(gene, check_three, same_vaccine, same_lineage, var1, var2)
}
