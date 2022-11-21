add_residuals <- function(df,x_var,y_var,out_var){
	genes = unique(df$gene)
	df_out <- lapply(seq_along(genes), function(i) { # adjust by genes, good
		gene_i <- genes[i]
		df_tmp <- df %>% filter(gene==gene_i) 
		df_tmp <- df_tmp[!is.na(df_tmp[[x_var]]),]
		df_tmp <- df_tmp[!is.na(df_tmp[[y_var]]),]
		lm_i <- lm(paste0(y_var,"~",x_var), data=df_tmp)
		df_tmp[[out_var]] <- lm_i$residuals
		df_tmp
	})
	bind_rows(df_out)
}

plot_reg <- function(data, x_var, y_var, x_lab=NA, y_lab=NA){
  label_x <- min(data[[x_var]][!is.infinite(data[[x_var]])],na.rm=T)
  label_y <- max(data[[y_var]][!is.infinite(data[[y_var]])],na.rm=T)
  tmp <- ggscatter(data, x = x_var, y = y_var,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    alpha=0.8,
    size=0.8,
    conf.int = TRUE # Add confidence interval
    )+
    stat_cor(method = "pearson", label.x = label_x, label.y = label_y)+
    NULL
  if(!is.na(x_lab)){tmp <- tmp+xlab(x_lab)}
  if(!is.na(y_lab)){tmp <- tmp+ylab(y_lab)}
  tmp
}