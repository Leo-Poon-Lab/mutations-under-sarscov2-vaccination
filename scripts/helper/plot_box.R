plot_box <- function(df_plot, x_var, y_var, color_var, x_lab="Group", y_lab) {
	p <- ggplot(df_plot, aes_string(x=x_var, y=y_var, color=color_var))+ # mutation rate
	geom_point(alpha=0.8, position=position_jitterdodge(jitter.width=0.1), size=0.8)+
	geom_boxplot(outlier.size=0, outlier.alpha=0, alpha=0.8)+
	ylab(y_lab)+
	xlab(x_lab)+
	theme_classic()+
	NULL

	p
}
