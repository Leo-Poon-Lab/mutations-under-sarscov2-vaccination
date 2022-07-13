library(tidyverse)
library(readxl)
library(ggsci)
library(patchwork)

colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine) <- c("Comirnaty", "CoronaVac", "Unvaccinated")

df_f5a_delta_cd8 <- read_excel("../draft/2022-06-28/From_ust/Data_Figure4A.xlsx", sheet=1)
df_f5a_delta_cd4 <- read_excel("../draft/2022-06-28/From_ust/Data_Figure4A.xlsx", sheet=2)
df_f5a_omicron_cd8 <- read_excel("../draft/2022-06-28/From_ust/Data_Figure4A.xlsx", sheet=3)
df_f5a_omicron_cd4 <- read_excel("../draft/2022-06-28/From_ust/Data_Figure4A.xlsx", sheet=4)
df_f5a_delta_cd8 <- df_f5a_delta_cd8 %>% mutate(type = "Delta (All epitopes)", T_cell="CD8+") %>% pivot_longer(`Non-vaccinated Delta`:`Sinovac Delta`)
df_f5a_delta_cd4 <- df_f5a_delta_cd4 %>% mutate(type = "Delta (All epitopes)", T_cell="CD4+") %>% pivot_longer(`Non-vaccinated Delta`:`Sinovac Delta`)
df_f5a_omicron_cd8 <- df_f5a_omicron_cd8 %>% mutate(type = "Omicron (All epitopes)", T_cell="CD8+") %>% pivot_longer(`Non-vaccinated Omicron`:`Sinovac Omicron`)
df_f5a_omicron_cd4 <- df_f5a_omicron_cd4 %>% mutate(type = "Omicron (All epitopes)", T_cell="CD4+") %>% pivot_longer(`Non-vaccinated Omicron`:`Sinovac Omicron`)

df_f5b_delta_cd8 <- read_excel("../draft/2022-06-28/From_ust/Data_Figure4B.xlsx", sheet=1)
df_f5b_delta_cd4 <- read_excel("../draft/2022-06-28/From_ust/Data_Figure4B.xlsx", sheet=2)
df_f5b_omicron_cd8 <- read_excel("../draft/2022-06-28/From_ust/Data_Figure4B.xlsx", sheet=3)
df_f5b_omicron_cd4 <- read_excel("../draft/2022-06-28/From_ust/Data_Figure4B.xlsx", sheet=4)
df_f5b_delta_cd8 <- df_f5b_delta_cd8 %>% mutate(type = "Delta (HK epitope-HLA pairs)", T_cell="CD8+") %>% pivot_longer(`Non-vaccinated Delta`:`Sinovac Delta`)
df_f5b_delta_cd4 <- df_f5b_delta_cd4 %>% mutate(type = "Delta (HK epitope-HLA pairs)", T_cell="CD4+") %>% pivot_longer(`Non-vaccinated Delta`:`Sinovac Delta`)
df_f5b_omicron_cd8 <- df_f5b_omicron_cd8 %>% mutate(type = "Omicron (HK epitope-HLA pairs)", T_cell="CD8+") %>% pivot_longer(`Non-vaccinated Omicron`:`Sinovac Omicron`)
df_f5b_omicron_cd4 <- df_f5b_omicron_cd4 %>% mutate(type = "Omicron (HK epitope-HLA pairs)", T_cell="CD4+") %>% pivot_longer(`Non-vaccinated Omicron`:`Sinovac Omicron`)

df_fs10a_delta_cd8 <- read_excel("../draft/2022-06-28/From_ust/Data_FigureS8A.xlsx", sheet=1)
df_fs10a_delta_cd4 <- read_excel("../draft/2022-06-28/From_ust/Data_FigureS8A.xlsx", sheet=2)
df_fs10a_omicron_cd8 <- read_excel("../draft/2022-06-28/From_ust/Data_FigureS8A.xlsx", sheet=3)
df_fs10a_omicron_cd4 <- read_excel("../draft/2022-06-28/From_ust/Data_FigureS8A.xlsx", sheet=4)
df_fs10a_delta_cd8 <- df_fs10a_delta_cd8 %>% mutate(type = "Delta (All epitopes)", T_cell="CD8+") %>% pivot_longer(`Non-vaccinated Delta`:`Sinovac Delta`)
df_fs10a_delta_cd4 <- df_fs10a_delta_cd4 %>% mutate(type = "Delta (All epitopes)", T_cell="CD4+") %>% pivot_longer(`Non-vaccinated Delta`:`Sinovac Delta`)
df_fs10a_omicron_cd8 <- df_fs10a_omicron_cd8 %>% mutate(type = "Omicron (All epitopes)", T_cell="CD8+") %>% pivot_longer(`Non-vaccinated Omicron`:`Sinovac Omicron`)
df_fs10a_omicron_cd4 <- df_fs10a_omicron_cd4 %>% mutate(type = "Omicron (All epitopes)", T_cell="CD4+") %>% pivot_longer(`Non-vaccinated Omicron`:`Sinovac Omicron`)

df_fs10b_delta_cd8 <- read_excel("../draft/2022-06-28/From_ust/Data_FigureS8B.xlsx", sheet=1)
df_fs10b_delta_cd4 <- read_excel("../draft/2022-06-28/From_ust/Data_FigureS8B.xlsx", sheet=2)
df_fs10b_omicron_cd8 <- read_excel("../draft/2022-06-28/From_ust/Data_FigureS8B.xlsx", sheet=3)
df_fs10b_omicron_cd4 <- read_excel("../draft/2022-06-28/From_ust/Data_FigureS8B.xlsx", sheet=4)
df_fs10b_delta_cd8 <- df_fs10b_delta_cd8 %>% mutate(type = "Delta (HK epitope-HLA pairs)", T_cell="CD8+") %>% pivot_longer(`Non-vaccinated Delta`:`Sinovac Delta`)
df_fs10b_delta_cd4 <- df_fs10b_delta_cd4 %>% mutate(type = "Delta (HK epitope-HLA pairs)", T_cell="CD4+") %>% pivot_longer(`Non-vaccinated Delta`:`Sinovac Delta`)
df_fs10b_omicron_cd8 <- df_fs10b_omicron_cd8 %>% mutate(type = "Omicron (HK epitope-HLA pairs)", T_cell="CD8+") %>% pivot_longer(`Non-vaccinated Omicron`:`Sinovac Omicron`)
df_fs10b_omicron_cd4 <- df_fs10b_omicron_cd4 %>% mutate(type = "Omicron (HK epitope-HLA pairs)", T_cell="CD4+") %>% pivot_longer(`Non-vaccinated Omicron`:`Sinovac Omicron`)


# Figure 5
df_f5a <- bind_rows(df_f5a_delta_cd8, df_f5a_delta_cd4, df_f5a_omicron_cd8, df_f5a_omicron_cd4)
df_f5a$name <- gsub("Non-vaccinated", "Unvaccinated", df_f5a$name)
df_f5a$name <- gsub("BioNTech", "Comirnaty", df_f5a$name)
df_f5a$name <- gsub("Sinovac", "CoronaVac", df_f5a$name)
df_f5a$vaccine <- gsub(" .+", "", df_f5a$name)
df_f5a$lineage <- gsub(".+ ", "", df_f5a$name)

p_5a <- ggplot(df_f5a) +
	geom_boxplot(aes(x=name, y=value, fill=vaccine), alpha=0.8)+
	facet_grid(rows=vars(T_cell), cols=vars(type), scales="free_x")+
	theme_classic()+
	scale_fill_manual(name="Vaccine",values=colors_vaccine)+
	theme(axis.text.x = element_text(angle = 15, hjust=1),
		axis.title.x = element_blank())+
	ylab("Overlapping T cell epitopes per mutation")+
	ggtitle("A")

df_f5b <- bind_rows(df_f5b_delta_cd8, df_f5b_delta_cd4, df_f5b_omicron_cd8, df_f5b_omicron_cd4)
df_f5b$name <- gsub("Non-vaccinated", "Unvaccinated", df_f5b$name)
df_f5b$name <- gsub("BioNTech", "Comirnaty", df_f5b$name)
df_f5b$name <- gsub("Sinovac", "CoronaVac", df_f5b$name)
df_f5b$vaccine <- gsub(" .+", "", df_f5b$name)
df_f5b$lineage <- gsub(".+ ", "", df_f5b$name)

p_5b <- ggplot(df_f5b) +
	geom_boxplot(aes(x=name, y=value, fill=vaccine), alpha=0.8)+
	facet_grid(rows=vars(T_cell), cols=vars(type), scales="free_x")+
	theme_classic()+
	scale_fill_manual(name="Vaccine", values=colors_vaccine)+
	theme(axis.text.x = element_text(angle = 15, hjust=1),
		axis.title.x = element_blank())+
	ylab("Overlapping T cell epitopes per mutation")+
	ggtitle("B")

p5 <- (p_5a / p_5b) + plot_layout(guides="collect") & theme(legend.position='right')
ggsave("../results/figure_5.pdf", width=8, height=8, plot=p5)

# Figure S10
df_fs10a <- bind_rows(df_fs10a_delta_cd8, df_fs10a_delta_cd4, df_fs10a_omicron_cd8, df_fs10a_omicron_cd4)
df_fs10a$name <- gsub("Non-vaccinated", "Unvaccinated", df_fs10a$name)
df_fs10a$name <- gsub("BioNTech", "Comirnaty", df_fs10a$name)
df_fs10a$name <- gsub("Sinovac", "CoronaVac", df_fs10a$name)
df_fs10a$vaccine <- gsub(" .+", "", df_fs10a$name)
df_fs10a$lineage <- gsub(".+ ", "", df_fs10a$name)

p_s10a <- ggplot(df_fs10a) +
	geom_boxplot(aes(x=name, y=value, fill=vaccine), alpha=0.8)+
	facet_grid(rows=vars(T_cell), cols=vars(type), scales="free_x")+
	theme_classic()+
	scale_fill_manual(name="Vaccine",values=colors_vaccine)+
	theme(axis.text.x = element_text(angle = 15, hjust=1),
		axis.title.x = element_blank())+
	ylab("Overlapping T cell epitopes per mutation")+
	ggtitle("A")

df_fs10b <- bind_rows(df_fs10b_delta_cd8, df_fs10b_delta_cd4, df_fs10b_omicron_cd8, df_fs10b_omicron_cd4)
df_fs10b$name <- gsub("Non-vaccinated", "Unvaccinated", df_fs10b$name)
df_fs10b$name <- gsub("BioNTech", "Comirnaty", df_fs10b$name)
df_fs10b$name <- gsub("Sinovac", "CoronaVac", df_fs10b$name)
df_fs10b$vaccine <- gsub(" .+", "", df_fs10b$name)
df_fs10b$lineage <- gsub(".+ ", "", df_fs10b$name)

p_s10b <- ggplot(df_fs10b) +
	geom_boxplot(aes(x=name, y=value, fill=vaccine), alpha=0.8)+
	facet_grid(rows=vars(T_cell), cols=vars(type), scales="free_x")+
	theme_classic()+
	scale_fill_manual(name="Vaccine", values=colors_vaccine)+
	theme(axis.text.x = element_text(angle = 15, hjust=1),
		axis.title.x = element_blank())+
	ylab("Overlapping T cell epitopes per mutation")+
	ggtitle("B")

ps10 <- (p_s10a / p_s10b) + plot_layout(guides="collect") & theme(legend.position='right')
ggsave("../results/figure_s10.pdf", width=8, height=8, plot=ps10)
