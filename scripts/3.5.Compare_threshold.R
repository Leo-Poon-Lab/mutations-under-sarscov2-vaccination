# We have ran a parallel analysis (scripts 3.0 to 3.4)
# using a MAF threshold of 0.05; In this script we try to 
# compare the difference between two sets of results.

# The comparison mainly focus on three aspects: 
# 1. Number of iSNVs (adjusted)
# 2. Nucleotide diversity
# 3. MAF
# 4. Selection pressure

library(tidyverse)
library(readxl)

# 1. Number of iSNVs (adjusted)
df_isnvs_test_0025 <- read_excel("../results/df_test_num_isnvs_by_vaccine_gene_more.xlsx")
df_isnvs_test_0050 <- read_excel("../results_005/df_test_num_isnvs_by_vaccine_gene_more_005.xlsx")

df_isnvs_test <- bind_rows(df_isnvs_test_0025 %>% mutate(MAF_threshold=0.025), df_isnvs_test_0050 %>% mutate(MAF_threshold=0.050))

df_isnvs_test %>% filter(check_three, gene=="Full genome") %>% mutate(pair=paste(var1, var2)) %>% group_by(pair) %>% filter(n()==1) %>% select(pair, MAF_threshold, same_vaccine, same_lineage, notation) %>% t() # significant results only appear in one group. The results show that, when using 0.05, two addtional significant results can be seen, the other 34 pairs are consistent.

# 2. Nucleotide diversity
df_diversity_test_0025 <- read_excel("../results/df_test_diversity_pi_by_vaccine_gene_more.xlsx")
df_diversity_test_0050 <- read_excel("../results_005/df_test_diversity_pi_by_vaccine_gene_more.xlsx")

df_diversity_test <- bind_rows(df_diversity_test_0025 %>% mutate(MAF_threshold=0.025), df_diversity_test_0050 %>% mutate(MAF_threshold=0.050))

df_diversity_test %>% filter(check_three, gene=="Full genome") %>% mutate(pair=paste(var1, var2)) %>% group_by(pair) %>% filter(n()==1) %>% select(pair, MAF_threshold, same_vaccine, same_lineage, notation) %>% t() # significant results only appear in one group. The results show that, when using 0.025, two addtional significant results can be seen, the other 32 pairs are consistent.

# 3. MAF
df_maf_test_0025 <- read_excel("../results/df_MAF_by_vaccine_gene_more.xlsx")
df_maf_test_0050 <- read_excel("../results_005/df_MAF_by_vaccine_gene_more.xlsx")

df_maf_test <- bind_rows(df_maf_test_0025 %>% mutate(MAF_threshold=0.025), df_maf_test_0050 %>% mutate(MAF_threshold=0.050))

df_maf_test %>% filter(check_three, gene=="Full genome") %>% mutate(pair=paste(var1, var2)) %>% group_by(pair) %>% filter(n()==1, same_vaccine) %>% select(pair, MAF_threshold, same_vaccine, same_lineage, notation) %>% t() # The results show that within same vaccination,no paris are consistent between 0.0025 and 0.050 thresholds. when using 0.025, three significant results can be seen. When using 0.05, two significant results can be seen
# 4. Selection pressure
df_maf_test %>% filter(check_three, gene=="Full genome") %>% mutate(pair=paste(var1, var2)) %>% group_by(pair) %>% filter(n()==1, same_lineage) %>% select(pair, MAF_threshold, same_vaccine, same_lineage, notation) %>% t() #The results show that within same lineages, when using 0.025, six addtional significant results can be seen (all related to BA5), the other 32 pairs are consistent.

# 4. Selection pressure
# TODO