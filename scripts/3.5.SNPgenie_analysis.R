# This part analysis the selection pressure and diversity of iSNVs in SARS-CoV-2 samples
##### Acknowledgement #####
# The results are based on output of SNPGenie from https://github.com/singing-scientist/snpgenie
# The analysis scripts are heavily based on https://github.com/krisp-kwazulu-natal/within-host-diversity-manuscript-analysis-code/blob/a276286680de3723e2b1e70f7a060750892cf8af/scripts/diversity_selection_analyses.R  
##### Acknowledgement #####

library(tidyverse)
library(lubridate)
library(boot)
library(writexl)
library(parallel)
library(RColorBrewer)
source("./helper/cal_nu_diveristy_pi.r")
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")

colors_lineage=c("#e41a1c", "#33a02c", "#1f78b4", "#ff7f00", "#f781bf", "#666666") 
names(colors_lineage) <- c("Alpha", "Delta", "Omicron", "B.1.36", "B.1.36.27", "B.1.1.63")
colors_vaccine=c("#a65628", "#7570b3", "#999999")
names(colors_vaccine)=c("BioNTech", "Sinovac", "Non-vaccinated")
df_orf_sim <- read_csv("../data/ORF_SCoV2_sim.csv")

df_meta <- read_csv("../results/df_samples_clean.csv", guess_max = 60000)
df_meta$lineage_sim <- factor(df_meta$lineage_sim, levels = names(colors_lineage))
df_meta$Vaccine <- factor(df_meta$Vaccine, levels = names(colors_vaccine))
df_meta$detection_lag <- as.numeric(dmy(df_meta$`Report date`) - dmy(df_meta$`Onset date`))

df_snvs_meta_add_qc <- read_csv("../results/df_snvs_meta_add_qc_bam.csv", guess_max=600000)
samples_analysed <- unique(df_snvs_meta_add_qc$sample)

# read codon results
files_tmp <- list.files("../results/snpgenie/", pattern="codon_results", recursive=T, full.names=T)
idx <- sapply(paste0("/", samples_analysed, "/"), function(x) {
	tmp <- grep(x, files_tmp)[1]
})

codon_results <- mclapply(files_tmp[idx], read_tsv, mc.cores=16)
codon_results <- bind_rows(codon_results)
codon_results$gene <- gsub("gene-", "", codon_results$product)
codon_results$gene[codon_results$gene=="orf1ab"] <- "ORF1ab"
codon_results$sample <- gsub(".vcf", "", codon_results$file)
stopifnot(all(unique(codon_results$gene) %in% df_orf_sim$sequence)) # normalize gene name
stopifnot(all(unique(codon_results$sample) %in% df_meta$sample)) # normalize sample name

df_all_codon <- codon_results %>% select(gene, site, codon) %>% unique()
df_samples_gene <- full_join(df_meta %>% select(sample, Vaccine, lineage_sim), df_all_codon, by = character())

codon_results <- left_join(df_samples_gene, codon_results %>% select(-product, -file)) # complete the codon results by adding samples without any iSNVs
codon_results <- codon_results %>% mutate_at(vars(num_overlap_ORF_nts:S_gdiv), function(x) {x[is.na(x)] <- 0;x})
stopifnot(sum(is.na(codon_results$N_diffs))==0)

# Rename product column
names(codon_results)[names(codon_results) == "gene"] <- "product_segment"

# Get unique product names
(uniq_product_segments <- unique(codon_results$product_segment))

# Number the codons for each product segment

# Obtain unique products and sites
(productSegment_uniqSites <- codon_results %>%
    group_by(product_segment, site) %>%
    summarise(
      count = n()
    )) # automatically sorts by product and site


# remove count column
productSegment_uniqSites <- dplyr::select(productSegment_uniqSites, -count)

# add codon_num column
productSegment_uniqSites$codon_num <- 0

# loop each product and site
for (this_product in unique(productSegment_uniqSites$product_segment)) {
  #this_product <- "ORF6"
  (num_uniq_sites <- nrow(productSegment_uniqSites[productSegment_uniqSites$product_segment == this_product, ]))
  productSegment_uniqSites[productSegment_uniqSites$product_segment == this_product, ]$codon_num <- 1:num_uniq_sites
}

# Examine
#View(productSegment_uniqSites)

# JOIN codon numbers
codon_results <- left_join(x = codon_results, y = productSegment_uniqSites, by = c("product_segment", "site"))
codon_results <- dplyr::select(codon_results, sample, product_segment, codon_num, everything())
#View(codon_results)

# check codon numbers
codon_results %>%
  group_by(product_segment) %>%
  summarise(
    highest_codon_num = max(codon_num)
  ) # confirmed

#product highest_codon_num
#    product_segment highest_codon_num
#    <chr>                       <dbl>
#  1 E                              76
#  2 M                             223
#  3 N                             420
#  4 ORF10                          39
#  5 ORF1ab                       7097
#  6 ORF3a                         276
#  7 ORF6                           62
#  8 ORF7a                         122
#  9 ORF8                          122
# 10 S                            1274

# Manually rename products to group different segments of same product together
codon_results$product <- codon_results$product_segment
codon_results$codon_num_ORF <- codon_results$codon_num
codon_results <- dplyr::select(codon_results, sample, product, product_segment, codon_num, codon_num_ORF, everything())


################################################################################
### piN/piS by SAMPLE (whole-genome), NOL regions only
(codon_results_NOL_bySample <- filter(codon_results, num_overlap_ORF_nts == 0) %>%
    group_by(sample) %>%
    summarise(
      N_diffs = sum(N_diffs),
      S_diffs = sum(S_diffs),
      N_sites = sum(N_sites),
      S_sites = sum(S_sites)
    ))

codon_results_NOL_bySample$piN <- (codon_results_NOL_bySample$N_diffs / codon_results_NOL_bySample$N_sites)
codon_results_NOL_bySample$piN[is.na(codon_results_NOL_bySample$piN)] <- 0
codon_results_NOL_bySample$piS <- (codon_results_NOL_bySample$S_diffs / codon_results_NOL_bySample$S_sites)
codon_results_NOL_bySample$piS[is.na(codon_results_NOL_bySample$piS)] <- 0
codon_results_NOL_bySample$piNpiS <- codon_results_NOL_bySample$piN / codon_results_NOL_bySample$piS

# test piN==piS
(codon_results_NOL_bySample_MEANPIN <- mean(codon_results_NOL_bySample$piN, na.rm = TRUE)) # 8.208488e-06
(codon_results_NOL_bySample_MEANPIS <- mean(codon_results_NOL_bySample$piS, na.rm = TRUE)) # 1.465037e-05
(codon_results_NOL_bySample_MEANPINPIS <- codon_results_NOL_bySample_MEANPIN / codon_results_NOL_bySample_MEANPIS) # 0.5602922
(codon_results_NOL_bySample_WILCOXP <- wilcox.test(codon_results_NOL_bySample$piN, codon_results_NOL_bySample$piS, paired = TRUE)) # p-value < 2.2e-16

###############################################################################
### IMPORT GROUP INFO
codon_results_NOL_bySample <- left_join(x = codon_results_NOL_bySample, y = df_meta %>% select(sample,Vaccine,lineage_sim), by = "sample")

# RESULTS BY OUTBREAK
(codon_results_NOL_bySample$pi <- codon_results_NOL_bySample$piN + codon_results_NOL_bySample$piS)
summary(codon_results_NOL_bySample$pi)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 0.000e+00 5.333e-06 2.286e-05 3.071e-05  3.970e-04 
write_csv(codon_results_NOL_bySample, "../results/codon_results_NOL_bySample.csv")

# SUMMARIZE BY OUTBREAK
(codon_results_NOL_bySample_bygroup <- codon_results_NOL_bySample %>%
  group_by(Vaccine, lineage_sim) %>%
  summarise(
    count = n(),
    
    min_pi = min(pi, na.rm = TRUE),
    median_pi = median(pi, na.rm = TRUE),
    mean_pi = mean(pi, na.rm = TRUE),
    max_pi = max(pi, na.rm = TRUE),
    sd_pi = sd(pi, na.rm = TRUE),
    
    min_piN = min(piN, na.rm = TRUE),
    median_piN = median(piN, na.rm = TRUE),
    mean_piN = mean(piN, na.rm = TRUE),
    max_piN = max(piN, na.rm = TRUE),
    sd_piN = sd(piN, na.rm = TRUE),
    
    min_piS = min(piS, na.rm = TRUE),
    median_piS = median(piS, na.rm = TRUE),
    mean_piS = mean(piS, na.rm = TRUE),
    max_piS = max(piS, na.rm = TRUE),
    sd_piS = sd(piS, na.rm = TRUE)
  ))

codon_results_NOL_bySample_bygroup$mean_piNpiS <- codon_results_NOL_bySample_bygroup$mean_piN / codon_results_NOL_bySample_bygroup$mean_piS
write_csv(codon_results_NOL_bySample_bygroup, "../results/codon_results_NOL_bySample_bygroup.csv")


# pi, between groups
source("./helper/df_test.R")
## already Done in Figure 2 (./3.4.Plots_on_isnvs.R, # plot for pi)

#piN, piS
df_wilc_test <- cal_wilc_test(codon_results_NOL_bySample %>% mutate(gene="Full genome"), "piN", genes=c("Full genome"))
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test, "../results/df_test_diversity_piN_by_vaccine_gene_more.xlsx")

df_wilc_test <- cal_wilc_test(codon_results_NOL_bySample %>% mutate(gene="Full genome"), "piS", genes=c("Full genome"))
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test, "../results/df_test_diversity_piS_by_vaccine_gene_more.xlsx")

df_wilc_test <- cal_wilc_test(codon_results_NOL_bySample %>% mutate(gene="Full genome"), "piNpiS", genes=c("Full genome"))
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test, "../results/df_test_diversity_piNpiS_by_vaccine_gene_more.xlsx")

codon_results_NOL_bySample$piN_m_piS <- codon_results_NOL_bySample$piN-codon_results_NOL_bySample$piS
df_wilc_test <- cal_wilc_test(codon_results_NOL_bySample %>% mutate(gene="Full genome"), "piN_m_piS", genes=c("Full genome"))
df_wilc_test <- highlight_diff(df_wilc_test)
write_xlsx(df_wilc_test, "../results/df_test_diversity_piN_m_piS_by_vaccine_gene_more.xlsx")

###############################################################################
### SUMMARIZE RESULTS BY CODON MEANS, NOL ONLY
unique(codon_results$product_segment) # 
(codon_results_NOL_byProductCodon <- filter(codon_results, num_overlap_ORF_nts == 0) %>%
    #group_by(product, codon_num) %>%
   group_by(product_segment, codon_num) %>%
    summarise(
      N_sites = mean(N_sites, na.rm = TRUE),
      S_sites = mean(S_sites, na.rm = TRUE),
      N_diffs = mean(N_diffs, na.rm = TRUE),
      S_diffs = mean(S_diffs, na.rm = TRUE)
    ))
names(codon_results_NOL_byProductCodon)[names(codon_results_NOL_byProductCodon) == "product_segment"] <- "product"


###############################################################################
# BOOTSTRAP PROCESS

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
codon_results_NOL_byProductCodon$num_defined_seqs <- 6 # seems this variable is for other purpose, so it is set fixed in the script, comment from Haogao

### SUMMARIZE RESULTS BY GENE
(codon_results_NOL_byProductCodon_summary <- codon_results_NOL_byProductCodon %>%
    group_by(product) %>%
    summarise(
      N_sites = sum(N_sites, na.rm = TRUE),
      S_sites = sum(S_sites, na.rm = TRUE),
      N_diffs = sum(N_diffs, na.rm = TRUE),
      S_diffs = sum(S_diffs, na.rm = TRUE)
    )
)

codon_results_NOL_byProductCodon_summary$dN <- codon_results_NOL_byProductCodon_summary$N_diffs / codon_results_NOL_byProductCodon_summary$N_sites
codon_results_NOL_byProductCodon_summary$dS <- codon_results_NOL_byProductCodon_summary$S_diffs / codon_results_NOL_byProductCodon_summary$S_sites
codon_results_NOL_byProductCodon_summary$dNdS <- codon_results_NOL_byProductCodon_summary$dN / codon_results_NOL_byProductCodon_summary$dS
#View(codon_results_NOL_byProductCodon_summary)

(codon_results_NOL_byProductCodon_summary_LONG <- codon_results_NOL_byProductCodon_summary %>%
    pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value"))


####################################################################################################
### BOOTSTRAP EACH STATISTIC

############################################################################################################
# *BASIC* BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT
dNdS_diff_boot_fun <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
  
  # Function for dN
  dN_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    return(dN)
  }
  
  # Function for dN
  dS_function <- function(D, indices) {
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    return(dS)
  }
  
  # Function for dN - dS
  dN_m_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_m_dS <- dN - dS
    return(dN_m_dS)
  }
  
  # Function for dN/dS
  dN_over_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_over_dS <- dN / dS
    return(dN_over_dS)
  }
  
  # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
  
  (dN <- sum(as.vector(codon_results[ , paste0(numerator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(numerator, "_sites")]), na.rm = TRUE))
  (dS <- sum(as.vector(codon_results[ , paste0(denominator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(denominator, "_sites")]), na.rm = TRUE))
  (dNdS <- dN / dS)
  
  # Run the BOOTSTRAPS
  # boot dN
  (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
  (dN <- boot_dN$t0)
  (boot_dN_SE <- sd(boot_dN$t))
  
  # boot dS
  (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dS <- boot_dS$t0)
  (boot_dS_SE <- sd(boot_dS$t))
  
  # boot dN - dS
  (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_m_dS <- boot_dN_m_dS$t0)
  (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
  (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
  (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
  
  # boot dN/dS
  (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_over_dS <- boot_dN_over_dS$t0)
  (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
  (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
  (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
  
  ### NEW: ASL (acheived significance level)
  boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0) # 345
  boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0) # 0
  boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0) # 655
  ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
               boot_dN_m_dS_SE, boot_dN_m_dS_P, 
               boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
               sep = "\t"))
}


############################################################################################################
### ANALYSIS VARIABLES: same as before
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 10000
NCPUS <- 16


############################################################################################################
### INITIALIZE DATA FRAME: intrahost
intrahost_results_bootstrap <- data.frame(gene_name = character(),
                                          num_bootstraps = integer(),
                                          min_defined_codons = integer(),
                                          num_codons = integer(),
                                          N_sites = numeric(),
                                          S_sites = numeric(),
                                          N_diffs = numeric(),
                                          S_diffs = numeric(),
                                          num_replicates = integer(),
                                          dN = numeric(),
                                          dS = numeric(),
                                          dNdS = numeric(),
                                          dN_m_dS = numeric(),
                                          boot_dN_SE = numeric(),
                                          boot_dS_SE = numeric(),
                                          boot_dN_over_dS_SE = numeric(),
                                          boot_dN_over_dS_P = numeric(),
                                          boot_dN_m_dS_SE = numeric(),
                                          P_value = numeric(),
                                          boot_dN_gt_dS_count = integer(), 
                                          boot_dN_eq_dS_count = integer(), 
                                          boot_dN_lt_dS_count = integer(), 
                                          ASL_dN_gt_dS_P = numeric(), 
                                          ASL_dN_lt_dS_P = numeric())


### LOOP EACH DATA SUBSET (4 minutes with 6 CPUs)
(genes_all <- c("Full genome", sort(unique(codon_results_NOL_byProductCodon$product))))

for (this_gene in genes_all) {
	print(this_gene)
    # this_gene = "Full genome"
    # Filter by gene, frame, and minimum number of defined codons
	if (this_gene == "Full genome") {
	    this_data <- filter(codon_results_NOL_byProductCodon, num_defined_seqs >= MIN_DEFINED_CODONS) 
	} else {
	    this_data <- filter(codon_results_NOL_byProductCodon, product == this_gene, num_defined_seqs >= MIN_DEFINED_CODONS) 
	}
    
    if(nrow(this_data) >= 2) {
      # LEADING SUMMARY COLUMNS:
      N_sites <- sum(this_data$N_sites, na.rm = T)
      S_sites <- sum(this_data$S_sites, na.rm = T)
      N_diffs <- sum(this_data$N_diffs, na.rm = T)
      S_diffs <- sum(this_data$S_diffs, na.rm = T)
      
      summary_data <- paste(nrow(this_data),
                            N_sites, S_sites, 
                            N_diffs, S_diffs, 
                            sep = "\t")
      
      # BOOTSTRAP THE ONE RATIO
      boot_dNdS <- dNdS_diff_boot_fun(this_data, 'N', 'S', NBOOTSTRAPS, NCPUS)
      
      # RECORD HEADER
      boot_vector_names <- c('num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                             'num_replicates', 
                             'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                             'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
      
      boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
      
      # Add names
      names(boot_dNdS_vector) <- boot_vector_names
      
      # Prepare additional rows
      intrahost_results_bootstrap_ADDITION <- data.frame(gene_name = this_gene,
                                                         num_bootstraps = NBOOTSTRAPS,
                                                         min_defined_codons = MIN_DEFINED_CODONS,
                                                         num_codons = as.integer(boot_dNdS_vector['num_codons']),
                                                         N_sites = as.numeric(boot_dNdS_vector['N_sites']),
                                                         S_sites = as.numeric(boot_dNdS_vector['S_sites']),
                                                         N_diffs = as.numeric(boot_dNdS_vector['N_diffs']),
                                                         S_diffs = as.numeric(boot_dNdS_vector['S_diffs']),
                                                         num_replicates = as.integer(boot_dNdS_vector['num_replicates']),
                                                         dN = as.numeric(boot_dNdS_vector['dN']),
                                                         dS = as.numeric(boot_dNdS_vector['dS']),
                                                         dNdS = as.numeric(boot_dNdS_vector['dNdS']),
                                                         dN_m_dS = as.numeric(boot_dNdS_vector['dN_m_dS']),
                                                         boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE']),
                                                         boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE']),
                                                         boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE']),
                                                         boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P']),
                                                         boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE']),
                                                         P_value = as.numeric(boot_dNdS_vector['P_value']),
                                                         boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count']),
                                                         boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count']),
                                                         boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count']),
                                                         ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P']),
                                                         ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P']))
      
      # Add the 4 new rows to results
      intrahost_results_bootstrap <- rbind(intrahost_results_bootstrap, intrahost_results_bootstrap_ADDITION)
    }
}
# warnings: NAs introduced by coercion (OKAY)

names(intrahost_results_bootstrap) <- c('gene_name',
                                        'num_bootstraps', 'min_defined_codons', 'num_codons', 
                                        'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                        'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                        'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                        'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')


#View(intrahost_results_bootstrap)

### Manual 2-sided ASL P-value
intrahost_results_bootstrap$P_ALS <- NA
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P < intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P < intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_gt_dS_P
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P > intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P > intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_lt_dS_P
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$P_ALS) & intrahost_results_bootstrap$P_ALS == 0, ]$P_ALS <- 1 / NBOOTSTRAPS

### FDR
intrahost_results_bootstrap$Q_ASL_BH <- p.adjust(intrahost_results_bootstrap$P_ALS, method = "BH")
intrahost_results_bootstrap$Q_Z_BH <- p.adjust(intrahost_results_bootstrap$P_value, method = "BH")

### SAVE
write_tsv(intrahost_results_bootstrap, "../results/intrahost_results_bootstrap_products.tsv")

#View(intrahost_results_bootstrap)


###############################################################################
### SUMMARIZE RESULTS BY CODON MEANS, NOL ONLY, NOW BY OUTBREAK!

# JOIN outbreak metadata
codon_results$outbreak <- paste(codon_results$Vaccine, codon_results$lineage_sim)
unique(codon_results$outbreak)

# SUMMARIZE
(codon_results_NOL_bygroupProductCodon <- filter(codon_results, num_overlap_ORF_nts == 0) %>%
    group_by(outbreak, product_segment, codon_num) %>%
    summarise(
      N_sites = mean(N_sites, na.rm = TRUE),
      S_sites = mean(S_sites, na.rm = TRUE),
      N_diffs = mean(N_diffs, na.rm = TRUE),
      S_diffs = mean(S_diffs, na.rm = TRUE)
    ))
names(codon_results_NOL_bygroupProductCodon)[names(codon_results_NOL_bygroupProductCodon) == "product_segment"] <- "product"


###############################################################################
# BOOTSTRAP PROCESS

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
codon_results_NOL_bygroupProductCodon$num_defined_seqs <- 6

### SUMMARIZE RESULTS BY OUTBREAK AND GENE
(codon_results_NOL_bygroupProductCodon_summary <- codon_results_NOL_bygroupProductCodon %>%
    group_by(outbreak, product) %>%
    summarise(
      N_sites = sum(N_sites, na.rm = TRUE),
      S_sites = sum(S_sites, na.rm = TRUE),
      N_diffs = sum(N_diffs, na.rm = TRUE),
      S_diffs = sum(S_diffs, na.rm = TRUE)
    )
)

codon_results_NOL_bygroupProductCodon_summary$dN <- codon_results_NOL_bygroupProductCodon_summary$N_diffs / codon_results_NOL_bygroupProductCodon_summary$N_sites
codon_results_NOL_bygroupProductCodon_summary$dS <- codon_results_NOL_bygroupProductCodon_summary$S_diffs / codon_results_NOL_bygroupProductCodon_summary$S_sites
codon_results_NOL_bygroupProductCodon_summary$dNdS <- codon_results_NOL_bygroupProductCodon_summary$dN / codon_results_NOL_bygroupProductCodon_summary$dS
#View(codon_results_NOL_bygroupProductCodon_summary)

(codon_results_NOL_bygroupProductCodon_summary_LONG <- codon_results_NOL_bygroupProductCodon_summary %>%
    pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value"))


####################################################################################################
### BOOTSTRAP EACH STATISTIC: *SAME* BOOTSTRAP FUNCTION (dN - dS) as above

############################################################################################################
### ANALYSIS VARIABLES: same as before
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 10000
NCPUS <- 16


############################################################################################################
### INITIALIZE DATA FRAME: intrahost
intrahost_bygroup_results_bootstrap <- data.frame(gene_name = character(),
                                          num_bootstraps = integer(),
                                          min_defined_codons = integer(),
                                          num_codons = integer(),
                                          N_sites = numeric(),
                                          S_sites = numeric(),
                                          N_diffs = numeric(),
                                          S_diffs = numeric(),
                                          num_replicates = integer(),
                                          dN = numeric(),
                                          dS = numeric(),
                                          dNdS = numeric(),
                                          dN_m_dS = numeric(),
                                          boot_dN_SE = numeric(),
                                          boot_dS_SE = numeric(),
                                          boot_dN_over_dS_SE = numeric(),
                                          boot_dN_over_dS_P = numeric(),
                                          boot_dN_m_dS_SE = numeric(),
                                          P_value = numeric(),
                                          boot_dN_gt_dS_count = integer(), 
                                          boot_dN_eq_dS_count = integer(), 
                                          boot_dN_lt_dS_count = integer(), 
                                          ASL_dN_gt_dS_P = numeric(), 
                                          ASL_dN_lt_dS_P = numeric(), 
                                          outbreak = character()) ## ADD OUTBREAK


### LOOP EACH DATA SUBSET 
for (this_outbreak in sort(unique(codon_results_NOL_bygroupProductCodon$outbreak))) {
  	print("--------------")	
	print(this_outbreak)
  	print("--------------")
  for (this_gene in genes_all) {
	print(this_gene)
    # Filter by gene, frame, and minimum number of defined codons
	if (this_gene == "Full genome") {
	    this_data <- filter(codon_results_NOL_bygroupProductCodon, outbreak == this_outbreak, num_defined_seqs >= MIN_DEFINED_CODONS) 
	} else {
	    this_data <- filter(codon_results_NOL_bygroupProductCodon, product == this_gene, outbreak == this_outbreak, num_defined_seqs >= MIN_DEFINED_CODONS) 
	}
    
    if(nrow(this_data) >= 2) {
      # LEADING SUMMARY COLUMNS:
      N_sites <- sum(this_data$N_sites, na.rm = T)
      S_sites <- sum(this_data$S_sites, na.rm = T)
      N_diffs <- sum(this_data$N_diffs, na.rm = T)
      S_diffs <- sum(this_data$S_diffs, na.rm = T)
      
      summary_data <- paste(nrow(this_data),
                            N_sites, S_sites, 
                            N_diffs, S_diffs, 
                            sep = "\t")
      
      # BOOTSTRAP THE ONE RATIO
      boot_dNdS <- dNdS_diff_boot_fun(this_data, 'N', 'S', NBOOTSTRAPS, NCPUS)
      
      # RECORD HEADER
      boot_vector_names <- c('num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                             'num_replicates', 
                             'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                             'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
      
      boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
      
      # Add names
      names(boot_dNdS_vector) <- boot_vector_names
      
      # Prepare additional rows
      intrahost_bygroup_results_bootstrap_ADDITION <- data.frame(gene_name = this_gene,
                                                                    num_bootstraps = NBOOTSTRAPS,
                                                                    min_defined_codons = MIN_DEFINED_CODONS,
                                                                    num_codons = as.integer(boot_dNdS_vector['num_codons']),
                                                                    N_sites = as.numeric(boot_dNdS_vector['N_sites']),
                                                                    S_sites = as.numeric(boot_dNdS_vector['S_sites']),
                                                                    N_diffs = as.numeric(boot_dNdS_vector['N_diffs']),
                                                                    S_diffs = as.numeric(boot_dNdS_vector['S_diffs']),
                                                                    num_replicates = as.integer(boot_dNdS_vector['num_replicates']),
                                                                    dN = as.numeric(boot_dNdS_vector['dN']),
                                                                    dS = as.numeric(boot_dNdS_vector['dS']),
                                                                    dNdS = as.numeric(boot_dNdS_vector['dNdS']),
                                                                    dN_m_dS = as.numeric(boot_dNdS_vector['dN_m_dS']),
                                                                    boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE']),
                                                                    boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE']),
                                                                    boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE']),
                                                                    boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P']),
                                                                    boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE']),
                                                                    P_value = as.numeric(boot_dNdS_vector['P_value']),
                                                                    boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count']),
                                                                    boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count']),
                                                                    boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count']),
                                                                    ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P']),
                                                                    ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P']),
                                                                    outbreak = this_outbreak)
      
      # Add the 4 new rows to results
      intrahost_bygroup_results_bootstrap <- rbind(intrahost_bygroup_results_bootstrap, intrahost_bygroup_results_bootstrap_ADDITION)
    }
  }
} # end this outbreak
# warnings: NAs introduced by coercion (OKAY)

names(intrahost_bygroup_results_bootstrap) <- c('gene_name',
                                        'num_bootstraps', 'min_defined_codons', 'num_codons', 
                                        'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                        'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                        'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                        'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P',
                                        "outbreak") ## ADDED THIS


#View(intrahost_bygroup_results_bootstrap)

### Manual 2-sided ASL P-value
intrahost_bygroup_results_bootstrap$P_ALS <- NA
intrahost_bygroup_results_bootstrap[! is.na(intrahost_bygroup_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_bygroup_results_bootstrap$ASL_dN_gt_dS_P < intrahost_bygroup_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_bygroup_results_bootstrap[! is.na(intrahost_bygroup_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_bygroup_results_bootstrap$ASL_dN_gt_dS_P < intrahost_bygroup_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_gt_dS_P
intrahost_bygroup_results_bootstrap[! is.na(intrahost_bygroup_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_bygroup_results_bootstrap$ASL_dN_gt_dS_P > intrahost_bygroup_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_bygroup_results_bootstrap[! is.na(intrahost_bygroup_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_bygroup_results_bootstrap$ASL_dN_gt_dS_P > intrahost_bygroup_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_lt_dS_P
intrahost_bygroup_results_bootstrap[! is.na(intrahost_bygroup_results_bootstrap$P_ALS) & intrahost_bygroup_results_bootstrap$P_ALS == 0, ]$P_ALS <- 1 / NBOOTSTRAPS

### FDR
intrahost_bygroup_results_bootstrap$Q_ASL_BH <- p.adjust(intrahost_bygroup_results_bootstrap$P_ALS, method = "BH")
intrahost_bygroup_results_bootstrap$Q_Z_BH <- p.adjust(intrahost_bygroup_results_bootstrap$P_value, method = "BH")

### SAVE
write_tsv(intrahost_bygroup_results_bootstrap, "../results/intrahost_bygroup_results_bootstrap_products.tsv")

#View(intrahost_bygroup_results_bootstrap)



####################################################################################################
# RELOAD ALL
intrahost_results_bootstrap <- read_tsv("../results/intrahost_results_bootstrap_products.tsv")
intrahost_results_bootstrap$outbreak <- "Combined" # add this now
intrahost_bygroup_results_bootstrap <- read_tsv("../results/intrahost_bygroup_results_bootstrap_products.tsv")

# COMBINE
intrahost_results_bootstrap <- rbind(intrahost_results_bootstrap, intrahost_bygroup_results_bootstrap)

### CONVERT TO LONG
intrahost_results_bootstrap_LONG <- intrahost_results_bootstrap %>%
  pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value")
# so the 2 measures of d (dN and dS) are on different lines

#View(intrahost_results_bootstrap_LONG)

####################################################################################################
### PLOT WHOLE GENES, VARIOUS RATIOS
intrahost_results_bootstrap_LONG$d_measure <- factor(intrahost_results_bootstrap_LONG$d_measure, levels = c('dN', 'dS'))

# Add error bars
intrahost_results_bootstrap_LONG$d_SE_min <- NaN
intrahost_results_bootstrap_LONG$d_SE_max <- NaN
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_min <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_value - intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_max <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_value + intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_min <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_value - intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_max <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_value + intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_SE_min < 0 & ! is.na(intrahost_results_bootstrap_LONG$d_SE_min), ]$d_SE_min <- 0 # OK if fail because none negative

# significance for Q-value
intrahost_results_bootstrap_LONG$significance <- ""
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$Q_Z_BH) & intrahost_results_bootstrap_LONG$Q_Z_BH < 0.05, ]$significance <- '*'
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$Q_Z_BH) & intrahost_results_bootstrap_LONG$Q_Z_BH < 0.01, ]$significance <- '**'
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$Q_Z_BH) & intrahost_results_bootstrap_LONG$Q_Z_BH < 0.001, ]$significance <- '***'
intrahost_results_bootstrap_LONG$significance[intrahost_results_bootstrap_LONG$significance==""] <- NA

# ORDER AND NAME THE GENES
unique(intrahost_results_bootstrap_LONG$gene_name)
unique(codon_results$product_segment)

#gene_ids_sorted <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
unique(intrahost_results_bootstrap_LONG$gene_name)
gene_ids_sorted <- c("Full genome","ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N", "ORF10")

#gene_names_sorted <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
gene_names_sorted <- gene_ids_sorted

intrahost_results_bootstrap_LONG$gene_name <- factor(intrahost_results_bootstrap_LONG$gene_name, 
                                                     levels = gene_ids_sorted, 
                                                     labels = gene_names_sorted)

unique(intrahost_results_bootstrap_LONG$outbreak)
names(colors_vaccine)
names(colors_lineage)
df_tmp <- df_meta %>% group_by(lineage_sim,Vaccine) %>% summarise(N=n()) %>% arrange(Vaccine, lineage_sim)
sum(df_tmp$N)
(labels <- paste0(df_tmp$Vaccine, " ", df_tmp$lineage_sim, "\n(N = ", df_tmp$N, ")"))
labels <- gsub("^Non-", "Un", labels)
intrahost_results_bootstrap_LONG$outbreak_ori <- intrahost_results_bootstrap_LONG$outbreak 
# intrahost_results_bootstrap_LONG$outbreak <- intrahost_results_bootstrap_LONG$outbreak_ori
intrahost_results_bootstrap_LONG$outbreak <- factor(intrahost_results_bootstrap_LONG$outbreak, levels = c("Combined", "BioNTech Delta", "BioNTech Omicron", "Sinovac Delta", "Sinovac Omicron", "Non-vaccinated Alpha", "Non-vaccinated Delta", "Non-vaccinated Omicron", "Non-vaccinated B.1.36", "Non-vaccinated B.1.36.27", "Non-vaccinated B.1.1.63"), labels=c("Combined\n(N=2053)", labels))


# For horizontal lines
horizontal_lines_data <- filter(intrahost_results_bootstrap_LONG, significance %in% c('*', '**', '***')) %>%
  group_by(gene_name) %>% # frame
  summarise(
    y_value = max(d_SE_max)
  )

intrahost_error_bar_colors <- dplyr::select(intrahost_results_bootstrap_LONG, outbreak, gene_name, d_measure, d_value, d_SE_min, d_SE_max) # <-- CHANGE AS NEEDED
intrahost_error_bar_colors$this_color <- ""
intrahost_error_bar_colors[intrahost_error_bar_colors$d_measure == 'dN', ]$this_color <- 'pink'
intrahost_error_bar_colors[intrahost_error_bar_colors$d_measure == 'dS', ]$this_color <- 'lightblue'
intrahost_error_bar_colors$d_measure <- factor(intrahost_error_bar_colors$d_measure, levels = c('dN', 'dS'))

# Boxes for all based on which π is higher


### PREPARE PLOT

### normalized dNdS based on difference
intrahost_results_bootstrap_LONG$dNdS_norm <- NaN
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$dNdS), ]$dNdS_norm <- 
  intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$dNdS), ]$dN_m_dS / 
  max(abs(intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$dNdS), ]$dN_m_dS))

### PLOT diversity chart ###
gene_names_allLevels <- gene_names_sorted # old: c('1ab', 'S', '3a', "3c", 'E', 'M', '6', '7a', '7b', '8', 'N', '9b', '9c', '10') 
pi_corr_factor <- 10000
max_dNdS <- max(intrahost_results_bootstrap_LONG$dNdS, na.rm = TRUE)
max_d_SE_max <- max(intrahost_results_bootstrap_LONG$d_SE_max, na.rm = TRUE)

### PLOT
intrahost_results_bootstrap_LONG$outbreak

n_color <- nrow(intrahost_results_bootstrap_LONG)/2
(intrahost_allProductsbygroup_PLOT <- ggplot(data = filter(intrahost_results_bootstrap_LONG, gene_name %in% gene_names_allLevels),
                                    mapping = aes(x = gene_name, y = d_value * pi_corr_factor, group = d_measure)) + #fill = d_measure)) +
    # Backdrop boxes based on which pi is higher
    geom_bar(data = filter(intrahost_results_bootstrap_LONG, d_measure == 'dN', gene_name %in% gene_names_allLevels), mapping = aes(y = Inf, fill = dNdS_norm), stat = 'identity') +
    
    geom_errorbar(data = filter(intrahost_error_bar_colors, gene_name %in% gene_names_allLevels), 
                  mapping =  aes(ymin = d_SE_min * pi_corr_factor, ymax = d_SE_max * pi_corr_factor), # , color = d_measure
                  color = rep(c('pink', 'lightblue'), n_color),
                  position = position_dodge(width = 0.5), width = 0, size = 1) +
    geom_point(data = filter(intrahost_error_bar_colors, gene_name %in% gene_names_allLevels), 
               mapping =  aes(y = d_SE_min * pi_corr_factor), # , color = d_measure
               color = rep(c('pink', 'lightblue'), n_color),
               position = position_dodge(width = 0.5), size = 0.29) + # size = 1
    geom_point(data = filter(intrahost_error_bar_colors, gene_name %in% gene_names_allLevels), 
               mapping =  aes(y = d_SE_max * pi_corr_factor), # color = d_measure
               color = rep(c('pink', 'lightblue'), n_color),
               position = position_dodge(width = 0.5), size = 0.29) +
    geom_point(stat = 'identity', position = position_dodge(width = 0.5), pch = 21, size = 2, stroke = 0, 
               fill = rep(c(brewer.pal(9, 'Set1')[1], brewer.pal(9, 'Set1')[2]), n_color), mapping=aes(color = d_measure)) + 
    
    facet_grid(outbreak ~ ., scales = 'free', space = 'free_x',) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
        #   legend.position = 'bottom',
          legend.title = element_blank(),
        #   axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), # angle = 90, hjust = 1, vjust = 1
          #axis.text.x = element_text(size = 9, face = "italic"), # angle = 90, hjust = 1, vjust = 1
          axis.text.y = element_text(size = 8),
          #axis.ticks.y.right = element_line(colour = '#DCB118'),
          #axis.text.y.right = element_text(color = '#DCB118'),
          axis.title.y = element_text(size = 8),
          panel.border = element_rect(),
          strip.text = element_text(size = 8),
        #   strip.text.y = element_blank(),
		  strip.text.y = element_text(angle = 0),
          #legend.background = element_rect(fill = 'white', color = 'black'),
          strip.background = element_blank()) +
    xlab("") + 
    ylab(bquote('Differences per site ('*'x 10'^'-4'*')')) + 
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    scale_y_continuous(breaks = scales::pretty_breaks(3), expand = expand_scale(mult = c(0, 0.1))) + # 0.1
    scale_fill_gradient2(low = brewer.pal(9, "Blues")[5], mid = 'white', high = brewer.pal(9, "Reds")[5], midpoint = 0, na.value = "white")) 
# ggsave("../results/helper.pdf")
# save_pptx("../results/helper.pptx")

# SAVE SOURCE
intrahost_results_bootstrap_LONG$significance_P <- ""
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$P_value) & intrahost_results_bootstrap_LONG$P_value < 0.1, ]$significance_P <- '*'
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$P_value) & intrahost_results_bootstrap_LONG$P_value < 0.05, ]$significance_P <- '**'
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$P_value) & intrahost_results_bootstrap_LONG$P_value < 0.01, ]$significance_P <- '***'
write_xlsx(intrahost_results_bootstrap_LONG, "../results/intrahost_allProductsbygroup_PLOT_SOURCE.xlsx")


df_significance_label <- intrahost_results_bootstrap_LONG %>% group_by(outbreak) %>% mutate(group_max=max(d_SE_max)) %>% group_by(gene_name, outbreak) %>% summarise(significance_P=significance_P[1], d_SE_max = group_max*0.9, d_measure=d_measure[1])
p_out <- intrahost_allProductsbygroup_PLOT + geom_text(aes(x=gene_name, y=d_SE_max * pi_corr_factor, label=significance_P), data=df_significance_label)
# SAVE PLOT
ggsave(filename = "../results/intrahost_allProductsbygroup_PLOT.pdf", width = 6, height = 6*sqrt(2)-1, plot=p_out)
save_pptx(file = "../results/intrahost_allProductsbygroup_PLOT.pptx", width = 6, height = 6*sqrt(2)-1,  plot=p_out)


# PIVOT WIDE AND SAVE
intrahost_results_bootstrap_LONG_WIDE <- pivot_wider(data = dplyr::select(intrahost_results_bootstrap_LONG, -significance, -dNdS_norm, -significance_P), 
                                                                names_from = d_measure, values_from = c(d_value, d_SE_min, d_SE_max))
write_xlsx(intrahost_results_bootstrap_LONG_WIDE, "../results/intrahost_results_bootstrap_LONG_WIDE.xlsx")
View(intrahost_results_bootstrap_LONG_WIDE)

# Both ratios above/below 1?
(ratio_direction_by_outbreak <- filter(intrahost_results_bootstrap, outbreak != "Combined") %>%
  group_by(gene_name) %>%
  summarise(
    pos = sum(dNdS > 1.3),
    neg = sum(dNdS < 0.7)
  ))
ratio_direction_by_outbreak$diff <- FALSE
ratio_direction_by_outbreak[! is.na(ratio_direction_by_outbreak$pos) & ratio_direction_by_outbreak$pos > 0 & 
                              ! is.na(ratio_direction_by_outbreak$neg) & ratio_direction_by_outbreak$neg > 0, ]$diff <- TRUE
#View(ratio_direction_by_outbreak)


###############################################################################
### SLIDING WINDOWS

# SAVE CODON-MEAN RESULTS for each PRODUCT; INCLUDE OL and NOL codons

###############################################################################
### Per SAMPLE per GENE: SUMMARIZE RESULTS BY CODON MEANS
(codon_results_byProductCodon <- filter(codon_results) %>%
   group_by(product_segment, codon_num) %>%
   summarise(
     N_sites = mean(N_sites, na.rm = TRUE),
     S_sites = mean(S_sites, na.rm = TRUE),
     N_diffs = mean(N_diffs, na.rm = TRUE),
     S_diffs = mean(S_diffs, na.rm = TRUE)
   ))

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
codon_results_byProductCodon$num_defined_seqs <- 6
codon_results_byProductCodon$product <- codon_results_byProductCodon$product_segment
# LOOP AND DO
dir.create("../results/intrahost_sliding_windows/")
system("chmod 755 ./SNPGenie_sliding_windows.R")
mclapply(seq(10, 50, 10), function(size) {
	print("############")
	print(size)
	print("############")
	for (this_product in unique(codon_results_byProductCodon$product_segment)) {
		print(this_product)
		file_tsv <- paste0("../results/intrahost_sliding_windows/codon_results_byProductCodon_", this_product, "_size_", size, ".tsv")
		file_out <- paste0("../results/intrahost_sliding_windows/SNPGenie_sliding_windows_", this_product, "_size_", size, ".out")
		write_tsv(filter(codon_results_byProductCodon, product_segment == this_product), file_tsv)
		cmd <- paste0("./SNPGenie_sliding_windows.R ", file_tsv, " N S ", size, " 1 1000 6 NONE 8 > ", file_out)
		system(cmd)
	}
},mc.cores=8)


###############################################################################
###############################################################################
### PLOT sliding windows
mclapply(seq(10, 50, 10), function(WINDOW_SIZE) {
	# WINDOW_SIZE <- 30
	# Combine files :
	# IN: ../results/intrahost_sliding_windows_size30
	# DO: cat codon_results_byProductCodon_*.tsv | grep -v product_segment > codon_results_byProductCodon_WINDOWS_POOLED.tsv
	files_product_window_rst <- list.files("../results/intrahost_sliding_windows/", pattern=paste0("size_", WINDOW_SIZE, "_WINDOWS"), full.names=T)
	codon_results_byProductCodon_WINDOWS_POOLED <- lapply(files_product_window_rst, read_tsv, 
		# col_names = c("product", "codon_num", "N_sites", "S_sites", "N_diffs", "S_diffs", "num_defined_seqs", "sw_ratio", 
        #                  "sw_start", "sw_center", "sw_end", "sw_num_replicates", "sw_N_diffs", "sw_S_diffs", "sw_N_sites", "sw_S_sites", 
        #                  "sw_dN", "sw_dS", "sw_dNdS", "sw_dN_m_dS", "sw_boot_dN_SE", "sw_boot_dS_SE", "sw_boot_dN_over_dS_SE", "sw_boot_dN_over_dS_P", 
        #                  "sw_boot_dN_m_dS_SE", "sw_boot_dN_m_dS_P", "sw_boot_dN_gt_dS_count", "sw_boot_dN_eq_dS_count", "sw_boot_dN_lt_dS_count", 
        #                  "sw_ASL_dN_gt_dS_P", "sw_ASL_dN_lt_dS_P", "sw_ASL_dNdS_P")
		)
	codon_results_byProductCodon_WINDOWS_POOLED <- do.call(rbind, codon_results_byProductCodon_WINDOWS_POOLED)

	# Renumber codons starting at 1 for each product
	#codon_results_byProductCodon_WINDOWS_POOLED$codon_num_ORF <- codon_results_byProductCodon_WINDOWS_POOLED$codon_num
	codon_results_byProductCodon_WINDOWS_POOLED$sw_center_ORF <- codon_results_byProductCodon_WINDOWS_POOLED$sw_center
	for (this_product in unique(codon_results_byProductCodon_WINDOWS_POOLED$product)) {
	#this_product <- "E"
	
	# codon
	codon_results_byProductCodon_WINDOWS_POOLED[codon_results_byProductCodon_WINDOWS_POOLED$product == this_product, ]$codon_num <- 
		codon_results_byProductCodon_WINDOWS_POOLED[codon_results_byProductCodon_WINDOWS_POOLED$product == this_product, ]$codon_num - 
		min(codon_results_byProductCodon_WINDOWS_POOLED[codon_results_byProductCodon_WINDOWS_POOLED$product == this_product, ]$codon_num) + 1
	}

	# Renumber sw_center
	codon_results_byProductCodon_WINDOWS_POOLED$sw_center <- codon_results_byProductCodon_WINDOWS_POOLED$codon_num + ((WINDOW_SIZE - 1) / 2)
	codon_results_byProductCodon_WINDOWS_POOLED[is.na(codon_results_byProductCodon_WINDOWS_POOLED$sw_start), ]$sw_center <- NA
	#View(codon_results_byProductCodon_WINDOWS_POOLED)

	# # Rename nsp12
	# codon_results_byProductCodon_WINDOWS_POOLED[codon_results_byProductCodon_WINDOWS_POOLED$product %in% c("nsp12_1", "nsp12_2"), ]$product <- "nsp12"

	# Make a long form for plotting dN and dS
	(codon_results_byProductCodon_WINDOWS_POOLED_LONG <- codon_results_byProductCodon_WINDOWS_POOLED %>%
		pivot_longer(cols = c('sw_dN', 'sw_dS'), names_to = "d_measure", values_to = "d_value"))

	### DEFINE products
	(uniq_products_LONG <- unique(codon_results_byProductCodon_WINDOWS_POOLED_LONG$product)) # nsp not included because only 8 codons!
	# "E", "M", "N", "ORF10", "ORF3a", "ORF3c", "ORF3d", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "ORF9c", "S", "SiORF1", "SiORF2", "nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", "nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", 
	gene_names_sorted # recall
	#"ORF1ab" "S"      "ORF3a"  "E"      "M"      "ORF6"   "ORF7a"  "ORF7b"  "ORF8"   "N"      "ORF10" 

	uniq_products_LONG_SORTED <- gene_names_sorted
	uniq_products_LONG_SORTED_KEEPERS <- gene_names_sorted

	# FACTOR products
	codon_results_byProductCodon_WINDOWS_POOLED_LONG$product <- factor(codon_results_byProductCodon_WINDOWS_POOLED_LONG$product, 
																	levels = uniq_products_LONG_SORTED)

	###############################################################################
	### SUMMARIZE piN and piS for each product (NOT LONG)
	(codon_results_byProductCodon_byProduct <- filter(codon_results_byProductCodon_WINDOWS_POOLED) %>%
	group_by(product) %>%
	summarise(
		N_sites = sum(N_sites, na.rm = TRUE),
		S_sites = sum(S_sites, na.rm = TRUE),
		N_diffs = sum(N_diffs, na.rm = TRUE),
		S_diffs = sum(S_diffs, na.rm = TRUE)
	))

	#View(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10"))
	#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10")$S_diffs) # 0.008516171
	#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10")$S_sites) # 26.33383
	#0.008516171/26.33383 # 0.0003233928

	#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10", codon_num >= 1, codon_num <= 30)$S_diffs) # 0.008516171
	#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10", codon_num >= 1, codon_num <= 30)$S_sites) # 21.0009
	#0.008516171/21.0009 # 0.0004055146

	#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10", codon_num >= 10, codon_num <= 39)$S_diffs) # 0.008516171
	#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10", codon_num >= 10, codon_num <= 39)$S_sites) # 20.66626
	#0.008516171/20.66626 # 0.0004120809

	# UNDERSTAND IT NOW! EXPLANATION:
	# I was confused that the mean piS for all of ORF10 could be LOWER than piS for EVERY WINDOW. However, the only S differences are observed in
	# the middle of the gene. Thus, because the window size is 30 codons, every window overlaps the only S differences that exist. As a result,
	# the only effect of calculating the total piS for the whole gene is to add S sites but now S differences, and the denominator but not
	# numerator grows. This is why the total value of piS for the gene (horitzontal dashed line) can be lower than that of any given window.

	# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
	codon_results_byProductCodon_byProduct$num_defined_seqs <- 6

	# ADD STATS
	codon_results_byProductCodon_byProduct$dN <- codon_results_byProductCodon_byProduct$N_diffs / codon_results_byProductCodon_byProduct$N_sites
	codon_results_byProductCodon_byProduct$dS <- codon_results_byProductCodon_byProduct$S_diffs / codon_results_byProductCodon_byProduct$S_sites
	codon_results_byProductCodon_byProduct$dNdS <- codon_results_byProductCodon_byProduct$dN / codon_results_byProductCodon_byProduct$dS
	#View(codon_results_byProductCodon_byProduct)

	# FACTOR MEANS
	codon_results_byProductCodon_byProduct$product <- factor(codon_results_byProductCodon_byProduct$product, 
															levels = uniq_products_LONG_SORTED)

	# Create STANDARD ERROR COLUMN
	#View(codon_results_byProductCodon_WINDOWS_POOLED_LONG)
	codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_SE <- codon_results_byProductCodon_WINDOWS_POOLED_LONG$sw_boot_dN_SE
	codon_results_byProductCodon_WINDOWS_POOLED_LONG[codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_measure == "sw_dS", ]$d_value_SE <- 
	codon_results_byProductCodon_WINDOWS_POOLED_LONG[codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_measure == "sw_dS", ]$sw_boot_dS_SE
	codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_CI_min <- 
	codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value - codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_SE
	codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_CI_max <- 
	codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value + codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_SE

	# reset negative mins to 0
	codon_results_byProductCodon_WINDOWS_POOLED_LONG[!is.na(codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_CI_min) & 
													codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_CI_min < 0, ]$d_value_CI_min <- 0

	# Overlapping genes OLGs to mask
	OLGs_to_mask <- data.frame(
	product = c("S", "ORF3a", "ORF3a", "ORF3a", "N", "N"), # "nsp12", 
	name = c("S-iORF", "ORF3c", "ORF3d", "ORF3b", "ORF9b", "ORF9c"), # "nsp11", 
	start_codon = c(61, 22+(WINDOW_SIZE-22), 44, 141, 4+(WINDOW_SIZE-4), 154), # 8, 
	end_codon = c(101, 64, 102, 164, 102, 228), # 14, 
	sw_center = rep(0, 6), # 
	d_value = rep(0, 6), # 
	d_measure = rep("sw_dN", 6) # 
	)

	# PLOT
	pi_corr_factor2 <- 1000

	(codon_results_byProductCodon_WINDOWS_POOLED_LONG_PLOT <- ggplot(data = filter(codon_results_byProductCodon_WINDOWS_POOLED_LONG, product %in% uniq_products_LONG_SORTED_KEEPERS), mapping = aes(x = sw_center, y = pi_corr_factor2 * d_value, color = d_measure)) + 
		
		geom_rect(data = OLGs_to_mask, mapping = aes(xmin = start_codon, xmax = end_codon, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "grey", color = NA, 
				size = .3) + 
		geom_line(size = .3) +
		#geom_vline(xintercept = 614) + # to see S 614
		geom_hline(data = filter(codon_results_byProductCodon_byProduct, product %in% uniq_products_LONG_SORTED_KEEPERS),
				mapping = aes(yintercept = pi_corr_factor2 * dS), linetype = "dashed", color = brewer.pal(9, 'Set1')[2], 
				size = .3) +
		geom_hline(data = filter(codon_results_byProductCodon_byProduct, product %in% uniq_products_LONG_SORTED_KEEPERS),
				mapping = aes(yintercept = pi_corr_factor2 * dN), linetype = "dashed", color = brewer.pal(9, 'Set1')[1], 
				size = .3) +
		geom_ribbon(mapping = aes(ymin = pi_corr_factor2 * d_value_CI_min, ymax = pi_corr_factor2 * d_value_CI_max, fill = d_measure), alpha = 0.25, linetype = 0, 
					size = .3) +
		theme_bw() +
		theme(strip.background = element_blank(), 
			plot.title = element_text(hjust = 0.5),
			axis.ticks = element_line(size = 0.3),
			legend.title = element_blank(),
			legend.position = 'right',
			axis.text = element_text(size = 8),
			strip.text = element_text(face = c("bold")),
		  	strip.text.y = element_text(angle = 0),
			panel.grid = element_blank()) +
		xlab("Codon") + 
		ylab(bquote('Nucleotide diversity ('*'x 10'^'-3'*')')) + 
		scale_x_continuous(expand = c(0, 0)) +
		scale_y_continuous(breaks = scales::pretty_breaks(3), expand = c(0, 0)) + 
		# coord_flip()+
		# facet_grid(rows=vars(factor(product, levels=gene_names_sorted)), scales = "free", space="free") +
		facet_wrap(vars(factor(product, levels=gene_names_sorted)), scales = "free", ncol=2) +
		scale_color_manual(values = c(brewer.pal(9, 'Set1')[1], brewer.pal(9, 'Set1')[2]))
	)

	file_out <- paste0("../results/", "codon_results_byProductCodon_WINDOWS_POOLED_LONG_size", WINDOW_SIZE, "_PLOT.pdf")
	ggsave(file_out, width=8, height=10)

})


### CANDIDATE WINDOWS for 30 codons

# PIVOT to SHORT again
#View(codon_results_byProductCodon_WINDOWS_POOLED_LONG)
codon_results_byProductCodon_WINDOWS_POOLED_WIDE <- pivot_wider(data = codon_results_byProductCodon_WINDOWS_POOLED_LONG, 
                                                                 names_from = d_measure, values_from = c(d_value, d_value_SE, d_value_CI_min, d_value_CI_max))
#View(codon_results_byProductCodon_WINDOWS_POOLED_WIDE)

# JOIN the PRODUCT'S values
product_dS_values <- dplyr::select(codon_results_byProductCodon_byProduct, product, dS)
names(product_dS_values) <- c("product", "product_dS")
codon_results_byProductCodon_WINDOWS_POOLED_WIDE <- left_join(x = codon_results_byProductCodon_WINDOWS_POOLED_WIDE, product_dS_values, by = "product")

# Get the candidates
CANDIDATE_WINDOWS_SIZE30 <- filter(codon_results_byProductCodon_WINDOWS_POOLED_WIDE, 
                                   d_value_CI_min_sw_dN > d_value_CI_max_sw_dS, 
                                   d_value_CI_min_sw_dN > product_dS)
#View(CANDIDATE_WINDOWS_SIZE30)

# SAVE CANDIDATES
write_tsv(CANDIDATE_WINDOWS_SIZE30, "../results/CANDIDATE_WINDOWS_SIZE30.tsv")


################################################################################
# MANUALLY DEFINED START TO END OF CANDIDATE WINDOWS IN CANDIDATE_WINDOWS_SIZE30.xlsx, saved in CANDIDATE_WINDOWS_SIZE30_table.tsv
# LOOP and define/test each window; list codons with nonsyn diffs; etc.
CANDIDATE_WINDOWS_SIZE30_table <- read_tsv("../results/CANDIDATE_WINDOWS_SIZE30_table.tsv")
CANDIDATE_WINDOWS_SIZE30_table$N_diffs <- NaN
CANDIDATE_WINDOWS_SIZE30_table$S_diffs <- NaN
CANDIDATE_WINDOWS_SIZE30_table$N_sites <- NaN
CANDIDATE_WINDOWS_SIZE30_table$S_sites <- NaN
CANDIDATE_WINDOWS_SIZE30_table$N_diff_codons <- ""
CANDIDATE_WINDOWS_SIZE30_table$piN_max_codon <- NaN

# Additional columns
CANDIDATE_WINDOWS_SIZE30_table$num_bootstraps <- NaN
CANDIDATE_WINDOWS_SIZE30_table$num_codons <- NaN
#CANDIDATE_WINDOWS_SIZE30_table$N_sites_check <- NA
#CANDIDATE_WINDOWS_SIZE30_table$S_sites_check <- NA
#CANDIDATE_WINDOWS_SIZE30_table$N_diffs_check <- NA
#CANDIDATE_WINDOWS_SIZE30_table$S_diffs_check <- NA
CANDIDATE_WINDOWS_SIZE30_table$num_replicates <- NaN
CANDIDATE_WINDOWS_SIZE30_table$dN_check <- NaN
CANDIDATE_WINDOWS_SIZE30_table$dS_check <- NaN
CANDIDATE_WINDOWS_SIZE30_table$dNdS_check <- NaN
CANDIDATE_WINDOWS_SIZE30_table$dN_m_dS_check <- NaN
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_SE <- NaN
CANDIDATE_WINDOWS_SIZE30_table$boot_dS_SE <- NaN
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_over_dS_SE <- NaN
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_over_dS_P <- NaN
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_m_dS_SE <- NaN
CANDIDATE_WINDOWS_SIZE30_table$P_value <- NaN
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_gt_dS_count <- NaN
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_eq_dS_count <- NaN
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_lt_dS_count <- NaN
CANDIDATE_WINDOWS_SIZE30_table$ASL_dN_gt_dS_P <- NaN
CANDIDATE_WINDOWS_SIZE30_table$ASL_dN_lt_dS_P <- NaN

# Bootstrap parameters
NBOOTSTRAPS <- 10000
NCPUS <- 8
# no min number codons because NA here

#View(CANDIDATE_WINDOWS_SIZE30_table)

# LOOP
for (i in 1:nrow(CANDIDATE_WINDOWS_SIZE30_table)) {
  #i <- 15 # 21 # 15
  cat(i, " ")
  this_row <- CANDIDATE_WINDOWS_SIZE30_table[i, ]
  this_product <- as.character(this_row$product)
  this_start <- as.integer(this_row$codon_start)
  this_end <- as.integer(this_row$codon_end)
  
  this_data_subset <- filter(codon_results_byProductCodon_WINDOWS_POOLED_WIDE, # codon_results_byProductCodon, 
                             product == this_product, # product_segment == this_product, # 
                             codon_num >= this_start,
                             codon_num <= this_end)
  
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$num_codons <- nrow(this_data_subset)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_diffs <- sum(this_data_subset$N_diffs, na.rm = TRUE)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$S_diffs <- sum(this_data_subset$S_diffs, na.rm = TRUE)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_sites <- sum(this_data_subset$N_sites, na.rm = TRUE)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$S_sites <- sum(this_data_subset$S_sites, na.rm = TRUE)
  
  # Find codons with N diffs
  #View(this_data_subset)
  this_data_subset_NDiffCodons <- filter(this_data_subset, N_diffs > 0)$codon_num
  this_data_subset_NDiffCodons_string <- paste(this_data_subset_NDiffCodons, collapse = ", ")
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_diff_codons <- this_data_subset_NDiffCodons_string
  
  # Find codon with MAX piN
  this_data_subset$piN <- this_data_subset$N_diffs / this_data_subset$N_sites
  #this_data_subset[is.nan(this_data_subset$piN), ]$piN <- NA
  this_data_subset_maxPiNCodon <- as.integer(this_data_subset[! is.na(this_data_subset$piN) & this_data_subset$piN == max(this_data_subset$piN, na.rm = TRUE), ]$codon_num)
  #View(this_data_subset)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$piN_max_codon <- this_data_subset_maxPiNCodon
  
  #View(this_data_subset)
  
  # BOOTSTRAP THIS FULL WINDOW
  boot_dNdS <- dNdS_diff_boot_fun(this_data_subset, 'N', 'S', NBOOTSTRAPS, NCPUS)
  
  # RECORD HEADER
  boot_vector_names <- c(#'num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                         'num_replicates', 
                         'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                         'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
  
  #boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
  boot_dNdS_vector <- unlist(str_split(string = paste(boot_dNdS, sep = "\t"), pattern = "\t"))
  
  # Add names
  names(boot_dNdS_vector) <- boot_vector_names
  
  # Populate additional rows
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$num_bootstraps = NBOOTSTRAPS
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$num_codons = as.integer(boot_dNdS_vector['num_codons'])
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_sites_check = as.numeric(boot_dNdS_vector['N_sites'])
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$S_sites_check = as.numeric(boot_dNdS_vector['S_sites'])
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_diffs_check = as.numeric(boot_dNdS_vector['N_diffs'])
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$S_diffs_check = as.numeric(boot_dNdS_vector['S_diffs'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$num_replicates = as.integer(boot_dNdS_vector['num_replicates'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$dN_check = as.numeric(boot_dNdS_vector['dN'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$dS_check = as.numeric(boot_dNdS_vector['dS'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$dNdS_check = as.numeric(boot_dNdS_vector['dNdS'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$dN_m_dS_check = as.numeric(boot_dNdS_vector['dN_m_dS'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$P_value = as.numeric(boot_dNdS_vector['P_value'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P'])
  
}

# Calculate basic estimates
CANDIDATE_WINDOWS_SIZE30_table$dN <- CANDIDATE_WINDOWS_SIZE30_table$N_diffs / CANDIDATE_WINDOWS_SIZE30_table$N_sites
CANDIDATE_WINDOWS_SIZE30_table$dS <- CANDIDATE_WINDOWS_SIZE30_table$S_diffs / CANDIDATE_WINDOWS_SIZE30_table$S_sites
CANDIDATE_WINDOWS_SIZE30_table$dNdS <- CANDIDATE_WINDOWS_SIZE30_table$dN / CANDIDATE_WINDOWS_SIZE30_table$dS

# SAVE
#write_tsv(CANDIDATE_WINDOWS_SIZE30_table, "../results/CANDIDATE_WINDOWS_SIZE30_table_FULL.tsv")
write_xlsx(CANDIDATE_WINDOWS_SIZE30_table, "../results/CANDIDATE_WINDOWS_SIZE30_table_FULL2.xlsx")





