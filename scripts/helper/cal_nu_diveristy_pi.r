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

