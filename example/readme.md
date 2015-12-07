
## 
ipython --pylab

# match SNPs to genes
from SNPs2Genes import *

time SNPs2Genes("snps.txt","gene.txt",5000)

- CPU times: user 1h 6min 18s, sys: 11.6 s, total: 1h 6min 30s
- Wall time: 1h 6min 33s

###################################
## Run the gene enrichment analysis

from gea import *
### calculate the gene scores
c = calculateGeneScores('snp2genes.txt','max')
	-- The gene scores were calculated
	-- Now correct for the gene length
	-- Choose a thresBin value smaller than 659

- CPU times: user 26.6 s, sys: 1.79 s, total: 28.4 s
- Wall time: 30.3 s



### correct for the gene legth
p = correctGeneLength(c,200), where 200 is the thresBin value. 

- CPU times: user 5.41 s, sys: 44.1 ms, total: 5.46 s
- Wall time: 5.47 s


### extract the outlier genes
outlierGenes(p,0.001)

- CPU times: user 34.9 ms, sys: 1.5 ms, total: 36.4 ms
- Wall time: 35.8 ms

## Examine if the distribution of the pathway lengths (10,25,50 and 100) follow normality

time ipython testScoreDistribution.py "finalGeneScores.txt" 10000

- real	0m9.819s
- user	0m7.299s
- sys	0m0.653s

######## Based on the figure PathwayDistribution.pdf, the distribution of the scores does not follow  normality, so the argument -perm is choosen for the GSEA analysis.

## Run the gene set enrichment analysis

time ipython run.py "pathway" "finalGeneScores.txt" 10 "perm" 1000 2 200 0.15 0.15
	-- Calculation of p-values and q-values before pruning were calculated and saved in the pathsBeforePruning.txt file
	-- Calculation of p-values after pruning were calculated and saved in the finalPathsAfterPruning.txt file
	-- the random pathways were created
	-- the pruning for each of the random list of pathways is done
	-- The thresholds were calculated and saved in the thresholds.txt file

- real    417m17.910s
- user    412m44.412s
- sys     0m57.587s

######## As it was expected no "significant" pathways were identified after pruning.
######## We note that the parameter values might be different for each dataset. In this example we used a small number of permutations (1000) and a small number of random pathway sets to calculate the expected distribution (2), just to illustrate how the tool is working. Larger values should be used to estimate significance robustly. 

