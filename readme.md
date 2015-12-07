Author: Alexandra Vatsiou (alex.vatsiou@gmail.com),
Created: 24/2/2015,
Last update: 6/12/2015

SEL-GSEA is a python tool that performs a Gene Enrichment Analysis (GEA) and a Gene Set Enrichment Analysis (GSEA) using SNP data and SNP scores.

###################
Run the GEA as follows:
    
    from gea import *
    c = calculateGeneScores('snp2genes.txt','max',10)
             - sys.argv[1] = the SNP file (format: SNP Position\t Score\t GeneID\t Symbol\t Start Position\t End Position
	     - sys.argv[2] = the way that the user chooses to calculate the gene scores
        			'max': take the maximum score of all the SNP scores in the gene
        			'mean': take the mean of all the SNP scores in the gene


    p = correctGeneLength(c,200)
              - sys.argv[1] = output from calculateGeneScores
              - sys.argv[2] = thresBin value
                              An approximate thresBin value is provided when you run calculateGeneScore
             
    outlierGenes(p,0.001), where 0.001 the cut-off level to extract the outlier genes
	      - sys.argv[1] = output from correctGeneLength
              - sys.argv[2] = cut-off threshold to extract the outlier genes


The output of the GEA are 3 files:
    - the geneScores_max.eps or geneScores_mean.eps that shows the distribution of the gene scores before correction
    - the bins.eps histogram that shows the frequency of the gene lengths
    - bins.txt file that describes how GEA splits the genes in bins to correct for the very long genes that have many SNPS.
          The format of the file is as follows:
          #snps_per_gene\t #genes_in_bin\t mean_of_the_bin\t std_of_the_bin\t median_of_the_bin\n
    - finalGeneScores.txt that has the final gene scores file that has the geneIDs with the corresponding scores after correction   
    - corrected_geneScores.eps that shows the distribution of the gene scores after correction
    - significantGenes.txt that has the outlier gene IDs and the corresponding gene scores


####################################################################################################################################################################################
Before you run the GSEA, run the testScoreDistribution.py as follows to examine whether the distribution of the pathway scores follow normality.

        python (or ipython) testScoreDistribution.py 'finalGeneScores.py 10000'
                - sys.argv[1] = the file with the gene IDs and gene scores produced with GEA.
                - sys.argv[2] = number of samples that are drawn for the pathway lentghs (10,25,50,and 100).
 
       The output of testScoreDistribution.py is:
                - PathwayDistribution.pdf that shows the overlap of the distributions of the pathways drawn under normality and with permutations.
                  If the two distributions overlap for all the pathway lengths, choose the option -norm for GSEA. Otherwise choose -perm. 
                  The p-value on the top of each graph indicates how different the two distributions are. If the p-value is <=0.05, the distributions could be significantly different.

#############
The Gene Set Enrichment Analysis takes 9 arguments: 

    -sys.argv[1] = the pathway file in the format of Biosystems (PathID\tGeneID\n). This file can also be produced by OLOPath.db, a database that integrates Biosystems, IntPath and GeneSetDB databases and performs automatically the overlapping procedure as the user wants and reduces redundancy. 
    -sys.argv[2] = the gene file which is a file with the geneIDs and the corresponding scores (geneID\tgeneScore\n). The user can either provide this file or produce it from the GEA.
    -sys.argv[3] = the minimun number of genes that one wants to have in each of the gene sets. Very small genes may spoil the statistics. A logical number is 10
    -sys.argv[4] = the “norm” or “perm” argument to calculate the p-values either based on the normal distribution or on permutations respectively.
    -sys.argv[5] = number of permutations to calculate the p-values and correct empirically. A logical number is >=100000. The larger this number is, the more robust the analysis is. However, after a certain number of permutations, the result doesn’t change.
    -sys.argv[6] = number of random pathway sets to be created for the randomization procedure that is used to created the expected distribution (e.g. 200, the whole enrichment analysis will be repeated 200 times with random gene sets to create the expected distribution)
    -sys.argv[7] = number of bins to split to create the map file p-values to q-values (e.g. 200 bins).
    -sys.argv[8] = qvalue threshold before pruning
    -sys.argv[9] = qvalue threshold after pruning


The output of the GSEA are files:

	- pathsBeforePruning.txt which is the pathway scores before pruning. It has the following format: LengthOfPathway\t Score\t PathID\t Pvalue\t Qvalues\n
        - pathsAfterPruning.txt which is the pathway scores after pruning. It has the following format: LengthOfPathway\t Path_Score\t PathID\t Pvalue\n
        - thresholds.txt which is the file that maps the p-values to q-values. 
        - singificantPathwaysBeforePruning.txt which is the significant pathways according to the q-value threshold that the user defines before pruning
        - singificantPathwaysAfterPruning.txt which is the significant pathways according to the q-value threshold that the user defines after pruning  


Run the GSEA as follows:

	python (ipython) run.py ‘pathway’ ‘finalGeneScores.txt’ 10 ‘perm’(or ‘norm’) 1000 50 1000 0.20 0.20

################
Additional tools:

	SNPs2Genes: match SNPs to genes 
                Inputs:
                	-sys.argv[1] = a SNP file with the SNP positions and the corresponding SNP scores. Any format is supported as far as the SNP positions are in the second column and scores in the third column (e.g. snp.txt).
    		        -sys.argv[2] = a gene file with chromosome, geneIDs, start, end position of the gene and symbol (e.g. gene.txt).
		        -sys.argv[3] = the distance (k) in base pairs that accounts for SNPs that are k bp upstream or downstream of the gene.  
                      
                   
                The output is a snp2genes.txt file that has the following format: 'SNP Position\t Score\t GeneID\t Symbol\t Start Position\t End Position\n'

		Run as follows:
		SNPs2Genes('snps.txt','gene.txt',5000)
		   
  

Dependencies:

To run the GEA and GSEA, you need python and pylab. They are easy and free to install. It is suggested to install Canopy Enthought for the corresponding system that automatically installs python and all the Python packages.


Please feel free to contact me to alex.vatsiou@gmail.com for any questions. 


