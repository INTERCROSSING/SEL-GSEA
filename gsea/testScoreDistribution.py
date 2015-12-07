"""
Examine if the distribution of random pathways of different length follow normality or not.

inputs:
    - finalGeneScores.txt; file produced by the GEA
    - number of samples that are drawn for the pathway lentghs (10,25,50,and 100). 

output:
    - PathwayDistribution.pdf that shows the overlap of the distributions of the pathways drawn under normality and with permutations.
      If the two distributions overlap for all the pathway lengths, choose the option -norm for GSEA. Otherwise choose -perm. The p-value on the top of each graph indicates how different the two distributions are. If the p-value is <=0.05, the distributions could be significantly different.
Run as follows:
python (or ipython) testScoreDistribution.py 'finalGeneScores.py 10000'

"""


import sys
import numpy
import math
from pylab import *
from random import *
from scipy import stats

genes=[]
scores=[]
geneScores=open(str(sys.argv[1]),"r")
geneScores=geneScores.readlines()
for i in range(1,len(geneScores)):
    geneScores[i]=geneScores[i].split('\t')
    genes.append(int(geneScores[i][0]))
    scores.append(float(geneScores[i][1]))

mean=numpy.mean(scores)
std=numpy.std(scores)
varianceGenes = std*std
lengths = [10,25,50,100]

for l in range(0,len(lengths)):
    all_randomSumstat=[]
    for k in range(0,int(sys.argv[2])):
        all_randomSumstat.append(sum(sample(scores,lengths[l])))

    subplot(2,2,l+1)           
    hist(all_randomSumstat,bins=50,alpha=0.5,label="permutations")

    variancePerm = varianceGenes*(lengths[l])
    stdPerm = math.sqrt(variancePerm)
    r = numpy.random.normal(loc=mean*lengths[l], scale=stdPerm, size=int(sys.argv[2]))
    hist(r,bins=50,alpha=0.5,label='Normal randomization')
    ### run RankSum test to check if the two distributions are the same.
    ### if the p-value is <=0.05, we can say that the two distributions are significantly different.
    z_stat, p_val = stats.ranksums(all_randomSumstat, r)

    title("pathway size = {c}, p-value= {p}".format(c=lengths[l],p=round(p_val,5)), fontsize=10)
    
    
legend(prop={'size':6})
    
savefig('PathwayDistribution.pdf')
print("If the two distributions overlap for all the pathway lengths, choose the option -norm for GSEA. Otherwise choose -perm. The p-value on the top of each graph indicates how different the two distributions are. If the p-value is <=0.05, the distributions could be significantly different.")
