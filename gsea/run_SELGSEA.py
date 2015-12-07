
"""
The GSEA takes 8 arguments 
    -sys.argv[1] = output file name
    -sys.argv[2] = geneFile
    -sys.argv[3] = the minimun number of genes that we want to have in each of the pathways
    -sys.argv[4] = the norm or perm option to calculate the p-values either based on the normal distribution or on permutations
    -sys.argv[5] = number of permutations
    -sys.argv[6] = number of random gene sets that are created to calculate significance empirically
    -sys.argv[7] = number of bins that defines the range of the map pvalue(p) to FDR(p)
    -sys.argv[8] = qvalue threshold to extract the significant pathways before pruning
    -sys.argv[9] = qvalue threshold to extract the significant pathways after pruning

You can run every single step individually by opening a terminal and copy-paste the commands or
you can run directly the run.py file as follows:

python (ipython) run.py pathway finalGeneScores.txt 10 10000 50 200 0.05 0.05

Before you run the GSEA, run the testScoreDistribution.py script to examine if the distribution of pathways of different lengths follow normality or not.
If they follow normality, choose the option -norm in the argv[4]
If they don't, the option -perm should be chosen

"""

from SnpGsea import *
from significance import *
from pylab import *
from pickle import *
from shutil import *
import sys
import os


## extract pathways in 2 .dat files (pathIds and pathGenes), which is easier to be used in python
extractPathways(sys.argv[1])

## calculate the sumstat score and the p-values
sumstat('PathScores.txt','PathIds.dat', 'PathGenes.dat', sys.argv[2],int(sys.argv[3]),str(sys.argv[4]),int(sys.argv[5]))
copyfile('currentPathIds.dat', 'PathIdsBeforePruning.dat')
copyfile('currentPathGenes.dat', 'PathGenesBeforePruning.dat')


## calculate the q-values 
qvalue('PathScores.txt')
print("Calculation of p-values and q-values before pruning were calculated and saved in the pathsBeforePruning.txt file")


## pruning
## we calculate the p-values with sumstat at every single step a pathway is removed.
## the pruned pathways with the corresponding new p_values are saved in the list l1

l1=[]
copyfile( sys.argv[2], "geneScores.txt")
while len(load(open("pvalues.dat")))>1:
    result = pruning('pvalues.dat','keys.dat', "geneScores.txt",'currentPathIds.dat', 'currentPathGenes.dat',int(sys.argv[3]))
    l1.append(result)
    if len(load(open('currentPathIds.dat')))>1:
        sumstat('pathAfterPruning.txt', 'currentPathIds.dat', 'currentPathGenes.dat', "geneScores.txt", int(sys.argv[3]), str(sys.argv[4]), int(sys.argv[5]))
    else:
        break

x = open('finalPathsAfterPruning.txt','w')
x.write('Length\tID\tPvalue\n')
for i in range(0,len(l1)):
    x.write(str(len(l1[i][1]))+"\t"+str(l1[i][0][1])+"\t"+str(l1[i][0][0])+'\n')
x.close()
dump((l1),open('l1.dat',"w"))
print("Calculation of p-values after pruning were calculated and saved in the finalPathsAfterPruning.txt file")



################################
# to calculate significance this time, we follow a randomization procedure
# we produce x random gene sets and repeat the pruning procedure in order to create a map with FDR

###  randomization to produce n lists of random gene sets
randomization('pathsBeforePruning.txt',sys.argv[2],str(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]))
print('the random pathways were created')


############################################
##pruning of each of the random gene set
allRandoml1 = []
## x lists of random p-values of random gene sets before pruning
allRandomPvalues=pickle.load(open('randomPvalues.dat'))
allRandomGeneSets=pickle.load(open('randomPathGenes.dat'))
## do the pruning for each of the lists of random gene sets
for i in range(0,len(allRandomPvalues)):
    print(i)
    # extract the current p-value list
    dump((allRandomPvalues[i]),open('pvalues.dat',"w"))
    dump((allRandomGeneSets[i]),open('currentPathGenes.dat',"w"))

    l1_random=[]
    copyfile(sys.argv[2], "geneScores.txt")
    copyfile("sorted_keys.dat", 'keys.dat')
    copyfile("PathIds.dat", 'currentPathIds.dat')
    while len(load(open("pvalues.dat")))>1:
        result1 = pruning('pvalues.dat','keys.dat', "geneScores.txt",'currentPathIds.dat', 'currentPathGenes.dat',int(sys.argv[3]))
        l1_random.append(result1)
        if len(load(open('currentPathIds.dat')))>1:
            sumstat('RandomPaths.txt', 'currentPathIds.dat', 'currentPathGenes.dat', "geneScores.txt", int(sys.argv[3]), str(sys.argv[4]),int(sys.argv[5]))
        else:
            break
    allRandoml1.append(l1_random)

dump((allRandoml1),open('allRandoml1.dat',"w"))

## save only the pvalues after pruning
all_random_pvalues_afterPruning =[]
for i in range(0,len(allRandoml1)):
    current = []
    for k in range(0,len(allRandoml1[i])):
        current.append(allRandoml1[i][k][0][0])
    all_random_pvalues_afterPruning.append(current)
print('the pruning for each of the random list of pathways is done')


# calculate the mean proportion of pathways in bins before pruning
vpBeforePruning = thresholds(allRandomPvalues,int(sys.argv[7]))

## calculate the mean proportion of pathways in bins after pruning
vpAfterPruning = thresholds(all_random_pvalues_afterPruning,int(sys.argv[7]))

## create the map of p-values and FDR values
FDR(l1,vpAfterPruning,int(sys.argv[7]))
print('The thresholds were calculated and saved in the thresholds.txt file')



# remove the files that were saved for the purpose of the GSEA and are not needed to examine the final result.
os.remove('currentPathIds.dat')
os.remove('currentPathGenes.dat')
os.remove("geneScores.txt")
os.remove('keys.dat')
os.remove('l1.dat')
os.remove("PathGenes.dat")
os.remove("PathGenesBeforePruning.dat")
os.remove("PathIds.dat")
os.remove("PathIdsBeforePruning.dat")
os.remove('pvalues.dat')
os.remove('RandomPaths.txt')
os.remove('randomPvalues.dat')
os.remove('randomPathGenes.dat')
os.remove('sorted_keys.dat')
os.remove('pathScores.txt')
os.remove('pathAfterPruning.txt')
os.remove('allRandoml1.dat')

                
################################################
## significance

significanceBeforePruning('pathsBeforePruning.txt',float(sys.argv[8]))
significanceAfterPruning('finalPathsAfterPruning.txt','thresholds.txt',float(sys.argv[9]))





