"""
SnpGsea has the following functions:
    - extractPathways: extract the pathways from the txt file in .dat files which is easier to manipulate
    - sumstat: calculate the pathway score
    - qvalue: to calculate the q-values for each pathway from the p-values
    - pruning: trimming to correct for significance due to overlapping among the pathways
    - randomization: create X random pathways to calculate empirically the significance of the trimmed pathways
    - thresholds: calculate the mean proportion of pathways in bins of p-values
    - FDR: create the map of p-values(p) and FDR(p).
"""



from pylab import *
import numpy
from scipy.stats import *
from pickle import *
from random import randrange
from random import *
from qvalue import *
from ast import literal_eval

def extractPathways(pathfile):
    x = open(pathfile,'r')
    x=x.readlines()
    length = len(x)
    PathIds = []
    all_paths = []
    p = []
    for i in range(0,length):
        x[i]=x[i].split('\t')
        if x[i][0][0] == '[':
            p.append((literal_eval(x[i][0]),int(x[i][1])))
        else:
            p.append((int(x[i][0]),int(x[i][1])))
    p = sorted(p)
    
    
    for i in range(0,len(p)):
        if i == 0:
            PathGenes = []
            PathIds.append(int(p[i][0]))
            PathGenes.append(int(p[i][1]))
        else:
            if (p[i][0]) == (p[i-1][0]):
                PathGenes.append(int(p[i][1]))
                if (i == len(p)-1):
                    all_paths.append(PathGenes)
            else:
                all_paths.append(PathGenes)
                PathGenes = []
                PathIds.append((p[i][0]))
                PathGenes.append((p[i][1]))

    pickle.dump((PathIds),open("PathIds.dat","w"))
    pickle.dump((all_paths),open("PathGenes.dat","w"))


def sumstat(pathfile,pathId,paths,genesfile,minGenes,distr,noPerm):
    """
      Sumstat function calculates the sumstat score and the p-values for each pathway
      It takes 7 arguments:
          - pathfile is the output file name
          - pathId is the dat file with all the pathway ids extracted from the extractPathways function 
          - paths is the dat file with all the genes of each pathway extracted from the extractPathways function
          - genesfile is the txt file with all the gene ids and the gene scores
          - minGenes is the minimum number of genes that the user wants in the pathway
          - distr defines the way that p-values are calculated. It has two options 'norm' and 'perm'. Run the testScoreDistribution module to choose the appropriate option.
          - noPerm is the number of permutations to calculate the p-values
    """
    minGenes = int(minGenes)
    noPerm = int(noPerm)
    # create the file to write down the pathIds and the score according to SUMSTAT
    pathScores=open(pathfile,'w')
    pathScores.write('Length\tScore\tID\tPvalue\n')

    # read the human path Ids and the corresponding genes that belong to each of them
    setId=load(open(pathId))
    gene_sets= load(open(paths))

    # read the gene scores 
    genes=[]
    scores=[]
    geneScores=open(genesfile,"r")
    geneScores=geneScores.readlines()
    for i in range(1,len(geneScores)):
        geneScores[i]=geneScores[i].split('\t')
        genes.append(int(geneScores[i][0]))
        scores.append(float(geneScores[i][1]))
    mean=numpy.mean(scores)
    std=numpy.std(scores)
    varianceGenes = std*std


    # calculate the sumstat for each pathway that belong to the final list
    all_sumstat=[] # save all the sumstat scores
    length=[]      # save all the lengths of the pathways that have a score
    all_ids=[]         # save the pathIds that have a sumstat score
    all_geneSets = []

    for sets in range(0,len(gene_sets)):
        stat = []         # save all the gene scores that exist in the pathway
        count_genes = 0   # count if the genes that exist in the pathway and have a score are more than minGenes
        for gene in range(0,len(gene_sets[sets])):
            if gene_sets[sets][gene] in genes:
                index=genes.index(gene_sets[sets][gene])
                stat.append(scores[index])
                count_genes+=1

        # if the pathway has more that minGenes gene scores assigned to it, we calculate the sumstat. minGenes is defined by the user
        if count_genes>minGenes:
            sumstat=sum(stat)
            length.append(count_genes)
            all_sumstat.append(sumstat)
            all_ids.append(setId[sets])
            all_geneSets.append(gene_sets[sets])

    pickle.dump((all_ids),open("currentPathIds.dat","w"))
    pickle.dump((all_geneSets),open("currentPathGenes.dat","w"))
  


##########################################
##  Calculate the p-values according to the distribution (-norm or perm) that the user defined



##  To optimize the procedure, we make the dictionary according to
##  the length of the sets, where
##  the length is the key and 
##  the value is a list of 2 lists where the first has all the scores of the pathways with the corresponding length and the second has the pathIDs

    dic={}
    # save the initial values 
    p=[]
    id=[]
    p.append(all_sumstat[0])
    id.append(all_ids[0])
    init=[] 
    init.append(p)
    init.append(id)

    dic[length[0]]=init
    for i in range(1,len(length)):
        if length[i] in dic.keys():
            dic[length[i]][0].append(all_sumstat[i]) 
            dic[length[i]][1].append(all_ids[i])
        else:
            p=[]
            id=[]
            p.append((all_sumstat[i]))
            id.append(all_ids[i])
            new_entry=[]
            new_entry.append(p)
            new_entry.append(id)
            dic[length[i]]=new_entry

    ### create the random distributions to calculate the significance according to the length of the gene sets
    all_pvalues = []
    all_keys = []
    for key,values in dic.items():

        if distr == 'norm':
            variancePerm = varianceGenes*int(key)
            stdPerm = math.sqrt(variancePerm)
            all_randomSumstat = numpy.random.normal(loc=mean*int(key), scale=stdPerm, size=noPerm)

        # make the randomization with noPerm permutations.      
        elif distr == 'perm':
            all_randomSumstat=[]
            for k in range(0,noPerm):
                all_randomSumstat.append(sum(sample(scores, int(key))))

        # calculate the p-value directly without saving the random values 
        for sc in range(0,len(values[0])):
            prop=0
            for rand in all_randomSumstat:
                if rand>=values[0][sc]:
                    prop+=1
            p_value=float(prop)/float(len(all_randomSumstat))
            all_pvalues.append(p_value)
            all_keys.append(str(values[1][sc]))
            pathScores.write(str(key)+'\t'+str(values[0][sc])+'\t'+str(values[1][sc])+'\t'+str(p_value)+"\n")
    pathScores.close()
    
   
    dump((all_pvalues),open('pvalues.dat',"w"))
    dump((all_keys),open('keys.dat',"w"))
    



def qvalue(pathfile):

    """
    Once we have the p-values for each pathway, we calculate the q-values according to Storey et al. 
    This function take as an input the pathfile, a txt file (e.g. 'pathScores.txt') that contains the length of the pathway, pathID, sumstat score and the p-value

      The p-values are extracted to calculate the q-values.
    - the output is the pathScoresQvalues, a txt file that contains all the previous info mentioned above plus the q-values for each path

    """
     #save all the data to lists (length of pathway, pathId, pathscore, pvalue)
    length = []
    pathscores = []
    path_ids = []
    pv = []
    y=open(pathfile,'r')
    y=y.readlines()
    for i in range(1,len(y)):
        y[i]=y[i].split('\t')        
        length.append(y[i][0])
        pathscores.append(y[i][1])
        path_ids.append(y[i][2])
        pv.append(float(y[i][3]))

    # estimate the q_values and save again the updated data
    all_pvalues = numpy.array(pv)
    qv = list(estimate(all_pvalues))

    qvFile = open('pathsBeforePruning.txt','w')
    qvFile.write('Length\tScore\tID\tPvalue\tQvalues\n')
    for i in range(0,len(length)):
        qvFile.write(str(length[i])+'\t'+str(pathscores[i])+'\t'+str(path_ids[i])+'\t'+str(pv[i])+'\t'+str(qv[i])+"\n")
    qvFile.close()



def pruning(pvalues,keys,genesfile,pathId,paths,minGenes):

    """

    Pruning function conducts the pruning to correct for significance due to overlapping pathways
 
    It takes 6 arguments:
        - pvalues has all the pathway p-values 
        - keys is a list of 2 lists where the first has all the scores of the pathways with the corresponding length and the second has the pathIDs
        - genesfile is the txt file with all the gene ids and the gene scores
        - pathId is the dat file with all the pathway ids extracted from the extractPathways function 
        - paths is the dat file with all the genes of each pathway extracted from the extractPathways function
        - minGenes is the minimum number of genes that the user wants in the pathway
    """

    
    
    #save all the data to lists (length of pathway, pathId, pathscore, pvalue)
    all_pvalues1=pickle.load(open(str(pvalues))) 
    all_pvalues = all_pvalues1[:]
    all_keys=pickle.load(open(str(keys)))       

    #read the human path Ids and the corresponding genes that belong to each of them
    setId=pickle.load(open(pathId))      

    # convert all to strings
    for i in range(0,len(setId)):
        setId[i] = str(setId[i])
    gene_sets= pickle.load(open(paths))  
  
    # read the gene scores 
    genes=[]
    scores=[]
    geneScores=open(genesfile,"r")
    geneScores=geneScores.readlines()
    for i in range(1,len(geneScores)):
        geneScores[i]=geneScores[i].split('\t')
        genes.append(int(geneScores[i][0]))
        scores.append(float(geneScores[i][1]))

       
    l=zip(all_pvalues,all_keys)
    l=sorted(l)

    #make a copy of the pathways because we need to remove elements.
    geneslistsCopy=[]
    for copy in gene_sets:
        copy=copy[:]
        geneslistsCopy.append(copy)
  
    # remove the first pathway with the lowest p-value    
    removedPath=l[0]
    l.remove(l[0])
 
    # find the id of this pathway in order to extract all the genes that are in this pathway
    index=setId.index(removedPath[1])
    genes2remove=list(gene_sets[index])

    # remove the genes from the gene list
    genesNext=genes[:]
    scoresNext=scores[:]
    for g1 in range(0,len(genes)):
        if genes[g1] in genes2remove:
            genesNext.remove(genes[g1])
            scoresNext.remove(scores[g1])
 
    genes=genesNext[:]
    scores=scoresNext[:]

    # create the new gene file, which will be used from sumstat to calculate the new p-values for the trimmed pathways
    genefileNext = open("geneScores.txt",'w')
    genefileNext.write("geneID\tfinal Normalized Score\n")
    for g1 in range(0,len(genes)):
        genefileNext.write(str(genes[g1])+"\t"+str(scores[g1])+'\n')
    genefileNext.close()       
    
    # remove from the other gene sets the genes of the removed pathway
    # here we remove the genes, we don't remove pathways, even if there is no gene, the list will be empty
    for glist in range(0,len(gene_sets)):
        for g in range(0,len(gene_sets[glist])):
            if gene_sets[glist][g] in genes2remove:
                geneslistsCopy[glist].remove(gene_sets[glist][g])
   
    # remove the genesets that are less than length minGenes 

    geneslistsNext=[]
    setIdNext=[]
    
    for glist in geneslistsCopy:
        if len(glist)>minGenes:
            index2keep = geneslistsCopy.index(glist)
            geneslistsNext.append(glist)
            setIdNext.append(setId[index2keep])
    dump((geneslistsNext),open('currentPathGenes.dat',"w"))
    dump((setIdNext),open('currentPathIds.dat',"w"))
    return removedPath,genes2remove
    

def randomization(pathfile,genesfile,distr,noPerm,randGeneSets):
    """
    Create X random pathways and calculate their p-values.
    It has 5 inputs:
        - pathfile is the name of file that has all the path info
        - genesfile is the txt file with all the gene ids and the gene scores
        - distr defines the way that p-values are calculated. It has two options 'norm' and 'perm'. Run the testScoreDistribution module to choose the appropriate option.
        - noPerm is the number of permutations to calculate the p-values
        - randGeneSets is the number of random sets of pathways generated. These are being used to calculate the significance of the trimmed pathways empirically

    """

    #save all the path data (length of pathway, pathId, pathscore, pvalue)
    length = []
    pathscores = []
    all_keys = []
    all_keys1 = []
    all_pvalues = []
    all_qvalues = []
    y=open(pathfile,'r')
    y=y.readlines()
    for i in range(1,len(y)):
        y[i]=y[i].split('\t')        
        length.append(int(y[i][0]))
        pathscores.append(y[i][1])
        all_keys1.append(str(y[i][2]))
        if y[i][2][0] == '[':
            all_keys.append(literal_eval(y[i][2]))
        else:
            all_keys.append(int(y[i][2]))
        all_pvalues.append(float(y[i][3]))
        all_qvalues.append(float(y[i][4]))


    # read the gene scores 
    genes=[]
    scores=[]
    geneScores=open(genesfile,"r")
    geneScores=geneScores.readlines()
    for i in range(1,len(geneScores)):
        geneScores[i]=geneScores[i].split('\t')
        genes.append(int(geneScores[i][0]))
        scores.append(float(geneScores[i][1]))
    mean=numpy.mean(scores)
    std=numpy.std(scores)
    varianceGenes = std*std

    # sort the path data according to the length
    tuples = zip(length,pathscores,all_keys)
    tuples = sorted(tuples)
    length,pathscores,all_keys = zip(*tuples)
 
    # we create random gene sets, each of which have len(pathscores) pathways
    all_permutations=[]
    all_permPathGenes = []
    for perm in range(0,randGeneSets):
        all_randomScores = []
        all_randomGenes = []
        for path in range(0,len(length)):
            randomScores = []
            randomGenes = []
            for i in range(length[path]):
                random_gene = randrange(1,len(scores),1)
                randomScores.append(scores[random_gene])
                randomGenes.append(genes[random_gene])
            all_randomScores.append(sum(randomScores))
            all_randomGenes.append(randomGenes)
        all_permutations.append(all_randomScores)
        all_permPathGenes.append(all_randomGenes)
    
    # calculate and save the p-values for each of the different gene sets
    # each of the list in the lists has a set of pathway p-values.
    # The p-value corresponds to the pathway after they have been sorted according to their length

    lists = [[] for i in range(randGeneSets)]
    lengthSet=list(set(length))

    for size in lengthSet:
        if distr == 'norm':
            variancePerm = varianceGenes*int(size)
            stdPerm = math.sqrt(variancePerm)
            all_randomSumstat = numpy.random.normal(loc=mean*int(size), scale=stdPerm, size=noPerm)
        # make the randomization with nomPer permutations. noPerm is defined by the user     
        elif distr == 'perm':
            all_randomSumstat=[]
            for k in range(0,noPerm):
                all_randomSumstat.append(sum(sample(scores, int(size))))

        for allSize in range(0,len(length)):
            if length[allSize] == size:
                for per in range(0,len(all_permutations)):
                    prop=0
                    for rand in all_randomSumstat:
                        if rand>=all_permutations[per][allSize]:
                            prop+=1
                    p_value=float(prop)/float(len(all_randomSumstat))
                    lists[per].append(p_value)

    dump(list(all_keys),open('PathIds.dat',"w"))                    
    dump((lists),open('randomPvalues.dat',"w"))
    dump((all_permPathGenes),open('randomPathGenes.dat',"w"))
    dump((list(all_keys1)),open('sorted_keys.dat',"w"))


def thresholds(all_pvalues,bins):
    """
    Calculate the mean proportion of pathways in bins of p-values considering all the random set of pathways generated by the function randomization
    It takes 2 inputs:
        - all_pvalues are all the p-values from all the random pathways
        - bins is the step between [0,1.02] which simply defines the range of p-values for which later on corresponding FDR value are calculated.
          A greater value for the bins parameter results in a more detailed map (p-value(p) to FDR(p)) 
        

"""
    thresholds=list(numpy.arange(0,1.02,1/float(bins)))
    mean=[]
    all_mean_prop=[]
    for thres in range(0,len(thresholds)):
        all_props=[]
        for i in range(0,len(all_pvalues)):
            prop=0
            for current_pvalue in range(len(all_pvalues[i])):
                if all_pvalues[i][current_pvalue]<=float(thresholds[thres]) and all_pvalues[i][current_pvalue]>=float(thresholds[thres-1]):
                    prop+=1
            all_props.append(prop)
        mean.append(numpy.mean(all_props))
        all_mean_prop.append((thresholds[thres],numpy.mean(all_props)))
    return all_mean_prop



def FDR(x,y,bins):

    """
    Empirical correction using a histogram based method which compares the distributions of the real and the random p-values
    This approach is based on (Mosig et al., 2001; Nettleton D, 2006; Bancroft, 2009) explained in detail in the SI

    It takes 3 inputs:
      x is the pvalues of the real data
      y is the mean proportion of the random pathways that were calculated by the threshold function
      bins is is the step between [0,1.02] which simply defines the range of p-values for which corresponding FDR value are calculated.

"""
                               
    # real values
    thresholds=list(numpy.arange(0,1.02,1/float(bins)))

    #extract the p-values from the real data (x)
    real_pvalues=[]
    all_reaLen=[]
    for thres in range(0,len(thresholds)):
        current=[]
        for i in range(0,len(x)):
           if (x[i][0][0])<float(thresholds[thres]) and (x[i][0][0])>=float(thresholds[thres-1]):
               current.append(x[i][0][0])
        all_reaLen.append(len(current))
        real_pvalues.append((thresholds[thres], len(current)))


    #expected distribution (y)
    all_randomLen=[]
    for i in range(0,len(y)):
        all_randomLen.append(y[i][1])

    #plot the real and expected distribution in a histogram 
    all_reaLen = all_reaLen[1:]
    all_randomLen = all_randomLen[1:]
    thresholds=thresholds[1:]
    bar(thresholds,all_reaLen,width=0.02,label="Real",color='green')
    plot(thresholds,all_randomLen,'ro',label="Expected",color='red',linewidth=4)
    xlim(0,1)
    xlabel("p-values")
    ylabel('frequency')
    legend()
    savefig('Real_vs_Expected_Distribution.pdf')


    # calculate the m0
    m=len(x)
    for i in range(0,len(all_reaLen)):
        if all_reaLen[i]<all_randomLen[i]:
             m0=sum(all_reaLen[i+1:])/(1-thresholds[i])
             m0=m0/float(m)
             break



    #calculate R(p*)
    all_r=[]
    for thres in thresholds:
        current_r=[]
        for i in range(0,len(x)):
           
            if x[i][0][0]<=thres:
                current_r.append(x[i])
        all_r.append(len(current_r))
        

    # calculate the FDR
    all_fdr=[]
    thresholds_new = []
    for thres in range(0,len(thresholds)):
        if (float(all_r[thres]))!=0.0:
            thresholds_new.append(thresholds[thres])
            all_fdr.append((m0*m*thresholds[thres])/float(all_r[thres]))

    
    if all_fdr!=[]:
    # create the map p-value /q-value
        qvalues=[]
        mapFile = open('thresholds.txt','w')
        mapFile.write('Pvalue\tFDR\n')
        # calculate the q-value
        for thres in range(0,len(thresholds_new)):
            mapFile.write(str(thresholds_new[thres])+'\t'+str(min(all_fdr[thres:]))+'\n')
            qvalues.append(min(all_fdr[thres:]))    
        mapFile.close()







    
        

