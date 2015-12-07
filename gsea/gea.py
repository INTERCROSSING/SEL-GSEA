""" Implementation of Gene Enrichment Analysis (GEA)
    Author: Alexandra Vatsiou
    Created: 24/2/2015
    Last update 6/12/2015

    ########################################################################################################################
    Description of all the three functions
    
    calculateGeneScores: calculates the gene scores
        
        inputs: 
            - Snp2genes.txt file with format: SNP Position\t Score\t GeneID\t Symbol\t Start Position\t End Position
              this file can be produced by the SNPs2Genes module
            - parameter to calculate the gene scores. Two options are provided:
                'max': takes the maximum score of all the SNP scores in the gene
                'mean': takes the mean of all the SNP scores in the gene

         outputs:
            - the geneScores_max.eps or geneScores_mean.eps that shows the distribution of the gene scores before correction
            - the bins.eps histogram that shows the frequency of the gene lengths           


    correctGeneLength: correct the gene scores for gene length bias
        inputs:
            - the output of calculateGeneScores
            - the thresBin value
        outputs:
            - bins.txt file that describes how GEA splits the genes in bins to correct for the very long genes that have many SNPS.
              The format of the file is as follows:
              #snps_per_gene\t #genes_in_bin\t mean_of_the_bin\t std_of_the_bin\t median_of_the_bin\n
            - finalGeneScores.txt that has the final gene scores file that has the geneIDs with the corresponding scores after correction
            - corrected_geneScores.eps that shows the distribution of the gene scores after correction


    outlierGenes (optional): extracts the outlier genes based on a cut-off thershold
        inputs:
            - the output of correctGeneLength 
            - the cut-off threshold to extract the outlier genes (e.g. 0.001)
        outputs:
            - significantGenes.txt that has the outlier gene IDs and the corresponding gene scores

###############################################################################################################################
    To run the gea.py open a terminal and run the following commands:

    from gea import *
    c = calculateGeneScores('snp2genes.txt','max')
    p = correctGeneLength(c,200), where 200 is the thresBin value. 
    outlierGenes(p,0.001), where 0.001 the cut-off level to extract the outlier genes

    
"""

import numpy
from pylab import *
import os
import pickle


def calculateGeneScores(snp2genes,metric):   
    # extract the geneID and the SNP scores from the snp2genes.txt file
    geneFile = open(snp2genes,"r")
    geneFile = geneFile.readlines()

    score = []
    genes = []
    for i in range(1,len(geneFile)):
        geneFile[i] = geneFile[i].split('\t')
        score.append(float(geneFile[i][1]))
        genes.append(int(geneFile[i][2]))

    # make a tuple with the geneID and scores and sort it
    tuple=zip(genes,score)
    tuple=sorted(tuple) 
    data =zip(*tuple)

    ##################################################################
    # calculate the gene scores with the parameter -max or -mean
    
    #initial conditions
    s = []             # save all the SNP scores that belong to each gene
    geneLength = 0    # save the length of each gene 
    l = []            # save all the lengths of the genes
    gene_scores = [] # save all the gene scores
    ids = []         # save the geneIDs

    for i in range(1,len(data[0])-1):
        if data[0][i] == data[0][i+1]:
            geneLength+=1
            s.append(data[1][i])
        else:
            if metric == 'max':
                gene_scores.append(max(s))
            if metric == 'mean':
                gene_scores.append(mean(s))
            l.append(geneLength)
            ids.append(data[0][i])
            geneLength = 0

            #re-initialize the parameters for the next geneID
            s = []
            geneLength+=1
            s.append(data[1][i])

    # save the histogram with the gene scores
    fig = figure()
    hist(gene_scores,bins=100)
    xlabel("Gene Scores ({c})".format(c=metric),fontsize=18)
    title("Distribution of Gene Scores",fontsize=18)
    fig.savefig('geneScores_{c}.eps'.format(c=metric))        
    print("The gene scores were calculated")
    print("Now correct for the gene length")
    
    #####################################################################################################################
    # binning of genes
    # GEA groups the genes in bins according to their length. The important parameter for this is the thresBin value.
    # to give an approximate idea of the thresbin value, GEA produces the histogram, below.
    # it is advised to take the 

    # make the histogram with bins equals to the maximum number of gene length
    fig = figure()
    hist(l,bins=max(l),color='blue')
    p,edges = histogram(l,bins=max(l))
    xlabel("Gene Length",fontsize=18)
    fig.savefig('bins.eps')
    print("It is advised to choose a thresBin value smaller than {c}".format(c=max(p[10:])))
    return l,gene_scores,ids


def correctGeneLength(c,thresBin):
    # p = Number of genes according to their length (bins)
    p,edges = histogram(c[0],bins=max(c[0]))
    p= list(p)

    # remove the bins that have no genes
    val=0
    while val in p:
        p.remove(val)

    #count how many genes have a unique gene length in the dataset
    count=0
    for i in range(0,len(p)):
        if p[i]==1:
            count=count+1
           
    # lset is the sorted list of the length of the genes 
    lset=(list(set(c[0])))

    # tuple1 is a tuple that has the length of the genes (in SNPs) and the frequency of the genes that have the specific length
    tuple1=zip(lset,p)
    tuple1=sorted(tuple1)
    data1 = zip(*tuple1)


    # if the number of genes in the bin exceed the thresBin value, GEA keeps the bin
    # otherwise GEA adds consecutive bins lengths till it reaches the thresBin threshold

    group=[] #the group of gene lengths that belong to each bin
    bins=[] # the number of genes that are in each bin
    k=0    
    while k<(len(data1[1])-int(thresBin)):
        binSnp=[]
        
        #if the number of the snps in the gene exceed the thresBin, keep the bin
        if data1[1][k]>=int(): 
            bins.append(data1[1][k]) 
            binSnp.append(data1[0][k])
            k=k+1
            
        # if the number of genes in the bin does not exceed the thresBin value ,
        # add consequtive bins till it reaches the thresBin value
        else:       
            binSnp.append(data1[0][k])
            j=k
            c=data1[1][k]
            if k<(len(data1[1])-int(thresBin)):
                while c<=int(thresBin):
                    c+=data1[1][j+1]
                    binSnp.append(data1[0][j+1])
                    j=j+1
                bins.append(c)
                k=j+1

            else:
                bins.append(sum(data1[1][k:]))           
                for y in range(k+1,len(data1[0])):
                    binSnp.append(data1[0][y])
                k=(len(data1[1])-int(thresBin))           
        group.append(binSnp)



    #################################################
    ### standardise the gene scores using the mean and standard deviation of the bin that each gene belongs to.
                      


    binGeneScores=[]         # save all the SNP scores of each bin of all genes before normalization to compare in the end 
    normalizedGeneScores=[]  # save the final gene scores after normalization
    normalizedGeneIDs=[]     # save the corresponding gene iDs

    for Nsnp in group:
        correct4GeneScore=[]  # save the scores for each bin to calculate the mean and std of each bin
        finalGeneIDs=[]       # save the geneIDs of each bin
        for snp in Nsnp:
            for length in range(0,len(c[0])):
                if c[0][length]==snp:               
                    correct4GeneScore.append(c[1][length])
                    finalGeneIDs.append(c[2][length])
        mean=float(numpy.mean(correct4GeneScore))
        std=float(numpy.std(correct4GeneScore))
        for score in range(0,len(correct4GeneScore)):
            if std!=0:
                new_geneScore=float((float(correct4GeneScore[score]-mean))/std)
                normalizedGeneScores.append(new_geneScore)
                normalizedGeneIDs.append(finalGeneIDs[score])
        binGeneScores.append(correct4GeneScore)
        





    #######################################################################################
    ### save the results of the bins and the final gene scores after the correction

    x = open('bins.txt','w')
    x.write('#SNPs per gene\t # genes in bin\t mean of the bin \t sd of the bin \n')
    for b in range(0,len(binGeneScores)):
        x.write(str(group[b])+'\t'+str(bins[b])+'\t'+str(numpy.mean(binGeneScores[b]))+"\t"+str(numpy.std(binGeneScores[b]))+"\n")
    x.close()

    y = open('finalGeneScores.txt','w')
    y.write('geneID\t final_Normalized_Score\n')
    for score in range(0,len(normalizedGeneScores)):
        y.write(str(normalizedGeneIDs[score])+'\t'+str(normalizedGeneScores[score])+'\n')
    y.close()
    tupleNorm = zip(normalizedGeneScores,normalizedGeneIDs)
    fig = figure()
    hist(normalizedGeneScores,bins=100)
    xlabel("Gene Scores",fontsize=18)
    title("Distribution of Gene Scores after correction",fontsize=18)
    fig.savefig('corrected_geneScores.eps')
    return tupleNorm




def outlierGenes(n,threshold):
    # extract the outlier genes based on the cut-off level that the user defines
    z = open('significantGenes.txt','w')
    length=len(n)
    cutoff = int(length*float(threshold))
    #t = zip(normalizedGeneScores, normalizedGeneIDs)
    t = sorted(n)
    significantGenes=t[length-cutoff:]
    z.write("Genes\tGene_Scores\n")
    for i in range(0,len(significantGenes)):
        z.write(str(significantGenes[i][1])+'\t'+str(significantGenes[i][0])+'\n')
    z.close()

       
