"""
Author: Alexandra Vatsiou
created: 24/2/2015
last update: 6/12/2015


Match Snps to genes if they are withing the gene transcript or k bp upstream or downstream

Three arguments are needed:
     -sys.argv[1] = a SNP file with the SNP positions and the corresponding SNP scores.
       Any format is supported as far as the SNP positions are in the second column and scores in the third column.
       e.g. snp.txt
     -sys.argv[2] = a gene file with chromosome, geneIDs, start, end position of the gene and symbol.
       e.g. gene.txt
     -sys.argv[3] = the distance (k) in base pairs that accounts for SNPs that are k bp upstream or downstream of the gene.  

The output is a snp2genes.txt file that has the following format:
'SNP Position\t Score\t GeneID\t Symbol\t Start Position\t End Position\n'

Run as follows:
SNPs2Genes('snps.txt','gene.txt',100)
"""


def SNPs2Genes(snpfile, genefile, distance):


    distance = int(distance)
    ### snps
    snps = open(snpfile,'r')
    snps = snps.readlines()
    ch = []
    pos = []
    score=[]
    for i in range(1,len(snps)):
        snps[i] = snps[i].split('\t')
        ch.append(int(snps[i][0]))
        pos.append(int(snps[i][1]))
        score.append(float(snps[i][2]))

    ### genes
    chromosome = []
    geneIds = []
    start = []
    end = []
    symbol = []

    genes = open(genefile,'r')
    genes = genes.readlines()
    for i in range(1,len(genes)):
        genes[i]=genes[i].split('\t')
        chromosome.append(int(genes[i][0]))
        geneIds.append(int(genes[i][1]))
        start.append(int(genes[i][2]))
        end.append(int(genes[i][3]))
        symbol.append(genes[i][4][0:-1])

    #### match the SNPs to the genes according to their position
    y = open('snp2genes.txt','w')
    y.write('SNP Position\tScore\tGeneID\tSymbol\tStart Position\tEnd Position\n')
    for p in range(0,len(pos)):
        for chd in range(0,len(chromosome)):
                if chromosome[chd]==ch[p]:
                    if pos[p]>=(start[chd]-distance) and pos[p]<=(end[chd]+distance):
                            y.write(str(pos[p])+'\t'+str(score[p])+'\t'+str(geneIds[chd])+'\t'+str(symbol[chd])+'\t'+str(start[chd])+'\t'+str(end[chd])+'\n')
    y.close()
                
                
    
