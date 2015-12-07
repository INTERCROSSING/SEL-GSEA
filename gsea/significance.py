"""
Calculate the significance of the pathways according to a q-value threshold that the user will give as parameter
"""



import sys

def significanceBeforePruning(paths,cutoff):
    x = open(paths,'r')
    x=x.readlines()

    pathId = []
    pvalue = []
    qvalue = []
    length = []
    for i in range(1,len(x)):
        x[i] = x[i].split('\t')
        pathId.append(x[i][2])
        pvalue.append(float(x[i][3]))
        qvalue.append(float(x[i][4]))
        length.append(int(x[i][0]))
    data = zip(qvalue,pvalue,pathId)
    data = sorted(data)
    qvalue,pvalue,pathId = zip(*data)
    y = open('singificantPathwaysBeforePruning.txt','w')
    y.write('PathID\tPathway_Length\tP_value\tQ_value\n')
    for i in range(0,len(pvalue)):
        if qvalue[i]<=cutoff:
            y.write(str(pathId[i])+'\t'+str(length[i])+'\t'+str(pvalue[i])+'\t'+str(qvalue[i])+'\n')
    y.close()




def significanceAfterPruning(paths,thresholds,cutoff):
    x = open(paths,'r')
    x=x.readlines()
    pathId = []
    pvalue = []
    length = []
    for i in range(1,len(x)):
        x[i] = x[i].split('\t')
        pathId.append(x[i][1])
        pvalue.append(float(x[i][2]))
        length.append(int(x[i][0]))



    thres = open(thresholds,'r')
    thres = thres.readlines()
    pv = []
    qv = []
    for i in range(1,len(thres)):
        thres[i] = thres[i].split('\t')
        pv.append(float(thres[i][0]))
        qv.append(float(thres[i][1]))

    y = open('singificantPathwaysAfterPruning.txt','w')
    y.write('PathID\tPathway_Length\tP_value\tQ_value\n')
    for i in range(0,len(qv)):
        if qv[i]<=cutoff:
            current_qv = qv[i]
            current_pv = pv[i]
            for k in range(0,len(pvalue)):
                if pvalue[k]<=current_pv:
                    y.write(str(pathId[k])+'\t'+str(length[k])+'\t'+str(pvalue[k])+'\t'+str(current_qv)+'\n')
       
    y.close()
        
    
