
# coding: utf-8

# In[ ]:


import networkx as nx
import time
from itertools import chain
from collections import Counter
import pandas as pd


# In[ ]:


def getASBS(wcm, mcores, AS_path, ucores_path):
    '''
    from WCM table get analog series
    '''

    G = nx.Graph()

    t0 = time.time()

    for row in wcm.iterrows(): # this operation can be parallelized by splitting data and composing the graphs after
        w,c,m = row[1]
        G.add_edge(m,w) #only using metacores (smaller network, same info)

    
    AS = [G.subgraph(c) for c in nx.connected_components(G)]
    
    iclus = 0

    with open(AS_path + "ASW.tsv", "w+") as asw, open(AS_path+"ASM.tsv", "w+") as asm:
        asw.write("ASID\tWID\n")
        asm.write("ASID\tMID\tMnmols\tASnmols\n")
    

    with open(AS_path + "ASW.tsv", "a") as asw, open(AS_path+"ASM.tsv", "a") as asm:
        for x in AS:
            iclus += 1
        
            tmp = Counter(chain.from_iterable(x.edges())) #count occurrences
            m = [(x,tmp[x]) for x in tmp.keys() if x.count("M") > 0] #get cores
            w = [x for x in tmp.keys() if x.count("W") > 0] #get mols
        
            [asw.write("AS" + str(iclus)+"\t"+iw+"\n") for iw in w]
            [asm.write("AS" + str(iclus)+"\t"+im[0] + "\t" + str(im[1]) + "\t"+ str(len(w)) +"\n") for im in m]

    t1 = time.time()
    print("Analog series generation: {} mins".format(round((t1-t0)/60,2)))
    
    #Now, find the metacores with coverage equal to the size of their AS. 
    #Finally, select the largest metacore to represent the AS.
    
    asm = pd.read_csv(AS_path+"ASM.tsv", sep="\t")
    
    asm = asm[asm.ASnmols > 1] #first requirement: AS size > 1
    
    mcores = mcores.set_index("MID")

    asbs = asm[(asm.Mnmols == asm.ASnmols)] #strict requirement: complete AS coverage; metaseries are taken out
    asbs = pd.merge(asbs, mcores, left_on="MID", right_index=True)

    largest = asbs.groupby("ASID")["nhatoms"].apply(max)
    asbs = pd.merge(asbs, pd.DataFrame({"largest": largest}), left_on="ASID", right_index=True)

    asbs = asbs[asbs.nhatoms == asbs.largest]

    asbs = asbs[["ASID", "MID", "metacore", "ASnmols"]]

    asbs = asbs[[not i for i in asbs.ASID.duplicated()]] #Warning: check first? Counter(asbs.ASID).most_common()
    
    return asbs


def getAS(wcm, mcores, prefix):
    '''
    from WCM table get analog series
    '''

    G = nx.Graph()

    t0 = time.time()

    for row in wcm.iterrows(): # this operation can be parallelized by splitting data and composing the graphs after
        w,c,m = row[1]
        G.add_edge(m,w) #only using metacores (smaller network, same info)

    
    AS = [G.subgraph(c) for c in nx.connected_components(G)]
    
    iclus = 0

    with open(prefix + "ASW.tsv", "w+") as asw, open(prefix+"ASM.tsv", "w+") as asm:
        asw.write("ASID\tWID\n")
        asm.write("ASID\tMID\tMnmols\tASnmols\n")
    

    with open(prefix + "ASW.tsv", "a") as asw, open(prefix+"ASM.tsv", "a") as asm:
        for x in AS:
            iclus += 1
        
            tmp = Counter(chain.from_iterable(x.edges())) #count occurrences
            m = [(x,tmp[x]) for x in tmp.keys() if x.count("M") > 0] #get cores
            w = [x for x in tmp.keys() if x.count("W") > 0] #get mols
        
            [asw.write("AS" + str(iclus)+"\t"+iw+"\n") for iw in w]
            [asm.write("AS" + str(iclus)+"\t"+im[0] + "\t" + str(im[1]) + "\t"+ str(len(w)) +"\n") for im in m]

    t1 = time.time()
    print("Analog series generation: {} mins".format(round((t1-t0)/60,2)))
    
