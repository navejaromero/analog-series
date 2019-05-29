
# coding: utf-8

# # Core processing

# In[ ]:


import time
import pandas as pd
import dask.dataframe as dd
import multiprocessing as mp
from scripts.wash import rmw, getha, fwash


# In[ ]:


def ncuts(smi):
    return smi.count("*")


def cor2mol(smi):
    return fwash(smi.replace("*", "[H]"))


# In[ ]:


def ucore(recap_path,ncpu):
    '''
    Input: a path to fragmented files
    Output: unique cores, metacores and mapping file to washed structures
    '''
    t0 = time.time()
    cores = dd.read_csv(recap_path+"*", sep ="\t", header=None) #read all fragmentation files

    ucores = pd.DataFrame({"cores": cores.loc[:,0].drop_duplicates().compute()}) #get unique cores
    tmp=len(ucores)
    ucores["CRID"] = [i+str(j)  for i,j in zip(["C"] * tmp, range(1,tmp+1))] #add a core ID 

    with mp.Pool(processes = ncpu) as pool: # get hcut canonical version and number of cuts
        ucores["hcut"] = pool.map(cor2mol,ucores.cores)
        ucores["ncuts"] = pool.map(ncuts, ucores.cores)
    
    tmp = ucores.hcut.drop_duplicates()
    mcores = pd.DataFrame({"metacore": tmp, "MID":["M" + str(i) for i in range(1,len(tmp)+1)]})
    
    ucores = pd.merge(ucores, mcores.set_index("metacore"), right_index=True, left_on="hcut")
    
    with mp.Pool(processes = ncpu) as pool: # Calculate number of heavy atoms and molecular weight
        mcores["nhatoms"] = pool.map(getha,mcores.metacore)
        mcores["mw"] = pool.map(rmw,mcores.metacore)
    
    t1 = time.time()
    print("Uniquify cores: {} mins".format(round((t1-t0)/60,2)))
    
    ## Annotate fragmented molecules to cores (and metacores)
    cwdf = dd.read_csv(recap_path+"*", header = None, sep="\t").rename(columns={0: "core", 1: "WID"})
    cwdf = dd.merge(ucores.set_index("cores")[["CRID", "MID"]], cwdf, left_index=True, right_on="core")[["WID", "CRID", "MID"]]
    cwdf = cwdf.compute()
    
    return ucores, mcores, cwdf

