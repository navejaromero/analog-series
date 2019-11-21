
# coding: utf-8

# ### Open Source Analog Series
# Steps:
# 1. Wash and get unique washed structures
# 2. Fragment molecules (get cores)
# 3. Metacores (H-filled cores)
# 4. Analog Series (clusters)
# 
# Minimal version for running inline.

# #### Preamble

# In[51]:


from scripts.wash import washdf
from scripts.fragment import parfragdf
from scripts.coreproc import ucore
from scripts.asbs import getAS

import argparse
import shutil
import os
import multiprocessing as mp
import pandas as pd


# #### Parser

# In[3]:


parser = argparse.ArgumentParser(description='To get the cores for a database of compounds')
parser.add_argument('-i','--infile', help='Input database', required=True)
parser.add_argument('-p','--prefix', help='Prefix for output file', default="out_")
parser.add_argument('-c','--coreprop', help='Minimum scaffold/molecule proportion', default=2/3)
parser.add_argument('-s','--sep', help='Separator in input file', default=",")
parser.add_argument('-t','--maxt', help='Maximum time (secs) per molecule for processing', default=60)
parser.add_argument('-smi','--smilescol', help='Name of column with SMILES', default="Molecule")
parser.add_argument('--ncpu', help='number of CPU to use (if positive) or to keep free (if negative)', default= -1)
args = vars(parser.parse_args())


# #### Paths calculated

# In[7]:


path_processed = "tmp/" #temporary folder will be deleted at the end

recap_path =  path_processed + "fragments/recap/recap_"

ucores_path = path_processed + "fragments/unique_cores/"


# In[ ]:


args["ncpu"] = int(args["ncpu"])

if args["ncpu"] < 0:
    args["ncpu"] = mp.cpu_count() + args["ncpu"]
    


# In[6]:


if not os.path.exists(args["infile"]):
    raise ValueError("ERROR: Input file " + args["infile"] + " does not exist") 


# In[8]:


if os.path.exists(path_processed) and os.path.isdir(path_processed):
    shutil.rmtree(path_processed)
    
[os.mkdir(path_processed+i) for i in ["", "/fragments", 
                                          "/fragments/recap",
                                      "fragments/unique_cores"] if not os.path.exists(path_processed+i)]


# # 1. Wash

# In[7]:


washed, uwashed = washdf(args["infile"], args["smilescol"],args["ncpu"], args["sep"]) #washed and unique washed


# # 2. Fragment

# In[8]:


parfragdf(uwashed, args["ncpu"], recap_path)


# # 3. Unique cores

# In[11]:


ucores, mcores, wcm = ucore(recap_path, args["ncpu"])


# # 4. Write output

# In[52]:


getAS(wcm, mcores, args["prefix"])


# In[47]:


uwashed.rename(index=str, columns={"wid": "WID", "nhatoms": "mol_nHeavyAtoms", "mw": "mol_MW"}, inplace=True)
wcm.drop(["CRID"],axis=1,inplace=True)


# In[48]:


mcores.rename(index=str, columns={"metacore": "core", "nhatoms": "core_nHeavyAtoms", "mw": "core_MW"}, inplace=True)


# In[50]:


final = pd.merge(pd.merge(pd.merge(washed, uwashed, on="washed"),wcm, on="WID"), mcores, on="MID")


# In[50]:


final.to_csv(args["prefix"] + "cores.tsv", sep="\t", index=False)
shutil.rmtree(path_processed)

