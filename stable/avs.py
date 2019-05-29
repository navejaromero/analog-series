
# coding: utf-8

# # Virtual screening of analogs
# After running the getcores script, virtual screening of analogs for a given compound can be carried out following the next steps:
# 
# 1. Enumeration of putative cores in the molecule
# 2. Search in cores dictionary
# 3. Map cores to molecules

# In[1]:


from scripts.fragment import fragment
from scripts.coreproc import cor2mol
from scripts.Rtables import Rcore
import pandas as pd
import argparse
import os


# In[26]:


parser = argparse.ArgumentParser(description='To perfore core virtual screening for a list of queries in a proprocessed database. Run get-cores.py first')
parser.add_argument('-q','--queries', help='Input queries', required=True)
parser.add_argument('-qid','--queryid', help='Queries ID column name', required=True)
parser.add_argument('-d','--database', help='Input database', required=True)
parser.add_argument('-did','--databaseid', help='ID column name in database', required=True)
parser.add_argument('-p','--prefix', help='Prefix for output file', default="vs_")
parser.add_argument('-c','--coreprop', help='Minimum scaffold/molecule proportion', default=2/3)
parser.add_argument('-s','--sep', help='Separator in input file', default=",")
parser.add_argument('-smi','--smilescol', help='Name of column with SMILES', default="Molecule")
parser.add_argument('-r','--rtables', help='If true, R tables are written down', default=True)
args = vars(parser.parse_args())


# args = dict()
# args["queries"] = "queries.txt"
# args["queryid"] = "mol_ID"
# args["database"] = "out_cores.tsv"
# args["databaseid"] = "mol_ID"
# args["prefix"] = "vs_"
# args["coreprop"] = 2/3
# args["sep"] = "\t"
# args["smilescol"] = "Smile"
# args["qid"] = "mol_ID"
# args["rtables"] = True

# In[3]:


insmi = pd.read_csv(args["queries"], sep=args["sep"], engine="python")
fcores = pd.read_csv(args["database"], sep="\t", engine="python")


# In[4]:


output = pd.DataFrame()


# In[5]:


for row in insmi.iterrows():
    cores = {cor2mol(i) for i in fragment(row[1][args["smilescol"]],c=args["coreprop"])}
    hitcores = fcores[fcores.core.isin(cores)]
    hitcores["queryid"] = row[1][args["queryid"]]
    hitcores["query"] = row[1][args["smilescol"]]
    output = pd.concat([output,hitcores])


# In[6]:


output.to_csv(args["prefix"] + args["queries"].split(".")[0] + "_" + args["database"].split(".")[0]+".tsv", index=None, sep="\t")


# In[7]:


print("{} total analogs found for {} queries".format(len(output), len(output.queryid.unique())))


# In[8]:


noqueries = output[~output.mol_ID.isin(output.queryid.unique())]


# In[9]:


noqueries.to_csv(args["prefix"] + "rmqueries_" + args["queries"].split(".")[0] + "_" + args["database"].split(".")[0]+".tsv", index=None, sep="\t")


# In[10]:


print("{} analogs (excluding queries) found for {} queries".format(len(noqueries), len(noqueries.queryid.unique())))


# In[15]:


if args["rtables"] == True:
    os.mkdir("./Rtables/")
    os.mkdir("./Rtables/0index/")
    j=1
    for core in noqueries.core.unique():
        dirname = "./Rtables/" + str(j)
        os.mkdir(dirname)
        analogs = noqueries[noqueries.core == core][["washed", args["databaseid"]]].drop_duplicates()
        analogs.index = analogs[args["databaseid"]]
        analogs = analogs["washed"]
        queries = noqueries[noqueries.core == core][["query", "queryid"]].drop_duplicates()
        queries.index = ["query:" + str(i) for i in queries["queryid"]]
        queries = queries["query"]
        analogs = queries.append(analogs)
        fig, rtab = Rcore(core, analogs, filename = dirname + "/Rtab")
        svg2png(bytestring=fig.data, write_to="./Rtables/0index/" + str(j))
        j+=1

