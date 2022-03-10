#!/usr/bin/env python
# coding: utf-8

# # Virtual screening of analogs
# After running the getcores script on the queries, virtual screening of analogs for a given compound can be carried out following the next steps:
# 
# 1. Enumeration of putative cores in the molecule
# 2. Search in cores dictionary
# 3. Map cores to molecules
# 4. Get R group tables
# 
# Most is done off-memory and in parallel

# In[8]:


from scripts.Rtables import Rcore
from scripts.vs import write_core_hits, find_line
from cairosvg import svg2png
import pandas as pd
import argparse
import os
import shutil
from functools import partial
import multiprocessing as mp
import numpy as np


# In[2]:


#Parser for running as .py
parser = argparse.ArgumentParser(description='To perform core virtual screening for a list of queries in a huge preprocessed database. Run big-data-cores.py first')
parser.add_argument('-q','--queries', help='Input queries file', required=True)
parser.add_argument('-qid','--queryid', help='Queries ID column name', required=True)
parser.add_argument('-qs','--qsep', help='Separator in queries file', default="\t")
parser.add_argument('-qsmi','--qsmilescol', help='Name of column with SMILES in queries file', default="washed")
parser.add_argument('-d','--database', help='Database folder (output from big-data...py)', required=True)
parser.add_argument('-o','--output', help='Prefix for output folder', default="vs_output/")
parser.add_argument('-p','--prefix', help='Prefix for output files', default="vs_")
#parser.add_argument('-qfrag','--fragqueries', help='Whether queries should be fragmented or taken as is', default=True)
parser.add_argument('-c','--coreprop', help='Minimum scaffold/molecule proportion', default=2/3)
parser.add_argument('-r','--rtables', help='If true, R tables are written down', default=True)
parser.add_argument('-n','--ncpu', help='Number of cpu to use. If negative, it is substracted from total', default=-1)
args = vars(parser.parse_args())


# #Manual setup when running as notebook
# args = dict()
# args["output"] = "vs_output/"
# args["path_VS"] = "/VS_tb/"
# args["queries"] = "input/tb/out_cores.tsv"
# args["queryid"] = "query_ID"
# args["qsep"] = "\t"
# args["qsmilescol"] = "core"
# args["database"] = "ZINC-frags/"
# args["prefix"] = "vs_"
# args["coreprop"] = 2/3
# args["rtables"] = True
# args["ncpu"] = -2
# #args["fragqueries"] = False

# In[10]:


ncpu = int(args["ncpu"])
if ncpu < 0:
    ncpu = mp.cpu_count() + ncpu


# In[12]:


path_VS = args["output"] + args["path_VS"]


# In[16]:


queries = pd.read_csv(args["queries"], sep=args["qsep"], engine="python")


# In[17]:


#queries


# In[22]:


if not os.path.exists(args["output"]):
    os.mkdir(args["output"])

if os.path.exists(path_VS):
    shutil.rmtree(path_VS)

[os.mkdir(path_VS + i) for i in ["", "hits", "Rtables", "hit_ids"]]


# # 1. Find core hits (only structures)

# In[23]:


def find_cores(file,cores,indir,outdir):
    with open(indir+file, "r") as dfile:
        for line in dfile:
            if any(core in line for core in cores):
                with open(outdir +file,"a+") as outfile:
                    outfile.write(line)
    with open(outdir +"/finished","a+") as outfile:
        outfile.write(file+"\n")


# In[24]:


queries = pd.read_csv(args["queries"], sep="\t")[["MID", "core"]].drop_duplicates()


# In[25]:


frags_dir = args["database"]+"frags/split/" 
cores = list(set(["|{}|".format(each) for each in queries.core]))

part_find_cores = partial(find_cores, cores =cores,
                               indir=args["database"]+"frags/split/", 
                               outdir=path_VS + "hits/")

files = os.listdir(frags_dir)


# In[26]:


with mp.Pool(processes = ncpu) as pool:  #apply wash function in parallel
    pool.map(part_find_cores,files)      # run ~5 hrs


# # 2. Add original IDs

# In[27]:


def get_hits(file, hits):
    
    df = pd.read_csv(file, sep="\t").set_index("washed")
    labs = df.index.intersection(hits)
    
    if len(labs) > 0:
        return df.loc[labs]


# In[28]:


#read hit files

data = pd.DataFrame()

for file in os.listdir(path_VS + "/hits"):
        
    try:
        if "txt" in file:
            data = pd.concat([data,pd.read_csv(path_VS + "/hits/" + file, sep="\t", header=None)])
    
    except:
        None


# In[29]:


hits =list()

for row in data.iterrows():
    hits.append(row[1].iloc[0].split(",")[0][2:-1])


# In[30]:


part_get_hits= partial(get_hits, hits=hits)


# In[31]:


files = [args["database"] + "uwashed/" + file for file in os.listdir(args["database"] + "uwashed/")]

with mp.Pool(processes = ncpu) as pool:  #apply function in parallel
    tmp = pool.map(part_get_hits,files)  #run ~5hrs


# In[34]:


hits= pd.concat(tmp)
#hits["washed"] = hits.index
hits = pd.DataFrame(hits.groupby("washed").zinc_id.apply(lambda x: ",".join(x)))


# In[36]:


datadf = list()
for row in data.iterrows():
    tmp = row[1].iloc[0][1:].split(",")
    hit_str = tmp[0][1:-1]
    zinc_id = hits.loc[hit_str].iloc[0]
        
    for i in tmp:
        datadf.append({"hit_str": hit_str, "zinc_id": zinc_id, "core": i[1:-1]})


# In[37]:


datadf=pd.DataFrame(datadf).merge(queries, on="core")


# In[38]:


queriesstr = pd.read_csv(args["queries"], sep="\t")


# In[39]:


final = queriesstr[[args["queryid"],"Molecule", "washed", "WID", "MID"]].merge(datadf, on="MID")


# In[40]:


final.to_csv(path_VS + "/hit_ids/final.tsv", index=None, sep="\t")


# In[41]:


final = pd.read_csv(path_VS + "/hit_ids/final.tsv", sep="\t")
#final.head()


# In[50]:


print(f"Found a total of {len(final)} hits matching {len(final.MID.unique())} cores.")


# # 3. Create R-group tables

# In[45]:


print("Creating R-group tables")

data=final
if args["rtables"] == True:
    rpath = path_VS + "/Rtables/"
    os.mkdir(rpath + "0index/")

    for mid in data.MID.unique():
        core = queries.loc[queries.MID == mid, "core"].iloc[0]
        dirname = rpath + mid
        os.mkdir(dirname)
        analogs = data[data.core == core][["hit_str", "zinc_id"]].drop_duplicates()
        analogs.index = analogs["zinc_id"]
        analogs = analogs["hit_str"]
        queriestmp = data[data.core == core][["washed", args["queryid"]]].drop_duplicates()
        queriestmp.index = ["query:" + str(i) for i in queriestmp[args["queryid"]]]
        queriestmp = queriestmp["washed"]
        analogs = queriestmp.append(analogs)
        fig, rtab = Rcore(core, analogs, filename = dirname + "/Rtab")
        svg2png(bytestring=fig.data, write_to=rpath + "0index/" + mid + ".png")
        with open(dirname+"/core.smi", "wt") as f:
            f.write(core)


# In[ ]:


print("DONE")

