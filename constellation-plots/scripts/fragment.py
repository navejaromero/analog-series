
# coding: utf-8

# # Fragment molecule
# ## Input: Molecule, core minimum proportion 
# ## Output: Cores or fragments

# In[1]:


from rdkit.Chem import Recap
from rdkit import Chem
from itertools import compress
import signal
import time
import multiprocessing as mp
from functools import partial


# In[17]:


def fragment(smi, core = True, c = 2/3):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    else:
        tot = mol.GetNumHeavyAtoms()
        o = sorted(Chem.Recap.RecapDecompose(mol).GetAllChildren().keys())
        if core:
            return list(compress(o, [Chem.MolFromSmiles(x).GetNumHeavyAtoms() >= c * tot for x in o])) + [smi]
        else:
            return list(o) + [smi]


# In[3]:


def handler(signum, frame):
    raise 

def timeout_fragment(smi, maxtime = 60,core = True, c = 2/3):
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(maxtime)
    try:
        res = fragment(smi, core, c)
        signal.alarm(0)
    except Exception:
        res = "too long"
    
    return res
    


# In[4]:


def err2none(string):
    '''
    Converts Wash function error to NoneType
    '''
    if string.lower().find("error") > -1:
        return None
    else:
        return string
    
def write_fragment(smi,wid,path, maxtime = 60, core = True, c = 2/3):
    '''
    Parsed fragments
    '''
    
    frags = timeout_fragment(smi, maxtime, core, c)
    
    if frags:
        
        if frags == "too long":
            with open(path.split("recap/recap_")[0] + "toolong.tsv", "a+") as f:
                f.write(smi + "\t" + wid + "\n")
    
        else:
            with open(path, "a+") as f:
                [f.write(i + "\t" + wid + "\n") for i in frags]
    
    
def fragdf(tup, path, maxtime = 60, core = True, c = 2/3):
    '''
    Fragment a dataframe
    '''
    
    idf, num = tup
    filename = path + str(num)
    [ write_fragment(washed, wid, filename, maxtime, core, c) for washed, wid in zip(idf.washed, idf.wid)] 


# In[ ]:


def parfragdf(uwashed, ncpu, recap_path):
    
    n = 10000  #chunk row size
    list_df = [uwashed[i:i+n] for i in range(0,uwashed.shape[0],n)] #split washed in order to parallelize

    start = 0

    fnums = range(start, start + len(list_df)) 
    
    with mp.Pool(processes = ncpu) as pool: # Fragment using RECAP rules. Parallel, off-memory
        t0 = time.time()
        pool.map(partial(fragdf, path=recap_path), zip(list_df,fnums))
        t1 = time.time()
    print("Fragmentation time: {} mins".format(round((t1-t0)/60,2)))

