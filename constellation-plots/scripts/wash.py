
# coding: utf-8

# # Wash function

# Input: SMILES
# Output: Canonical SMILES
# 
# Script for standardisation of molecules. Includes removing salts and neutralizing the charges in the molecules. The largest fragments are kept.

# In[4]:


from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import time 
import multiprocessing as mp
from rdkit.Chem.Descriptors import MolWt as mw


# In[ ]:


def rmw(smi):
    '''
    Gets Rounded Average Molecular Weight for a string
    '''
    return round(mw(Chem.MolFromSmiles(smi.replace("*", "[H]"))))

def getha(string):
    '''
    Counts number of heavy atoms
    '''
    
    return Chem.MolFromSmiles(string).GetNumHeavyAtoms()


# The order of execution I propose is as follows:
# 1. Get the largest fragment in the molecule (that was originally thought to remove salts and seems like it works).
# 
# I do not know how good is to use such strict conditions for the selection of the largest fragment, because I think that sometimes, when the size difference between fragments is not "big enough", it is difficult to decide which one to keep (I am thinking about which one is the responsible for biological activity).
# 
# Another reason to start here is because if you start by filering elements you could reject molecules with metals even when such metals are part of the salts, for example.
# 
# 2. Filter elements.
# 
# In this part I found you use sets instead of lists to improve the speed, I supose. However, order is important when looking for "oops" in the string because by doing the search the way you write it, two-symbol elements are not considered at all. That is why I opted to delete the allowed symbols, and in such case you have to start by two-symbol elements to avoid confussions (and that is why the order matters).
# 
# 3. Neutralise Charges
# 
# I asked you by email if the neutralisation step is necessary because charge may play a crucial role in biological activity. If the algorithm you are going to implement is sensible to molecules' charges, of course it is necessary. 
# 
# The neutralisation function seems to work well but there are some things that should be noted. For example, the following groups:
# 
# Non-neutralisable:
# 
# 1.- Boron compounds (e.g. Tetraphenylborates)
# 
# 2.- Quaternary ammonium cations
# 
# 3.- Sulfonium ions
# 
# Neutralizable but not supported by the function:
# 
# 1.- Negatively charged amines
# 
# 2.- Hydroxylamines
# 
# Neutralizable, not supported by the function and questionable:
# 
# 1.- Carbanions
# 
# Personally I would not perform the neutralisation step because all of the processes involved consist of adding or removing hydrogens, which consierably modify the behaviour and the "chemistry" of compounds. But again, if it is needed, I can add the missing cases.

# In[2]:


def getlargestfrag(smi):
    
    if not isinstance(smi, str):
        return "Error_getlargestfrag. Input should be a string, not " + str(type(smi))
    if smi is None:
        return "Error_getlargestfrag. Molecule is None"
    
    largest = "[*]"
    
    for i in smi.split('.'):
        try:
            if Chem.MolFromSmiles(i).GetNumHeavyAtoms() > Chem.MolFromSmiles(largest).GetNumHeavyAtoms():
                largest = i
        except:
            return "Error_getlargestfrag. Molecule is None"
            
    if largest == "[*]":
        return None
    else:
        return Chem.MolToSmiles(Chem.MolFromSmiles(largest))


# In[3]:


def ElementFilter(smi, AllowFrags=False, AddTwoSymElems=None, AddOneSymElems=None, AddNotAllowed=None, Allowed=None):
    '''
    Input:  smi -> string, AllowFrags -> whether it should allow fragments or not 
            optional: AddAllowed -> a list of additionally allowed elements, 
                      Allowed -> a list of allowed elements to override default list
    Output: a detailed error when string contains unsupported characters, 
            and the original string otherwise
    '''
    if not isinstance(smi, str):
        return "Error_ElementFilter. Input should be a string, not " + str(type(smi))
    
    if Allowed is None:
        OneSymElems = ["H","B","C","N","O","F","P","S","I","c", "n", "o", "p", "s"]
        if AddOneSymElems:
            OneSymElems = OneSymElems + AddOneSymElems
        
        TwoSymElems = ["Si","Se","Br","Cl","se"]
        if AddTwoSymElems:
            TwoSymElems = TwoSymElems + AddTwoSymElems
        
        AllElems = TwoSymElems + OneSymElems
        chars = ["(",")","[","]","=","/","\\","+","@","-","#","%"]
        numbers = [str(x) for x in list(range(10))]
        Allowed = AllElems + chars + numbers
        if AllowFrags is True:
            Allowed = Allowed + ["*", "R"]
            
    NotAllowed = ["[Sc","[Co","[Cn","[In","[Po","[Rn","[Ho","[Np","[No"]
    if AddNotAllowed:
        NotAllowed = NotAllowed + AddNotAllowed
        
    oops = [i for i in NotAllowed if i in smi]
    
    if oops:
        return "Error_ElementFilter. Molecule contains unsupported characters: " + ", ".join(oops) + ". Original mol: " + smi
    else:
        oops = smi
        for i in Allowed:
            oops = oops.replace(i,"")    
    
        if oops:
            return "Error_ElementFilter. Molecule contains unsupported characters: " + ", " + oops + ". Original mol: " + smi
        else:
            return smi


# In[4]:


""" starts contribution from Hans de Winter """

def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None

def NeutraliseCharges(mol, reactions=None):
    
    
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    newmol = mol
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while newmol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(newmol, reactant, product)
            newmol = rms[0]
    if replaced:
        return newmol
    else:
        return mol

""" ends contribution from Hans de Winter """


# In[19]:


def Wash(smi, NonStereo = True,AllowFrags=False, AddTwoSymElems=None, AddOneSymElems=None, AddNotAllowed=None, Allowed=None):
    
    smi = getlargestfrag(smi)
    if smi.lower().find("error") >= 0:
        return smi, None
    
    smi =ElementFilter(smi, AllowFrags, AddTwoSymElems, AddOneSymElems, AddNotAllowed, Allowed)
    if smi.lower().find("error") >= 0:
        return smi, None
    
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return "Error. Cannot parse molecule: " + smi, None
    
    if smi.find("+") + smi.find("-") > -2:
        mol = NeutraliseCharges(mol)
    
    if NonStereo is True:
        AllChem.RemoveStereochemistry(mol)
    
    postsmi = Chem.MolToSmiles(mol)
    return postsmi, mol


# In[ ]:


def err2none(string):
    '''
    Converts Wash function error to NoneType
    '''
    if string.lower().find("error") > -1:
        return None
    else:
        return string
    
def fwash(string):
    '''
    To get only washed smiles from Wash function
    '''
    return err2none(Wash(string)[0])


# In[ ]:


def uniqwash(df, ncpu):
    '''
    get unique smiles with descriptors
    '''
    
    washed = pd.DataFrame({"washed": df.washed.drop_duplicates().dropna()}) #get unique washed
    washed["wid"] = ["W" + str(x) for x in range(len(washed))] #add washed ID: WID
    with mp.Pool(processes = ncpu) as pool: # Calculate number of heavy atoms and molecular weight
        t0 = time.time()
        washed["nhatoms"] = pool.map(getha,washed.washed)
        washed["mw"] = pool.map(rmw,washed.washed)
        t1 = time.time()
    print("Uniquify washed: {} mins".format(round((t1-t0)/60,2)))
    
    return washed


# In[ ]:


def washdf(infile, smilescol, ncpu, sep=","):
    '''
    final function to read a dataframe and add a washed column
    '''

    df = pd.read_csv(infile, sep=sep, engine='python') #read compounds
    df = df.rename(columns = {smilescol: "Molecule"}) #change smiles name
    #print(df.head())

    with mp.Pool(processes = ncpu) as pool:  #apply wash function in parallel
        t0 = time.time()
        df["washed"] = pool.map(fwash,df.Molecule)
        t1 = time.time()
        print("Wash done in: {} mins".format(round((t1-t0)/60,2)))
        
    return df, uniqwash(df,ncpu)

