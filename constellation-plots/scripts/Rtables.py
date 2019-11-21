
# coding: utf-8

# # R group tables 
# ## This script contains relevant functions for generating and visualizing R group tables

# In[1]:


from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults
from IPython.display import SVG
import re
import pandas as pd
from cairosvg import svg2png

# In[2]:


def getFrags(analog, core):
    '''
    For a given pair analog-core, it gets the fragments with their attachment site
    '''
    wid = analog.index[0]
    analog = Chem.MolFromSmiles(analog.values[0])
    core = Chem.MolFromSmiles(core)
    
    frags = Chem.MolToSmiles(Chem.ReplaceCore(analog, core, analog.GetSubstructMatch(core),labelByIndex=True), isomericSmiles=True).split(".")
    
    subs = {"WID": wid}
    
    for j in frags:
        try:
            subs[int(re.search("\[(.*?)\*\]",j).group(1))] = re.sub("\[(.*?)\*\]","[*]",j)
        except:
            subs[0] = re.sub("\*","[*]",j)
    
    return subs #, set(subs.keys()) - {"WID"}

def Rtab(analogs, core):
    '''
    input: labeled list of analogs and core
    output: R group table
    '''
    frags = list()
    for k in range(len(analogs)):
        kfrags = getFrags(analogs[[k]], core)
        frags.append(kfrags)
    
    res = pd.DataFrame(frags).set_index("WID").fillna("[H]")
    res[res == ""] = "[H]"
    
    keep = (res.apply(lambda x: (x == "[H]").mean()) != 1)
    keep = keep[keep].index
    res = res[keep]
    
    return res

def mol_with_atom_index( mol ):
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx() ) )
    return mol

def Rcore(core, analogs, filename= None):
    Rtable = Rtab(analogs,core)
    mol = mol_with_atom_index(Chem.MolFromSmiles(core))
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(400,200)
    drawer.DrawMolecule(mol, highlightAtoms=list(Rtable.columns),highlightBonds=[])
    drawer.FinishDrawing()
    fig = SVG(drawer.GetDrawingText().replace('svg:',''))
    
    if filename:        
        #figfile = filename + ".svg"
        #with open(figfile, "w") as f:
        #    f.write(fig.data)
        svg2png(bytestring=fig.data, write_to=filename +".png")
        
        Rfile = filename + ".csv"
        Rtable.to_csv(Rfile)
            
    return fig, Rtable

def moltosvg(mol,molSize=(450,150),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:','')
DrawingOptions.bondLineWidth=1.8

