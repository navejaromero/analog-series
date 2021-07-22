# analog-series
## José J. Naveja
## May 29th 2019
## naveja@comunidad.unam.mx

Scripts for breaking down a collection of molecules into analog series, getting R group tables and performing virtual screening. Off-memory and parallel computing approaches are implemented.

Requirements:

-Linux system
-dask
-networkx

To annotate a library with metacores and AS, run get-cores.py. Usage: 

usage: get-cores.py [-h] -i INFILE [-p PREFIX] [-c COREPROP] [-s SEP]
                    [-t MAXT] [-smi SMILESCOL] [--ncpu NCPU]

To get the CCR scaffolds for a database of compounds

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input database
  -p PREFIX, --prefix PREFIX
                        Prefix for output file
  -c COREPROP, --coreprop COREPROP
                        Minimum scaffold/molecule proportion
  -s SEP, --sep SEP     Separator in input file
  -t MAXT, --maxt MAXT  Maximum time (secs) per molecule for processing
  -smi SMILESCOL, --smilescol SMILESCOL
                        Name of column with SMILES
  --ncpu NCPU           number of CPU to use (if positive) or to keep free (if
                        negative)
