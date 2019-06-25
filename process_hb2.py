#Return the protein hbonds formed with each nucleotide in the RNA in given .pdb, OR,
#Return tne RNA hbonds formed with each residue in the protein in given .pdb.

import numpy as np
import sys
import os
import getopt

def usage():
    print('''python process_hb2.py -i input.hb2 [...] -p file.pdb -R -P -S -d bond_length -o output.txt
             Options:
             -h, --help

             REQUIRED
             -i, --input     The input .hb2 file from running hbplus on .pdb file.

             -p, --pdb       Original .pdb file. Default: file has same root as .hb2 file + .pdb extension.


             OPTIONAL
             -P, --protein   [Default] Output hydrogen bond information for each nucleotide.

             -R, --RNA       Output hydrogen bond information for each nucleotide in RNA in structure.

             -S, --summary   Output summary statistics for protein-RNA hydrogen bonds in structure.

             -d, --dadist    Specify a threshold donor-acceptor distance above which to ignore reported hydrogen bonds.
                            [hbplus has its own default for dadist, but you may want a seconary cutoff.]

             -o, --output    Specify an output file name. Default is stdout.

          ''')
