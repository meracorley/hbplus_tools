# This process_hb2.py was written to process the information in .hb2 files.
# .hb2 files are generated with HBPLUS analysis of hydrogen bonds in PDB files.
# Acquire HBPLUS here: https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/ 
# Generate .hb2 file like so:
$ echo file.pdb | clean 
$ hbplus -d 3.5 -h 2.7 file.new file.pdb

$ python process_hb2.py -h
python process_hb2.py -i input.hb2 -p file.pdb -R -P pr|aa|b -S pr|aa|b -d bond_length -o output.txt
             Options:
             -h, --help

             REQUIRED
             -i, --input     The input .hb2 file from running hbplus on PDB file. Not required if using -S option.
             
             OPTIONAL
             -p, --pdb       Original PDB file. Default: file has same root as .hb2 file + ".pdb" extension.
            
             -d, --dadist    Specify a threshold donor-acceptor distance above which to ignore reported hydrogen bonds.
                            [hbplus has its own default for dadist, but you may want a secondary cutoff.]

             -o, --output    Specify an output file name. Default is "hb2_statistics".
            
             Pick one of:
             -P, --protein   Summary statistics for the given protein-RNA structure, either RNA (pr), amino acids (aa), or bases (b).

             -R, --RNA       List each nucleotide in the RNA and the hbonds formed with backbone, sugar, and base.
             
             -H, --hbonds    Output information for each protein-RNA hbond in .hb2 file.

             -S, --summary   [Default] Summary statistics--protein-RNA (pr), amino acids (aa), or bases (b)--by domain type.
                             Based on provided PDB files categorized by domain (RRM, KH, dsRBD, ZnF, YTH, PUF, DEAD, CSD).
