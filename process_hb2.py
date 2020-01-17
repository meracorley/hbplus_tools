#Return the protein hbonds formed with each nucleotide in the RNA in given .pdb, OR,
#Return tne RNA hbonds formed with each residue in the protein in given .pdb.

import numpy as np
import sys
import os
import getopt
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
import operator
from scipy import stats as scipystats
import copy

def getChainsFromPDB2(pdbfile): #returns dict containing info about each chain in PDB file
    #more correct handling of hetatms in chain
    RNA = {'A':1,'C':1,'U':1,'G':1}
    protein = {"ALA":1,"ARG":1,"ASN":1,"ASP":1,"CYS":1,"GLU":1,"GLN":1,"GLY":1,"HIS":1,"ILE":1,
               "LEU":1,"LYS":1,"MET":1,"PHE":1,"PRO":1,"SER":1,"THR":1,"TRP":1,"TYR":1,"VAL":1}
    file = open(pdbfile,'r') #'5m73.hb2'
    lines = file.read().splitlines()
    file.close()
    chains = {} #code for later if I want to return sequence of whole chains...
    for i in range(0,len(lines)):
        thisline = lines[i]
        if thisline[0:4]!="ATOM" or thisline[0:4]!="HETA":
            continue
        chain = thisline[21:22]
        if chain in chains:
            continue
        resnum = int(thisline[22:26].strip(' '))
        restype = thisline[17:20].strip(' ')
        chains
        TER = False
        while not TER:
            i+=1
            thisline = lines[i]
            if thisline[0:3]=="TER":
                TER = True
        ##not finished...
    
def getChainsFromPDB(pdbfile): #returns dict containing info about each chain in PDB file
    #format: {chainLetter:['RNA|protein',length]}
    #example: {'A':['RNA',27]}
    RNA = {'A':1,'C':1,'U':1,'G':1}
    protein = {"ALA":1,"ARG":1,"ASN":1,"ASP":1,"CYS":1,"GLU":1,"GLN":1,"GLY":1,"HIS":1,"ILE":1,
               "LEU":1,"LYS":1,"MET":1,"PHE":1,"PRO":1,"SER":1,"THR":1,"TRP":1,"TYR":1,"VAL":1} 
    file = open(pdbfile,'r') #'5m73.hb2'
    lines = file.read().splitlines()
    file.close()
    chains = {} #code for later if I want to return sequence of whole chains...
    for line in lines:
        #figure out chains and how many bases/res in each
        thisline = line #.split() #each column occupies an exact number of char, not always white space between them
        #see: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        #print(thisline)
        if thisline[0:4]!="ATOM":
            continue
        chain = thisline[21:22]
        resnum = int(thisline[22:26].strip(' '))
        restype = thisline[17:20].strip(' ')
        if chain in chains:
            if resnum > chains[chain][3]: #have reached a new residue. each entry in pdb file is for an atom, many atom per residue
                chains[chain][1]+=(resnum-chains[chain][3])
            chains[chain][3] = resnum
            if restype in RNA:
                chains[chain][0] = "RNA"
            elif restype in protein:
                chains[chain][0] = "protein"
            continue
        if restype in protein:
            chains[chain] = ["protein",1,resnum,resnum]
        elif restype in RNA:
            chains[chain] = ["RNA",1,resnum,resnum]
        else:
            chains[chain] = ["other",1,resnum,resnum]
    #print(chains)
    return chains #format: {chainLetter: [type,num_res_with_data,res_number_of_first_res,num_last_res]}
                                #(sometimes there isnt data on some residues)

def getSeqFromPDB(pdbfile):
    #making a dictionary of the residues in the PDF file such that any
    #residue number in the .hb2 file (in format A0248 for ex.) can be used to fetch its base/amino acid
    #then I can iterate through all bases in RNA(s) and print out any interactions
    file = open(pdbfile,'r') #'5m73.hb2'
    lines = file.read().splitlines()
    file.close()
    residues = {}

    for line in lines:
        thisline = line #line.split()
        if thisline[0:4]!="ATOM" and thisline[0:4]!="HETA": #modified bases, like 5BU are labeled as HETATM in .pdbs
            continue
        chain = thisline[21:22]
        restype = thisline[17:20].strip(' ')
        resnum = thisline[22:26].strip(' ')
        while len(resnum)<4: #hb2 resnums are of standard length 4...they dont expect more than 9999 residues i guess
            resnum = "0"+resnum
        residues[chain+resnum] = restype
    return residues

def getDonor(line,chainDict): #returns [residuenum, amino/HOH/base, atomtype] of donor as dict
    resnum = line[0:5]
    resnumInt = int(resnum[1:])
    name = line[6:9].lstrip(' ')
    atomtype = line[9:13].strip(' ')
    if resnum[0] in chainDict and resnumInt<=chainDict[resnum[0]][3] and resnumInt>=chainDict[resnum[0]][2]:
        chaintype = chainDict[resnum[0]][0]
    else:
        chaintype="other"
    return {'residue':resnum,'molecule':name,'atom':atomtype,'type':chaintype}
def getAcceptor(line,chainDict): #returns [residuenum, amino/HOH/base, atomtype] of acceptor
    thisline = line[14:]
    resnum = thisline[0:5]
    resnumInt = int(resnum[1:])
    name = thisline[6:9].lstrip(' ')
    atomtype = thisline[9:13].strip(' ')
    #chaintype = chainDict[resnum[0]][0]
    if resnum[0] in chainDict and resnumInt<=chainDict[resnum[0]][3] and resnumInt>=chainDict[resnum[0]][2]:
        chaintype = chainDict[resnum[0]][0]
    else:
        chaintype="other"
    return {'residue':resnum,'molecule':name,'atom':atomtype,'type':chaintype}

def getHbonds(hb2file,pdbfile,giveRNAinteractions=False,convertMolecule=False):
    #pdbfile = hb2file.split('.')[0]+".pdb"
    #Grab everything below starting line9 (is this true of every hb2 file?)
    waterdict = {} #dict to store water mol and things hbonding with them
    waterbondsdict = {} #storing hbondnums or halengths for water bonds
    waterbondsdict2 = {}
    #using to find base-protein pairs that both interact with the same water molecule
    #and thus are coordinated by it
    interactions = [] #for each interaction: [RNAres,RNAbase,RNAmoeity,aminores,amino,aminoMoeity,waterCoordinated?,H-bondnum]
    #hbondnum refers to the number assigned to this interaction in the hb2 file
    RNA = {'A':1,'C':1,'U':1,'G':1}
    protein = {"ALA":1,"ARG":1,"ASN":1,"ASP":1,"CYS":1,"GLU":1,"GLN":1,"GLY":1,"HIS":1,"ILE":1,
               "LEU":1,"LYS":1,"MET":1,"PHE":1,"PRO":1,"SER":1,"THR":1,"TRP":1,"TYR":1,"VAL":1} 
    rnaAtom = {"O3'":"backbone","O5'":"backbone","P":"backbone","OP1":"backbone","OP2":"backbone","OP3":"backbone",
               "O2'":"sugar","O4'":"sugar"} #anything not in this list is base atom
    protAtom = {"O":"mainchain","N":"mainchain"} #anything not in this list is amino sidechain atom
    file = open(hb2file,'r') #'5m73.hb2'
    lines = file.read().splitlines()[8:]
    file.close()
    for line in lines:
        #M is mainchain amino, S is sidechain amino atom. H is more het molecule (is not amino)
        #H can refer to RNA, but also water and other small molecules
        thistype = line[33:35] #looking for MH/HM or SH/HS or HH
        if "H" not in thistype:
            continue
        chains = getChainsFromPDB(pdbfile)
        acc = getAcceptor(line,chains)
        donor = getDonor(line,chains)
        #can determine residue types with RNA and prot dicts, but this won't pick up modified bases like 5BU
        #and I don't know all the modified bases there may be. So, also trying to determine residue type by
        #by chain identity, hence accChainType and donorChainType. acc['type'] will give type (RNA or protein or other)
        dadist = line[28:32].strip(' ') 
        hadist = line[53:57].strip(' ')
        hbondnum = line[70:75].strip(' ')
        waterCoord = "False"
        
        if thistype=="HM" or thistype=="HS": #RNA or HOH is donor, protein is acceptor
            if donor['molecule']=="HOH":
                waterres = donor['residue']
                if waterres in waterdict and waterdict[waterres]['type']=="RNA": #gives the residue coordinated by water
                    donor = waterdict[waterres] #if coordinated residue is RNA...
                    waterCoord = "True"
                    hbondnum += ","+waterbondsdict[waterres]
                    dadist= str(round(float(dadist)+float(waterbondsdict[waterres]),2))
                    hadist+=","+waterbondsdict2[waterres]
                else:
                    waterdict[waterres] = acc #[resnum,aminoacid,atomtype] of interacting prtein
                    waterbondsdict[waterres] = dadist #hbondnum
                    waterbondsdict2[waterres] = hadist
                    continue
            if donor['type']=="RNA":
            #else the H atom is a base (M/S is mainchain or sidechain of protein)
                interactions.append([donor['residue'],donor['molecule'],donor['atom'],
                                    acc['residue'],acc['molecule'],acc['atom'],waterCoord,dadist,hadist])

        elif thistype=="HH":
            if donor['type']=="RNA" and acc['molecule']=='HOH':
                rna = donor
                water = acc
            elif acc['type']=="RNA" and donor['molecule']=='HOH':
                rna = acc
                water = donor
                #'''
            #This includes rna base pairs in hbond output
            elif acc['type']=="RNA" and donor['type']=="RNA" and giveRNAinteractions: #RNA basepair
                interactions.append([donor['residue'],donor['molecule'],donor['atom'],
                                     acc['residue'],acc['molecule'],acc['atom'],waterCoord,dadist,hadist])
                interactions.append([acc['residue'],acc['molecule'],acc['atom'],
                                     donor['residue'],donor['molecule'],donor['atom'],waterCoord,dadist,hadist])
                continue
                #'''
            else: #type HH can be two waters pairing or other
                #but Im only interested in RNA and water under HH type
                continue
            waterres = water['residue']
            if waterres in waterdict and waterdict[waterres]['type']=="protein":
                waterCoord = "True"
                hbondnum += ","+waterbondsdict[waterres]
                dadist=str(round(float(dadist)+float(waterbondsdict[waterres]),2))
                hadist+=","+waterbondsdict2[waterres]
                prot = waterdict[waterres]
                interactions.append([rna['residue'],rna['molecule'],rna['atom'],
                                     prot['residue'],prot['molecule'],prot['atom'],waterCoord,dadist,hadist])
            else:
                waterdict[waterres] = rna
                waterbondsdict[waterres] = dadist #hbondnum
                waterbondsdict2[waterres] = hadist
                continue
                
        else: #thistype== "MH or SH; protein is the donor, RNA or water the acceptor"
            if acc['molecule']=="HOH":
                waterres = acc['residue']
                if waterres in waterdict and waterdict[waterres]['type']=="RNA":
                    acc = waterdict[waterres]
                    waterCoord = "True"
                    hbondnum += ","+waterbondsdict[waterres]
                    dadist=str(round(float(dadist)+float(waterbondsdict[waterres]),2))
                    hadist+=","+waterbondsdict2[waterres]
                else:
                    waterdict[waterres] = donor
                    waterbondsdict[waterres] = dadist #hbondnum
                    waterbondsdict2[waterres] = hadist
                    continue
            if acc['type']=="RNA":
                interactions.append([acc['residue'],acc['molecule'],acc['atom'],
                                     donor['residue'],donor['molecule'],donor['atom'],waterCoord,dadist,hadist])
    #Replace RNAatoms in column 2 and proteinAtoms in col 5 with their type (backbone, base, etc)
    if convertMolecule:
        for i in range(0,len(interactions)):
            if interactions[i][2] in rnaAtom:
                interactions[i][2] = rnaAtom[interactions[i][2]]
            else:
                interactions[i][2] = "base"
            if chains[interactions[i][3][0]][0]=="RNA": #RNA basepair
                if interactions[i][5] in rnaAtom:
                    interactions[i][5] = rnaAtom[interactions[i][5]]
                else:
                    interactions[i][5] = "base"
                continue
            if interactions[i][5] in protAtom:
                interactions[i][5] = "main"
            else:
                interactions[i][5] = "side"
            
    return interactions

def listInteractions(hb2file,pdbfile,dadist=10):
    interactions = getHbonds(hb2file,pdbfile,False,False)
    #First print out RNA sequence. Find longest RNA chain, go thru each residue and fetch nucleotide.
    chains = getChainsFromPDB(pdbfile)
    #print(chains)
    residues = getSeqFromPDB(pdbfile)
    #print(residues)
    bestChain = ""
    bestChainCount = 0
    for letter in chains: #sometimes structures will have multiple chains
        #of the same RNA. pick the one with most complete data
        chain = chains[letter]
        chainCount = 0
        if chain[0]=='RNA': #residues in chainA for ex, will range from A001 to A00N, where N is chain length
            #chain length is stored in chains[chain][1]
            if chain[1]>bestChainCount:
                bestChainCount = chain[1]
                bestChain = letter
    #print("bestChain",bestChain)
    chain = chains[bestChain]

    #Go through interaction list and collapse on nucleotide, make comma sep list of columns 3-8
    if len(interactions)==0:
        print("No interactions with RNA bases.")
        return
    nucList = {}
    uniqueInteractions = []
    inter = np.array(interactions,dtype='object') #dtype 'str' only allows up to 8 chars, use object instead
    #interactions = inter[inter[:,6]!='True'] #ignore water-coordinated hbonds
    interactions = inter
    
    for i in range(0,len(interactions)):
        thisnuc = interactions[i][0]
        bondLen = float(interactions[i][7])
        if bondLen >= dadist:
            continue
        #if interactions[i][6]=="True": #trying out ignoring water coordinated bonds
        #    continue
        if thisnuc in nucList:
            #add columns 2-9 to position in uniqueInteractions inidicated by nucList
            pos = nucList[thisnuc]
            for j in range(2,8):
                uniqueInteractions[pos][j]+= ","+interactions[i][j] 
        else:
            nucList[thisnuc] = len(uniqueInteractions)
            uniqueInteractions.append(np.copy(interactions[i,0:8]))
            
    return uniqueInteractions

def getFeatures(chain,bestChain,nucList,uniqueInteractions,residues): 
    #requires that getHbonds returned backbone/sugar/base types instead of atoms
    #returns the number and lengths of Hbonds for each nucleotide falling into 6 categories:
    #backbone-protein, sugar-protein, base-protein, backbone-RNA, sugar-RNA, base-RNA, 
    #rather, returns the lengths of each, 0 if none of a given category, comma separated if multiple
    #example: A 0 2.3 1.1,4.2 0 0 0
    features = np.zeros((chain[3]+1-chain[2],7),dtype='object')
    positionDict = {"backbone-protein":1,"sugar-protein":2,"base-protein":3,
                    "backbone-RNA":4,"sugar-RNA":5,"base-RNA":6}
    rnaOrProt = {"base":"RNA","backbone":"RNA","sugar":"RNA","side":"protein","main":"protein",}
    counter = 0
    for i in range(chain[2],chain[3]+1): #now print out all the interactions for each base/res in chain
        resnum = str(i)
        while len(resnum)<4: #hb2 resnums are of standard length 4...pad with zeros
            resnum = "0"+resnum
        resnum = bestChain+resnum
        if resnum in nucList:
            pos = nucList[resnum]
            buildStr = ""
            Hbonds = uniqueInteractions[pos] 
            thisPart = Hbonds[2].split(',')
            thatPart = Hbonds[5].split(',')
            lengths = Hbonds[7].split(',')
            for j in range(0,len(thisPart)):
                bondType = thisPart[j]+"-"+rnaOrProt[thatPart[j]]
                listPos = positionDict[bondType]
                if features[counter][listPos] == 0:
                    features[counter][listPos] = lengths[j]
                else:
                    features[counter][listPos]+= ","+lengths[j]

            #format example: E0148 C base,sugar C0237,E0099 LYS,G side,base False,False 3.1,2.83
        if resnum in residues:
            features[counter][0] = residues[resnum] #nucleotide identity
            #print(residues[resnum]+" 0 0 0 0 0 0")
        else: #some bases(resnum) are missing from the pdb structure. print out N
            features[counter][0] = "N"
            #print("N 0 0 0 0 0 0")
        counter+=1

    return features

def interactionsPerBase(hb2file,pdbfile,dadist=10):
    interactions = getHbonds(hb2file,pdbfile,True,True)
    #First print out RNA sequence. Find longest RNA chain, go thru each residue and fetch nucleotide.
    chains = getChainsFromPDB(pdbfile)
    #print(chains)
    residues = getSeqFromPDB(pdbfile)
    #print(residues)
    bestChain = ""
    bestChainCount = 0
    for letter in chains: #sometimes structures will have multiple chains
        #of the same RNA. pick the one with most complete data
        chain = chains[letter]
        chainCount = 0
        if chain[0]=='RNA': #residues in chainA for ex, will range from A001 to A00N, where N is chain length
            #chain length is stored in chains[chain][1]
            if chain[1]>bestChainCount:
                bestChainCount = chain[1]
                bestChain = letter
    #print("bestChain",bestChain)
    chain = chains[bestChain]
    fullSeq = ""
    for i in range(chain[2],chain[3]+1): #get all nucleotides
        resnum = str(i)
        while len(resnum)<4: #hb2 resnums are of standard length 4...pad with zeros
            resnum = "0"+resnum
        resnum = bestChain+resnum
        if resnum in residues:
            fullSeq+= residues[resnum]
        else:
            fullSeq+= "N" #sometimes residues are missing in the pdb file. Ex res 70-100 in 5AOX.pdb 
    #print(fullSeq)
    #Go through interaction list and collapse on nucleotide, make comma sep lits of columns 3-8
    if len(interactions)==0:
        print("No interactions with RNA bases.")
        return
    nucList = {}
    uniqueInteractions = []
    inter = np.array(interactions,dtype='object') #dtype 'str' only allows up to 8 chars, use object instead
    #interactions = inter[inter[:,6]!='True'] #ignore water-coordinated hbonds
    interactions = inter
    
    for i in range(0,len(interactions)):
        thisnuc = interactions[i][0]
        bondLen = float(interactions[i][7])
        if bondLen >= dadist:
            continue
        #if interactions[i][6]=="True": #trying out ignoring water coordinated bonds
        #    continue
        if thisnuc in nucList:
            #add columns 2-9 to position in uniqueInteractions inidicated by nucList
            pos = nucList[thisnuc]
            for j in range(2,8):
                uniqueInteractions[pos][j]+= ","+interactions[i][j] 
        else:
            nucList[thisnuc] = len(uniqueInteractions)
            uniqueInteractions.append(np.copy(interactions[i,0:8]))
    
    #Now for each in RNA, return Hbond features.
    features = getFeatures(chain,bestChain,nucList,uniqueInteractions,residues)
    return features #returns list of each nucleotide in format:
    #Base bonds_between_protein_backbone bonds_between_protein_sugar bonds_between_protein_base bonds_between_rna_backbone
    #bonds_between_RNA_sugar bonds_between_RNA_base

def interactionsPerProtein(hb2file,pdbfile,metric,dadist=10,notrnabound=False):
    #metric should be one of: "pr", "aa", or "b"
    #dadist and notrnabound tell this whether to ignore hbonds with dadist>dadist and
    #notrnabound=True tells it to ignore rna-protein hbonds if the given RNA moiety interacts 
    #with RNA (with a bond length < dadist)
    #interactions format:
    #[['A0003', 'A', 'N1', 'C0178', 'TYR', 'OH', 'False', '2.90', '2.04'],...,]
    #amino acid dict to count aa frequencies among each domain's interacting residues
    AA = {"ALA":0,"ARG":0,"ASN":0,"ASP":0,"CYS":0,"GLU":0,"GLN":0,"GLY":0,"HIS":0,"ILE":0,
               "LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0}
    AAlist = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
               "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
    baseList = ['A','C','U','G']
    base = {'A':0,'C':0,'U':0,'G':0} #count base freqs among residue-base interactions as well?
    metricLengths = {"aa":20, "pr":8, "b":4}
    
    rnaAtom = {"O3'":"backbone","O5'":"backbone","P":"backbone","OP1":"backbone","OP2":"backbone","OP3":"backbone",
               "O2'":"sugar","O4'":"sugar"} #anything not in this list is base atom
    protAtom = {"O":"mainchain","N":"mainchain"} #anything not in this list is amino sidechain atom
    
    #pick the best (most complete) RNA chain. Only look at interactions in this one.
    chains = getChainsFromPDB(pdbfile) #dict of each chain and its type (protein | RNA)
    #residues = getSeqFromPDB(pdbfile) #dict of each residue and its base or amino acid
    '''
    bestChain = ""
    bestChainCount = 0
    for letter in chains: #sometimes structures will have multiple chains
        #of the same RNA or protein. pick the one with most complete data
        chain = chains[letter]
        chainCount = 0
        if chain[0]=='protein': #residues in chainA for ex, will range from A001 to A00N, where N is chain length
            #chain length is stored in chains[chain][1]
            if chain[1]>bestChainCount:
                bestChainCount = chain[1]
                bestChain = letter
    '''
    #print("bestchain",bestChain)
    if notrnabound:
        #print("Considering RNA-RNA bonds")
        interactions2 = []
        rnaRna = {} #store rna moieties that pair with rna
        interactions = getHbonds(hb2file,pdbfile,True)
        if len(interactions)==0:
            arr = np.zeros(metricLengths[metric])
            arr.fill(np.nan)
            return arr
        for hbond in interactions:
            rnakey = hbond[0]+"_"+hbond[2] #rnakey[0] gives the chain letter, 'A' for example
            key2 = hbond[3]+"_"+hbond[5]
            #if rnakey[0]!=bestChain: #only looking at interwction with the best RNA chain
            #    continue
            if chains[rnakey[0]][0]=="RNA" and chains[key2[0]][0]=="RNA" and float(hbond[7])<dadist:
                rnaRna[rnakey] = 1
                rnaRna[key2] = 1
        for hbond in interactions:
            if hbond[2] in rnaAtom: #ignoring hbonds that are not base-protein
                continue
            rnakey = hbond[0]+"_"+hbond[2]
            protkey = hbond[3]+"_"+hbond[5]
            if rnakey not in rnaRna and chains[protkey[0]][0]!="RNA":
                interactions2.append(hbond)
        #print(interactions2)
        interactions = interactions2

        if len(interactions)==0:
            return np.zeros(metricLengths[metric],dtype='float')
        #then go through interactions and only save RNA-protein interactions for which
        #that rna moiety is not paired with RNA
        
    else:
        #print("not considering RNA-RNA bonds")
        interactions = getHbonds(hb2file,pdbfile,False,False)
        #print(interactions)
        if len(interactions)==0:
            arr = np.zeros(metricLengths[metric])
            arr.fill(np.nan)
            return arr
    
    #Go through interactions and add residues+moeity as key to dictionary (both aminos and bases)
    protDict = {} #{A0002_N1:"main|side"}
    rnaDict = {} #{C0178_OP1:"backbone|sugar|base"}
    
    for hbond in interactions: #['B0162', 'U', 'OP1', 'D0135', 'ASN', 'ND2', 'False', '2.95', '1.97']
        rnakey = hbond[0]+"_"+hbond[2] #rnakey[0] gives the chain letter, 'A' for example
        protkey = hbond[3]+"_"+hbond[5]
        
        #if protkey[0]!=bestChain or float(hbond[7])>dadist: #only looking at interaction with the best RNA chain
        #    continue #some of the structures are with 2strands of RNA...then picking only one RNA chain doesnt make sense
        #trying a protein "bestChain" instead of RNA
        if float(hbond[7])>dadist:
            continue
            
        if rnakey not in rnaDict:
            if hbond[1] in base: #stats on freq of nucleotides used
                base[hbond[1]]+=1
            if hbond[2] in rnaAtom: #Get type (base|back|sugar) of this uniq RNA moiety
                rnaDict[rnakey] = rnaAtom[hbond[2]]
            else:
                rnaDict[rnakey] = "base"
                
        if protkey not in protDict:
            if hbond[4] in AA: #stats on freq of aa used
                AA[hbond[4]]+=1
            if hbond[5] in protAtom: #Get type (side | mainchain) of this uniq Protein moiety
                protDict[protkey] = protAtom[hbond[5]]
            else:
                protDict[protkey] = "sidechain"

    if metric=="aa":
        s = sum(AA.values())
        AAvals = np.zeros(len(AAlist),dtype='float')
        x = 0
        for aa in AAlist:
            AAvals[x] = AA[aa]*1./s*100
            x+=1
        return AAvals
    
    if metric=="b":
        s = sum(base.values())
        baseVals = np.zeros(len(baseList),dtype='float')
        if s==0:
            baseVals.fill(np.nan)
        else:
            x=0
            for b in baseList:
                baseVals[x] = base[b]*1./s*100
                x+=1
        return baseVals
                
    protChainsDict = {}
    for key in protDict:
        if key[0] in protChainsDict:
            protChainsDict[key[0]]+=1
        else:
            protChainsDict[key[0]] = 1
    #Get the prot chain with max number of RNA hbonds
    maxProtChain = max(protChainsDict.items(), key=operator.itemgetter(1))[0]
    aveInteractingAtomsPerProtein = len(protDict)*1./(len(protChainsDict))
    interactingAtomsPerProtein = protChainsDict[maxProtChain] #I might prefer this metric to averaging interactions over chains
    #sometimes a structure have multiple protein units and one unit has only one Hbond and brings the average down
    
    if len(protChainsDict)==0:
        return np.zeros(metricLengths[metric], dtype='float') #returning zero *detectable* hbonds
    
    protResDict = {}
    interactingProtRes = 0
    for key in protDict:
        if key.split('_')[0] not in protResDict:
            protResDict[key.split('_')[0]] = 1
            if key[0]==maxProtChain: #the number of interacting residues in the maxProtChain, rather than average over all chains
                interactingProtRes+=1
    aveInteractingProtRes = len(protResDict)*1./len(protChainsDict)
    
    rnaResDict = {} #Count rna hbonds and residues that occur with maxProtChain
    rnaDict2 = {}
    for hbond in interactions: #['B0162', 'U', 'OP1', 'D0135', 'ASN', 'ND2', 'False', '2.95', '1.97']
        rnakey = hbond[0]+"_"+hbond[2] #rnakey[0] gives the chain letter, 'A' for example
        protkey = hbond[3]+"_"+hbond[5]
        if protkey[0]!=maxProtChain or float(hbond[7])>dadist: #only looking at interactions with the picked protein chain
            continue
        if rnakey not in rnaDict2:
            if hbond[2] in rnaAtom: #Get type (base|back|sugar) of this uniq RNA moiety
                rnaDict2[rnakey] = rnaAtom[hbond[2]]
            else:
                rnaDict2[rnakey] = "base"
        if rnakey.split('_')[0] not in rnaResDict:
            rnaResDict[rnakey.split('_')[0]] = 1
    #for key in rnaDict:
    #    if key.split('_')[0] not in rnaResDict:
    #        rnaResDict[key.split('_')[0]] = 1
    interactingRNAres = len(rnaResDict)
    rnaDict = rnaDict2
    
    #now count the % interacting RNA atoms that are backbone/sugar/base
    #and % interacting protein atoms that are sidechain/mainchain
    backboneSum = 0.
    sugarSum = 0.
    sidechainSum = 0.
    mainchainSum=0.
    for atom in rnaDict:
        if rnaDict[atom]=="backbone":
            backboneSum+=1
        elif rnaDict[atom]=="sugar":
            sugarSum+=1
           
    for atom in protDict:
        if protDict[atom]=="sidechain" and atom[0]==maxProtChain:
            sidechainSum+=1
        elif atom[0]==maxProtChain:
            mainchainSum+=1
    
    percentBackbone = backboneSum/len(rnaDict)*100
    percentSugar = sugarSum/len(rnaDict)*100
    percentBase = max(0,100-percentBackbone-percentSugar)
    
    percentSidechain = sidechainSum/(sidechainSum+mainchainSum)*100 #len(protDict)*100
    percentMainchain = 100-percentSidechain

    return np.array([interactingAtomsPerProtein, interactingProtRes, len(rnaDict2), interactingRNAres, 
                     percentBackbone,percentSugar,percentBase,percentSidechain], dtype='float')

def domainStats(domaintype,metric,dadist=10,notrnabound=False):
    #metric is one of: "pr", "aa" or "b" to ananlyze either protein stats, amino acids stats, or base stats
    metricLength = {"pr":8, "aa":20,"b":4}
    hb2files = []
    for file in os.listdir(domaintype):
        if file.endswith(".hb2"):
            hb2files.append(file)
    #print(len(hb2files))
    allstats = np.zeros((len(hb2files),metricLength[metric]),dtype='float') #8 col for features, 20 for aa, 4 for bases
    counter = 0
    for hb2file in hb2files:
        #print("file is", hb2file)
        #interactions = getHbonds(root+thistype+"/"+hb2file)
        #if len(interactions)==0:
        #    continue
        pdbfile = domaintype+"/"+hb2file.split('.')[0]+".pdb"
        stats = interactionsPerProtein(domaintype+"/"+hb2file,pdbfile,metric,dadist,notrnabound)
        #print(stats)
        allstats[counter] = stats
        counter+=1
    #print(thistype)
    return allstats #,np.mean(allstats[:,6]),np.mean(allstats[:,7]))

def writeToFile(array,filename):
    if len(array)==0 or len(filename)==0:
        return
    outfile = open(filename,'w')
    for line in array:
        thisline = ""
        for i in line:
            thisline+=str(i)+"\t"
        thisline = thisline.rstrip('\t')
        outfile.write(thisline+"\n")
    outfile.close()
    return

def violinPlotByDomain(allStats): #for plotting -S pr summary stats as a violin plot for each domain
    sns.set(font_scale = 1.5)
    types = ["kh","dsRBD","rrm","znf","puf","dead","yth","csd"]
    names = ["KH","dsRBD","RRM","ZnF",'PUF','DEAD','YTH','CSD']
    colors = ["#c3d21eff","#44a5d8ff","#40a0a0ff","#296399ff","#b39e9eff","#d89a44ff","#a04094ff","#f4570bff"]
    metrics = ["Number of Protein Hbonds","Number of Interacting Residues","Number of RNA Hbonds","Number of Interacting bases"
               ,"% Hbonds with RNA Backbone","% Hbonds with RNA Sugar","% Hbonds with RNA Base", "% Hbonds with Protein Sidechain"]
    maxpdbs = 0
    for domaintype in types:
        pdbcount = len([f for f in os.listdir(domaintype) if f.endswith('.pdb')])
        if pdbcount>maxpdbs: #count which domain has the most pdb files to process
            maxpdbs = pdbcount
            
    if TTEST: #t-test every domain against all others, one domain per row, one metric per column
        allTtests = np.zeros((len(types)+1,len(metrics)),dtype='object') #nxn matrix of ttest results
        allTtests.fill(np.nan)
        allTtests[0] = metrics
        counter = 0
        for domain in types:
            counter+=1
            Ns = np.count_nonzero(np.isfinite(allStats[domain]),axis=0)
            var = np.nanvar(allStats[domain],axis=0,ddof=1)
            varOverN = var/Ns
            sqrdSums = var*(Ns-1)
            mean = np.nanmean(allStats[domain],axis=0)

            b = []
            for domain2 in types:
                if domain2!=domain and len(b)==0:
                    b = copy.deepcopy(allStats[domain2])
                elif domain2!=domain:
                    b = np.append(b,allStats[domain2],axis=0)
            bNs = np.count_nonzero(np.isfinite(b),axis=0)
            bvar = np.nanvar(b,axis=0,ddof=1)
            bvarOverN = bvar/bNs
            bsqrdSums = bvar*(bNs-1)
            bmean = np.nanmean(b,axis=0)

            for i in range(len(metrics)):
                    df = (varOverN[i]+bvarOverN[i])**2/(varOverN[i]**2/(Ns[i]-1)+bvarOverN[i]**2/(bNs[i]-1))
                    s2 = (sqrdSums[i]+bsqrdSums[i])/(Ns[i]+bNs[i]-2)
                    t = (mean[i]-bmean[i])/(s2*(1/Ns[i]+1/bNs[i]))**0.5
                    p = 1 - scipystats.t.cdf(abs(t),df=df)
                    allTtests[counter][i] = str(2*p)+"_"+domain #two-tailed t-test, multiply pval by 2
        writeToFile(allTtests, "pr.ttests.txt")  
            
    for metric in range(0,8):
        allBase = np.zeros((maxpdbs,len(types)),dtype='float') #at most maxpdbs pdb files, 8 domain types
        allBase.fill(np.nan)
        colN = 0
        for domain in types:
            thisstats = allStats[domain]
            for i in range(0,len(thisstats)):
                allBase[i][colN] = thisstats[i][metric]
            #allBase[0:len(thisstats[:,0]),colN] = thisstats[:,0]
            colN+=1
        #calculate means and sd of each column (domain) if TTEST = True
        '''
        if TTEST: #pairwise t-tests comparing ave of each domain to ave of each other domain
            Ns = np.count_nonzero(np.isfinite(allBase),axis=0)
            var = np.nanvar(allBase,axis=0,ddof=1)
            varOverN = var/Ns
            sqrdSums = var*(Ns-1)
            mean = np.nanmean(allBase,axis=0)
            #t-test every pairwise domain comparison
            allTtests = np.zeros((len(types),len(types)),dtype='object') #nxn matrix of ttest results
            allTtests.fill(np.nan)
            for i in range(0,len(types)):
                for j in range(i+1,len(types)):
                    df = (varOverN[i]+varOverN[j])**2/(varOverN[i]**2/(Ns[i]-1)+varOverN[j]**2/(Ns[j]-1))
                    s2 = (sqrdSums[i]+sqrdSums[j])/(Ns[i]+Ns[j]-2)
                    t = (mean[i]-mean[j])/(s2*(1/Ns[i]+1/Ns[j]))**0.5
                    p = 1 - scipystats.t.cdf(abs(t),df=df)
                    allTtests[i,j] = str(2*p)+"_"+names[i]+"_"+names[j] #two-tailed t-test, multiply pval by 2
            writeToFile(allTtests, metrics[metric].replace(" ","_")+"pairwise.ttests.txt")
        '''
        #PLOT 
        allBase = pd.DataFrame(allBase, columns=names)
        #plot stripcharts/violinplots for each domain
        #sns.set(font_scale=.25)
        sns.set_style('ticks')
        #put alldatasets together, each domain in one column
        plot = sns.violinplot(data=allBase, inner="box", palette=colors) #inner='point'
        plot = sns.stripplot(data=allBase, color="black", jitter=True)
        plot.set_ylabel(metrics[metric], size="small")
        plot.set_xticklabels(plot.get_xticklabels(), rotation=45)
        plot.set_xlabel("RBP domains", size="small")
        means = np.mean(allBase)
        top = np.max(np.max(allBase))
        for i in range(0,len(means)):
            plt.text(i, top+.2*top, str(round(means[i],1)), horizontalalignment='center', size='small', color='black', weight='semibold')
        #allBase
        plt.savefig(metrics[metric].replace(" ","_")+".protein_summary_by_domain.pdf")
        plt.clf()    
    
def plotSummaryStats(metric,domainValues): #for plotting AA or base summary statistics for each domain and aggregate
    if metric=="b":
        statlist = ["A","C","U","G"]
        xlab="Base"
        xlabAngle = 0
        outname="_base_summary.pdf"
        ymax = 100
        textPlace = [0.5,90]

    else:
        statlist = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
                    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
        statlist = ["A","R","N","D","C","E","Q","G","H","I",
                    "L","K","M","F","P","S","T","W","Y","V"]
        xlab="Amino Acid"
        xlabAngle = 45
        outname = "_aa_summary.pdf"
        ymax = 80
        textPlace = [0,70]
    domainNames = {'kh':"KH",'dsRBD':"dsRBD",'rrm':"RRM",'znf':"ZnF",'puf':"PUF",'dead':"DEAD",'yth':"YTH",'csd':"CSD"}
    domainColors = {'kh': "#c3d21eff",'dsRBD':"#44a5d8ff",'rrm':"#40a0a0ff","znf":"#296399ff",
                    'puf': "#b39e9eff",'dead':"#d89a44ff",'yth':"#a04094ff","csd":"#f4570bff"}
    
    allTtests = np.zeros((len(domainNames),len(statlist)),dtype='object') #nxn matrix of ttest results
    allTtests.fill(np.nan)
    counter=-1
    for domain in domainNames:
        counter+=1
        #TTEST comparing each domains ave to every other domain for each base or aa
        if TTEST:
            Ns = np.count_nonzero(np.isfinite(domainValues[domain]),axis=0)
            var = np.nanvar(domainValues[domain],axis=0,ddof=1)
            varOverN = var/Ns
            sqrdSums = var*(Ns-1)
            mean = np.nanmean(domainValues[domain],axis=0)

            #means/var of all other domains at each base/aa
            b = []
            for domain2 in domainNames:
                if domain2!=domain and len(b)==0:
                    b = copy.deepcopy(domainValues[domain2])
                elif domain2!=domain:
                    b = np.append(b,domainValues[domain2],axis=0)
            bNs = np.count_nonzero(np.isfinite(b),axis=0)
            bvar = np.nanvar(b,axis=0,ddof=1)
            bvarOverN = bvar/bNs
            bsqrdSums = bvar*(bNs-1)
            bmean = np.nanmean(b,axis=0)
                          
            for i in range(len(statlist)):
                    df = (varOverN[i]+bvarOverN[i])**2/(varOverN[i]**2/(Ns[i]-1)+bvarOverN[i]**2/(bNs[i]-1))
                    s2 = (sqrdSums[i]+bsqrdSums[i])/(Ns[i]+bNs[i]-2)
                    t = (mean[i]-bmean[i])/(s2*(1/Ns[i]+1/bNs[i]))**0.5
                    p = 1 - scipystats.t.cdf(abs(t),df=df)
                    allTtests[counter][i] = str(2*p)+"_"+domain+"_"+statlist[i] #two-tailed t-test, multiply pval by 2
            writeToFile(allTtests, metric+".ttests.txt")

        plt.clf()
        plt.ylim(0,ymax)
        for sample in domainValues[domain]:
            plt.plot(sample,'-',color=domainColors[domain])
        plt.plot(np.nanmean(domainValues[domain],axis=0),'--',color='black')
        plt.xticks(np.arange(len(statlist)), statlist, rotation=0)
        plt.xlabel(xlab,size="large")
        plt.ylabel("% Frequency in RNA interactions",size="large")
        plt.annotate(domainNames[domain],textPlace,weight="bold",size="large")
        plt.savefig(domain+outname)
        
    #Now plotting each domain's average in one plot:
    plt.clf()
    b = np.zeros((len(domainNames),len(statlist)),dtype='float') #overall average
    counter=-1
    for domain in domainNames:
        counter+=1
        b[counter] = np.nanmean(domainValues[domain],axis=0)
        plt.plot(np.nanmean(domainValues[domain],axis=0),'--',color=domainColors[domain])
    statlist = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
                    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
    if metric=="b":
        statlist = ["A","C","U","G"]
        plt.plot([28.9,16.5,34.3,20.2],'--',color='black') #percent of each base in sequence motifs from
        #dominquez et al, Mol Cell, 2018
        print(np.corrcoef(np.array([28.9,16.5,34.3,20.2]),np.nanmean(b,axis=0)))
    plt.xticks(np.arange(len(statlist)), statlist, rotation=xlabAngle)
    plt.xlabel(xlab,size="large")
    plt.ylabel("% Frequency in RNA interactions",size="large")
    #legend:
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color="#c3d21eff", lw=4),
                    Line2D([0], [0], color="#44a5d8ff", lw=4),
                    Line2D([0], [0], color="#40a0a0ff", lw=4),
                    Line2D([0], [0], color="#296399ff", lw=4),
                    Line2D([0], [0], color="#b39e9eff", lw=4),
                    Line2D([0], [0], color="#d89a44ff", lw=4),
                    Line2D([0], [0], color="#a04094ff", lw=4),
                    Line2D([0], [0], color="#f4570bff", lw=4)]
    if metric=="b":
        custom_lines.append(Line2D([0], [0], color="#000000", lw=4))
        plt.legend(custom_lines, ['KH','dsRBD','RRM','ZnF','PUF','DEAD','YTH','CSD','Dominguez'])
    else:
        plt.legend(custom_lines,['KH','dsRBD','RRM','ZnF','PUF','DEAD','YTH','CSD'])
    plt.savefig("all"+ outname)

def usage():
    print('''python process_hb2.py -i input.hb2 -p file.pdb -R -P pr|aa|b -S pr|aa|b -d bond_length -o output.txt
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
                             
             -t, --ttest     Perform t-tests between all domains for the summary=pr option.

          ''')

if __name__ == "__main__":
    #######Command line options
    HB2 = ""
    PDB = ""
    PROTEIN = False
    RNA = False
    SUMMARY = False
    OUTPUT = "hb2_statistics.txt"
    HBOND = False
    DADIST = 10
    TTEST = False
    
    argv = sys.argv[1:] #grabs all the arguments
    if len(argv)==0:
        usage()
        sys.exit()
    initialArgLen = len(argv)
    #print(argv)
    try:
        opts, args = getopt.getopt(argv, "hi:p:d:o:P:RHS:t", ["help","input=", 
                                                             "pdb=", "dadist=","output=","protein=","RNA",
                                                             "hbonds","summary=","ttest"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()

        elif opt in ("-i", "--input"):
            HB2 = arg
            if not os.path.isfile(HB2):
                print("File "+ HB2 + " does not exist.")
                sys.exit()

        elif opt in ("-p", "--pdb"):
            PDB = arg
            if not os.path.isfile(PDB):
                print("File "+ PDB + " does not exist.")
                sys.exit()

        elif opt in ("-d", "--dadist"):
            try:
                DADIST = float(arg)
            except ValueError:
                print(arg, "not a float, ignoring.")
                DADIST = 10
            
        elif opt in ("-o", "--output"):
            OUTPUT = arg

        elif opt in ("-P", "--protein"):
            PROTEIN = arg.lower()
            if PROTEIN!="pr" and PROTEIN!="aa" and PROTEIN!="b":
                print("-S option must be one of: pr | aa | b")
                sys.exit()
            
        elif opt in ("-H", "--hbonds"):
            HBOND = True

        elif opt in ("-R", "--RNA"):
            RNA = True   
            
        elif opt in ("-S", "--summary"):
            SUMMARY = arg.lower()
            if SUMMARY!="pr" and SUMMARY!="aa" and SUMMARY!="b":
                print("-S option must be one of: pr | aa | b")
                sys.exit()
        elif opt in ("-t", "--ttest"):
            TTEST = True
            
    if len(args)>0 and len(args)<initialArgLen:
        print("WARNING: Unused options", args)
    
    if bool(PROTEIN) + RNA + HBOND + bool(SUMMARY) > 1:
        print("Must choose one of -P | -R | -S  options.")
        sys.exit()
    if bool(PROTEIN) + RNA + HBOND + bool(SUMMARY) < 1:
        SUMMARY = True
        
    if PDB=="":
        PDB = HB2.rstrip('hb2')+"pdb"
    if not os.path.isfile(PDB) and not bool(SUMMARY):
        print("PDB file "+ PDB + " does not exist. PDB file matching .hb2 file required. Please provide a PDB file.")
    ########END OPTIONS
        
    #Analysis done based on whether P, R, or S is True
    stats = []
    if RNA:
        stats = interactionsPerBase(HB2,PDB,DADIST)
        header = np.array([["Base","Bond_len_prot_to_backbone","Bond_len_prot_to_sugar","Bond_len_prot_to_base",
                           "Bond_len_RNA_to_backbone","Bond_len_RNA_to_sugar","Bond_len_RNA_to_base"]], dtype='object')
        stats = np.append(header,stats,axis=0)
        writeToFile(stats,OUTPUT)
    
    elif HBOND:
        stats = listInteractions(HB2,PDB,DADIST)
        header = np.array([["RNA_num","Base","Base_atom","Prot_num","A.A.","Prot_atom","Water_coordinated?","Bond_length(s)"]],
                          dtype='object')
        stats = np.append(header,stats,axis=0)
        writeToFile(stats,OUTPUT)
     
    elif PROTEIN: #returns list of each protein residue interacting with RNA
        headers = {"pr":np.array([["Num_prot_hbonds","Num_interacting_residues","Num_RNA_hbonds","Num_interacting_bases"
                   ,"%hbonds_with_RNA_backbone","%hbonds_with_RNA_sugar","%hbonds_with_RNA_base",
                    "%hbonds_with_prot_sidechain"]],dtype=object),
                   "aa":np.array([["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
                    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]],dtype='object'),
                   "b":np.array([["A","C","U","G"]],dtype='object')}
        stats = interactionsPerProtein(HB2,PDB,PROTEIN,DADIST,False)
        stats = np.append(headers[PROTEIN],[stats],axis=0)
        writeToFile(stats,OUTPUT)
    
    else: #SUMMARY==True, iterate through domain type, for each structure run interactionsPerProtein
        #SUMMARY should be "pr", "aa", or "b"
        domains = ["kh","dsRBD","rrm","znf","puf","dead","yth","csd"]
        names = ["KH","dsRBD","RRM","ZnF",'PUF','DEAD','YTH','CSD']
        allStats = {} #store all metric stats for each domain type
        for domain in domains: #dictionary of domains
            thisstats = domainStats(domain,SUMMARY,DADIST,False)
            allStats[domain] = np.copy(thisstats)
        if SUMMARY=="aa" or SUMMARY=="b":
            plotSummaryStats(SUMMARY,allStats)
        else:
            violinPlotByDomain(allStats)
