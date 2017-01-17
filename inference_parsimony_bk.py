
import subprocess
import glob
import sys
import os
import csv
from pygraph.classes.graph import graph
from pygraph.algorithms.accessibility import connected_components
from itertools import combinations, chain


### [1]. This is just parsing the fasta file.  

## VARIABLES HERE: fastaFile- is the fasta file, specified in [5]
########### proteinNames - the names of each protein (ex: HPRR1950679)
def getAllProteinNames(fastaFile):
    proteinNames = []
    with open(fastaFile,"r") as lines:
        for line in lines:
            if line[0]=='>':
                line = line[1:]
                line = line.rstrip()
                words = line.split(' ')
                proteinNames.append(words[0])
    return proteinNames


# [2] parsimony.py - takes percolator xml output file and outputs parsimonous
#                protein identifications with peptide matches
# Outputs the identified proteins identified using parsimony. If a pair of 
# equally sized sets of proteins both cover all peptides in their group, the 
# set that matched first is returned.



def parse_percolator_txt(filename, q_value_cutoff= 0.00001): # Banu : This q_value_cutoff needs to be changed in each run of the program from 0.01 upto 0.1 ?. 
    with open(filename, 'r') as f:
        f= f.read().splitlines()
        peptides= {}
        for line in f[1:]:
            line = line.split('\t')
            q_value_str = line[9]
            pep_seq = line[11]
            protein_str = line[13]
            try:
                q_value = float(q_value_str)
            except ValueError:
                raise ValueError("Invalid q value {}".format(q_value_str))
            if pep_seq in peptides:
                raise RuntimeError("Duplicate peptide in percolator xml output")
            if q_value < q_value_cutoff:                               # Banu : Changed this from pass back to continue as in old script
                protein_ids = protein_str.split(',')
                peptides[pep_seq] = protein_ids
        return peptides

def parsimonous_protein_identification(peptides):
    """
    parsimonous_protein_identification - takes a dict of the form
        {<peptide_seq>: <protein_name>, [<protein_na,e> ...] } and returns the 
        proteins identified using parsimony.
    """
    detected_proteins = {}
    protein2peptides = {}

    # start with the uniquely determined proteins
    for peptide, proteins in peptides.items():
        if len(proteins) == 1:
            detected_proteins[proteins[0]] = [peptide]
            peptides.pop(peptide)
        else:
            for p in proteins:
                if not p in protein2peptides:
                    protein2peptides[p] = [peptide]
                else:
                    protein2peptides[p].append(peptide)

    # remaining peptides have multiple potential proteins, use parsimony
    g = graph()

    # identify protein clusters
    for peptide, proteins in peptides.items():
        for protein in proteins:
            if not g.has_node(protein):
                g.add_node(protein)
        for p1, p2 in combinations(proteins, 2):
            if not g.has_edge((p1, p2)):
                g.add_edge((p1, p2))
    connected = connected_components(g).items()
    clusters = {subgraph: set() for protein, subgraph in connected}
    for protein, subgraph in connected:
        clusters[subgraph] = clusters[subgraph].union(set((protein,)))

    def find_covering(proteins):
        peptides = set(chain(*(tuple(protein2peptides[p]) for p in proteins)))
        for k in range(1, len(proteins) + 1):
            for covering in combinations(proteins, k):
                covered = set(chain(*(tuple(protein2peptides[p]) for p in 
                    covering)))
                if len(covered) == len(peptides):
                    return covering
        return None

    # find the minimal protein covering of each cluster
    for cluster in clusters.values():
        covering = find_covering(cluster)
        if covering is None:
            print "Error, failed to cover " + str(subgraph)
            sys.exit(1)
        else:
            for protein in covering:
                detected_proteins[protein] = protein2peptides[protein]

    return detected_proteins


# [3] function to call the parsimony functions

def runParsimony(filename):
    peptides = parse_percolator_txt(filename)
    proteins = parsimonous_protein_identification(peptides)
    ## proteins puts out name, sequence 
    
    return proteins  # Banu : removed the part of making protein_list
    
    
 
# [4] Calculated FDR for each protein after parsimony ??
 

def getProteinWithFDR(targetFile,decoyFile):
    prot2q ={}
    targetproteins = runParsimony(targetFile)     # Banu: consider what the for loop with t does in Lukas' script below.
    decoyproteins = runParsimony(decoyFile)
    fdr = float(len(decoyproteins))/len(targetproteins)
    a = open('fdr_0.00001.txt', 'a') ######################################################################################################################################
    a.write(str(fdr) + '\t') ######################################################################################################################################
    a.close()    ######################################################################################################################################
    for prot in targetproteins:
     	if prot not in prot2q:
            prot2q[prot]=fdr
    return prot2q

															    
															# psuedo script from Lukas to calcluate FDR    
															#for t=0:tiny_step:1.0:
															 #  tagetProteins=parsimony(targetPeprtides_with_q_less_than_t)
															  # decoyProteins=parsimony(decoyPeprtides_with_q_less_than_t)
															   #fdr = len(decoyProteins)/len(targetProteins)
															   #for prot in targetProteins:
															    #  if prot not in prot2q:
															     #    prot2q[prot]=fdr





# MAIN PROGRAM


#[5]
fastaFile = "iPRG2016.fasta"  
spec_files = glob.glob("mz_ml_files/*.mzML")
setNames=[os.path.splitext(fileName)[0] for fileName in spec_files]

print spec_files, fastaFile
 

#[6]
'''
subprocess.call(["crux", "tide-index", 
      fastaFile, "prest-index"])

for setN in setNames:
    options = "--output-dir %s-output --overwrite T --compute-sp T --overwrite T --precursor-window 5.0 --precursor-window-type ppm" % setN  
    tide_call = "crux tide-search %s %s.mzML prest-index" % (options,setN) 
    print "Matching %s.mzML against the amino acid sequences."%(setN)
    subprocess.check_output(tide_call,  shell=True)




## [7]. Here we call percolator, also from crux. At this point we have several subdirectories that will be spit out from crux output. I think it is "prest-index", "crux-output", and "mz_ml_files/*_output" 
##So far we are only working with the mz_ml_files output.

for setN in setNames:
    print "Post processing %s\'s results."%(setN)
    options = "--output-dir %s-output --overwrite T" % setN  
    percolator_call = "crux percolator %s ./%s-output/tide-search.target.txt" % (options,setN)

    subprocess.check_output(percolator_call, shell=True)

'''


## [8] Only report the PrEST sequences

proteinRunDict = dict()
for prot in getAllProteinNames(fastaFile):
    
    if prot[:4]=='HPRR':
        proteinRunDict[prot] = []


for setN in setNames:
    print "Parsing %s\'s results."%(setN)
    protDict = getProteinWithFDR("%s-output/percolator.target.peptides.txt"%(setN),  # using [3] should just return a set of proteins. 
                    "%s-output/percolator.decoy.peptides.txt"%(setN))           ##### Banu : Not sure how this dataframe will hold protein list when the program runs from A1-D3. Please check!!!
    for prot in proteinRunDict:
        if prot in protDict:
            fdr = protDict[prot]
        else:
            fdr = float(1.0)
        proteinRunDict[prot].append(fdr)    
    
    
### Banu : We could use the function [4] to calculate protein level FDR, after we decide which q value threshold to use from examining the ROC curve.  In that case in [8] we call the getProteinWithFDR
### ... and use the script below to write the final result file.    


## [9]. Just write the results and add a description row on top 

with open("my_resultFile_bk_0.00001.txt","w") as outFile:
    csvWriter = csv.writer(outFile, delimiter = '\t',quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    csvWriter.writerow(["FDR"]+setNames)
    for prot in proteinRunDict:
        csvWriter.writerow([prot]+proteinRunDict[prot])
    outFile.close()





