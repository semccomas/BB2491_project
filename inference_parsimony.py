## copied from iprg2016_two_peptide_rule.ipyn
### subsequently copied from changed_example_soln.py and adding parsimony.py. Using this script to make actual changes and mess some things up 

### all imports copied here.

import subprocess
import glob
import sys
import os
import csv
from pygraph.classes.graph import graph
from pygraph.algorithms.accessibility import connected_components
from itertools import combinations, chain


### in [1]. This is just parsing the fasta file. We will call this function in [8] below. 

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





## in [2]. This is the two peptide rule. We call this function in [3] below. Even though we are using parsimony, it does not hurt to leave this function here, ...
### ... I never end up calling it at any point. It's good to be able to compare the output so that we know we didn't fuck up along the way. 

## VARIABLES HERE: fileName = the name of the file in the percolator output files. Specified in [5] EX: for me it is study_package/mz_ml_files/A1-output/percolator.target.peptides.txt.
###### isDecoy = either True or False depending on if its A1-output/percolator.decoy.peptides.txt or percolator.target.peptides.txt. 
###### proteinList = a list of 3 items. A score from percolator (NOT the Q val, something else but idk what, it's a column made in the percolator.target.peptides.txt), either True or False depending on if it is the decoy file, and the protein name...
#### ... so the total output looks like :('-0.65778631', True, 'CAQ31710'), ('0.17604649', True, 'HPRR4000237'), ('-0.37825835', True, 'HPRR3780010') etc.

## so this we will replace with the parsimony. I THINK the only thing we need to change is the FDR, otherwise they act the same (ie. both take as input only one file, either the target or the decoy)...
## ... and it will need to have the same output, proteinList (explained above)
## I will explain in [8] how they loop through the percolator.target and percolator.decoy files to get these outputs. 
def parsePeptieFileUsingTwoPeptidesRule(fileName, isDecoy):
    with open(fileName, "r") as fileIter:
        tabReader = csv.reader(fileIter,doublequote = False, delimiter = '\t')
        colNames = tabReader.next()
        proteinsCol=colNames.index("protein id")
        scoreCol=colNames.index("percolator score")
        lowNumber = -0.4711e10
        proteinDict = dict()
        for row in tabReader:
            proteins = row[proteinsCol]
            score = row[scoreCol]
            for protein in proteins.split(','):
                if protein in proteinDict:
                    if proteinDict[protein] == lowNumber:
                        # assign the score of the second best scoring peptide to the protein
                        proteinDict[protein] = score
                else:
                    proteinDict[protein] = lowNumber
    proteinList = []
    for prot in proteinDict:
        if proteinDict[prot] != lowNumber:
            proteinList.append((proteinDict[prot],isDecoy,prot))
    return proteinList


######## NOW WE ARE ADDING PARSIMONY HERE!!!!!!!!!!!!!!!!!! 


# parsimony.py - takes percolator xml output file and outputs parsimonous
#                protein identifications with peptide matches
#



# Outputs the identified proteins identified using parsimony. If a pair of 
# equally sized sets of proteins both cover all peptides in their group, the 
# set that matched first is returned.
#
# Tab-delimited output format: 
# <protein name>\t<peptide 1>\t...\t<peptide N>


### so now this will parse the file the same way as the two peptide rule. Filename here is setN in setNames, (see [8]) aka the percolator.decoy or .target. This should be changed to be in [8] as setN and the %s thing. This is the ...
### ... first part of the parsePeptideFileUsingTwoPeptideRule function (attributing values to variables and reading the file). 
### VARIABLES: peptides. A dictionary of the sequence attached to the protein name for each file. I think maybe later we will change this to the q value score so that it is matching ...
### .. the script with two peptide rule. 

def parse_percolator_txt(filename, q_value_cutoff= 0.01):
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
            if q_value > q_value_cutoff:
                pass
            protein_ids = protein_str.split(',')
            peptides[pep_seq] = protein_ids
        return peptides


### this is the second part of what is replacing the parsePeptideFileUsingTwoPeptideRule function. I don't understand really anything of this script's inner ...
### ... workings but I think it's just about parsimony. 

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


#### so what we are doing here is making a function to run the actual 2 functions from parsimony.
#### this we will put into getProteinWithFDR. 
### returns the same thing as the last bit in the parsimony.py script, just now in function form with added True and False. I think we won't ...
### ... end up needing the True or False in the file, since we will change how to FDR works but i am not sure yet.... 
def runParsimony(filename, isDecoy):
    peptides = parse_percolator_txt(filename)
    proteins = parsimonous_protein_identification(peptides)
    protein_list = [] 
    for name, seq in proteins.iteritems():
        if isDecoy:
            tup = (name, True)
            protein_list.append(tup)
        else:
            tup = (name, False)
            protein_list.append(tup)
    return protein_list


  



##in [3]. Will call this function in 8 below but we should also change this to take away the FDR part and also adapt it to the parsimony.

### ok so here we take the decoy and the target files from each pool (A1, A2, etc). 
### VARIABLES: proteins, is the output from parsePeptieFileUsingTwoPeptidesRule(therefore proteinList).
############ fdr = this is the number of decoys divided by the number of targets. IF i understand this correctly, it is recalculated at each line as it loops through "proteins" ...
############ ... i.e. proteinList. You also add this value so it is paired to the protein name under the output, protFdrDict. One also can indicate the threshold for the fdr and the script will report ...
############ ... back a number of proteins found at that FDR. 
############ protFdrDict = This is a dictionary in python, which is essentially a list of values that are paired to each other. Here, the paired values are the protein names and the fdr for that protein....
########### ... example:   'CAQ33735': 0.0,   'CAQ33734': 0.16252390057361377,    'HPRR4340063': 0.3311897106109325,    'HPRR1450141': 0.0,    
############ order and q are protein names and some number in the function, I don't totally understand their purpose. I'm thinking something about a q value threshold.



### ideas 18 December- compare the two outputs from parsePeptieFileUsingTwoPeptidesRule and runParsimony. Where do they differ? We are getting same ...
### ... FDR calc for all vals on runParsimony so something is off on the FDR calc. Number of peptides listed is pretty normal, we just are missing some form of calculation here
### or maybe even in the setN in setNames part. But problem is we don't have a way to calculate FDR otherwise. How does parsimony do this? Are we supposed to be adding the q value
### back in? 

#### UPDATE 20 December. As of right now, this is still definitely the part in the script that I think needs fixing. The problem is that we get the ...
##### ... same FDR number for each set (i.e. in row A2, they are all 0.02323534345 or something.) So obviously that needs fixing. Maybe we are supposed to incorporate that indiv q value

def getProteinWithFDR(targetFile,decoyFile):
   # proteins = parsePeptieFileUsingTwoPeptidesRule(targetFile,False)
   # proteins += parsePeptieFileUsingTwoPeptidesRule(decoyFile,True)
    proteins = runParsimony(targetFile, False)
    proteins += runParsimony(decoyFile, True)
    proteins.sort(reverse=True, key=lambda prot_tupple: prot_tupple[0])
    protFdrDict = dict()
    targets,decoys = 0,0
    order = []
    search = True
    for prot, isDecoy in proteins:
        if isDecoy:
            decoys += 1
        else:
            targets += 1
            fdr = decoys/float(targets)
            protFdrDict[prot] = fdr
            order.append(prot)
            if (search and fdr>0.01):
                print "Found %d proteins at FDR %f"%(targets,fdr)
                search = False
    order.reverse()
    # Set the q value to be the minimal FDR that includes the current protein
    q=1.0
    for prot in order:
        q = min(q,protFdrDict[prot])
        protFdrDict[prot] = q
    return protFdrDict

## skip in [4] because it was converting to mzml format. Already done.


### in [5] the script is just looking in your directory that you are running the script from to see where the mzml files are. In my directory ...
### ... they are in a subdirectory called mz_ml_files so that is why I specify that here. If you want to run this in your own computer, make sure that you change this name to the name of the directory where your mzml files are stored. 

fastaFile = "iPRG2016.fasta"  
spec_files = glob.glob("mz_ml_files/*.mzML")
setNames=[os.path.splitext(fileName)[0] for fileName in spec_files]

print spec_files, fastaFile




## in [6]. This requires that crux is installed and also in your path of executable files. For this, I write " export PATH=$PATH:crux-3.0.Darwin.i386/bin " in the command ...
### ... line each time I open a terminal. If you want it to be automatic, add it to your ~/.bashrc . I have crux in the same directory as this script but I don't think it's necessary to do that.
### UPDATE 20 December. I have these two parts commented out right now because I have already run crux once so I don't need to wait for it to run while I already can work with the output files
### if you want to run it with crux, take away the    '''   below and after the subprocess.call bits.   

'''
subprocess.call(["crux", "tide-index", 
      fastaFile, "prest-index"])

for setN in setNames:
    options = "--output-dir %s-output --overwrite T --compute-sp T --overwrite T --precursor-window 5.0 --precursor-window-type ppm" % setN  
    tide_call = "crux tide-search %s %s.mzML prest-index" % (options,setN) 
    print "Matching %s.mzML against the amino acid sequences."%(setN)
    subprocess.check_output(tide_call,  shell=True)




## in [7]. Here we call percolator, also from crux. At this point we have several subdirectories that will be spit out from crux output. I think it is "prest-index", "crux-output", and "mz_ml_files/*_output" 
##So far we are only working with the mz_ml_files output (or whatever you named your directory with the mzml files in it).

for setN in setNames:
    print "Post processing %s\'s results."%(setN)
    options = "--output-dir %s-output --overwrite T" % setN  
    percolator_call = "crux percolator %s ./%s-output/tide-search.target.txt" % (options,setN)

    subprocess.check_output(percolator_call, shell=True)

'''




## in [8]    as it says below in the comment, only report back prEST frags from the fasta file. This is where you call the function from [1]
proteinRunDict = dict()
for prot in getAllProteinNames(fastaFile):
    # Only report the PrEST sequences
    if prot[:4]=='HPRR':
        proteinRunDict[prot] = []

#proteinRunDict is all the prest frags from the fasta file. Not to be confused with protFdrDict. 


#### Ok so here is where the functions in [2] and [3] are called. the function parsePeptieFileUsingTwoPeptidesRule is called in getProteinWithFDR so that's why we don't see it below
#### This will loop through the files A1, A2, A3, B1, B2.... etc -output/percolator. I.e. setN will be A1 , and then when it loops through this for loop again it will be A2, then A3, etc. 
#### It then calls the function getProteinWithFDR and specifies that we only want A1-output/percolator.target(or .decoy).peptides.txt. None of the other things in A1-output.

##### VARIABLES: protFdrDict = the output of getProteinWithFDR, which was the protein and the fdr number (HPRR1951262': 0.6462459695992631)
##### proteinRunDict comes from the part above, it is the prEST fragment names. Now what this part of the script will do is see if the fragment was also in the protFdrDict, and 
##### if it is, assign the fdr number that was already attached to it (i.e if you use the example above, HPRR1951262 it will get the number 0.6462459695992631)
###### if it is NOT, give it the fdr number of 1. So this means our assumption before of the 0's and the 1's in the result file is backwards. A value of 0 means it's there, a value of 1 means ...
###### ... that is is not there, or that the fdr is 100% 


for setN in setNames:
    print "Parsing %s\'s results."%(setN)
    protFdrDict = getProteinWithFDR("%s-output/percolator.target.peptides.txt"%(setN),
                    "%s-output/percolator.decoy.peptides.txt"%(setN))
    for prot in proteinRunDict:
        if prot in protFdrDict:
            fdr = protFdrDict[prot]
        else:
            fdr = float(1.0)
        proteinRunDict[prot].append(fdr)




## in [9]. Just write the results and add a description row on top 

with open("my_resultFile2.txt","w") as outFile:
    csvWriter = csv.writer(outFile, delimiter = '\t',quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    csvWriter.writerow(["FDR"]+setNames)
    for prot in proteinRunDict:
        csvWriter.writerow([prot]+proteinRunDict[prot])
    outFile.close()





