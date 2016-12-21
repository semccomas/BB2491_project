## this is to figure out what to add to parsimony file and also to compare to the two peptide rule 


import sys
import numpy as np
#f = open(sys.argv[1]).read().splitlines()

# line[9] is q val
# line[11] is peptide sequence
# line[14] is peptide ID

'''

for line in f:
	line = line.split('\t')
	raw_input()
	print line[14]


with open(sys.argv[1], 'r') as f:
    f= f.read().splitlines()
    for line in f:
        line = line.split('\t')
        q_value_str = line[9]
        pep_seq = line[11]
        pep = line[14]
        print pep, pep_seq, q_value_str
#        raw_input()


'''


import subprocess
import glob
import sys
import os
import csv

## copied from iprg2016_two_peptide_rule.ipyn


### all imports copied here.

import subprocess
import glob
import sys
import os
import csv



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

#import csv







## in [2]. This is the two peptide rule. We call this function in [3] below.

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






### ok so here we take the decoy and the target files from each pool (A1, A2, etc). 
### VARIABLES: proteins, is the output from parsePeptieFileUsingTwoPeptidesRule(therefore proteinList).
############ fdr = this is the number of decoys divided by the number of targets. IF i understand this correctly, it is recalculated at each line as it loops through "proteins" ...
############ ... i.e. proteinList. You also add this value so it is paired to the protein name under the output, protFdrDict. One also can indicate the threshold for the fdr and the script will report ...
############ ... back a number of proteins found at that FDR. 
############ protFdrDict = This is a dictionary in python, which is essentially a list of values that are paired to each other. Here, the paired values are the protein names and the fdr for that protein....
########### ... example:   'CAQ33735': 0.0,   'CAQ33734': 0.16252390057361377,    'HPRR4340063': 0.3311897106109325,    'HPRR1450141': 0.0,    
############ order and q are protein names and some number in the function, I don't totally understand their purpose. I'm thinking something about a q value threshold.

##in [3]. Will call this in 8 below but we should also change this to take away the FDR part and also adapt it to the parsimony. 
def getProteinWithFDR(targetFile,decoyFile):
    proteins = parsePeptieFileUsingTwoPeptidesRule(targetFile,False)
    proteins += parsePeptieFileUsingTwoPeptidesRule(decoyFile,True)
    proteins.sort(reverse=True, key=lambda prot_tupple: prot_tupple[0])
    protFdrDict = dict()
    targets,decoys = 0,0
    order = []
    search = True
    for score,isDecoy,prot in proteins:
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

##### DONT FORGET THAT HERE WE WONT NEED FDR. SO WE DONT NEED THIS FDR THING HERE 
a = protFdrDict
b = proteinRunDict

'''


## in [9]. Just write the results and add a top column. 


with open("my_resultFile.txt","w") as outFile:
    csvWriter = csv.writer(outFile, delimiter = '\t',quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    csvWriter.writerow(["FDR"]+setNames)
    for prot in proteinRunDict:
        csvWriter.writerow([prot]+proteinRunDict[prot])
    outFile.close()

'''
