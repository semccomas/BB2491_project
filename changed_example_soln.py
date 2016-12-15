## copied from iprg2016_two_peptide_rule.ipyn


### all imports copied here

import subprocess
import glob
import sys
import os
import csv



### in [1]
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


## in [2]
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


##in [3]
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

##in [5]
fastaFile = "iPRG2016.fasta"  
spec_files = glob.glob("*.mzML")
setNames=[os.path.splitext(fileName)[0] for fileName in spec_files]

print spec_files, fastaFile




## in [6]


subprocess.call(["crux", "tide-index", 
      fastaFile, "prest-index"])

for setN in setNames:
    options = "--output-dir %s-output --overwrite T --compute-sp T --overwrite T --precursor-window 5.0 --precursor-window-type ppm" % setN  
    tide_call = "crux tide-search %s %s.mzML prest-index" % (options,setN) 
    print "Matching %s.mzML against the amino acid sequences."%(setN)
    subprocess.check_output(tide_call,  shell=True)





## in [7]
for setN in setNames:
    print "Post processing %s\'s results."%(setN)
    options = "--output-dir %s-output --overwrite T" % setN  
    percolator_call = "crux percolator %s ./%s-output/tide-search.target.txt" % (options,setN)

    subprocess.check_output(percolator_call, shell=True)






## in [8]
proteinRunDict = dict()
for prot in getAllProteinNames(fastaFile):
    # Only report the PrEST sequences
    if prot[:4]=='HPRR':
        proteinRunDict[prot] = []

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






## in [9]


with open("my_resultFile.txt","w") as outFile:
    csvWriter = csv.writer(outFile, delimiter = '\t',quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    csvWriter.writerow(["FDR"]+setNames)
    for prot in proteinRunDict:
        csvWriter.writerow([prot]+proteinRunDict[prot])
    outFile.close()


