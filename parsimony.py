#!/usr/bin/env python

# parsimony.py - takes percolator xml output file and outputs parsimonous
#                protein identifications with peptide matches
#


from pygraph.classes.graph import graph
from pygraph.algorithms.accessibility import connected_components
from itertools import combinations, chain
import sys

# Outputs the identified proteins identified using parsimony. If a pair of 
# equally sized sets of proteins both cover all peptides in their group, the 
# set that matched first is returned.
#
# This program assumes that the input percolator results XML validates.
#
# Tab-delimited output format: 
# <protein name>\t<peptide 1>\t...\t<peptide N>


fastaFile = "iPRG2016.fasta"  

### so now this will parse the file the same way as the two peptide rule. Filename here is sys. This should be changed to be in [8] as setN and the %s thing. This is the ...
### ... first part of the parsePeptideFileUsingTwoPeptideRule function (attributing values to variables and reading the file). 
### VARIABLES: peptides. A dictionary of the sequence attached to the protein name for each file. I think maybe later we will change this to the q value score so that it is matching ...
### .. the script with two peptide rule. 

def parse_percolator_txt(filename, q_value_cutoff= 0.01):
    f = open(filename).read().splitlines()
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
            pass # i have NO idea why but when continue is here, it does not get along with the decoy sets. Dont change this or nothing will work for decoy
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

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: parsimony.py <name of percolator xml output file>"
        sys.exit(1)
    per_filename = sys.argv[1]
    peptides = parse_percolator_txt(per_filename)
    proteins = parsimonous_protein_identification(peptides)
    for protein, peptides in proteins.items():
        print "{}\t{}".format(protein, "\t".join(peptides))



def runParsimony(filename, isDecoy):
    peptides = parse_percolator_txt(filename)
    proteins = parsimonous_protein_identification(peptides)
    protein_list = [] 
    #proteins = list(proteins)
    for name, sequence in proteins.iteritems():
        if isDecoy == True:
            tup = (name, True)
            protein_list.append(tup)
        else:
            tup = (name, False)
            protein_list.append(tup)
    return protein_list
 
a= runParsimony(sys.argv[1], True) #and make this F for target

       