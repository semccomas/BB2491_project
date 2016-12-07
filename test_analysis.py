## testing to see if we can compare the results from example_results to the actual information in the fasta file


import sys 
import numpy as np
import matplotlib.pyplot as plt

probs = open(sys.argv[1]).read().splitlines()     # study_package/example_results/my_iPRG2016_results.txt 
fasta_a = open(sys.argv[2]).read().splitlines()    # study_package/fasta/prest_pool_a.fasta


fasta_a_list = [ ] 
for line in fasta_a:
	if '>' in line:
		fasta_a_list.append(line[1:])
probs_names = [ ] 
def calculations (probs_names, number):
	probs_names = [ ] 
	for line in probs[1:]:
		line = line.split()
		if number == 3:
			if float(line[1]) + float(line[2]) + float(line[3]) == 3:
				probs_names.append(line[0][1:-1])
		else:
			if float(line[1]) + float(line[2]) + float(line[3]) > number:
				probs_names.append(line[0][1:-1])

	count = 0
	for line in fasta_a_list:
		if line in probs_names:
			count = count +1
	FPR = (len(probs_names) - count) / float(len(probs_names))
	TPR = count / float(len(probs_names))
	FNR = (len(fasta_a_list) - count) / float(len(probs_names))
	return FPR, TPR, FNR, len(probs_names), len(fasta_a_list), count

x_axis = []
y_axis = []
for num in xrange(4):
	print num
	x_axis.append(calculations(probs_names, num)[0])   #list of false positive rates at certain threshold, will be x axis on roc
	y_axis.append(calculations(probs_names, num)[1])   # list of true positive rates at certain threshold, will be y axis on roc
	print calculations(probs_names, num)
plt.plot(x_axis, y_axis)
plt.xlim(0,1)
plt.ylim(0,1)

plt.show()
 




