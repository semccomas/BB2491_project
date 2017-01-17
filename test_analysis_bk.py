## testing to see if we can compare the results from example_results to the actual information in the fasta file
 
import numpy as np
import matplotlib.pyplot as plt



## want ROC curve for many different q values.... 
## and then we want ROC curves for each 4 pools (A1-3, B1-3, C1-3, D1-3) as well as all the pools total. 


#result_00 = open("my_resultFile_0.0.txt").read().splitlines()

result_001 = open("my_resultFile_bk_0.01.txt").read().splitlines()
result_002 = open("my_resultFile_bk_0.02.txt").read().splitlines()
result_003 = open("my_resultFile_bk_0.03.txt").read().splitlines()
result_004 = open("my_resultFile_bk_0.04.txt").read().splitlines()
result_005 = open("my_resultFile_bk_0.05.txt").read().splitlines()
result_006 = open("my_resultFile_bk_0.06.txt").read().splitlines()
result_007 = open("my_resultFile_bk_0.07.txt").read().splitlines()
result_008 = open("my_resultFile_bk_0.08.txt").read().splitlines()
result_009 = open("my_resultFile_bk_0.09.txt").read().splitlines()
result_01 = open("my_resultFile_bk_0.1.txt").read().splitlines()
result_02 = open("my_resultFile_bk_0.2.txt").read().splitlines()
result_05 = open("my_resultFile_bk_0.5.txt").read().splitlines()
result_06 = open("my_resultFile_bk_0.6.txt").read().splitlines()
result_08 = open("my_resultFile_bk_0.8.txt").read().splitlines()
result_10 = open("my_resultFile_bk_1.0.txt").read().splitlines()
result_00001 = open("my_resultFile_bk_0.00001.txt").read().splitlines()

results = [result_00001, result_001, result_002, result_003, result_004, result_005, result_006, result_007, result_008, result_009, result_01, result_02, result_05, result_06, result_08, result_10]


'''
result_01 = open("my_resultFile_NEW_0.1.txt").read().splitlines()
result_02 = open("my_resultFile_NEW_0.2.txt").read().splitlines()
result_03 = open("my_resultFile_NEW_0.3.txt").read().splitlines()
result_04 = open("my_resultFile_NEW_0.4.txt").read().splitlines()
result_05 = open("my_resultFile_NEW_0.5.txt").read().splitlines()
result_06 = open("my_resultFile_NEW_0.6.txt").read().splitlines()
result_07 = open("my_resultFile_NEW_0.7.txt").read().splitlines()
result_08 = open("my_resultFile_NEW_0.8.txt").read().splitlines()
result_09 = open("my_resultFile_NEW_0.9.txt").read().splitlines()
result_10 = open("my_resultFile_NEW_1.0.txt").read().splitlines()

result_0_025 = open("my_resultFile_0.025.txt").read().splitlines()
result_0_05 = open("my_resultFile_0.05.txt").read().splitlines()
result_0_075 = open("my_resultFile_0.075.txt").read().splitlines()
result_0_125 = open("my_resultFile_0.125.txt").read().splitlines()
result_0_15 = open("my_resultFile_0.15.txt").read().splitlines()
result_0_175 = open("my_resultFile_0.175.txt").read().splitlines()
result_0_225 = open("my_resultFile_0.225.txt").read().splitlines()
result_0_25 = open("my_resultFile_0.25.txt").read().splitlines()
result_0_275 = open("my_resultFile_0.275.txt").read().splitlines()
result_0_325 = open("my_resultFile_0.275.txt").read().splitlines()
result_0_35 = open("my_resultFile_0.35.txt").read().splitlines()
result_0_375 = open("my_resultFile_0.375.txt").read().splitlines()
result_0_425 = open("my_resultFile_0.425.txt").read().splitlines()
result_0_45 = open("my_resultFile_0.45.txt").read().splitlines()
result_0_475 = open("my_resultFile_0.475.txt").read().splitlines()
'''
#results = [result_01, result_02, result_03, result_04, result_05, result_06, result_07, result_08, result_09, result_10]

#results = [result_0_025, result_0_05, result_0_075, result_01, result_0_125, result_0_15, result_0_175, result_02, result_0_225, result_0_25, result_0_275, result_03, result_0_325, result_0_35, result_0_375, result_04, result_0_425, result_0_45, result_0_475, result_05, result_06, result_07, result_08, result_09, result_10]


fasta_pool_a = open('actual_fasta/prest_pool_a.fasta').read().splitlines()
fasta_pool_b = open('actual_fasta/prest_pool_b.fasta').read().splitlines()
fasta_pool_c = open('actual_fasta/prest_1000_random.fasta').read().splitlines()  #Banu: changing this to 1000_random which will be used to calculate FPR in line 77.


print 'lens'
fasta_len_a = len(fasta_pool_a)/ 2
fasta_len_b = len(fasta_pool_b)/ 2
fasta_len_c = len(fasta_pool_c)/ 2


### this is for lukas's thing

fdr_0_00001 = open('fdr_0.00001.txt').read().splitlines()
fdr_0_01 = open('fdr_0.01.txt').read().splitlines()
fdr_0_02 = open('fdr_0.02.txt').read().splitlines()
fdr_0_03 = open('fdr_0.03.txt').read().splitlines()
fdr_0_04 = open('fdr_0.04.txt').read().splitlines()
fdr_0_05 = open('fdr_0.05.txt').read().splitlines()
fdr_0_06 = open('fdr_0.06.txt').read().splitlines()
fdr_0_07 = open('fdr_0.07.txt').read().splitlines()
fdr_0_08 = open('fdr_0.08.txt').read().splitlines()
fdr_0_09 = open('fdr_0.09.txt').read().splitlines()
fdr_0_1 = open('fdr_0.1.txt').read().splitlines()
fdr_0_2 = open('fdr_0.2.txt').read().splitlines()
fdr_0_5 = open('fdr_0.5.txt').read().splitlines()
fdr_0_6 = open('fdr_0.6.txt').read().splitlines()
fdr_0_8 = open('fdr_0.8.txt').read().splitlines()
fdr_1_0 = open('fdr_1.0.txt').read().splitlines()

fdrs = [fdr_0_00001, fdr_0_01, fdr_0_02, fdr_0_03, fdr_0_04, fdr_0_05, fdr_0_06, fdr_0_07, fdr_0_08, fdr_0_09, fdr_0_1, fdr_0_2, fdr_0_5, fdr_0_6, fdr_0_8, fdr_1_0]

a_fdr = [] 
b_fdr = []
c_fdr = []
for fdr in fdrs:
	for line in fdr:
		line = line.split()
		list2 = [float(x) for x in line]
		#print line
		a_fdr.append(sum(list2[:3])/ 3.0)
		b_fdr.append(sum(list2[3:6])/ 3.0)
		c_fdr.append(sum(list2[6:]) / 3.0)
print a_fdr, b_fdr, c_fdr



def compare(result_file, fasta_file, index):
	fasta_list = [ ] 
	for line in fasta_file:
		if '>' in line:
			fasta_list.append(line[1:])
	result_namesA1 = [ ]
	for line in result_file[1:]:
		line = line.split()
		FDR = line[int(index)]
		if float(FDR) != 1.0:
			result_namesA1.append(line[0][1:-1])
	count = 0 
	for line in fasta_list:
		if line in result_namesA1:
			count = count + 1
	#for line in result_namesA1:
	#	count = count + 1
	float(count)   
	return count #Banu



a_FPR_all_qval = [ ] 
b_FPR_all_qval = [ ] 
c_FPR_all_qval = [ ]
# d_FPR_all_qval = [ ] Not sure how to calculate TPR for this coz in principal it shouldn't contain any of A/B and using random_1000 seems meaningless. 
a_TPR_all_qval = [ ] 
b_TPR_all_qval = [ ] 
c_TPR_all_qval = [ ]
# d_TPR_all_qval = [ ]

for result in results:	
	a_FPR = [ ]
	b_FPR = [ ]
	a_FPR1 = [ ] #Banu introducing two varaibales coz FPR is calculated by comparing result to two fasta files.
	b_FPR1 = [ ]
	a_FPR2 = [ ]
	b_FPR2 = [ ]
	c_FPR= [ ] 
	# d_FPR = [ ] 
	a_TPR = [ ]
	b_TPR = [ ]
	c_TPR = [ ]
	c_TPR1 = [ ]
	c_TPR2 = [ ]
	# d_TPR = [ ] 


	for x in xrange(1, 13):
		if x < 4:
			a_TPR.append(float(compare(result, fasta_pool_a, x))) #Banu
			a_FPR.append(float(compare(result, fasta_pool_b, x)) + (float(compare(result, fasta_pool_c, x))))

		if x == 4 or x ==5 or x== 6:
			b_TPR.append(float(compare(result, fasta_pool_b, x)))
			b_FPR.append(float(compare(result, fasta_pool_a, x)) + (float(compare(result, fasta_pool_c, x))))
		if x == 7 or x ==8 or x == 9:
			c_TPR.append(float(compare(result, fasta_pool_a, x)) + (float(compare(result, fasta_pool_b, x))))
			c_FPR.append(float(compare(result, fasta_pool_c, x)))
			
			
		# if x == 10 or x == 11 or x == 12:
		#	d_FPR.append(float(compare(result, fasta_pool_d, x)))
		#	d_TPR.append(float(compare(result, fasta_pool_d, x)[1]))


	a_FPR = (sum(a_FPR) / 3) / (fasta_len_b + fasta_len_c)
	b_FPR = (sum(b_FPR) / 3) / (fasta_len_a + fasta_len_c)
	c_FPR = (sum(c_FPR) / 3) / (fasta_len_c) 
	#c_FPR = c_FPR * 2 
	# d_FPR = sum(d_FPR) / 3

	a_FPR_all_qval.append(a_FPR) 
	b_FPR_all_qval.append(b_FPR)
	c_FPR_all_qval.append(c_FPR)
	# d_FPR_all_qval.append(d_FPR)


	a_TPR = (sum(a_TPR) / 3) / (fasta_len_a)
	b_TPR = (sum(b_TPR) / 3) / (fasta_len_b)
	c_TPR = (sum(c_TPR) / 3) / (fasta_len_a + fasta_len_b)
	#c_TPR = c_TPR /2
	# d_TPR = sum(d_TPR) / 3

	a_TPR_all_qval.append(a_TPR) 
	b_TPR_all_qval.append(b_TPR)
	c_TPR_all_qval.append(c_TPR)
	# d_TPR_all_qval.append(d_TPR)



print 'FPR:'
print 'pool a', a_FPR_all_qval
print 'pool b', b_FPR_all_qval
print 'pool c', c_FPR_all_qval
#print 'pool d', d_FPR_all_qval
print
print
print 'TPR'
print 'pool a', a_TPR_all_qval
print 'pool b', b_TPR_all_qval
print 'pool c', c_TPR_all_qval
#print 'pool d', d_TPR_all_qval




### #FDR calcuations thing for lukas
a_fdr_real = [ ]
for F, T in zip(a_FPR_all_qval, a_TPR_all_qval):
	a_fdr_real.append(F/(F+ T))

b_fdr_real = [ ]
for F, T in zip(b_FPR_all_qval, b_TPR_all_qval):
	b_fdr_real.append(F/(F+ T))

c_fdr_real = [ ]
for F, T in zip(c_FPR_all_qval, c_TPR_all_qval):
	c_fdr_real.append(F/(F+ T))

print a_fdr_real, b_fdr_real, c_fdr_real




a = plt.scatter(a_FPR_all_qval, a_TPR_all_qval, color = 'blue')
b = plt.scatter(b_FPR_all_qval, b_TPR_all_qval, color = 'green')
c = plt.scatter(c_FPR_all_qval, c_TPR_all_qval, color = 'red')
a = plt.plot(a_FPR_all_qval, a_TPR_all_qval, label = 'Sample A', color = 'blue')
b = plt.plot(b_FPR_all_qval, b_TPR_all_qval, label = "Sample B", color = 'green')
c = plt.plot(c_FPR_all_qval, c_TPR_all_qval, label = 'Sample C', color = 'red')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend()
plt.show()
#plt.plot(d_FPR_all_qval, d_TPR_all_qval)
plt.title("ROC curve")
plt.savefig('plot_out.png')

'''s
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
 


'''

ranger = []
rangery= [] 

for x in range(0,10):
	ranger.append(x/ 10.0)
	rangery.append(x/ 10.0)

print ranger, rangery

a_comparison = plt.scatter(a_fdr, a_fdr_real, color = 'blue')
b_comparison = plt.scatter(b_fdr, b_fdr_real, color = 'green')
c_comparison = plt.scatter(c_fdr, c_fdr_real, color = 'red')
a_comparison = plt.plot(a_fdr, a_fdr_real, color = 'blue', label = "Pool A")
b_comparison = plt.plot(b_fdr, b_fdr_real, color = 'green', label = "Pool B")
c_comparison = plt.plot(c_fdr, c_fdr_real, color = 'red', label = 'Pool C')
real = plt.plot(ranger, rangery, color = 'grey')
plt.xlabel('estimated FDR')
plt.ylabel('real FDR')
plt.xlim([0,1])
plt.ylim([0,1])
plt.title("estimated vs actual FDR")
plt.legend()
plt.savefig('plot_estimated.png')














