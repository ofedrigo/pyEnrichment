#! /usr/bin/env python

import getopt,math,sys,sets,os,time
sys.path.append(os.path.realpath(os.path.dirname(sys.argv[0]))+"/tools/")
from mapNames import *

def diff(list1,list2):
	listDiff=sets.Set(list1).difference(sets.Set(list2))
	return list(listDiff)

def outputParameters(fileName,listOfvalues):
	fileHandle = open(fileName, 'w')
	for value in listOfvalues:
		fileHandle.write(str(value)+"\n")
	fileHandle.close()

unclass=None
inputFileName=""
outputFileName=""
unclassifiedCats = ["no phenotypic analysis","Biological process unclassified", "Molecular function unclassified", "Pathway unclassified", "Cellular component unclassified"]
suffix=str(time.time())
vals_in_cat_FILE="/tmp/Wilcoxonvals_in_cat_"+suffix
vals_not_in_cat_FILE="/tmp/Wilcoxonvals_not_in_cat_"+suffix
R_inp_FILE=	"/tmp/R_inpWilcoxon_"+suffix
R_outp_FILE="/tmp/R_outpWilconxon_"+suffix

# Process options.
try:
	opts, args = getopt.getopt(sys.argv[1:], "i:o:k:", ['inputFileName','outputFileName','keep unclassified'])
except getopt.GetoptError:
	print "Usage Error"
	sys.exit(2)
for opt,arg in opts:
	if opt in ("-i", "--inputFileName"):
		inputFileName=arg
	if opt in ("-o", "--outputFileName"):
		outputFileName=arg
	if opt in ("-k", "--keep unclassified"):
		if arg in ["yes","no"]:
			if arg=="yes":
				unclass=True
			else:
				unclass=False

if inputFileName!="" and outputFileName!="" and unclass!=None:

	inpData_toSort,inpData_values=getInputMappedFile(inputFileName) # get input data
	listOfGenes=sortKeysByVals(inpData_toSort) #let's sort it, in case it was not
	genes_per_cat,vals_per_gene,cats_per_gene,cats_per_geneDict=filterData(listOfGenes,inpData_values,unclass)  # here we eliminate unclassified genes

	# calculate it with R
	results={}
	genes=vals_per_gene.keys()
	pValPerCat={}
	for cat in genes_per_cat.keys():
		genes_in_cat= genes_per_cat[cat]
		vals_in_cat=[vals_per_gene[gene] for gene in genes_in_cat]
		outputParameters(vals_in_cat_FILE,vals_in_cat)
		vals_not_in_cat=[vals_per_gene[gene] for gene in (diff(genes,genes_in_cat))]
		outputParameters(vals_not_in_cat_FILE,vals_not_in_cat)
		r_inp_handle = open(R_inp_FILE, 'w')
		r_inp_handle.write("in_cat <- scan('"+vals_in_cat_FILE+"')\n")
		r_inp_handle.write("not_in_cat <- scan('"+vals_not_in_cat_FILE+"')\n")
		r_inp_handle.write("wilcox.test(in_cat, not_in_cat, alternative = 'less')\n") # Alternative hypothesis: genes within category are more differentially expressed than genes outside category (pvalue smaller).
		r_inp_handle.close()
		exitCode = os.system("(R --vanilla < "+R_inp_FILE+"  > "+R_outp_FILE+")")
		pvalue= openPameters(R_outp_FILE)
		results[cat]=pvalue
		print cat,pvalue
	os.remove(R_inp_FILE)
	os.remove(R_outp_FILE)
	os.remove(vals_in_cat_FILE)
	os.remove(vals_not_in_cat_FILE)
	# output it
	outputHandle= open(outputFileName, 'w')
	cats=sortKeysByVals(results)
	outputHandle.write("cat\tp-val\ttot occ\n")
	for cat in cats: outputHandle.write(cat+"\t"+str(results[cat])+"\t"+str(len(genes_per_cat[cat]))+"\n")
	outputHandle.close()