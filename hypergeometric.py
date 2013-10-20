#! /usr/bin/env python

import sys,os,time
sys.path.append(os.path.realpath(os.path.dirname(sys.argv[0]))+"/tools/")
from mapNames import *

unclass=None
inputFileName,outputFileName,propTop="","",-1
suffix=str(time.time())
R_inp_FILE=	"/tmp/R_inpFisher_"+suffix
R_outp_FILE="/tmp/R_outpFisher_"+suffix

# Process options.
try:
	opts, args = getopt.getopt(sys.argv[1:], "i:p:o:k:", ['inputFileName',"propTop",'outputFileName','keep unclassified'])
except getopt.GetoptError:
	print "Usage Error"
	sys.exit(2)
for opt,arg in opts:
	if opt in ["-p", "--propTop"]:
		propTop = float(arg)
		if propTop < 0.0 or propTop > 1.0: raise ValueError("BAD PROP_TOP "+opt[1]+", SHOULD BE BETWEEN 0 AND 1")
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

if inputFileName!="" and outputFileName!="" and propTop!=-1 and unclass!=None:

	inpData_toSort,inpData_values=getInputMappedFile(inputFileName) # get input data
	listOfGenes=sortKeysByVals(inpData_toSort) #let's sort it, in case it was not
	genes_per_cat,vals_per_gene,cats_per_gene,cats_per_geneDict=filterData(listOfGenes,inpData_values,unclass)  # here we eliminate unclassified genes

	results,results2={},{}
	listOfCats=sets.Set([])
	N = len(cats_per_gene) # total number of genes
	n = int(math.floor(propTop*N))	# total number of genes in the top propTop %
	topOccPerCat,bottomOccPerCat= {},{}
	for i in range(N):
		for cat in cats_per_gene[i]:
			if not topOccPerCat.has_key(cat): topOccPerCat[cat] = 0
			if not bottomOccPerCat.has_key(cat): bottomOccPerCat[cat] = 0
			if i < n:
				topOccPerCat[cat] += 1
			else:
				bottomOccPerCat[cat] +=1
			listOfCats.add(cat)
				
	for cat in list(listOfCats):
		if topOccPerCat[cat]>0:
			lineForR="ontology<-matrix(c("+",".join([str(topOccPerCat[cat]),str(bottomOccPerCat[cat]),str(n-topOccPerCat[cat]),str((N-n)-bottomOccPerCat[cat])])+"),nr=2,dimnames=list(c('In "+cat+"','Not in "+cat+"'),c('top"+str(propTop*100)+"%','bottom"+str(100-propTop*100)+"%')))"
			r_inp_handle = open(R_inp_FILE, 'w')
			r_inp_handle.write(lineForR+"\n")
			r_inp_handle.write("fisher.test(ontology,alternative ='greater')\n")
			r_inp_handle.close()
			exitCode = os.system("(R --vanilla < "+R_inp_FILE+" > "+R_outp_FILE+")")
			pvalue= openPameters(R_outp_FILE)
			results[cat]=pvalue
			print cat,pvalue
		else:
			results[cat]=1.0
	os.remove(R_inp_FILE)
	os.remove(R_outp_FILE)
	cats=sortKeysByVals(results)
	outputHandle = open(outputFileName, 'w')
	outputHandle.write("cat\tp-val\ttop "+str(propTop*100)+"% occ\ttot occ\n")
	for cat in cats: outputHandle.write(cat+"\t"+str(results[cat])+"\t"+str(topOccPerCat[cat])+"\t"+str(topOccPerCat[cat]+bottomOccPerCat[cat])+"\n")
	outputHandle.close()