#! /usr/bin/python

import getopt,sys,sets,time,os

suffix=str(time.time())
matrixFile="/tmp/matrixTree_"+suffix
R_inp_FILE="/tmp/R_inpTree_"+suffix
R_outp_FILE="/tmp/R_outpTree_"+suffix
inputFileName=""
outputFileName=""
mappedFile=""
thresh=1.0

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:o:t:m:", ['inputFileName','outputFileName',"mappedFile",'thresh'])
except getopt.GetoptError:
	print "Usage Error"
	sys.exit(2)
for opt,arg in opts:
	if opt in ("-i", "--inputFileName"):
		inputFileName=arg
	if opt in ("-o", "--outputFileName"):
		outputFileName=arg
	if opt in ("-t", "--thresh"):
		thresh=min(float(arg),thresh)
	if opt in ("-m", "--mappedFile"):
		mappedFile=arg

if inputFileName!="" and outputFileName!="" and mappedFile!="":
	categories=sets.Set([])
	fileHandle=open(inputFileName,"r")
	data=fileHandle.read().splitlines()
	fileHandle.close()
	for line in data[1:]:
		items=line.split("\t")
		cat=items[0]
		pvalue=float(items[1])
		if pvalue<=thresh: categories.add(cat)
	results={}

	fileHandle=open(mappedFile,"r")
	data=fileHandle.read().splitlines()
	fileHandle.close()
	for line in data[1:]:
		geneId,geneNames,pval,cats=line.split("\t")
		for cat in cats.split(";"):
			if cat in categories: results.setdefault(cat,sets.Set([])).add(geneId)
	listOfCat=results.keys()
	
	if listOfCat>=2:
		output=open(matrixFile,"w")
		output.write("\t".join(listOfCat))
		for i in range(len(listOfCat)):
			line=""
			for j in range(len(listOfCat)):
				dxy=1-(float(len(list(results[listOfCat[i]].intersection(results[listOfCat[j]]))))/float(min([len(list(results[listOfCat[i]])),len(list(results[listOfCat[j]]))])))
				if line!="": 
					line=line+"\t"+str(dxy)
				else:
					line=str(dxy)
			output.write("\n"+line)
		output.write("\n")
		output.close()
		
		r_inp_handle = open(R_inp_FILE, 'w')
		r_inp_handle.write('x=read.csv("'+matrixFile+'", header=TRUE,sep="\t")'+'\n')
		r_inp_handle.write('m<-as.dist(x)\n')
		r_inp_handle.write('h<-hclust(m,method="average")\n')
		r_inp_handle.write('pdf("'+outputFileName+'")\n')
		r_inp_handle.write('plot(h)\n')
		r_inp_handle.write('dev.off()')
		r_inp_handle.close()
		exitCode = os.system("(R --vanilla < "+R_inp_FILE+"> "+R_outp_FILE+")")
		os.remove(R_inp_FILE)
		os.remove(R_outp_FILE)
		os.remove(matrixFile)
	else:
		print "less than 2 category analyzed - I can't make a tree"
