#! /usr/bin/python

import getopt,sys

inputfileName,outputfileName="",""

try:
	opts,args=getopt.getopt(sys.argv[1:],"i:o:",['inputfileName','outputfileName'])
except getopt.GetoptError:
	print "Usage Error"
	sys.exit(2)
for opt,arg in opts:
	if opt in ("-i", "--inputfileName"):
		inputfileName=arg
	if opt in ("-o", "--outputfileName"):
		outputfileName=arg


if inputfileName!="" and outputfileName!="":
	fileHandle=open(inputfileName,"r")
	data=fileHandle.read().splitlines()
	fileHandle.close()
	results={}
	legends=data[0].split("\t")
	for line in data[1:]:
		geneId,geneNames,pval,cats=line.split("\t")
		for cat in cats.split(";"):
			results.setdefault(cat,[]).append(geneId)
	output=open(outputfileName,"w")
	output.write(legends[2]+"\tgeneNames\n")
	for cat in results:
		output.write(cat+"\t"+";".join(results[cat])+"\n")
	output.close()
		



