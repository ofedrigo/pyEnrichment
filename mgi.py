#! /usr/bin/env python

import sys,os
sys.path.append(os.path.realpath(os.path.dirname(sys.argv[0]))+"/tools/")
dirToMGI=os.path.realpath(os.path.dirname(sys.argv[0]))+"/mgi/"

from mapNames import *

cat_file=dirToMGI+'mgi_acc_to_name.yml'
uniprot_file=dirToMGI+'uniprot_id_to_mgi_accs.yml'

geneNameType="any"
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
	mgi,mgiDict=findCategories(uniprot_file,cat_file,"MGI") # get mgi names database and mgi database (only uniprot based)
	dictionary=findNames(geneNameType,mgi.keys()) # get dictionary to uniprot Ids
	allValues,allNames=getInputFile(inputfileName) # read input file name
	catsPerGene=assignCategories(allNames,mgi,mgiDict,dictionary,"MGI") #assign
	listOfGenes=sortKeysByVals(allValues) #sort list of genes from smaller to largest values
	saveMappedOntology(outputfileName,listOfGenes,allValues,allNames,catsPerGene,"MGI")