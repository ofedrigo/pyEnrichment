#! /usr/bin/env python

import sys,os
sys.path.append(os.path.realpath(os.path.dirname(sys.argv[0]))+"/tools/")
dirToPanther=os.path.realpath(os.path.dirname(sys.argv[0]))+"/panther/"

from mapNames import *

cat_files={"BP":dirToPanther+'panther_biol_proc_acc_to_name.yml','PW':dirToPanther+'panther_pathw_acc_to_name.yml','MF':dirToPanther+'panther_mol_funct_acc_to_name.yml'}
uniprot_files={"BP":dirToPanther+'uniprot_id_to_panther_biol_proc_accs.yml',"PW":dirToPanther+'uniprot_id_to_panther_pathw_accs.yml',"MF":dirToPanther+'uniprot_id_to_panther_mol_funct_accs.yml'}

geneNameType="any"
inputfileName,outputfileName,category="","",""
try:
	opts,args=getopt.getopt(sys.argv[1:],"i:o:c:",['inputfileName','outputfileName','category'])
except getopt.GetoptError:
	print "Usage Error"
	sys.exit(2)
for opt,arg in opts:
	if opt in ("-i", "--inputfileName"):
		inputfileName=arg
	if opt in ("-o", "--outputfileName"):
		outputfileName=arg
	if opt in ("-c", "--category"):
		if arg in ["BP","PW","MF"]:
			category=arg

if inputfileName!="" and outputfileName!="" and category!="":
	panther,pantherDict=findCategories(uniprot_files[category],cat_files[category],category) # get panther names database and panther database
	dictionary=findNames(geneNameType,panther.keys()) # get dictionary to uniprot Ids	
	allValues,allNames=getInputFile(inputfileName) # read input file name
	catsPerGene=assignCategories(allNames,panther,pantherDict,dictionary,category) #assign
	listOfGenes=sortKeysByVals(allValues) #sort list of genes from smaller to largest values
	saveMappedOntology(outputfileName,listOfGenes,allValues,allNames,catsPerGene,category)