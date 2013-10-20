#! /usr/bin/env python

import sys,os
sys.path.append(os.path.realpath(os.path.dirname(sys.argv[0]))+"/tools/")
dirToGo=os.path.realpath(os.path.dirname(sys.argv[0]))+"/go/"

from mapNames import *

cat_files={"BP":dirToGo+'go_biol_proc_id_to_name.yml','CC':dirToGo+'go_cell_comp_id_to_name.yml','MF':dirToGo+'go_mol_funct_id_to_name.yml'}
uniprot_files={"BP":dirToGo+'uniprot_id_to_go_biol_proc_ids.yml',"CC":dirToGo+'uniprot_id_to_go_cell_comp_ids.yml',"MF":dirToGo+'uniprot_id_to_go_mol_funct_ids.yml'}

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
		if arg in ["BP","CC","MF"]:
			category=arg

if inputfileName!="" and outputfileName!="" and category!="":

	go, goDict=findCategories(uniprot_files[category],cat_files[category],category) # get go names database and go database
	dictionary=findNames(geneNameType,go.keys()) # get dictionary to uniprot Ids
	allValues,allNames=getInputFile(inputfileName) # read input file name
	catsPerGene=assignCategories(allNames,go,goDict,dictionary,category) #assign
	listOfGenes=sortKeysByVals(allValues) #sort list of genes from smaller to largest values
	saveMappedOntology(outputfileName,listOfGenes,allValues,allNames,catsPerGene,category)