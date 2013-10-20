#! /usr/bin/env python

import math,getopt,sys,sets
import syck,os

dirToUniprot=os.path.realpath(os.path.dirname(sys.argv[0]))+"/uniprot/"
dirToRefSeq=os.path.realpath(os.path.dirname(sys.argv[0]))+"/refseq/"

unclassifiedCats = {"MGI":"no phenotypic analysis","BP":"Biological process unclassified","MF":"Molecular function unclassified","PW":"Pathway unclassified","CC":"Cellular component unclassified"}
geneTypes={'hgnc':dirToUniprot+'hgnc_symb_to_uniprot_ids.yml','refseq_Peptide':dirToUniprot+'refseq_pept_acc_to_uniprot_ids.yml','refseq_mRNA':dirToRefSeq+'refseq_mrna_acc_to_pept_acc.yml','entrez':dirToUniprot+'entrez_gene_id_to_uniprot_ids.yml',"ensembl":dirToUniprot+'ensembl_gene_id_to_uniprot_ids.yml',"gene_Symbol":dirToUniprot+'gene_name_to_uniprot_ids.yml',"uniprot_Name":dirToUniprot+'uniprot_prim_acc_to_id.yml'}
ontologyName={"MGI":"Mouse Genome Informatics Phenotype Database","BP":"biological process(es)","MF":"molecular function(s)","PW":"pathway(s)","CC":"cellular component(s)"}
allGeneTypes=["hgnc","entrez","ensembl","gene_Symbol","refseq_Peptide","refseq_mRNA","uniprot_Name","uniprot_ID","any"]

def openPameters(fileName):
	pvalue=-1
	fileHandle = open(fileName, 'r')
	inputLine=fileHandle.read().splitlines()
	fileHandle.close()
	for line in inputLine:
		index=line.find("p-value = ")
		if index!=-1:
			pvalue=float(line[index+len("p-value = "):])
	if pvalue==-1: 
		for line in inputLine:
			index=line.find("p-value < ")
			if index!=-1:
				pvalue=0.0 #float(line[index+len("p-value < "):])
	return pvalue

def get_yaml(fileName):
	fileHandle = open(fileName, 'r')
	fileData=fileHandle.read()
	results=syck.load(fileData)
	fileHandle.close()
	#print results
	return results

def concatenateDict(dict1,dict2):
	for keyName in dict2.keys():
		if dict1.has_key(keyName)==False:
			dict1[keyName]=dict2[keyName]
		else:
			list1= dict1[keyName]
			list2= dict2[keyName]			
			dict1[keyName]=list(sets.Set(list1) | sets.Set(list2))
	return dict1

def findCategories(fileName_ontology,fileName_ontologyDict,category):
	results1=get_yaml(fileName_ontology)
	results2=get_yaml(fileName_ontologyDict)
	print "Database for "+ontologyName[category]+" loaded"
	return results1,results2

def getVariantDictionary(fileName):
	dictionary={}
	dictionary_temp=get_yaml(geneTypes[geneNameType])
	for dictKey in dictionary_temp: dictionary[dictKey]=[dictionary_temp[dictKey]]
	return dictionary

def findThisName(geneNameType,uniprotList):
	dictionary={}
	if geneNameType in allGeneTypes[0:5]:  # get directly dictionary ["hgnc","entrez","ensembl","gene_Symbol","refseq_Peptide"] to uniprot Ids (provided by panther)
		dictionary=get_yaml(geneTypes[geneNameType])
		# ex:
		# 0 "hgnc" KRTAP13-1: KR131_HUMAN Q14D20_HUMAN
		# 1 "entrez" "140258": KR131_HUMAN
		# 2 "ensembl" ENSG00000198390: KR131_HUMAN;Q14D20_HUMANN
		# 3 "gene name" KRTAP13-1:  KR131_HUMAN,Q14D20_HUMAN
		# 4 "refseq_Peptide" NP_853630: KR131_HUMAN		
	elif geneNameType in allGeneTypes[5:8]: # get the other name types 5"refseq_mrNA",6"uniprot_Name",7"uniprot_ID" with uniprot and refseq tables
		if geneNameType==allGeneTypes[7]: # ex: 7 uniprot_ID Q5JYW3_HUMAN: Q5JYW3_HUMAN
			for uniprotID in uniprotList: dictionary[uniprotID]=[uniprotID]
		if geneNameType==allGeneTypes[6]: # ex: 6 uniprotName Q5JYW3: Q5JYW3_HUMAN
			dictionary={}
			dictionary_temp=get_yaml(geneTypes[geneNameType])
			for dictKey in dictionary_temp: dictionary[dictKey]=[dictionary_temp[dictKey]]
		if geneNameType==allGeneTypes[5]: # refseq_mRNA=>refseq_peptide=>uniprot
			dictionary_mRNA=get_yaml(geneTypes[geneNameType])
			dictionary_Peptide=get_yaml(geneTypes[allGeneTypes[4]])
			for mRNA in dictionary_mRNA.keys():
				peptide=dictionary_mRNA[mRNA]
				if dictionary_Peptide.has_key(peptide): dictionary[mRNA]=dictionary_Peptide[peptide]
	return dictionary

def findNames(geneNameType,uniprotList):
	dictionary={}
	# allGeneTypes=[0"hgnc",1"entrez",2"ensembl",3"gene_Symbol",4"refseq_Peptide",5"refseq_mrNA",6"uniprot_Name",7"uniprot_ID",8"any"]
	if geneNameType in allGeneTypes[0:8]: #not any, just one name type
		dictionary=findThisName(geneNameType,uniprotList)
	elif geneNameType=="any": # IT IS DONE SEQUENCTIALLY!
		for geneNameType_loop in allGeneTypes[0:8]: # loop through all name types and concatenate dictionnaries (key=uniprot)
			dictionary=concatenateDict(dictionary,findThisName(geneNameType_loop,uniprotList))
	print "Uniprot database loaded for "+geneNameType+" names"
	return dictionary

def assignCategories(allGenes,ontology,ontologyDict,uniprotDict,category):
	#assign 
	results={}
	for geneId in allGenes.keys():
		geneName=allGenes[geneId]
		uniprots=sets.Set([])
		for subGeneName in geneName.split(","):				
			if uniprotDict.has_key(subGeneName.strip(" ").strip('"')):
				for subName in uniprotDict[subGeneName.strip(" ")]: uniprots.add(subName)
		#print geneName,",".join(list(uniprots))
		if len(list(uniprots))==0:
			results[geneId]=unclassifiedCats[category]
		else:
			listOfOntologies=sets.Set([])
			for uniprot in list(uniprots):
				if ontology.has_key(uniprot):
					for thisCat in ontology[uniprot]: listOfOntologies.add(thisCat)
			if len(listOfOntologies)>0:
				results[geneId]=";".join([ontologyDict[thisCat] for thisCat in list(listOfOntologies)])
			else:
				results[geneId]=unclassifiedCats[category]
	print "Categories  "+ontologyName[category]+" assigned to genes"
	return results


def getInputFile(inputfileName): # read input file name
	results_vals,results_names={},{}
	fileHandle=open(inputfileName,'r')
	inputData=fileHandle.read().splitlines()
	fileHandle.close()
	for line in inputData[1:]:
		items=line.split("\t")
		geneId=items[0]
		geneName=items[1].strip('"')
		pvalue=float(items[2])
		results_names[geneId]=geneName
		results_vals[geneId]=pvalue
	print "Input file "+inputfileName+" loaded"
	return results_vals,results_names

def getInputMappedFile(inputFileName):
	inpData_toSort,inpData_values={},{}	
	fileHandle = open(inputFileName, 'r')
	data = fileHandle.read().splitlines()
	fileHandle.close()
	for line in data[1:]:
		items=line.split("\t")
		geneId=items[0]
		value=float(items[2])
		inpData_toSort[geneId]=value
		inpData_values[geneId]="\t".join(line.split("\t"))
	return inpData_toSort,inpData_values
		
def sortKeysByVals(d): # Returns the keys of dictionary d sorted by the values.
	pairs = [[v, k] for [k, v] in d.items()]
	pairs.sort()
	print "List of genes for input values has been ascendant sorted"
	return [pairs[i][1] for i in range(len(pairs))]

def sortKeysByValsReverse(d): # Returns the keys of dictionary d sorted by the values.
	pairs = [[v, k] for [k, v] in d.items()]
	pairs.sort(reverse=True)
	print "List of genes for input values has been descendant sorted"
	return [pairs[i][1] for i in range(len(pairs))]

def saveMappedOntology(outputfileName,listOfGenes,allValues,allNames,catsPerGene,category):
	outputFile = open(outputfileName, 'w')
	outputFile.write("geneId\tgeneName\tvalue\t"+ontologyName[category])
	for geneName in listOfGenes:
		outputFile.write("\n"+geneName+"\t"+allNames[geneName]+"\t"+str(allValues[geneName])+"\t"+catsPerGene[geneName])
	outputFile.close()
	print "Done"

def filterData(listOfGenes,inpData_values,unclass):
	genes_per_cat,vals_per_gene,cats_per_gene,cats_per_geneDict = {},{},[],{}
	#print unclassifiedCats
	for geneId in listOfGenes:
		line=inpData_values[geneId]
		items=line.split("\t")
		catsStr = items[3]
		pvalue=items[2]
		if (catsStr not in unclassifiedCats.values() and unclass==False) or unclass==True:  # here we eliminate unclassified genes if required
			for cat in catsStr.split(";"):
				genes_per_cat.setdefault(cat,[]).append(geneId)
				cats_per_geneDict.setdefault(geneId,[]).append(cat)
			vals_per_gene[geneId]=float(pvalue)
			cats_per_gene.append(catsStr.split(";"))
	return genes_per_cat,vals_per_gene,cats_per_gene,cats_per_geneDict