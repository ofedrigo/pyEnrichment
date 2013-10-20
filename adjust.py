import getopt,sys,sets,time,os

suffix=str(time.time())
R_inp_FILE="/tmp/R_inp_adjust_"+suffix
R_outp_FILE1="/tmp/R_outp_adjust1_"+suffix
R_outp_FILE2="/tmp/R_outp_adjust2_"+suffix
methods=["holm","hochberg","hommel","bonferroni","BH","BY","fdr"]
inputFileName=""
outputFileName=""
method=""
occ=1

# Process options.
try:
	opts, args = getopt.getopt(sys.argv[1:], "i:o:m:c:", ['inputFileName','outputFileName','method',"occupancy"])
except getopt.GetoptError:
	print "Usage Error"
	sys.exit(2)
for opt,arg in opts:
	if opt in ("-i", "--inputFileName"):
		inputFileName=arg
	if opt in ("-o", "--outputFileName"):
		outputFileName=arg
	if opt in ("-m", "--method"):
		method=arg
	if opt in ("-c", "--occupancy"):
		occ=int(arg)

if inputFileName!="" and outputFileName!="" and method!="":
	if method not in methods:
		print "correction method "+method+" unkown - please use holm, hochberg, hommel, bonferroni, BH, BY, or fdr"
	else:
		fileHandle=open(inputFileName,"r")
		data=fileHandle.read().splitlines()
		fileHandle.close()
		categories=[]
		pvals=[]
		info=[]
		legend="\t".join((data[0].split("\t"))[1:])
		occIndex=(data[0].split("\t")).index("tot occ")
		for line in data[1:]:
			items=line.split("\t")
			cat=items[0]
			pvalue=items[1]
			occupancy=int(items[occIndex])
			if occupancy>=occ:
				info.append("\t".join(items[1:]))
				categories.append(cat)
				pvals.append(pvalue)

		r_inp_handle = open(R_inp_FILE, 'w')
		r_inp_handle.write('summary_vector<-p.adjust(c('+",".join(pvals)+'),"'+method+'")\n')
		r_inp_handle.write('write(summary_vector, file= "'+R_outp_FILE1+'",sep = "\n")')
		r_inp_handle.close()
		exitCode = os.system("(R --vanilla < "+R_inp_FILE+"> "+R_outp_FILE2+")")
		os.remove(R_inp_FILE)
		os.remove(R_outp_FILE2)
		adjusted=[]
		fileHandle = open(R_outp_FILE1, 'r')
		adjusted=fileHandle.read().splitlines()
		fileHandle.close()
		os.remove(R_outp_FILE1)
		
		output = open(outputFileName, 'w')
		output.write("cat\t"+method+"\t"+legend+"\n")
		for i in range(len(categories)):
			cat=categories[i]
			infos=info[i]
			adj=adjusted[i]
			output.write(cat+"\t"+adj+"\t"+infos+"\n")
		output.close()
