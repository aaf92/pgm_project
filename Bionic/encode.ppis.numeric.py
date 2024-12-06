import pandas as pd 
import sys 

# commandline args for the input and the output files

infile = str(sys.argv[1]) # the file with the PPIs (co-complex and direct ppis union)
infile_gene_order = str(sys.argv[2]) # the order of the genes 
outfile = str(sys.argv[3]) # the name of the output file where the encoded PPIs need to be stored 

# read in the dataframe 

data = pd.read_csv(infile, header=None, sep="\t")
data_order = pd.read_csv(infile_gene_order,sep="\t",header=None)

# make the dictionaries 

genes = list(data_order[0])

genes2encode = {}

for idx,item in enumerate(genes):
	genes2encode[item] = idx+1

ppis = list(zip(data[0],data[1]))
count  = 0

handle = open(outfile,'w')

for p1,p2 in ppis:

	try:
		handle.write(str(genes2encode[p1])+"\t"+str(genes2encode[p2])+"\n")

	except KeyError:
		count = count +1 		
		continue 


handle.close()

print (count)

print("done!")
