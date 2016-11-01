#!/usr/bin/python

####
#This script parses a vcf file and prints out tumor minor allele frequencies
####
import sys


def usage():
	print ("Usage:")
	print ("vcf_concatenate.py vcfFile")
	sys.exit()


def parsing(vcfFile):
	"This script parses a vcf file"
	sample = vcfFile[0:17]
	#parsed lists
	header = []
	snps = []
	hetSNPs = []
	
	vcf = open(vcfFile)
	for line in vcf:
		line = line.rstrip("\n")
		line = line.split("\t")
		#save the column header list
		if (line[0][0] == "#"):
			if (line[0] == "#CHROM"):
				header = line
			continue
		snps.append(line)
	vcf.close()
	
	#find column numbers
	normal = header.index("NORMAL")
	tumor = header.index("TUMOR")
	CHROM = header.index("#CHROM")
	POS = header.index("POS")
	ID = header.index("ID")
	REF = header.index("REF")
	ALT = header.index("ALT")
	ID = header.index("ID")
	FILTER = header.index("FILTER")
	FORMAT = header.index("FORMAT")
	INFO = header.index("INFO")

	hetSNPs.append(["sample", "chromosome","position", "id", "ref", "alt", "filter", "MAF", "INFO"])
	for line in snps:
		#parse the format column
		snpformat = line[FORMAT].split(":")
		GT = snpformat.index("GT")
		AD = snpformat.index("AD")
		#parse the data columns
		normalData = line[normal].split(":")
		tumorData = line[tumor].split(":")
		#***BEGIN FILTERS***
		if (line[tumor] == "./."):
			continue
		if (tumorData[AD] == "."):
			continue			
		#***END FILTERS***
		tRef = tumorData[AD].split(",")[0]
		tAlt = tumorData[AD].split(",")[1]
		tTotal = int(tAlt) + int(tRef)
		if (tTotal < 10):
			continue
		if (int(tAlt) <= int(tRef)):
			tMAF = 100 * float(int(tAlt)) / tTotal
		else:
			tMAF = 100 * float(int(tRef)) / tTotal
		hetSNPs.append([str(sample), line[CHROM], line[POS], line[ID], line[REF], line[ALT], line[FILTER], str(tMAF), line[INFO]])

	#print the output
	
	file = open(sample + "concatenate.txt", "w")
	for line in hetSNPs:
		file.write("\t".join(line)+'\n')
	file.close()
	
if (len(sys.argv) != 2):
	usage()
else:
	parsing(sys.argv[1])

