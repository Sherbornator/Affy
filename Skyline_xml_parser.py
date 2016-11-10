#!/usr/bin/python

####
#This script parses a bunch of Skyline .xml reports and makes a data table from them
####

import sys
import re
import glob
import time

def usage():
	print ('Usage:')
	print ('Skyline_xml_parser.py path/to/dataFolder')
	sys.exit()
	
def parsing(dataFolder):
	print ('Parsing this location:' + dataFolder + ' recursively looking for Report***.xml files')
	
	#set up a table with the column headings for each variable
	table = []
	#NB - the order of the column headings is very important, otherwise the wrong data could be shown for the wrong values
	table.append(['pID', 'SKY92', 't(4;14)', 't(11;14)', 't(14;16) or t(14;20)', 'gain(1q)', 'hyperdiploid MM-gain(9)', 'del(13q)', 'del(17p)', 'MScluster', 'MFcluster', 'NFKB', 'CD2', 'CD38', 'CD38 min of range', 'CD38 max of range', 'CD138', 'CD138 min of range', 'CD138 max of range', 'CRBN', 'CRBN cut point', 'CRBN min of range', 'CRBN max of range', '% Present', 'R^2 hyb controls', 'Scale factor'])

	#find all the Report.xmls in the folder containing Skyline output folders looking recursively
	for xmlFile in glob.glob(dataFolder + '/**/Report*.xml'):

		#open the xml file and cut off end of line
		#add each line to dataframe 'report'
		report = []
		with open(xmlFile) as xml:
			for line in xml:
				line = line.rstrip('\n')
				report.append(line)
			xml.close()

		#define each desired value in the report. NB requires that report line numbers have not changed
		pID = re.search(r'<patientid>(\S+)</patientid>', report[4])
		SKY92 = re.search(r'<SKY92>(\S+) *(\S+)*<', report[37])
		#SKY92 (and a lot of the following fields) can be one word, or two words separated with a single whitespace.
		#This loop checks if there is only one word (group(2) == None), and if so sets the output to the first word as the 'None' value can not be treated as a string.
		#If there are two words (the else part of the loop), it joins them with a whitespace in-between.
		if SKY92.group(2) == None:
			SKY92o = SKY92.group(1)
		else:
			SKY92o = (SKY92.group(1) + ' ' + SKY92.group(2))
		#	
		T_4_14 = re.search(r'<T_4_14>(\S+) *(\S+)*<', report[38])
		if T_4_14.group(2) == None:
			T_4_14o = T_4_14.group(1)
		else:
			T_4_14o = (T_4_14.group(1) + ' ' + T_4_14.group(2))
		#	
		T_11_14 = re.search(r'<T_11_14>(\S+) *(\S+)*<', report[39])
		if T_11_14.group(2) == None:
			T_11_14o = T_11_14.group(1)
		else:
			T_11_14o = (T_11_14.group(1) + ' ' + T_11_14.group(2))
		#
		T_14_16_20 = re.search(r'<T_14_16_20>(\S+) *(\S+)*<', report[40])
		if T_14_16_20.group(2) == None:
			T_14_16_20o = T_14_16_20.group(1)
		else:
			T_14_16_20o = (T_14_16_20.group(1) + ' ' + T_14_16_20.group(2))
		#
		ADD1Q = re.search(r'<ADD1Q>(\S+) *(\S+)*<', report[41])
		if ADD1Q.group(2) == None:
			ADD1Qo = ADD1Q.group(1)
		else:
			ADD1Qo = (ADD1Q.group(1) + ' ' + ADD1Q.group(2))
		#
		ADD9 = re.search(r'<ADD9_H_MM>(\S+) *(\S+)*<', report[42])
		if ADD9.group(2) == None:
			ADD9o = ADD9.group(1)
		else:
			ADD9o = (ADD9.group(1) + ' ' + ADD9.group(2))
		#
		DEL13 = re.search(r'<DEL13>(\S+) *(\S+)*<', report[43])
		if DEL13.group(2) == None:
			DEL13o = DEL13.group(1)
		else:
			DEL13o = (DEL13.group(1) + ' ' + DEL13.group(2))
		#
		DEL17 = re.search(r'<DEL17P>(\S+) *(\S+)*<', report[44])
		if DEL17.group(2) == None:
			DEL17o = DEL17.group(1)
		else:
			DEL17o = (DEL17.group(1) + ' ' + DEL17.group(2))
		#	
		MScluster = re.search(r'<MScluster>(\S+) *(\S+)*<', report[45])
		if MScluster.group(2) == None:
			MSclustero = MScluster.group(1)
		else:
			MSclustero = (MScluster.group(1) + ' ' + MScluster.group(2))
		#	
		MFcluster = re.search(r'<MFcluster>(\S+) *(\S+)*<', report[46])
		if MFcluster.group(2) == None:
			MFclustero = MFcluster.group(1)
		else:
			MFclustero = (MFcluster.group(1) + ' ' + MFcluster.group(2))
		#	
		NFKB = re.search(r'<NFKB>(\S+) *(\S+)*<', report[47])
		if NFKB.group(2) == None:
			NFKBo = NFKB.group(1)
		else:
			NFKBo = (NFKB.group(1) + ' ' + NFKB.group(2))
		#	
		CD2 = re.search(r'<CD2>(\S+) *(\S+)*<', report[48])
		if CD2.group(2) == None:
			CD2o = CD2.group(1)
		else:
			CD2o = (CD2.group(1) + ' ' + CD2.group(2))
		#	
		CD38 = re.search(r'<CD38>(\S+)<', report[49])
		CD38MIN = re.search(r'<CD38MIN>(\S+)<', report[51])
		CD38MAX = re.search(r'<CD38MAX>(\S+)<', report[52])
		CD138 = re.search(r'<CD138>(\S+)<', report[53])
		CD138MIN = re.search(r'<CD138MIN>(\S+)<', report[55])
		CD138MAX = re.search(r'<CD138MAX>(\S+)<', report[56])
		CRBN = re.search(r'<CRBN>(\S+)<', report[57])
		CRBNCP = re.search(r'<CRBNCP>(\S+)<', report[58])
		CRBNMIN = re.search(r'<CRBNMIN>(\S+)<', report[59])
		CRBNMAX = re.search(r'<CRBNMAX>(\S+)<', report[60])
		QC_PC_PRESENT = re.search(r'<QC_PC_PRESENT>(\S+) *(\S+)*<', report[29])
		QC_RSQ_HYB = re.search(r'<QC_RSQ_HYB>(\S+)<', report[26])
		QC_SCALE_FACTOR = re.search(r'<QC_SCALE_FACTOR>(\S+)<', report[32])

		#order here is very important too
		table.append([pID.group(1), SKY92o, T_4_14o, T_11_14o, T_14_16_20o, ADD1Qo, ADD9o, DEL13o, DEL17o, MSclustero, MFclustero, NFKBo, CD2o, CD38.group(1), CD38MIN.group(1), CD38MAX.group(1), CD138.group(1), CD138MIN.group(1), CD138MAX.group(1), CRBN.group(1), CRBNCP.group(1), CRBNMIN.group(1), CRBNMAX.group(1), QC_PC_PRESENT.group(1), QC_RSQ_HYB.group(1), QC_SCALE_FACTOR.group(1)])

	#print the output file
	fname = (time.strftime('%Y%m%d') + '_Skyline_report.txt')
	print ('Done! Writing table to ' + fname)

	file = open(fname, 'w')
	for line in table:
		file.write('\t'.join(line)+'\n')
	file.close()
	
if (len(sys.argv) != 2):
	usage()
else:
	parsing(sys.argv[1])

