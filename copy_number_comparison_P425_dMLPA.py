#!/usr/bin/python

####
# This script parses a digital MLPA file and an array output and compares values for multiple samples
# All files must have digital Probe names in a column header "Probes"
# NB. YOU CAN'T HAVE plus/minus symbols in the files, it breaks the csv reader.
# NB. There are genes where there are more probes in the dMLPA than the P425.
# In these instances, the equivalent P425 probe was duplicated so that there was data available for each dMLPA probe:
#    ± CDKN2C probe 14652-L16304  01-051.212271	1p32.3	T-CDKN2C-1(S011456)-
#	 ± CDKN2C probe 14652-L16304  01-051.212271	1p32.3	T-CDKN2C-2(S013428)-
#	 ± CDKN2C probe 14652-L16304  01-051.212271	1p32.3	T-CDKN2C-3(S011457)-
#	FAM46C probe 18949-L24912  01-117.966975	1p12	T-FAM46C-01B(S011448)-
#	FAM46C probe 18949-L24912  01-117.966975	1p12	T-FAM46C-02C(S011386)-
#	RP11-541J2 probe 16229-L18482  01-162.573278	1q23.3	T-RP11-541J2-541J2(S011465)-
#	RP11-541J2 probe 16229-L18482  01-162.573278	1q23.3	T-RP11-541J2-541J2(S011464)-
#	RP11-541J2 probe 16229-L18482  01-162.573278	1q23.3	T-RP11-541J2-541J2(S011466)-
#	NCAPD2 probe 16235-L18488  12-006.510814	12p13.31	T-NCAPD2-5(S011429)-
#	NCAPD2 probe 16235-L18488  12-006.510814	12p13.31	T-NCAPD2-23(S011427)-
#	TP53 probe 02376-L24733  17-007.519217	17p13.1	T-TP53-4b(S010584)-
#	TP53 probe 02376-L24733  17-007.519217	17p13.1	T-TP53-4b(S010585)-
# NB. There is no equivalent to the PPAP2B/PLPP3 probe in the dMLPA (PPAP2B probe 18798-L24446  01-056.775291).
####

import sys
import csv
import math
import pandas as pd
import numpy as np

def usage():
	print ('Usage:')
	print ('copy_number_comparison.py data-one-name path/to/data-one.csv data-two-name path/to/data-two.csv')
	sys.exit()
	
def parsing(data1name, path1, data2name, path2):
    #Use pandas to import tables and set probe names as row names
    #Maintain 'data1/2name' as text variables so it can be called on later
    d1 = data1name
    d2 = data2name
    #d3 = data3name
    d1 = pd.read_csv(path1)
    d2 = pd.read_csv(path2)
    #d3 = pd.read_csv(path3)
    d1 = d1.set_index('Probes')
    d2 = d2.set_index('Probes')
    #d3 = d3.set_index('Probes')

    #Take forward only samples that are in all three data sets
    d1samp = list(d1.columns.values)
    d2samp = list(d2.columns.values)
    #d3samp = list(d3.columns.values)    
    samples = set(d1samp) & set(d2samp)
    print('There are ' + str(len(d1samp)) + ' samples in the ' + data1name + ' dataset, '\
    + str(len(d2samp)) + ' samples in the ' + data2name +\
    ' dataset, giving an overlap of ' + str(len(samples)) + ' samples in total.')
    
    #Take forward only probes that are in both data sets
    d1probes = list(d1.index.values)
    d2probes = list(d2.index.values)
    #d3probes = list(d3.index.values)
    probes = set(d1probes) & set(d2probes)
    print('There are ' + str(len(d1probes)) + ' probes in the ' + data1name + ' dataset, '\
    + str(len(d2probes)) + ' probes in the ' + data2name +\
    ' dataset, giving an overlap of ' + str(len(probes)) + ' probes in total.')

    d1sub = pd.DataFrame()
    for i in list(probes):
        data = d1.loc[i, :]
        d1sub = d1sub.append(data)

    d2sub = pd.DataFrame()
    for i in list(probes):
        data = d2.loc[i, :]
        d2sub = d2sub.append(data)

    #d3sub = pd.DataFrame()
    #for i in list(probes):
    #    data = d3.loc[i, :]
    #    d3sub = d3sub.append(data)

    d1sub = d1sub[list(samples)]
    d2sub = d2sub[list(samples)]
    #d3sub = d3sub[list(samples)]
        
    #Stack each table with multi-hierarchical indexing and raw data for printing later
    d1subsplit = d1sub.stack()
    d2subsplit = d2sub.stack()
    #d3subsplit = d3sub.stack()

    #Prompt user for the values of loss and gain in each data set
    d1loss = pd.DataFrame()
    d1lossvalue = float(input('Enter value below which the ' + data1name + ' is considered a loss: '))
    for i in list(probes):
        data = d1sub.loc[i] < d1lossvalue
        d1loss = d1loss.append(data)
    #Set index and column names here - these will be used when creating the final data frame
    d1loss.index.name = 'Probes'
    d1loss.columns.name = 'Samples'

    d1gain = pd.DataFrame()
    d1gainvalue = float(input('Enter value above which the ' + data1name + ' is considered a gain: '))
    for i in list(probes):
        data = d1sub.loc[i] > d1gainvalue
        d1gain = d1gain.append(data)
    
    d2loss = pd.DataFrame()
    d2lossvalue = float(input('Enter value below which the ' + data2name + ' is considered a loss: '))
    for i in list(probes):
        data = d2sub.loc[i] < d2lossvalue
        d2loss = d2loss.append(data)

    d2gain = pd.DataFrame()
    d2gainvalue = float(input('Enter value above which the ' + data2name + ' is considered a gain: '))
    for i in list(probes):
        data = d2sub.loc[i] > d2gainvalue
        d2gain = d2gain.append(data)
    
    #d3loss = pd.DataFrame()
    #d3lossvalue = float(input('Enter value below which the ' + data3name + ' is considered a loss: '))
    #for i in list(probes):
    #    data = d3sub.loc[i] < d3lossvalue
    #    d3loss = d3loss.append(data)

    #d3gain = pd.DataFrame()
    #d3gainvalue = float(input('Enter value above which the ' + data3name + ' is considered a gain: '))
    #for i in list(probes):
    #    data = d3sub.loc[i] > d3gainvalue
    #    d3gain = d3gain.append(data)
        
    #Set names for the table headers
    d1lossname = data1name + ' loss - ' + str(d1lossvalue)
    d1gainname = data1name + ' gain - ' + str(d1gainvalue)
    d2lossname = data2name + ' loss - ' + str(d2lossvalue)
    d2gainname = data2name + ' gain - ' + str(d2gainvalue)
    #d3lossname = data3name + ' loss - ' + str(d3lossvalue)
    #d3gainname = data3name + ' gain - ' + str(d3gainvalue)
    
    #Write to file 
    d1subsplit.to_csv(data1name + '_loss-' + str(d1lossvalue) + '_gain-' + str(d1gainvalue) + '_subgroupraw.csv')
    d2subsplit.to_csv(data2name + '_loss-' + str(d2lossvalue) + '_gain-' + str(d2gainvalue) + '_subgroupraw.csv')
    #d3subsplit.to_csv(data3name + 'subgroupraw.csv')

    #Make a table will all the binary data
    binary = pd.concat([d1loss.stack(), d1gain.stack(), d2loss.stack(), d2gain.stack()],\
    keys = [d1lossname, d1gainname, d2lossname, d2gainname],\
    axis = 1, join_axes = [d1loss.stack().index])
        
    #Code when the two data types agree and disagree
    binary['code1'] = np.where((binary[d1lossname] == 0) & (binary[d1gainname] == 0) & (binary[d2lossname] == 0) & (binary[d2gainname] == 0), 'norm', '')
    binary['code2'] = np.where((binary[d1lossname] == 1) & (binary[d1gainname] == 0) & (binary[d2lossname] == 0) & (binary[d2gainname] == 0), data1name + ' loss', '')
    binary['code3'] = np.where((binary[d1lossname] == 0) & (binary[d1gainname] == 1) & (binary[d2lossname] == 0) & (binary[d2gainname] == 0), data1name + ' gain', '')
    binary['code4'] = np.where((binary[d1lossname] == 0) & (binary[d1gainname] == 0) & (binary[d2lossname] == 1) & (binary[d2gainname] == 0), data2name + ' loss', '')
    binary['code5'] = np.where((binary[d1lossname] == 0) & (binary[d1gainname] == 0) & (binary[d2lossname] == 0) & (binary[d2gainname] == 1), data2name + ' gain', '')
    binary['code6'] = np.where((binary[d1lossname] == 1) & (binary[d1gainname] == 0) & (binary[d2lossname] == 1) & (binary[d2gainname] == 0), 'both loss', '')
    binary['code7'] = np.where((binary[d1lossname] == 0) & (binary[d1gainname] == 1) & (binary[d2lossname] == 0) & (binary[d2gainname] == 1), 'both gain', '')
    binary['code8'] = np.where((binary[d1lossname] == 1) & (binary[d1gainname] == 0) & (binary[d2lossname] == 0) & (binary[d2gainname] == 1), data1name + ' loss, ' + data2name + ' gain', '')
    binary['code9'] = np.where((binary[d1lossname] == 0) & (binary[d1gainname] == 1) & (binary[d2lossname] == 1) & (binary[d2gainname] == 0), data1name + ' gain, ' + data2name + ' loss', '')

    
    #Concat each column into one and delete
    binary['code'] = binary[['code1', 'code2', 'code3', 'code4', 'code5', 'code6', 'code7', 'code8', 'code9']].fillna('').sum(axis=1)
    binary.drop(['code1', 'code2', 'code3', 'code4', 'code5', 'code6', 'code7', 'code8', 'code9'], axis=1, inplace=True)
    
    #Write to file 
    binary.to_csv(data1name + '_' + data2name + '_binary.csv')


if (len(sys.argv) != 5):
	usage()
else:
	parsing(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
