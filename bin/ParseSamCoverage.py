#!/usr/bin/env python

from __future__ import division
import csv
import re
import argparse

parser = argparse.ArgumentParser(description='Count Coverage')

parser.add_argument('-i','--inputfile', help='Input txt file (Required)', required=True)
parser.add_argument('-o','--outputfile', help='Output txt file (Required)', required=True)


args = vars(parser.parse_args())

ifile  = open(args['inputfile'], "rb")
reader = csv.reader(ifile, delimiter='\t')
ofile  = open(args['outputfile'], "wb")

meanCoverage = []
percentage10 = []
percentage50 = []
percentage100 = []
percentage500 = []
percentage1000 = []

for row in reader:
	if not row[0].startswith("#"):
		meanCoverage.append(float(row[5]))
		percentage10.append(float(row[6]))
		percentage50.append(float(row[7]))
		percentage100.append(float(row[8]))
		percentage500.append(float(row[9]))
		percentage1000.append(float(row[10]))
		sample = row[11]

meanCov = str(sum(meanCoverage)/float(len(meanCoverage)))
mean10 = str(sum(percentage10)/float(len(percentage10)))
mean50 = str(sum(percentage50)/float(len(percentage50)))
mean100 = str(sum(percentage100)/float(len(percentage100)))
mean500 = str(sum(percentage500)/float(len(percentage500)))
mean1000 = str(sum(percentage1000)/float(len(percentage1000)))

ofile.write(sample+"\t"+meanCov+"\t"+mean10+"\t"+mean50+"\t"+mean100+"\t"+mean500+"\t"+mean1000+"\n")