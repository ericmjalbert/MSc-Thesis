#!/bin/python

import sys
from itertools import izip

if len(sys.argv) != 3: 
    print "[!] arg format: > calcSigmat.py input1 input2"
    exit()
fp1 = open(str(sys.argv[2]), "r")
fp2 = open(str(sys.argv[2]), "r")

sumVal = 0
for line1 in fp1.readlines():
    line2 = fp2.readline()
    if line1.strip():
        line2 = fp2.readline()
        print line1
#        sumVal += float(line1.split()[1]) - float(line2.split()[1])

fp2.close()
print str(sys.argv[1]) + " " + str(sumVal)

