#! /usr/bin/env python

import numpy

f = open('stability.dat', 'r')

lines = f.readlines()
f.close()

wave_y = []
check_x = []

lines = lines[::-1][1:-1] # Reverse list and ignore starting and ending '\n'
for line in lines:
    if len(line.split()) <= 1: # skip newline lines
        continue
    if float(line.split()[2]) >= 0.001 and float(line.split()[0]) not in check_x:
        wave_y.append(float(line.split()[1])) # Append y location
        check_x.append(float(line.split()[0])) # Append x location so you don't count it twice
        
if len(wave_y) in [33, 65, 129, 257]:
    print numpy.asarray(wave_y).mean(), numpy.asarray(wave_y).std(ddof = 1)
else:
    print numpy.asarray(wave_y).sum()/257, -1

