#! /usr/bin/env python

f = open('output.dat','r')
outName = "out1.dat"
out = open(outName, 'w')
count = 1
i = 0

lines = f.readlines()
f.close()


for line in lines: 
  if line == '\n':
    count += 1
    out.close()
    outName = "out" + str(count) + ".dat"
    out = open(outName,'w')
  out.write(line)
  #print '{0}, {1}'.format(count,line)
