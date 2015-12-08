#!/bin/bash

list=$(echo *Output/out* |\
  fmt -1 |\
  sort -n -t t -k 5)

rm wave_stddev.dat 2> /dev/null

for i in $list
do
  cp $i ./stability.dat
  python stability_stddev.py >> wave_stddev.dat
done

mv wave_stddev.dat *Output/
rm stability.dat
