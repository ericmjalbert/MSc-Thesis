#!/bin/bash
TEMP=$1
TEMP2=`echo $TEMP | cut -d '.' -f 1` # Takes off the '.txt'
DIRNAME=${TEMP2}ParalOutput  

gfortran -O3 -fopenmp -fdefault-real-8 PDE-ODEsolver.f90 CGdiag.f90 init_random_seed.f90
echo $1 | ./a.out

echo Seperating output files by timestep
./seperateOutput.py

# The 2> /dev/null ignores error messages 
# (since it complains if the directory already exists)
echo Moving datafiles into $DIRNAME
rm -r $DIRNAME
mkdir $DIRNAME 2> /dev/null 
mv total.dat $DIRNAME/
mv statReport.dat $DIRNAME/
cp $1 $DIRNAME/
rm output.dat
mv COprod.dat ./$DIRNAME/
mv travCheck.dat ./$DIRNAME/

find *.dat | wc > temp    # Counts files, lines, and words from input (the pipe "|")

read fileName fileCount wordCount < temp  #saves files lines and words from wc
echo Graphing $fileCount times
COUNTER=1
while [ $COUNTER -lt $(($fileCount+1)) ]    # -lt is actually <, Less Then
do
  n=$((1000+$COUNTER));  n=${n#1}    #adds leading zero by removing the 1 from 1XX 
  
  #Opens gnuplot and << EOF write the commands to gnuplot 
  gnuplot << EOF
  set terminal post enhanced color eps
  set view map
  unset surf
  set pm3d
  set xr[0:1]
  set yr[0:1]
  set zr[0.00001:1]
  set size square
  set xlabel "x"
  set ylabel "y"
  set output "biot$n.eps"
  splot "./out$COUNTER.dat" using 1:2:3 with lines
EOF
#    set view 90,90 is FOR TRAV WAVE VIEW; y-vs.-z
  mv out$COUNTER.dat $DIRNAME/
  mv biot$n.eps $DIRNAME/
  let COUNTER=COUNTER+1   #increment counter
done

#Combines all .eps files into a .pdf
echo Combining .eps files into biomassGraph.pdf
gs -sDEVICE=pdfwrite -dDEVICEWIDTHPOINTS=400 -dDEVICEHEIGHTPOINTS=300 -dFIXEDMEDIAswitch -dNOPAUSE -dBATCH -dSAFER -sOutputFile=$DIRNAME/biomassGraphs.pdf $DIRNAME/biot* > /dev/null
rm temp
rm $DIRNAME/biot*
echo Everything Complete!
