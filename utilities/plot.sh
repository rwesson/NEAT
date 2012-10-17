#!/bin/bash

prefix=$1

for filename in $prefix*binned
do
  echo "set terminal png enhanced" > tmp.prg
  echo "set output \""$filename".png\"" >> tmp.prg
#  echo "f(x)=a*exp((-(x-mu)**2)/(2*sigma**2))" >> tmp.prg
#  head -n 1 $filename | awk ' { print "mu="$1" ; sigma="($3-$2)/2 } ' >> tmp.prg
#  tail -n +3 $filename | sort -k 2 -g | tail -n 1 | awk ' { print "a="$2 } ' >> tmp.prg 
  echo "plot '"$filename"' every :::1::1 w l, '"$filename"' every :::0::0 using 1:(800):2:3 w xerr" >> tmp.prg
  gnuplot tmp.prg
done

rm tmp.prg
