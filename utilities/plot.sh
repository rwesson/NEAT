#!/bin/bash

prefix=$1

for filename in $prefix*binned
do
  echo "set terminal png enhanced" > tmp.prg
  echo "set output \""$filename".png\"" >> tmp.prg
  echo "plot '"$filename"' every :::1::1 w l, '"$filename"' every :::0::0 using 1:(200):2:3 w xerr" >> tmp.prg
  gnuplot tmp.prg
done

rm tmp.prg
