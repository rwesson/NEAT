#!/bin/bash

prefix=$1

for filename in $prefix*
do 
 step=`more $filename | sort -g | awk ' NR==500 { low=$1 } NR==9500 { high=$1 } END { print (high-low)/30 } '` 
  export step 
  if [ $step != 0 ]; then 
    more $filename | sort -g | awk -v step=$step ' step>0 && NR>50 && NR<9950 { print step*int($1/step) } ' | uniq -c > binned_$filename
    echo "set terminal postscript enhanced" > tmp.prg
    echo "set output \"binned_"$filename".ps" >> tmp.prg
    echo "plot 'binned_"$filename"' using 2:1 w l" >> tmp.prg
    gnuplot tmp.prg
  fi
done
mogrify -density 150 -format jpg -flatten *.ps
rm tmp.prg
