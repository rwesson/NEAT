#!/bin/bash

prefix=$1

for filename in $prefix*dens*
do
  calc=`uniq $filename | wc -l`
  if [ $calc != 1 ]; then 
    step=50
    export step
    more $filename | awk -v step=$step ' step>0 && NR>100 && NR<9900 { print step*int($1/step) } ' | uniq -c > binned_$filename
    echo "set terminal postscript enhanced" > tmp.prg
    echo "set output \"binned_"$filename".ps" >> tmp.prg
    echo "plot 'binned_"$filename"' using 2:1 w l" >> tmp.prg
    gnuplot tmp.prg
  fi
done

for filename in $prefix*temp*
do
  calc=`uniq $filename | wc -l`
  if [ $calc != 1 ]; then 
    step=50
    export step
    more $filename | awk -v step=$step ' step>0 && NR>100 && NR<9900 { print step*int($1/step) } ' | uniq -c > binned_$filename
    echo "set terminal postscript enhanced" > tmp.prg
    echo "set output \"binned_"$filename".ps" >> tmp.prg
    echo "plot 'binned_"$filename"' using 2:1 w l" >> tmp.prg
    gnuplot tmp.prg
  fi
done

for filename in $prefix*abund* $prefix*icf* $prefix*adf* *cHb
do 
 step=`more $filename | sort -g | awk ' NR==500 { low=$1 } NR==9500 { high=$1 } END { print (high-low)/20 } '` 
  export step 
  if [ $step != 0 ]; then 
    more $filename | awk -v step=$step ' step>0 && NR>100 && NR<9900 { print step*int($1/step) } ' | uniq -c > binned_$filename
    echo "set terminal postscript enhanced" > tmp.prg
    echo "set output \"binned_"$filename".ps" >> tmp.prg
    echo "plot 'binned_"$filename"' using 2:1 w l" >> tmp.prg
    gnuplot tmp.prg
  fi
done

mogrify -density 150 -format png -rotate 90 -flatten *.ps
rm tmp.prg *.ps
