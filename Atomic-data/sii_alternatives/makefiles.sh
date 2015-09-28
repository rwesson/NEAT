# 8 levels, number of Te points depends on file

for coll in RBS96 TZ10
do

  for A in PKW09 TZ10-PKW09 VVF96-KHOC93 VVF96-M82a
  do

    tpoints=`head -n 2 s_ii_coll_${coll}.dat | tail -n 1 | awk ' { print NF-2 } '`
    nlevels=8

    echo "2" > sii_${coll}_${A}.dat
    echo "%collision strengths: $coll" >> sii_${coll}_${A}.dat
    echo "%transition probabilities: $A" >> sii_${coll}_${A}.dat
    echo $nlevels $tpoints >> sii_${coll}_${A}.dat
    echo " 1  3s2.3p3 4S3/2" >> sii_${coll}_${A}.dat
    echo " 2  3s2.3p3 2D3/2" >> sii_${coll}_${A}.dat
    echo " 3  3s2.3p3 2D5/2" >> sii_${coll}_${A}.dat
    echo " 4  3s2.3p3 2P1/2" >> sii_${coll}_${A}.dat
    echo " 5  3s2.3p3 2P3/2" >> sii_${coll}_${A}.dat
    echo " 6  3s2.3p4 4P5/2" >> sii_${coll}_${A}.dat
    echo " 7  3s2.3p4 4P3/2" >> sii_${coll}_${A}.dat
    echo " 8  3s2.3p4 4P1/2" >> sii_${coll}_${A}.dat

    head -n 2 s_ii_coll_${coll}.dat | tail -n 1 | awk ' { print substr($0,5,999) } ' | awk ' gsub(" ","\n") { print $0 } ' | awk ' { print 10**$1 } ' >> sii_${coll}_${A}.dat

    echo "0" >> sii_${coll}_${A}.dat

    #echo level numbers and collision strengths for each temperature

    head -n 29 s_ii_coll_${coll}.dat | tail -n +3 | awk ' { print $1,$2,$3 } { for (i=4; i<=NF; i++) print "0 0 "$i } ' >> sii_${coll}_${A}.dat

    echo "0 0 0" >> sii_${coll}_${A}.dat

    #transition probabilities

    for i in {2..8}; do
      for j in `seq $i 8`; do
        head -n 10 s_ii_atom_${A}.dat | tail -n +3 | awk -v i=$i -v j=$j ' NR==j { print i-1,j,$(i-1) } ' >> sii_${coll}_${A}.dat
      done
    done

    #wavenumbers and statistical weights

    echo "  1  4        0.0" >> sii_${coll}_${A}.dat
    echo "  2  4    14852.9" >> sii_${coll}_${A}.dat
    echo "  3  6    14884.7" >> sii_${coll}_${A}.dat
    echo "  4  2    24524.8" >> sii_${coll}_${A}.dat
    echo "  5  4    24571.5" >> sii_${coll}_${A}.dat
    echo "  6  6    79395.4" >> sii_${coll}_${A}.dat
    echo "  7  4    79756.8" >> sii_${coll}_${A}.dat
    echo "  8  2    79962.6" >> sii_${coll}_${A}.dat

  done

done
