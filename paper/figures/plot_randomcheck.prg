set xlabel "Flux (arbitrary units)"
set ylabel "Frequency"

set terminal png enhanced size 960,720
set output "gaussian_test.png"

e=2.718281828
f(x) = n*e**-(((x-mu)**2)/(2*sigma**2)) #gaussian distribution
n=3284.
mu=303.2
sigma=6.064

fit f(x) './ha_6543_test_sorted' using 2:1 via n,mu,sigma

plot './ha_6543_test_sorted' using 2:1 notitle,f(x) lw 2 lt -1 notitle
set terminal wxt
