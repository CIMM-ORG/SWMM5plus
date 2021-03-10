# plot.plt

set xlabel "Time (s)"
set ylabel "Flowrate (m^3/s)"
m="data.txt"
set terminal x11 0
set nokey
set grid
set title "Q vs. T for AC (E1)"
set yrange [0:2]
set term png  
set output "Flowrate_AC.png"  
plot m using 1:2 with linespoints

