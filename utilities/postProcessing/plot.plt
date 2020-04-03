# plot.plt

set xlabel "Time (s)"
set ylabel "Depth (m)"
m="data.txt"
set terminal x11 0
set nokey
set grid
#set title "Downstream Boundary Condition"
#set yrange [0:1]
#set term png  
#set output "Downstream Boundary Condition_flowrate.png"  
plot m using 1:2 with linespoints

