# plot.plt
set xlabel "Time"
set ylabel "Depth"
m="data.txt"
set terminal x11 0
set nokey
set grid
set title "Depth vs. Time"
plot m using 1:2 with linespoints