# plot.plt
<<<<<<< HEAD
set xlabel "Time (s)"
set ylabel "Depth (m)"
=======
set xlabel "Time"
set ylabel "Depth"
>>>>>>> master
m="data.txt"
set terminal x11 0
set nokey
set grid
<<<<<<< HEAD
set title "Downstream Boundary Condition"
#set yrange [0:1]
set term png  
set output "Downstream Boundary Condition_flowrate.png"  
plot m using 1:2 with linespoints

=======
set title "Depth vs. Time"
plot m using 1:2 with linespoints
>>>>>>> master
