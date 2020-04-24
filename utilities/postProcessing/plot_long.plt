# plot_long.plt
set xlabel "Length (m)"
set ylabel "Depth (m)"
m = "data.txt"
n = "long_profile.txt"
o = "weir.txt"
p = "weir2.txt"
set terminal x11 0
set nokey
set grid
set title "Time: 0s"
set yrange [0:6]
set xrange [0:2000]
set terminal gif animate delay 8
set output "test_orif_eta.gif"
i = 1
j = 19
plot n using 1:2 with lines linecolor rgb "black" linewidth 2
replot o using 1:2 with lines linecolor rgb "red" linewidth 2
replot p using 1:2 with lines linecolor rgb "red" linewidth 2
do for [k = 1:301]{
	replot m every ::i::j using 1:2 with lines linecolor rgb "blue" linewidth 1
	i = i + 19
	j = j + 19
	set title sprintf("Time: %f s",k*100)
}

 