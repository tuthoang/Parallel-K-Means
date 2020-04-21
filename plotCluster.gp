set key left bottom left
set key font ",8"
set key spacing 1

set style increment user
set style line 1 lc rgb 'red'
set style line 2 lc rgb 'green'
set style line 3 lc rgb 'blue'
set style line 4 lc rgb 'yellow'
set style line 5 lc rgb 'purple'
set style line 6 lc rgb 'orange'
set style line 7 lc rgb 'green'
set style line 8 lc rgb 'brown'
set style data points


# Plot centroids, then plot border to make it pop out more
set term png
set output "clusters.png"
plot 'data_colors.txt' using 1:2:3 linecolor variable pt 7 ps 2 t ''
replot for [i=0:8] 'centroids.txt' every ::i::i using 1:2:3 linecolor variable pt 7 ps 2 t 'Cluster '.i,\
        for [i=0:8] 'centroids.txt' every ::i::i using 1:2 pt 8 ps 2 lc rgb "black" lw 1 t ''

exit