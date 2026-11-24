set terminal pdfcairo enhanced color size 6.8,3.2 font "Helvetica,8.5"
set output "fig_ef_vs_3fdd_tpfa.pdf"
set multiplot layout 1,2

set logscale xy
set format x "%g"
set format y "10^{%T}"
set grid
set xlabel "Nombre de nœuds N"
set ylabel "Erreur relative L^2"

set title "Ex. 2 (K discontinu x=1) : 4 méthodes sur domaine en L"
set key bottom left font ",7"
plot "ex2_history.dat"                       using 2:4 with linespoints pt 4 lw 1.8 lc rgb "#0077BB" title "EF P1 (SPD)", \
     "../04_VF_TPFA/ex2_history.dat"         using 2:4 with linespoints pt 6 lw 1.4 lc rgb "#888888" title "VF TPFA (stagne)", \
     "../01_3FDD_Rectangles/ex2_history.dat" using 3:5 with linespoints pt 5 lw 1.7 lc rgb "#009988" title "3FDD-Rectangles", \
     "../02_3FDD_Voronoi/ex2_history.dat"    using 3:5 with linespoints pt 7 lw 1.7 lc rgb "#882255" title "3FDD-Voronoï"

set title "Ex. 3 (Aniso 50:1, coin) : 4 méthodes sur domaine en L"
set key bottom left font ",7"
plot "ex3_history.dat"                       using 2:4 with linespoints pt 4 lw 1.8 lc rgb "#0077BB" title "EF P1 (SPD)", \
     "../04_VF_TPFA/ex3_history.dat"         using 2:4 with linespoints pt 6 lw 1.4 lc rgb "#888888" title "VF TPFA (stagne)", \
     "../01_3FDD_Rectangles/ex3_history.dat" using 3:5 with linespoints pt 5 lw 1.7 lc rgb "#009988" title "3FDD-Rectangles", \
     "../02_3FDD_Voronoi/ex3_history.dat"    using 3:5 with linespoints pt 7 lw 1.7 lc rgb "#882255" title "3FDD-Voronoï"

unset multiplot
