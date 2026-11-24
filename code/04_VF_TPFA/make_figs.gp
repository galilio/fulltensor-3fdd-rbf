set terminal pdfcairo enhanced color size 6.8,3.2 font "Helvetica,8.5"
set output "fig_tpfa_echec_vs_3fdd.pdf"
set multiplot layout 1,2

set logscale xy
set format x "%g"
set format y "10^{%T}"
set grid
set xlabel "Nombre de volumes / nœuds N"
set ylabel "Erreur relative L^2"

set title "VF TPFA : isotrope (Ex. 0) vs tenseurs pleins (Ex. 1–3)"
set key bottom left font ",7"
plot "ex0_history.dat" using 2:4 with linespoints pt 4 lw 1.9 lc rgb "#009988" title "Ex. 0 : K = I (isotrope, converge)", \
     "ex1_history.dat" using 2:4 with linespoints pt 6 lw 1.6 lc rgb "#0077BB" title "Ex. 1 : K plein const. (stagne)", \
     "ex2_history.dat" using 2:4 with linespoints pt 8 lw 1.6 lc rgb "#EE7733" title "Ex. 2 : K discontinu (stagne)", \
     "ex3_history.dat" using 2:4 with linespoints pt 9 lw 1.7 lc rgb "#CC3311" title "Ex. 3 : Aniso 50:1 (stagne)"

set title "Ex. 3 (Aniso 50:1, coin) : VF TPFA vs 3FDD (Rect & Voronoï)"
set key bottom left font ",7"
plot "ex3_history.dat"                       using 2:4 with linespoints pt 6 lw 1.6 lc rgb "#888888" title "VF TPFA (non consistant)", \
     "../01_3FDD_Rectangles/ex3_history.dat" using 3:5 with linespoints pt 4 lw 1.8 lc rgb "#0077BB" title "3FDD-Rectangles (Méthode 1)", \
     "../02_3FDD_Voronoi/ex3_history.dat"    using 3:5 with linespoints pt 7 lw 1.8 lc rgb "#882255" title "3FDD-Voronoï (Méthode 2)"

unset multiplot
