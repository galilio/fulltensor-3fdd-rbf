# 1. Domaine non convexe en L : initial (2 rectangles) vs Ex 2 final vs Ex 3 final
set terminal pdfcairo enhanced color size 7.4,2.8 font "Helvetica,8"
set output "fig1_lshape_adapt.pdf"
set multiplot layout 1,3 title "Domaine non convexe en L : décomposition initiale (M=2) et raffinement adaptatif (rectangles + points ARS)"

set xrange [-0.05:2.05]
set yrange [-0.05:2.05]
set size ratio -1
set xlabel "x"
set ylabel "y"
unset key

set title "Initial : M = 2 rectangles, N_u = 102"
plot "ex3_init_rects.dat" with lines lc rgb "#111111" lw 1.4, \
     "ex3_init_nodes.dat" using 1:2 with points pt 7 ps 0.35 lc rgb "#CC3311"

set title "Ex. 2 (K discontinu) : M = 24, N_u = 929"
plot "ex2_final_rects.dat" with lines lc rgb "#111111" lw 0.9, \
     "ex2_final_nodes.dat" using 1:2 with points pt 7 ps 0.25 lc rgb "#0077BB"

set title "Ex. 3 (Aniso 50:1, coin) : M = 36, N_u = 1634"
plot "ex3_final_rects.dat" with lines lc rgb "#111111" lw 0.85, \
     "ex3_final_nodes.dat" using 1:2 with points pt 7 ps 0.22 lc rgb "#008866"

unset multiplot

# 2. Courbes de convergence (e_L2 et e_ener vs N_u) pour les 3 exemples
set terminal pdfcairo enhanced color size 6.8,3.2 font "Helvetica,8.5"
set output "fig2_convergence_3ex.pdf"
set multiplot layout 1,2

set logscale xy
set format x "%g"
set format y "10^{%T}"
set grid
set xlabel "Nombre total de points N_u"

set title "Erreur relative L^2 (3 exemples)"
set ylabel "e_{L^2}"
set key bottom left font ",7.5"
plot "ex1_history.dat" using 3:5 with linespoints pt 4 lw 1.7 lc rgb "#0077BB" title "Ex. 1 : K plein constant", \
     "ex2_history.dat" using 3:5 with linespoints pt 7 lw 1.7 lc rgb "#009988" title "Ex. 2 : K discontinu (x=1)", \
     "ex3_history.dat" using 3:5 with linespoints pt 9 lw 1.9 lc rgb "#CC3311" title "Ex. 3 : Aniso 50:1 + coin (1,1)"

set title "Erreur relative énergie (3 exemples)"
set ylabel "e_{énergie}"
set key bottom left font ",7.5"
plot "ex1_history.dat" using 3:7 with linespoints pt 4 lw 1.7 lc rgb "#0077BB" title "Ex. 1 : K plein constant", \
     "ex2_history.dat" using 3:7 with linespoints pt 7 lw 1.7 lc rgb "#009988" title "Ex. 2 : K discontinu (x=1)", \
     "ex3_history.dat" using 3:7 with linespoints pt 9 lw 1.9 lc rgb "#CC3311" title "Ex. 3 : Aniso 50:1 + coin (1,1)"

unset multiplot
