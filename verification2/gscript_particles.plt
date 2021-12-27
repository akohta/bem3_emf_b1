# Gnuplot script file
# Usage : gnuplot -e 'file="particles_filename"' gscript_particles.plt

if ( !exists("file") ) file="v2.particles"

set terminal postscript eps color enhanced "Arial" 20 size 3in,7in
set output "particles.eps"

set view equal xyz
set view 75,32,1,1
set xtics -0.08, 0.04, 0.08 
set ytics -0.08, 0.04, 0.08
set ticslevel 0
unset colorbox
set grid x y z vertical
set pointsize 0.3

set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}" 
set zlabel "{/Arial-Italic z}"

splot file using 1:2:3:4 pt 7 lc 'green' title ""
