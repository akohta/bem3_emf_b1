# Gnuplot script file
# Usage : gnuplot -e 'file="particles_filename"' gscript_particles.plt

if ( !exists("file") ) file="v2.particles"

set terminal postscript eps color enhanced "Arial" 20 size 3in,7in
set output "particles.eps"

#set palette model RGB rgbformulae 35,13,10
set palette model RGB rgbformulae 22,13,-31
set view equal xyz
set view 75,32,1,1
set xtics -0.08, 0.04, 0.08 
set ytics -0.08, 0.04, 0.08
set ticslevel 0
set cbtics 1
set cblabel "Sphere ID"
set grid x y z vertical

set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}" 
set zlabel "{/Arial-Italic z}"

splot file using 1:2:3:4 w d palette title ""
