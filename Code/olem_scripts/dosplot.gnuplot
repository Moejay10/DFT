
reset
set terminal X11 enhanced font "arial,18" 
set output
set xlabel "E - E_{Fermi} (eV)"
set ylabel "DOS (states/eV)"

plot 'DOS0' using 1:2 with lines title "Total DOS"
pause -1

set term post eps color enhanced "Arial" 18
set output "Total_DOS.eps"
replot

#set term png
#set output "Total_DOS.png"
#replot
