
reset
set terminal X11 enhanced
set style line 1 lt 1 lw 1.5
set style line 2 lt 2 lw 1.5
set style line 3 lt 3 lw 1.5

set xlabel "E - E_{Fermi} (eV)"
set ylabel "DOS (states/atom*eV)"

plot 'DOS$1' using 1:2 with lines ls 1 title "LDOS $1 s",\
     'DOS$1' using 1:3 with lines ls 2 title "LDOS $1 p",\
     'DOS$1' using 1:4 with lines ls 3 title "LDOS $1 d"
pause -1

set term post eps color enhanced
set output "LDOS$1.eps"
replot

set term png
set output "LDOS$1.png"
replot

