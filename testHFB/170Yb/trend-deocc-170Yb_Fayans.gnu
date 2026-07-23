#set terminal postscript eps enhanced color "TimesNewRoman,30"
set terminal epslatex color colortext font "Helvetica,14" solid round standalone 
set colorsequence classic
#set output "trend-deocc-170Yb_Fayans.eps"
set output "trend-deocc-170Yb_Fayans.tex"
set border 15 lw 1
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set size 1.03,1.55
set style line 1 lt 1 lc rgb "#FF0000" lw 4 pt 7 ps 0.8
set style line 2 lt 1 lc rgb "#00AA00" lw 4 pt 7 ps 0.8
set style line 3 lt 1 lc rgb "#0000FF" lw 4 pt 7 ps 0.8
set style line 4 lt 1 lc rgb "#AA00AA" lw 4 pt 7 ps 0.8
set style line 5 lt 1 lc rgb "#009977" lw 4 pt 7 ps 0.8
set style line 6 lt 1 lc rgb "#000000" lw 4
set style line 7 lt 0 lc rgb "#FF0000" lw 3
set style line 8 lt 0 lc rgb "#000000" lw 3
set multiplot
#
set size 0.8,0.6
set origin 0.2,0.2
set xrange [6.5:19]
set xlabel "$N_{\\mathrm{states}}^{1/3}$" offset 0
set yrange [-1410:-1373]
set ylabel "total energy  (MeV)" offset 1
set form x 
set key top left
p  "trend-deocc-170Yb_Fayans.dat"  i 2 u ($3**0.3333):4 t "$\\varepsilon_{cut}$ = 20 MeV" w lp ls 1,\
  "trend-deocc-170Yb_Fayans.dat"  i 3 u ($3**0.3333):4 t "$\\varepsilon_{cut}$ = 40 MeV" w lp ls 2,\
  "trend-deocc-170Yb_Fayans.dat"  i 4 u ($3**0.3333):4 t "$\\varepsilon_{cut}$ = 60 MeV" w lp ls 3
#
set size 0.8,0.6
set origin 0.2,0.8
set yrange [5.26:5.55]
set ylabel "charge radius  (fm)" offset 0
set xlabel "" offset 0
set form x ""
set key bottom right
set title "$^{170}$Yb, FyDr (Fayans), HFB" offset 0,-0.3
p  "trend-deocc-170Yb_Fayans.dat"  i 2 u ($3**0.3333):12 t "$\\varepsilon_{cut}$ = 20 MeV" w lp ls 1,\
  "trend-deocc-170Yb_Fayans.dat"  i 3 u ($3**0.3333):12 t "$\\varepsilon_{cut}$ = 40 MeV" w lp ls 2,\
  "trend-deocc-170Yb_Fayans.dat"  i 4 u ($3**0.3333):12 t "$\\varepsilon_{cut}$ = 60 MeV" w lp ls 3
#
#!eps2pdf -f trend-deocc-170Yb_Fayans.eps
#!rm -f trend-deocc-170Yb_Fayans.eps
#!evince trend-deocc-170Yb_Fayans.pdf
#!open -a preview trend-deocc-170Yb_Fayans.pdf
q
