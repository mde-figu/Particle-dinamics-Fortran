 set terminal postscript eps enhanced color font "Helvetica,10"
 set output "potencial.eps"
 set xlabel "tempo(s)"
 set ylabel "potencial"
 m = "saida.dat"
 set encoding utf8
 set terminal x11 0
 set title "Gráfico Potencial"
 plot "saida.dat" u 1:2 with lines lc rgb "red", \
 	"saida.dat" u 1:3 with lines lc rgb "green" ,\
