 set terminal postscript eps enhanced color font "Helvetica,10"
 set output "centrodemassa.eps"
 set xlabel "tempo(s)"
 set ylabel "Velocidade(m/s)"
 m = "saida.dat"
 set encoding utf8
 set terminal x11 0
 set title "Gráfico Velocidade centro de massa"
 plot "saida.dat" u 1:2 with lines lc rgb "red", \
 	"saida.dat" u 1:3 with lines lc rgb "green" ,\
