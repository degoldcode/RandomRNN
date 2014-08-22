reset
set terminal postscript enhanced color eps "Helvetica, 12" size 5,1
START = 0
#set size ratio 0.25
set offsets graph 0, 0, 0.1, 0.1

## Ranges
START = 0
END = 5000
set xrange [-0.5+START:-0.5+END+START]
#set yrange [-0.5:599.5]
#set cbrange [0.:1.]
set key outside

## Tics
#set xtics 2000
#set ytics 100
#set cbtics 1.

set output "../plots/1dplot.eps"
plot "../results/my_data.mat" u 0:1 w l lt 1 lc rgb "red"  t "Input", "../results/my_data_test.mat" u 0:1 w l lt 1 lc rgb "violet"  t "Input test", "../results/my_teacher.mat" u 0:1 w l lt 1 lc rgb "green"  t "Teacher" 
set output
set output "../plots/1dplotout.eps"
plot "../results/my_out_test.mat" u 0:1 w l lt 1 lc rgb "green"  t "Output test", "../results/my_out.mat" u 0:1 w l lt 1 lc rgb "blue"  t "Output"
#, "../results/my_data_test.mat" u 0:1 w l lt 1 lc rgb "violet"  t "Input test"
set output
