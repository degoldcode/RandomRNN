reset
set terminal postscript enhanced color eps "Helvetica, 12"

set size square
set offsets graph 0, 0, 0.1, 0.1

## Ranges
rows = "`cat ../results/init_W.mat | wc -l`"
columns = "`head ../results/init_W.mat -n1 | wc -w`"
set xrange [-0.5:columns+0.5]
set yrange [-0.5:rows+0.5]
#set cbrange [0.:1.]

## Tics
set xtics columns/5
set ytics rows/5
set cbtics 0.2

## Labels
set xlabel "Unit index {/Helvetica-Oblique j}"
set ylabel "Unit index {/Helvetica-Oblique i}"
set cblabel "Weight"

##### Weights

set output "../plots/init_w.eps"
plot "../results/init_W.mat" matrix with image

set output "../plots/end_w.eps"
plot "../results/end_W.mat" matrix with image

set cbrange []
set cbtics autofreq
set output "../plots/delta_w.eps"
plot "../results/delta_W.mat" matrix with image
