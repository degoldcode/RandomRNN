reset
set terminal postscript enhanced color eps "Helvetica, 12" size 5,1
START = 0
#set size ratio 0.25
set offsets graph 0, 0, 0.1, 0.1

## Ranges
set xrange [-0.5:19999.5]
set yrange [-0.5:599.5]
set cbrange [0.:1.]

## Tics
set xtics 2000
set ytics 100
set cbtics 1.

## Labels
set xlabel "Time {/Helvetica-Oblique t} [ts]"
set ylabel "Unit index"
set cblabel "State"

#set pal gray
##### Reservoir activations
rows = "`cat ../results/res_st_train.mat | wc -l`"
columns = "`head ../results/res_st_train.mat -n1 | wc -w`"
set xrange [-0.5+START:columns+0.5]
set yrange [-0.5:rows+0.5]
set xtics autofreq
set output "../plots/res_st_train.eps"
plot "../results/res_st_train.mat" matrix with image

rows = "`cat ../results/res_st_test.mat | wc -l`"
columns = "`head ../results/res_st_test.mat -n1 | wc -w`"
set xrange [-0.5+START:columns+0.5]
set yrange [-0.5:rows+0.5]
set output "../plots/res_st_test.eps"
plot "../results/res_st_test.mat" matrix with image

set output "../plots/res_pst_test.eps"
plot "../results/res_pst_test.mat" matrix with image

rows = "`cat ../results/my_data.mat | wc -l`"
columns = "`head ../results/my_data.mat -n1 | wc -w`"
set xrange [-0.5+START:columns-0.5]
set yrange [-0.5:rows-0.5]
set ytics rows/5
## Labels
set xlabel "Time {/Helvetica-Oblique t} [ts]"
set ylabel "Unit index"
set cblabel "State"

set output "../plots/input.eps"
plot "../results/my_data.mat" matrix with image
set output "../plots/inputtest.eps"
plot "../results/my_data_test.mat" matrix with image


rows = "`cat ../results/my_out.mat | wc -l`"
columns = "`head ../results/my_out.mat -n1 | wc -w`"
set xrange [-0.5+START:columns-0.5]
set yrange [-0.5:rows-0.5]

set output "../plots/output.eps"
plot "../results/my_out.mat" matrix with image
set output "../plots/outputtest.eps"
plot "../results/my_out_test.mat" matrix with image
set output "../plots/teacher.eps"
plot "../results/my_teacher.mat" matrix with image
set cbrange []
#set xtics out
set palette gray negative
#set xtics 250,500,5500
set output "../plots/dout.eps"
plot "../results/dout.mat" matrix with image
