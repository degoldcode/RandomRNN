reset
set terminal postscript enhanced color eps "Helvetica, 12" size 5,1

#set size ratio 0.25
set offsets graph 0, 0, 0.1, 0.1

## Ranges
set xrange [-0.5:19999.5]
set yrange [-0.5:599.5]
set cbrange [0.:1.]

## Tics
set xtics 5000
set ytics 200
set cbtics 1.

## Labels
set xlabel "Time {/Helvetica-Oblique t} [ts]"
set ylabel "Unit index"
set cblabel "State"

#set pal gray
##### Reservoir activations

set output "../plots/res_st_train.eps"
plot "../results/res_st_train.mat" matrix with image

set output "../plots/res_pst_train.eps"
plot "../results/res_pst_train.mat" matrix with image

set output "../plots/res_dst_train.eps"
plot "../results/res_dst_train.mat" matrix with image

set xrange [-0.5:4999.5]
set output "../plots/res_st_test.eps"
plot "../results/res_st_test.mat" matrix with image

set output "../plots/res_pst_test.eps"
plot "../results/res_pst_test.mat" matrix with image

set output "../plots/res_dst_test.eps"
plot "../results/res_dst_test.mat" matrix with image

set xtics 1000
set ytics 20
set xrange [-0.5:4999.5]
set yrange [-0.5:99.5]

## Labels
set xlabel "Time {/Helvetica-Oblique t} [ts]"
set ylabel "Unit index"
set cblabel "State"

set output "../plots/input.eps"
plot "../results/my_data.mat" matrix with image
set output "../plots/output.eps"
plot "../results/my_out.mat" matrix with image
set cbrange [-1.:1.]
set output "../plots/dout.eps"
plot "../results/dout.mat" matrix with image
set cbrange [0.:1.]
#set output "../plots/inputdata.eps"
#plot "../input/inputtrain.mat" matrix with image
