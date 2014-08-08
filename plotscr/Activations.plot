reset
set terminal postscript enhanced color eps "Helvetica, 12" size 5,1

#set size ratio 0.25
set offsets graph 0, 0, 0.05, 0.05

## Ranges
set xrange [-0.5:24999.5]
set yrange [-0.5:399.5]
set cbrange [0.:1.]

## Tics
set xtics 1000
set ytics 200
set cbtics 1.

## Labels
set xlabel "Time {/Helvetica-Oblique t} [ts]"
set ylabel "Unit index"
set cblabel "State"

#set pal gray
##### Reservoir activations

set output "../plots/res_st.eps"
plot "../results/res_st.mat" matrix with image

set output "../plots/res_pst.eps"
plot "../results/res_pst.mat" matrix with image

set output "../plots/res_dst.eps"
plot "../results/res_dst.mat" matrix with image

set yrange [-0.5:4.5]

## Labels
set xlabel "Time {/Helvetica-Oblique t} [ts]"
set ylabel "Unit index"
set cblabel "State"

set output "../plots/my_data.eps"
plot "../results/my_data.mat" matrix with image
