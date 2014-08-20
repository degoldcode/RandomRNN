g++ main.cpp SORN.cpp ESN.cpp -std=c++11 -o start -O1 -larmadillo
./start
cd plotscr
gnuplot *.plot
