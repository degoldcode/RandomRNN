g++ main.cpp SORN.cpp -std=c++11 -o start -O1 -larmadillo
./start
cd plotscr
gnuplot *.plot
