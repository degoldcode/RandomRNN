/*
 * main.cpp
 *
 *  Created on: Aug 6, 2014
 *      Author: degoldschmidt
 */

#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include "SORN.h"
#include "timer.h"
using namespace std;
using namespace arma;

const int trainLen = 22000;				// time interval for training (in ts)
const int testLen = 5500;				// time interval for testing (in ts)
const int initLen = 0;					// time interval for initiating (in ts)
const int missingLen = 50;				// time interval of missing data (in ts)

const int inSize = 100; 				// Number of input neurons (input dim)
const int outSize = 100;				// Number of input neurons (input dim)
const int resSize = 500;				// Number of reservoir neurons (resSize = 0.8*Ne + 0.2*Ni)
const double a = 1.0; 					// leaking rate

SORN * my_sorn;
mat my_data;
mat my_data_test;
mat my_output;
mat my_output_test;

// if data needs to be written to file
ofstream my_write;

void load_data(string file){
	// load input data
	my_data.load(file.c_str(), raw_ascii);
	// transpose
	my_data = my_data.t();
	// repeat gestures, so it has testLen-initLen cols
	my_data = repmat(my_data, 1, (testLen-initLen)/my_data.n_cols);
	printf("Input dim: %u x %u\n\n", my_data.n_rows, my_data.n_cols);
	// save input data
	my_data.save("./results/my_data.mat", raw_ascii);
}

void load_test_data(int missingts, bool opt_randnoise){
	my_data_test = my_data;
	// time step, where it starts insert errors
	int insertionts = 250;
	// time length of one gesture
	int gesture_len = 500;
	// number of gestures in data
	int num_gestures = my_data_test.n_cols/gesture_len;
	// matrix to insert
	mat insert;
	if(opt_randnoise)
		insert = randu<mat>(my_data_test.n_rows, missingts);
	else
		insert = zeros<mat>(my_data_test.n_rows,missingts);
	for(int i = 0; i < num_gestures; i++){
		int idx_start = i*gesture_len+insertionts;
		my_data_test.cols(idx_start, idx_start + missingts-1) = insert;
	}
	// save test input data
	my_data_test.save("./results/my_data_test.mat", raw_ascii);
}

int main(){
	// create timer to measure runtime of your simulation
	Timer timer(true);

	// create SORN
	my_sorn = new SORN(resSize, a);

	// loading input
	load_data("./input/Input10d.txt");
	// loading input w/ missing time steps
	load_test_data(missingLen/*ts*/,false/*true=RANDOM NOISE; false=ZEROS*/);

	// plastic network
	my_sorn->train(my_data, trainLen);
	// static network with readout training
	my_output = my_sorn->test(my_data, my_data, testLen, true);
	// testing static network with independent data (missing ts)
	my_output_test = my_sorn->test(my_data_test, my_data, testLen, false);

	// save outputs
	my_output.save("./results/my_out.mat", raw_ascii);
	my_output_test.save("./results/my_out_test.mat", raw_ascii);
	mat dout = (my_output_test-my_output);
	dout.save("./results/dout.mat", raw_ascii);

	// delete SORN
	delete my_sorn;

	// timer prints runtime
	auto elapsed_secs = timer.Elapsed();
	printf("Executed. Runtime = %4.3f s.\n", elapsed_secs.count()/1000.);
}
