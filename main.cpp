/*
 * main.cpp
 *
 *  Created on: Aug 6, 2014
 *      Author: degoldschmidt
 */

#include <armadillo>
#include <iostream>
#include <fstream>
#include "SORN.h"
#include "timer.h"
using namespace std;
using namespace arma;

const int trainLen = 20000;				// time interval for training (in ts)
const int testLen = 5000;				// time interval for testing (in ts)
const int initLen = 1000;				// time interval for initiating (in ts) NOT USED ATM

const int inSize = 100; 				// Number of input neurons (input dim)
const int outSize = 100;				// Number of input neurons (input dim)
const int resSize = 600;				// Number of reservoir neurons (resSize = 0.8*Ne + 0.2*Ni)
const double a = 1.0; 					// leaking rate
const double sR = 1.0;					// spectral radius

SORN * my_sorn;
mat my_data;
mat my_output;

ofstream my_write;						// if data needs to be written

int main(){
	Timer timer(true);												// create timer to measure runtime of your simulation

	my_sorn = new SORN(resSize, sR);								// create SORN
	my_data.load("./input/inputtrain.mat", raw_ascii);				// load input data
	my_data.resize(100,500);										// only take first gesture
	my_data = repmat(my_data, 1, (testLen-initLen)/my_data.n_cols);	// repeat gesture, so it has testLen-initLen cols
	my_data.save("./results/my_data.mat", raw_ascii);				// save input data

	my_sorn->train(my_data, trainLen);								// plastic network
	my_output = my_sorn->test(my_data, my_data, testLen, true);		// static network with readout training
	my_output.save("./results/my_out.mat", raw_ascii);				// save outputs

	delete my_sorn;

	auto elapsed_secs = timer.Elapsed();
	printf("Executed. Runtime = %4.3f s.\n", elapsed_secs.count()/1000.);
}
