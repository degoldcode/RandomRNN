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

const int trainLen = 25000;				// time interval for training (in ts)
const int testLen = 5000;				// time interval for testing (in ts)
const int initLen = 100;				// time interval for initiating (in ts)

const int inSize = 1; 					// Number of input neurons (input dim)
const int outSize = 1;					// Number of input neurons (input dim)
const int resSize = 400;				// Number of reservoir neurons (resSize = 0.8*Ne + 0.2*Ni)
const double a = 1.0; 					// leaking rate
const double sR = 1.0;					// spectral radius

SORN * my_sorn;
mat my_data;

ofstream my_write;

int main(){
	Timer timer(true);							// create timer to measure runtime of your simulation
	my_write.open("./results/spectral_radius.dat");

	my_sorn = new SORN(resSize, sR);					// create SORN
	my_data = randu<mat>(inSize, trainLen);				// create training input data (random uniform [0,1])
	//my_data.zeros(inSize, testLen);					// no input data
	my_data.save("./results/my_data.mat", raw_ascii);	//save input data

	my_sorn->train(my_data, trainLen);

	delete my_sorn;
	my_write.close();

	auto elapsed_secs = timer.Elapsed();
	printf("Executed. Runtime = %4.3f s.\n", elapsed_secs.count()/1000.);
}
