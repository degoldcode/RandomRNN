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
#include "ESN.h"
#include "SORN.h"
#include "timer.h"
using namespace std;
using namespace arma;

const int trainLen = 10;				// time interval for training (in ts)
const int testLen = 5000;				// time interval for testing (in ts)
const int initLen = 0;					// time interval for initiating (in ts)
const int missingLen = 0;				// time interval of missing data (in ts)
const int numGestures = 11;				// number of gestures
const int gesture_len = 500;			// time length of each gesture
const int numTrials = 1;

const int inSize = 1; 				// Number of input neurons (input dim)
const int outSize = 1;				// Number of input neurons (input dim)
const int resSize = 200;				// Number of reservoir neurons (resSize = 0.8*Ne + 0.2*Ni)
const double a = 1.0; 					// leaking rate
const double sR = 0.95; 					// spectral radius

SORN * my_sorn;
ESN * my_esn;
cube data_cube;							// cube holding data for each gesture
mat noisa;								// noise data
mat my_data;							// input for training
mat my_datat;							// input for training
mat my_data_test;						// input for testing
mat my_data_testt;						// input for testing
mat my_teacher;							// teacher for readout training
mat my_output;
mat my_outputt;
mat my_output_test;
mat my_output_testt;

// if data needs to be written to file
ofstream my_write;

enum{sine,sin2cos,sqrf};

void load_data(string file){
	// load input data
	my_data.load(file.c_str(), raw_ascii);
	// transpose
	//my_data = my_data.t();
	// load gesture cube
	data_cube.zeros(my_data.n_rows, gesture_len, numGestures);
	for(int i = 0; i < numGestures; i++)
		data_cube.slice(i) = my_data.cols(i*500, i*500+499);
	//my_data = join_rows(data_cube.slice(0),data_cube.slice(1));
	//my_data = join_rows(my_data,data_cube.slice(4));
	//my_data.resize(10,500);										//only first gesture
	// repeat gestures, so it has testLen-initLen cols
	my_data = repmat(my_data, 1, (testLen-initLen)/my_data.n_cols);
	printf("Input dim: %u x %u\n\n", my_data.n_rows, my_data.n_cols);
	// save input data
	my_data.save("./results/my_data.mat", raw_ascii);
}

void load_data(int opt){
	// load input data
	my_data = zeros(1,500);
	for(int tstep = 0; tstep < my_data.n_cols; tstep++){
		if(opt == sine)
			my_data(0, tstep) = 0.5*(sin((2*M_PI/500)*tstep) + 1.0);
		if(opt == sin2cos)
			my_data(0, tstep) = sin((4*M_PI/500)*tstep)*cos((2*M_PI/500)*tstep);
		if(opt == sqrf)
			my_data(0, tstep) = sin((2*M_PI/500)*tstep);
	}
	if(opt == sqrf){
		my_data = 0.5*(sign(my_data) + 1.0);
	}
	my_data = repmat(my_data, 1, (testLen-initLen)/my_data.n_cols);
	my_datat = my_data.t();
	printf("Input dim: %u x %u\n\n", my_data.n_rows, my_data.n_cols);
	// save input data
	my_datat.save("./results/my_data.mat", raw_ascii);
}

void load_test_data(int missingts, bool opt_randnoise, bool just_last){
	my_data_test = my_data;
	// time step, where it starts insert errors
	int insertionts = 100;
	// number of gestures in data
	int num_gestures = my_data_test.n_cols/gesture_len;
	// matrix to insert
	mat insert;
	if(opt_randnoise)
		insert = randu<mat>(my_data_test.n_rows, missingts);
	else
		insert = zeros<mat>(my_data_test.n_rows, missingts);
	if(!just_last){
		for(int i = 0; i < num_gestures; i++){
			int idx_start = i*gesture_len+insertionts;
			if(missingts>0)
				my_data_test.cols(idx_start, idx_start + missingts-1) = insert;
		}
	}
	else{
		insert = zeros<mat>(my_data_test.n_rows,missingts);
		if(missingts>0)
			my_data_test.cols(my_data_test.n_cols-missingts, my_data_test.n_cols-1) = insert;
	}
	my_data_testt = my_data_test.t();
	// save test input data
	my_data_testt.save("./results/my_data_test.mat", raw_ascii);
}

void load_teacher(mat in, int numgestures, int len){
	my_teacher.zeros(numgestures, in.n_cols);
	for(int i = 0; i < in.n_cols; i++)
		my_teacher((i/500)%numgestures, i) = 1.0;
	my_teacher.save("./results/my_teacher.mat", raw_ascii);
}

int main(){
	// create timer to measure runtime of your simulation
	Timer timer(true);

	// loading input
	load_data(sqrf);
	noisa.randu(my_data.n_rows, trainLen);
	// loading input w/ missing time steps
	load_test_data(missingLen/*ts*/, true/*true=RANDOM NOISE; false=ZEROS*/, false);
	// loading teacher data
	load_teacher(my_data, 3, testLen);

	for(int trialIdx = 0; trialIdx < numTrials; trialIdx++){
		// create SORN
		//my_sorn = new SORN(resSize, a);
		my_esn = new ESN(resSize, a, sR);

		// plastic network
		//my_sorn->train(my_data, trainLen);
		// static network with readout training
		my_output = my_esn->test(my_data, my_data, testLen, true);
		my_outputt = my_output.t();
		// testing static network with independent data (missing ts)
		my_output_test = my_esn->test(my_data_test, my_data, testLen, false);
		my_output_testt = my_output_test.t();
	}

	// save outputs
	my_outputt.save("./results/my_out.mat", raw_ascii);
	my_output_testt.save("./results/my_out_test.mat", raw_ascii);
	mat dout = (my_output_testt-my_outputt);
	dout.save("./results/dout.mat", raw_ascii);

	// delete SORN
	delete my_sorn;

	// timer prints runtime
	auto elapsed_secs = timer.Elapsed();
	printf("Executed. Runtime = %4.3f s.\n", elapsed_secs.count()/1000.);
}
