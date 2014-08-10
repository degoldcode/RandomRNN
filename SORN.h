/*
 * SORN.h
 *
 *  Created on: Aug 6, 2014
 *      Author: degoldschmidt
 */

#ifndef SORN_H_
#define SORN_H_

#include <armadillo>
#include <string>
#include <sstream>
using namespace std;
using namespace arma;

class SORN {
public:
	SORN(int mN, double des_sR);
	~SORN();

	double max_eigenval(mat & in);
	mat max_of(mat A, mat B);
	mat min_of(mat A, mat B);
	void save_matrices(string mode);
	void set_input_con(mat in);
	void set_matrix(int time);
	void train(mat data, int time);
	mat test(mat data, mat teacher, int time, bool trainOut);
	void update(vec in);

	vec x;				// state of reservoir excitatory units
	vec R;				// recurrent drive
	vec xp;				// pseudo-state of reservoir excitatory units (no input drive)
	vec y;				// state of reservoir inhibitory units
	vec the;			// excitatory thresholds
	vec thi;			// inhibitory thresholds
	double te_max;		// maximum excitatory thresholds
	double ti_max;		// maximum inhibitory thresholds

	mat H;				// reservoir state matrix
	mat Hp;				// reservoir pseudo-state matrix
	mat dH;

	mat W;				// excitatory weights of reservoir
	mat Wie;			// excitatory-to-inhibitory weights of reservoir
	mat Wei;			// inhibitory-to-excitatory weights of reservoir
	mat Wu;				// input weights
	mat init_W;			// initial reservoir weights
	mat end_W;			// final reservoir weights
	mat delta_W;		// changes of reservoir weights
	mat Wout;			// readout weights
	uvec idx_W;			// indices of positive weights
	double mean_dw;		// mean change in weights (running sum)

	int N;				// total number of reservoir units
	int Ne;				// number of reservoir excitatory units (0.8*N)
	int Ni;				// number of reservoir inhibitory units (0.2*N)
	int Nu;				// number of reservoir units directly driven by one of the inputs
	int No;				// number of readout units

	double w_sp;		// sparseness of reservoir weight matrix
	//double w_sR;		// spectral radius of reservoir weight matrix

	double leak_rate;	// leaking rate (1 => "nothing is leaking")

	double mu_STDP;		// learning rate STDP
	double mu_IP;		// learning rate IP
	vec target_rate;	// target rate for IP
	bool STDP;			// Boolean (if STDP is on)
	bool SN;			// Boolean (if STDP is on)
	bool IP;			// Boolean (if STDP is on)
	int T_mat;			// maximum time interval saved in matrices
};

#endif /* SORN_H_ */
