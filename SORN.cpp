/*
 * SORN.cpp
 *
 *  Created on: Aug 6, 2014
 *      Author: degoldschmidt
 */

#include <cstdlib>
#include <cstdio>
#include "SORN.h"

SORN::SORN(int mN, double des_sR){
	arma_rng::set_seed_random();
	printf("\nASSEMBLE A SELF-ORGANIZING RECURRENT NETWORK\n\n");

	N = mN;
	Ne = int(0.8*N);
	Ni = int(0.2*N);
	Nu = 10;

	w_sp = 10./Ne;
	leak_rate = 1.0;
	target_rate = (2.*Nu/Ne)*ones<vec>(Ne);
	mu_IP = 0.001;
	mu_STDP = 0.001;

	//Random initial weights
	printf("Generating weight matrices (EE, IE, EI) with random connectivity (EE sparseness = %1.3f)...\n", w_sp);
	sp_mat w_mat = sprandu<sp_mat>(Ne, Ne, w_sp);
	W = w_mat;
	Wie = randu<mat>(Ni, Ne);
	Wei = randu<mat>(Ne, Ni);
	//W.diag() = zeros<vec>(W.n_rows);
	//Wie.diag() = zeros<vec>(Wie.n_rows);
	//Wei.diag() = zeros<vec>(Wei.n_cols);
	//printf("Diagonale has norm: %f (EE), %f (IE), %f (EI)\n", norm(W.diag()), norm(Wie.diag()), norm(Wei.diag()));
	printf("Done.\n\n");
	printf("Index actual connections (W(i,j) > 0.)...\n");
	idx_W = find(W>0);
	printf("Found %u connections.\n", idx_W.n_elem);
	printf("Average degree per unit = %f .\n\n", double(idx_W.n_elem/Ne));
	//printf("Computing spectral radius...\n");
	//w_sR = max_eigenval(W);
	//printf("sR = %4.3f\n", w_sR);
	//printf("Set spectral radius to %f...\n", des_sR);
	printf("Normalize weights row-wise...\n");
	W = normalise(W,2,1);
	Wei = normalise(Wei,2,1);
	Wie = normalise(Wie,2,1);
	printf("Done.\n\n");
	//W  *= (des_sR/w_sR);
	//w_sR = max_eigenval(W);								//Just to be sure...
	//printf("sR = %4.3f\n\n", w_sR);



	//Random initial states
	printf("Initialize states and thresholds...\n");
	x = randu<vec>(Ne)-0.5;
	R = x;
	x = sign(sign(x)+1.);
	xp = x;
	y = randu<vec>(Ni)-0.5;
	y = sign(sign(y)+1.);
	te_max = 0.5;
	ti_max = 1.0;
	the = te_max * randu<vec>(Ne);
	thi = ti_max * randu<vec>(Ni);

	STDP = true;
	SN = true;
	IP = true;

	printf("Done.\n\n");
}


SORN::~SORN(){

}


vec SORN::connect(vec & in, int in_degr, int connect_dim){
	vec out = zeros(connect_dim);
	imat idx = randi<imat>(in_degr, in.n_elem, distr_param(0, connect_dim));
	for(int j = 0; j < idx.n_cols; j++){
		for(int i = 0; i < idx.n_rows; i++)
			out(idx(i,j)) += 1.0;
	}
	return out;
}


/*
 * Returns the maximum absolute eigenvalue of a given matrix
 */
double SORN::max_eigenval(mat & in){
	cx_colvec eigval;
	cx_mat eigvec;
	eig_gen (eigval, eigvec, in);

	colvec abs_eigval = abs(eigval);
	colvec abs_eigval_sorted = sort(abs_eigval, 1);
	return abs_eigval_sorted(0);
}


/*
 * Returns a matrix the same size as A and B with the largest elements
 * taken from A or B. The dimensions of A and B must match.
 */
mat SORN::max_of(mat A, mat B){
	mat out(A.n_rows,A.n_cols);
	for(int i = 0; i < A.n_rows; i++){
		for(int j = 0; j < A.n_cols; j++){
			if(A(i,j)>B(i,j))
				out(i,j) = A(i,j);
			else
				out(i,j) = B(i,j);
		}
	}
	return out;
}


/*
 * Returns a matrix the same size as A and B with the smallest elements
 * taken from A or B. The dimensions of A and B must match.
 */
mat SORN::min_of(mat A, mat B){
	mat out(A.n_rows,A.n_cols);
	for(int i = 0; i < A.n_rows; i++){
		for(int j = 0; j < A.n_cols; j++){
			if(A(i,j)<B(i,j))
				out(i,j) = A(i,j);
			else
				out(i,j) = B(i,j);
		}
	}
	return out;
}

void SORN::update(vec in){
    vec lastx  = x;

	// Input
    //printf("---input\n");
	vec u = Wu*in;

	// Network activation
	//printf("---states\n");
	R = W * x - Wei * y - the;
	xp = sign(sign(R) + 1.);
	x = sign(sign(R + u) + 1.);
	x = (1.-leak_rate)*lastx + leak_rate*x;					// leaky-integrated discrete-time continuous-value units
	y = sign(sign(Wie * x - thi) + 1.);

    // Plasticity
    if(IP){
    	//printf("---IP\n");
        the += mu_IP * (x-target_rate);  					// intrinsic plasticity
    }
    if (STDP){
    	//printf("---STDP\n");
		mat A = lastx * x.t();
		mat delta = mu_STDP * (A.t() - A);
		W(idx_W) += delta(idx_W);  							//additive STDP
		//W *= 2.;

		W = min_of(W, ones<mat>(Ne, Ne));   				// clip weights to [0,1]
		W = max_of(W, zeros<mat>(Ne, Ne));

		if(SN){
			//printf("---Norm\n");
			W = normalise(W,2,1);							// synaptic normalization
		}
    }
}

/*
 * Trains the SORN for a given time interval (int time) with given input data (mat data)
 */
void SORN::train(mat data, int time){
	STDP = true;
	SN = true;
	IP = true;
	printf("START TRAINING PHASE FOR %u TIMESTEPS\nSTDP = %s\nSN = %s\nIP = %s\n\n", time, STDP ? "ON" : "OFF", SN ? "ON" : "OFF", IP ? "ON" : "OFF");

	//Set up input connections
	double input_conn = (1.0*data.n_rows*Nu)/(1.*Ne);
	printf("Set up input connections to the reservoir with connectivity = %f\n\n", input_conn);
	sp_mat A = sprandu<sp_mat>(Ne, data.n_rows, input_conn);
	sp_mat B = spones(A);
	Wu = B;
	uvec Wu_idx = find(Wu > 0.0);
	cout << Wu_idx.n_elem << endl;

	init_W = W;
	H.zeros(Ne, time);
	Hp.zeros(Ne, time);

	for(int t = 0; t < time; t++){			// test run
		H.col(t) = x;
		Hp.col(t) = xp;
		update(data.col(t));
		if(t%(time/10)==0)
			cout << 100*t/time << "% done.\n";
	}
	end_W = W;
	delta_W = end_W - init_W;
	W.save("./results/end_W.mat", raw_ascii);

	//save matrices
	init_W.save("./results/init_W.mat", raw_ascii);
	delta_W.save("./results/delta_W.mat", raw_ascii);
	end_W.save("./results/end_W.mat", raw_ascii);
	H.save("./results/res_st.mat", raw_ascii);
	Hp.save("./results/res_pst.mat", raw_ascii);
	dH = H-Hp;
	dH.save("./results/res_dst.mat", raw_ascii);
}


void SORN::test(mat data, int time){

}
