/*
 * SORN.cpp
 *
 *  Created on: Aug 6, 2014
 *      Author: degoldschmidt
 */

#include <cstdlib>
#include <cstdio>
#include "SORN.h"


/*
 *
 * Constructor of the SORN class. Creates a random graph of excitatory and inhibitory units.
 *
 * PARAM:	(int)		$mN 	= number of excitatory units
 * 			(double)	$des_sR	= desired spectral radius (for ESN)
 *
 * RETURN:	class object of SORN
 *
 */
SORN::SORN(int mN, double des_sR){
	arma_rng::set_seed_random();
	printf("\nASSEMBLE A SELF-ORGANIZING RECURRENT NETWORK WITH %u EXCITATORY UNITS\n\n", mN);

	//Number of units
	Ne = mN;
	Ni = int(0.2*Ne);
	N = Ne+Ni;
	Nu = 100;

	//Random initial weights
	w_sp = 10./Ne;
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
	printf("Normalize weights row-wise...\n");
	W = normalise(W,2,1);
	Wei = normalise(Wei,2,1);
	Wie = normalise(Wie,2,1);
	// Save initial weight matrix
	init_W = W;
	init_W.save("./results/init_W.mat", raw_ascii);
	printf("Done.\n\n");
	mean_dw = 0.0;


	//Random initial states
	printf("Initialize states & thresholds...\n");
	leak_rate = 1.;
	x.randu(Ne);
	y.randu(Ni);
	xp.randu(Ne);
	te_max = 0.5;
	ti_max = 1.;
	the = te_max * randu<vec>(Ne);
	thi = ti_max * randu<vec>(Ni);
	printf("Leaking rate = %1.3f\nMaximum excitatory threshold = %1.3f\nMaximum inhibitory threshold = %1.3f\n", leak_rate, te_max, ti_max);
	printf("Done.\n\n");

	printf("Initialize learning parameters...\n");

	double tr = 2.*(1.*Nu)/(1.*Ne);
	target_rate = tr*ones<vec>(Ne);//(2.*Nu/Ne)*ones<vec>(Ne);
	mu_IP = 0.001;//0.5;
	mu_STDP = 0.001;
	printf("Target firing rate = %1.6f\nLearning rate STDP = %1.3f\nLearning rate IP = %1.3f\n", tr, mu_STDP, mu_IP);

	STDP = true;
	SN = true;
	IP = true;
	T_mat = 25000;

	printf("Done.\n\n");
}


/*
 *
 * Destructor of the SORN class. Resets state histories and weight matrices.
 *
 * PARAM: (void)
 *
 * RETURN: (void)
 *
 */
SORN::~SORN(){
	Hp.reset();
	H.reset();
	W.reset();
}


/*
 *
 * Returns the maximum absolute eigenvalue of a given matrix.
 *
 * PARAM:	(mat)		$in	= matrix of which you want to determine the maximum absolute eigenvalue
 *
 * RETURN:	(double)	maximum absolute eigenvalue of $in
 *
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
 *
 * Returns a matrix the same size as A and B with the largest elements
 * taken from A or B. The dimensions of A and B must match.
 *
 * PARAM:	(mat) $A	= first matrix to compare
 * 			(mat) $B	= second matrix to compare
 *
 * RETURN:	(mat) matrix with elements being the maximum of respective elements of $A and $B
 *
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
 *
 * Returns a matrix the same size as A and B with the smallest elements
 * taken from A or B. The dimensions of A and B must match.
 *
 * PARAM:	(mat) $A	= first matrix to compare
 * 			(mat) $B	= second matrix to compare
 *
 * RETURN:	(mat) matrix with elements being the minimum of respective elements of $A and $B
 *
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


/*
 *
 * Saves state history matrices of reservoir states, pseudostates, and their difference.
 *
 * PARAM:	(string)	$mode = string (example: "bla") for defining a suffix of the file name
 *
 * RETURN:	(void)
 *
 */
void SORN::save_matrices(string mode){
	stringstream ss;
	ss << "./results/res_st_" << mode << ".mat";
	H.save(ss.str().c_str(), raw_ascii);
	ss.str( string() );
	ss.clear();
	ss << "./results/res_pst_" << mode << ".mat";
	Hp.save(ss.str().c_str(), raw_ascii);
	dH = H-Hp;
	ss.str( string() );
	ss.clear();
	ss << "./results/res_dst_" << mode << ".mat";
	dH.save(ss.str().c_str(), raw_ascii);
}


/*
 *
 * Sets up the input connections $Wu to the reservoir. $Nu units of $Ne get randomly chosen
 * and evenly connected with the input dimensions.
 *
 * PARAM:	(mat)	in	= input data to SORN
 *
 * RETURN:	(void)
 *
 */
void SORN::set_input_con(mat in){
	double input_conn = (1.0*Nu)/(1.*Ne*in.n_rows);
	if(Nu < in.n_rows)
		printf("!!! Warning: Not every input dimension is connected to the reservoir.\n\n");
	if(input_conn > 0.5)
		printf("!!! Warning: Reservoir is too input-driven. Increase reservoir size for respective input dimension.\n\n");
	printf("Set up input connections to the reservoir with connectivity = %f\n\n", input_conn);
	Wu = zeros<mat>(Ne, in.n_rows);
	ivec rand_idx = randi<ivec>(Nu, distr_param(0,Ne-1));
	for(int idx_col = 0; idx_col < Nu; idx_col++){
		if(idx_col%in.n_rows == 0 && Nu-idx_col<in.n_rows){
			ivec rand_idx2 = randi<ivec>(Nu-idx_col, distr_param(0,in.n_rows-1));
			Wu(rand_idx(idx_col), rand_idx2(idx_col%in.n_rows)) = 1.0;
		}
		else
			Wu(rand_idx(idx_col), idx_col%in.n_rows) = 1.0;
	}
	printf("Reservoir units driven by input = %u\n\n", int(accu(Wu)));
}


/*
 *
 * Sets up state history matrices for given time. $T_mat is the maximum size of the
 * matrices (memory issue).
 *
 * PARAM:	(int)	$time	= number of columns for state history matrices (usually time interval of run)
 *
 * RETURN:	(void)
 *
 */
void SORN::set_matrix(int time){
	if(time < T_mat){
		H.zeros(Ne, time);
		Hp.zeros(Ne, time);
	}
	else{
		H.zeros(Ne, T_mat);
		Hp.zeros(Ne, T_mat);
	}
}


/*
 *
 * Runs the plastic network for a given time interval with given input data.
 *
 * PARAM:	(mat)	$data	= input data for training run
 * 			(int)	$time	= time interval to run
 *
 * RETURN:	(void)
 *
 */
void SORN::train(mat data, int time){
	STDP = true;
	SN = true;
	IP = true;
	printf("START TRAINING PHASE FOR %u TIMESTEPS\nSTDP = %s\nSN = %s\nIP = %s\n\n", time, STDP ? "ON" : "OFF", SN ? "ON" : "OFF", IP ? "ON" : "OFF");

	set_input_con(data);
	set_matrix(time);

	for(int t = 0; t < time; t++){
		if(t%T_mat==0)
			save_matrices("train");
		H.col(t%T_mat) = x;
		Hp.col(t%T_mat) = xp;

		mat dW = -W;
		update(data.col(t%data.n_cols));
		dW += W;
		mean_dw += accu(abs(dW))/dW.n_elem;


		if(t%5000==0)
			printf("%6u ts: <dw> = %f\n", t, mean_dw/t);
		if(t==time-1 && t%T_mat>100)
			save_matrices("train");
	}
	end_W = W;
	end_W.save("./results/end_W.mat", raw_ascii);
	delta_W = end_W - init_W;
	delta_W.save("./results/delta_W.mat", raw_ascii);
}


/*
 *
 * Runs the static network for a given time with given input and teacher data.
 *
 * PARAM:	(mat)	$data		= input data for testing run
 * 			(mat)	$teacher	= teacher data for testing run
 * 			(int)	$time		= time interval to run
 * 			(bool)	$trainOut	= true, if output weights will be trained after run
 *
 * RETURN:	(mat)	output matrix of states (Wout*H)
 *
 */
mat SORN::test(mat data, mat teacher, int time, bool trainOut){
	STDP = false;
	SN = false;
	IP = false;
	printf("START TESTING PHASE FOR %u TIMESTEPS\n\n", time);

	//set_input_con(data);
	set_matrix(time);

	for(int t = 0; t < time; t++){			// test run
		H.col(t) = x;
		Hp.col(t) = xp;
		update(data.col(t%data.n_cols));
		if(t%1000==0)
			printf("%6u ts.\n", t);
	}
	save_matrices("test");
	No = data.n_rows;
	Wout.zeros(Ne, No);

	if(trainOut){
		//	Resize state history to match teacher data (cropping the first time steps -> initLen in main.cpp)
		mat shortHp = Hp;
		mat shortH = H;
		shortHp.resize(Hp.n_rows, teacher.n_cols);
		if(accu(abs(teacher)) != 0.0){
			Wout = teacher * pinv(shortHp);												// direct method: Pseudo-inverse
			//Wout = teacher * H.t() * pinv(H*H.t());									// normal functions
			//alpha = 0.7;																// regularization factor
			//Wout = teacher * H.t() * pinv(H*H.t() + alpha*alpha*ones<mat>(Ne,Ne));	// Tikhonov regularization
		}
	}

	mat A = (Wout*shortHp)-teacher;
	A.save("./results/dout.mat", raw_ascii);
	return Wout*H;
}


/*
 *
 * Updates the network for one time step.
 *
 * PARAM:	(vec)	$in	= current input vector
 *
 * RETURN:	(void)
 *
 */
void SORN::update(vec in){
    vec lastx  = x;
    vec lastxp  = xp;

	// Input
    //printf("---input\n");
	vec u = Wu*in;

	// Network activation
	//printf("---states\n");
	R = W * x - Wei * y - the;
	xp = sign(sign(R) + 1.);
	xp = (1.-leak_rate)*lastxp + leak_rate*xp;
	x = sign(sign(R + u) + 1.);
	x = (1.-leak_rate)*lastx + leak_rate*x;		// leaky-integrated discrete-time continuous-value units
	y = sign(sign(Wie * x - thi) + 1.);

    // Plasticity
    if(IP){
    	//printf("---IP\n");
        the += mu_IP * (x-target_rate);  		// intrinsic plasticity
    }
    if (STDP){
    	//printf("---STDP\n");
		mat A = lastx * x.t();
		mat delta = mu_STDP * (A.t() - A);
		W(idx_W) += delta(idx_W);  				// additive STDP
		//W += delta;

		W = min_of(W, ones<mat>(Ne, Ne));  		// clip weights to [0,1]
		W = max_of(W, zeros<mat>(Ne, Ne));

		if(SN){
			//printf("---Norm\n");
			W = normalise(W,2,1);				// synaptic normalization
		}
    }
}
