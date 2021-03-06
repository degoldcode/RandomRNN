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
 * PARAM:	(int)		$mN 			= number of excitatory units
 * 			(double)	$des_sR			= desired spectral radius (for ESN)
 * 			(bool)		$opt_verbose	= if trial information shall be printed in console
 *
 * RETURN:	class object of SORN
 *
 */
SORN::SORN(int mN, double leak, bool opt_verbose, double mte, double mti){
	arma_rng::set_seed_random();
	VERBOSE = opt_verbose;
	(VERBOSE)?printf("\nASSEMBLE A SELF-ORGANIZING RECURRENT NETWORK WITH %u EXCITATORY UNITS\n\n", mN):VERBOSE;

	//#include "network.conf"

	//Number of units
	Ne = mN;
	Ni = int(0.2*Ne);
	N = Ne+Ni;
	Nu = 10;//int(0.05*Ne);

	//Random initial weights
	w_sp = 10./Ne;
	(VERBOSE)?printf("Generating weight matrices (EE, IE, EI) with random connectivity (EE sparseness = %1.3f)...\n", w_sp):VERBOSE;
	sp_mat w_mat = sprandu<sp_mat>(Ne, Ne, w_sp);
	W = w_mat;
	Wie = randu<mat>(Ni, Ne);
	Wei = randu<mat>(Ne, Ni);
	W.diag() = zeros<vec>(W.n_rows);
	Wie.diag() = zeros<vec>(Wie.n_rows);
	Wei.diag() = zeros<vec>(Wei.n_cols);
	//(VERBOSE)?printf("Diagonale has norm: %f (EE), %f (IE), %f (EI)\n", norm(W.diag()), norm(Wie.diag()), norm(Wei.diag()));
	(VERBOSE)?printf("Done.\n\n"):VERBOSE;
	(VERBOSE)?printf("Index actual connections (W(i,j) > 0.)...\n"):VERBOSE;
	idx_W = find(W>0.);
	(VERBOSE)?printf("Found %u connections.\n", idx_W.n_elem):VERBOSE;
	(VERBOSE)?printf("Average degree per unit = %f .\n\n", 1.0*idx_W.n_elem/Ne):VERBOSE;
	//(VERBOSE)?printf("Computing spectral radius...\n");
	//w_sR = max_eigenval(W);
	//(VERBOSE)?printf("sR = %4.3f\n", w_sR);
	(VERBOSE)?printf("Normalize weights row-wise...\n"):VERBOSE;
	W = normalize(W);
	Wei = normalize(Wei);
	Wie = normalize(Wie);
	// Save initial weight matrix
	init_W = W;
	init_W.save("./results/init_W.mat", raw_ascii);
	(VERBOSE)?printf("Done.\n\n"):VERBOSE;
	mean_dw = 0.0;


	//Random initial states
	(VERBOSE)?printf("Initialize states & thresholds...\n"):VERBOSE;
	leak_rate = leak;
	x.zeros(Ne);
	y.zeros(Ni);
	xp.zeros(Ne);
	te_max = mte;
	ti_max = mti;
	the = te_max * randu<vec>(Ne);
	thi = ti_max * randu<vec>(Ni);
	(VERBOSE)?printf("Leaking rate = %1.3f\nMaximum excitatory threshold = %1.3f\nMaximum inhibitory threshold = %1.3f\n", leak_rate, te_max, ti_max):VERBOSE;
	(VERBOSE)?printf("Done.\n\n"):VERBOSE;

	(VERBOSE)?printf("Initialize learning parameters...\n"):VERBOSE;

	double tr = 0.2;//2.*(1.*Nu)/(1.*Ne);
	target_rate = tr*ones<vec>(Ne);//(2.*Nu/Ne)*ones<vec>(Ne);
	mu_IP = 0.001;//0.5;
	mu_STDP = 0.001;
	(VERBOSE)?printf("Target firing rate = %1.6f\nLearning rate STDP = %1.3f\nLearning rate IP = %1.3f\n", tr, mu_STDP, mu_IP):VERBOSE;

	STDP = true;
	SN = true;
	IP = true;
	sampling_rate = 1;
	T_mat = int(25000/sampling_rate);

	(VERBOSE)?printf("Done.\n\n"):VERBOSE;
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


mat SORN::normalize(mat in){
	double sum;
	for(int i = 0; i < in.n_rows; i++){
		sum = 0.0;
		for(int j = 0; j < in.n_cols; j++)
			sum += in(i,j);
		in.row(i) /= sum;
	}
	return in;
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
	//dH = H-Hp;
	//ss.str( string() );
	//ss.clear();
	//ss << "./results/res_dst_" << mode << ".mat";
	//dH.save(ss.str().c_str(), raw_ascii);
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
	Nu = 10;//in.n_rows;
	double input_conn = 10./Ne;///10*(1.0*Nu)/(1.*Ne*in.n_rows);
	if(Nu < in.n_rows)
		(VERBOSE)?printf("!!! Warning: Not every input dimension is connected to the reservoir.\n"):VERBOSE;
	if(input_conn > 0.5)
		(VERBOSE)?printf("!!! Warning: Reservoir is too input-driven. Increase reservoir size for respective input dimension.\n\n"):VERBOSE;
	(VERBOSE)?printf("Set up input connections to the reservoir with connectivity = %f\n\n", input_conn):VERBOSE;
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
	(VERBOSE)?printf("%u x %u: Reservoir units driven by input = %u\n", Wu.n_rows, Wu.n_cols, int(accu(Wu))):VERBOSE;
	(VERBOSE)?printf("Done.\n\n"):VERBOSE;
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
	(VERBOSE)?printf("Set state history matrices...\n"):VERBOSE;
	if(time < T_mat){
		H.zeros(Ne, time);
		Hp.zeros(Ne, time);
	}
	else{
		H.zeros(Ne, T_mat);
		Hp.zeros(Ne, T_mat);
	}
	(VERBOSE)?printf("Done.\n\n"):VERBOSE;
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
	(VERBOSE)?printf("START TRAINING PHASE FOR %u TIMESTEPS\nSTDP = %s\nSN = %s\nIP = %s\n\n", time, STDP ? "ON" : "OFF", SN ? "ON" : "OFF", IP ? "ON" : "OFF"):VERBOSE;
	(VERBOSE)?printf("Input dim: %u x %u\n\n", data.n_rows, data.n_cols):VERBOSE;
	set_input_con(data);
	set_matrix(time);

	for(int t = 0; t < time; t++){
		if(t%sampling_rate==0){
			H.col(t%T_mat) = x;
			//Hp.col(t%T_mat) = xp;
//		if(t%1000==0){
//			(VERBOSE)?printf("Computing spectral radius...\n");
//			w_sR = max_eigenval(W);
//			(VERBOSE)?printf("sR = %4.3f\n", w_sR);
//		}
		}
		mat dW = -W;
		update(data.col(t%data.n_cols));
		dW += W;
		mean_dw += accu(abs(dW))/dW.n_elem;


		if(t%1000==0)
			(VERBOSE)?printf("%6u ts: <dw> = %e, Input = %3.6f, X = %3.6f, W = %3.6f\n", t, mean_dw/t, accu(Wu*data.col(t%data.n_cols))/Wu.n_rows, accu(x)/x.n_elem, accu(W)/W.n_elem):VERBOSE;
		if(t==time-1)
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
	(VERBOSE)?printf("START TESTING PHASE FOR %u TIMESTEPS\n\n", time):VERBOSE;

	x.zeros(Ne);
	y.zeros(Ni);
	xp.zeros(Ne);
	dist.reset();

	//if(trainOut)
		//set_input_con(data);
	set_matrix(time);

	for(int t = 0; t < time; t++){			// test run
		H.col(t%T_mat) = x;
		Hp.col(t%T_mat) = xp;
		update(data.col(t%data.n_cols));
		if(t%1000==0){
			(VERBOSE)?printf("%6u ts: Input = %3.6f, X = %3.6f, W = %3.6f\n", t, accu(Wu*data.col(t%data.n_cols))/Wu.n_rows, accu(x)/x.n_elem, accu(W)/W.n_elem):VERBOSE;
		}
	}
	nonzero_x = find(H > 0.);
	Hzero = 1.0*nonzero_x.n_elem/(1.0*Ne*time);
	umat inactive = any(H, 1);
	//cout << inactive << endl;
	Izero = (1.0*(Ne-as_scalar(accu(inactive))))/(1.*Ne);

	save_matrices("test");
	No = teacher.n_rows;

	if(trainOut)
		Wout.zeros(No, Ne);
	mat shortHp = Hp;
	mat shortH = H;
	//	Resize state history to match teacher data (cropping the first time steps -> initLen in main.cpp)
	if((Hp.n_cols - teacher.n_cols) > 0)
		shortHp.shed_cols(0, Hp.n_cols - teacher.n_cols - 1);

	if(trainOut){
		if(accu(abs(teacher)) != 0.0){
			//Wout = teacher * pinv(shortHp);												// direct method: Pseudo-inverse
			//Wout = teacher * H.t() * pinv(H*H.t());									// normal functions
			double alpha = 0.7;																// regularization factor
			Wout = teacher * shortHp.t() * pinv(shortHp*shortHp.t() + alpha*alpha*ones<mat>(Ne,Ne));	// Tikhonov regularization
		}
	}
	mat sub = (Wout*shortHp);
	mat A = sub-teacher;
	A.save("./results/dout.mat", raw_ascii);
	pred_error = accu(abs(A))/(A.n_elem);
	//if(!trainOut)
		//cout << "NRMSE = " << pred_error << endl;
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
    //(VERBOSE)?printf("---input\n");
	vec u = Wu*in;

	// Network activation
	//(VERBOSE)?printf("---states\n");

	imat A = randi<imat>(1, 1, distr_param(0, x.n_elem-1));
	x_ = x;
	if(x_(A(0,0)) == 0.)
		x_(A(0,0)) = 1.0;
	else
		x_(A(0,0)) = 0.0;

	R = W * x - Wei * y - the;
	R_ = W * x_ - Wei * y - the;
	xp = sign(sign(R) + 1.);
	//xp = (1.-leak_rate)*lastxp + leak_rate*xp;

	//cout << mean(x) << "\t - " << mean(y) << "\t - " << mean(the)  << endl;
	x = sign(sign(R + u) + 1.);
	//x = (1.-leak_rate)*lastx + leak_rate*x;		// leaky-integrated discrete-time continuous-value units
	x_ = sign(sign(R_ + u) + 1.);
	//cout << "\t" << sum(abs(x_-x)) << endl;
	y = sign(sign(Wie * x - thi) + 1.);

	dist(sum(abs(x_-x)));

    // Plasticity
    if(IP){
    	//(VERBOSE)?printf("---IP\n");
        the += mu_IP * (x-target_rate);  		// intrinsic plasticity
    }
    if (STDP){
    	//(VERBOSE)?printf("---STDP\n");
		mat A = lastx * x.t();
		mat delta = mu_STDP * (A.t() - A);
		W(idx_W) += delta(idx_W);  				// additive STDP
		//W += delta;

		W = min_of(W, ones<mat>(Ne, Ne));  		// clip weights to [0,1]
		W = max_of(W, zeros<mat>(Ne, Ne));

		if(SN){
			//(VERBOSE)?printf("---Norm\n");
			W = normalize(W);				// synaptic normalization
		}
    }
}
