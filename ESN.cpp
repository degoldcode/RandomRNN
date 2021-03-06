/*
 * ESN.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: degoldschmidt
 */


#include "ESN.h"

/*
 *
 * Constructor of the ESN class. Creates a random graph of reservoir units.
 *
 * PARAM:	(int)		$mN 	= number of excitatory units
 * 			(double)	$leak	= leaking rate (leaky-integrator)
 * 			(double)	$sR	= desired spectral radius (for ESN)
 *
 * RETURN:	class object of SORN
 *
 */
ESN::ESN(int mN, double leak, double sR){
	arma_rng::set_seed_random();
	printf("\nASSEMBLE AN ECHO STATE NETWORK WITH %u UNITS\n\n", mN);

	//#include "network.conf"

	//Number of units
	N = mN;
	Nu = 10;//int(0.05*Ne);

	//Random initial weights
	Wu = zeros(1,1);
	w_sp = 5./N;
	printf("Generating weight matrices (EE, IE, EI) with random connectivity (EE sparseness = %1.3f)...\n", w_sp);
	sp_mat w_mat = sprandn<sp_mat>(N, N, w_sp);
	W = w_mat;
	W.diag() = zeros<vec>(W.n_rows);
	printf("Done.\n\n");
	printf("Index actual connections (W(i,j) > 0.)...\n");
	idx_W = find(abs(W)>0);
	printf("Found %u connections.\n", idx_W.n_elem);
	printf("Average degree per unit = %f .\n\n", 1.0*idx_W.n_elem/N);
	printf("Computing spectral radius...\n");
	w_sR = max_eigenval(W);
	printf("sR = %4.3f\n", w_sR);
	W *= sR/w_sR;
	//printf("Normalize weights row-wise...\n");
	//W = normalize(W);
	// Save initial weight matrix
	init_W = W;
	init_W.save("./results/init_W.mat", raw_ascii);
	printf("Done.\n\n");
	mean_dw = 0.0;


	//Random initial states
	printf("Initialize states & thresholds...\n");
	leak_rate = leak;
	x.randu(N);
	xp.randu(N);
	b.randu(N);
	printf("Leaking rate = %1.3f\nSpectral radius = %1.3f\n", leak_rate, sR);
	printf("Done.\n\n");

	printf("Initialize learning parameters...\n");

	double tr = 0.1;//2.*(1.*Nu)/(1.*Ne);
	target_rate = tr*ones<vec>(N);//(2.*Nu/Ne)*ones<vec>(Ne);
	mu_IP = 0.001;//0.5;
	mu_STDP = 0.001;
	printf("Target firing rate = %1.6f\nLearning rate STDP = %1.3f\nLearning rate IP = %1.3f\n", tr, mu_STDP, mu_IP);

	STDP = false;
	SN = false;
	IP = false;
	sampling_rate = 1;
	T_mat = int(25000/sampling_rate);

	printf("Done.\n\n");
}


/*
 *
 * Destructor of the ESN class. Resets state histories and weight matrices.
 *
 * PARAM: (void)
 *
 * RETURN: (void)
 *
 */
ESN::~ESN(){
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
double ESN::max_eigenval(mat & in){
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
mat ESN::max_of(mat A, mat B){
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
mat ESN::min_of(mat A, mat B){
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


mat ESN::normalize(mat in){
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
void ESN::save_matrices(string mode){
	stringstream ss;
	ss << "./results/esn_res_st_" << mode << ".mat";
	H.save(ss.str().c_str(), raw_ascii);
	ss.str( string() );
	ss.clear();
	ss << "./results/esn_res_pst_" << mode << ".mat";
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
void ESN::set_input_con(mat in){
	Nu = in.n_rows;
	double input_conn = (1.0*Nu)/(1.*N*in.n_rows);
	if(Nu < in.n_rows)
		printf("!!! Warning: Not every input dimension is connected to the reservoir.\n");
	if(input_conn > 0.5)
		printf("!!! Warning: Reservoir is too input-driven. Increase reservoir size for respective input dimension.\n\n");
	printf("Set up input connections to the reservoir with connectivity = %f\n\n", input_conn);
	Wu = zeros<mat>(N, in.n_rows);
	ivec rand_idx = randi<ivec>(Nu, distr_param(0,N-1));
	for(int idx_col = 0; idx_col < Nu; idx_col++){
		if(idx_col%in.n_rows == 0 && Nu-idx_col<in.n_rows){
			ivec rand_idx2 = randi<ivec>(Nu-idx_col, distr_param(0,in.n_rows-1));
			Wu(rand_idx(idx_col), rand_idx2(idx_col%in.n_rows)) = 1.0;
		}
		else
			Wu(rand_idx(idx_col), idx_col%in.n_rows) = 1.0;
	}
	printf("Reservoir units driven by input = %u\n", int(accu(Wu)));
	printf("Done.\n\n");
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
void ESN::set_matrix(int time){
	printf("Set state history matrices...\n");
	if(time < T_mat){
		H.zeros(N, time);
		Hp.zeros(N, time);
	}
	else{
		H.zeros(N, T_mat);
		Hp.zeros(N, T_mat);
	}
	printf("Done.\n\n");
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
mat ESN::test(mat data, mat teacher, int time, bool trainOut){
	printf("START TESTING PHASE FOR %u TIMESTEPS\n\n", time);

	if(Wu.n_rows != N)
		set_input_con(data);
	set_matrix(time);

	for(int t = 0; t < time; t++){			// test run
		H.col(t%T_mat) = x;
		Hp.col(t%T_mat) = xp;
		update(data.col(t%data.n_cols));
		if(t%1000==0){
			printf("%6u ts: Input = %3.6f, X = %3.6f, W = %3.6f\n", t, accu(Wu*data.col(t%data.n_cols))/Wu.n_rows, accu(x)/x.n_elem, accu(W)/W.n_elem);
		}
	}
	if(trainOut)
		save_matrices("test_out");
	else
		save_matrices("test");
	No = teacher.n_rows;

	if(trainOut){
		Wout.zeros(No, N);
		mat shortHp = Hp;
		mat shortH = H;
		//	Resize state history to match teacher data (cropping the first time steps -> initLen in main.cpp)
		if((Hp.n_cols - teacher.n_cols) > 0)
			shortHp.shed_cols(0, Hp.n_cols - teacher.n_cols - 1);
		if(accu(abs(teacher)) != 0.0){
			Wout = teacher * pinv(shortHp);												// direct method: Pseudo-inverse
			//Wout = teacher * H.t() * pinv(H*H.t());									// normal functions
			//alpha = 0.7;																// regularization factor
			//Wout = teacher * H.t() * pinv(H*H.t() + alpha*alpha*ones<mat>(Ne,Ne));	// Tikhonov regularization
		}
		mat sub = (Wout*shortHp);
		mat A = sub-teacher;
		A.save("./results/dout.mat", raw_ascii);
	}
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
void ESN::update(vec in){
    vec lastx  = x;
    vec lastxp  = xp;

	// Input
    //printf("---input\n");
	vec u = Wu*in;

	// Network activation
	//printf("---states\n");
	R = W * x + b;
	xp = sign(sign(R) + 1.);
	xp = (1.-leak_rate)*lastxp + leak_rate*xp;
	x = sign(sign(R + u) + 1.);
	x = (1.-leak_rate)*lastx + leak_rate*x;		// leaky-integrated discrete-time continuous-value units

    // Plasticity
    if(IP){
    	//printf("---IP\n");
        b += mu_IP * (x-target_rate);  		// intrinsic plasticity
    }
    if (STDP){
    	//printf("---STDP\n");
		mat A = lastx * x.t();
		mat delta = mu_STDP * (A.t() - A);
		W(idx_W) += delta(idx_W);  				// additive STDP
		//W += delta;

		W = min_of(W, ones<mat>(N, N));  		// clip weights to [0,1]
		W = max_of(W, zeros<mat>(N, N));

		if(SN){
			//printf("---Norm\n");
			W = normalize(W);				// synaptic normalization
		}
    }
}
