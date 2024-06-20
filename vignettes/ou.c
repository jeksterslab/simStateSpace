#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


/**
 * This function takes a double and gives back a double
 * It computes the logistic function (i.e. the inverse of the logit link function)
 * @param x, the double value e.g. a normally distributed number
 * @return logistic(x), the double value e.g. a number between 0 and 1
 */
double mathfunction_logistic(const double x){
	double value = 1.0/(1.0 + exp(-x));
	return value;
}

/**
 * This function takes a gsl_vector and modifies its second argument (another gsl_vector)
 * It computes the softmax function (e.g. for multinomial logistic regression)
 * @param x, vector of double values e.g. a vector of normally distributed numbers
 * @param result, softmax(x), e.g. a vector of numbers between 0 and 1 that sum to 1
 */
void mathfunction_softmax(const gsl_vector *x, gsl_vector *result){
	/* Elementwise exponentiation */
	size_t index=0;
	for(index=0; index < x->size; index++){
		gsl_vector_set(result, index, exp(gsl_vector_get(x, index)));
	}
	
	/* Sum for the scaling coeficient */
	double scale = 0.0;
	for(index=0; index < x->size; index++){
		scale += gsl_vector_get(result, index);
	}
	
	/* Multiply all elements of result by 1/scale */
	gsl_blas_dscal(1/scale, result);
}


void function_dx_dt(double t, size_t regime, const gsl_vector *x, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){
	switch (regime) {
		case 0:
			gsl_vector_set(F_dx_dt, 0, (param[3] * (gsl_vector_get(x, 0) - param[0])) + (param[4] * (gsl_vector_get(x, 1) - param[1])) + (param[5] * (gsl_vector_get(x, 2) - param[2])));
			gsl_vector_set(F_dx_dt, 1, (param[6] * (gsl_vector_get(x, 0) - param[0])) + (param[7] * (gsl_vector_get(x, 1) - param[1])) + (param[8] * (gsl_vector_get(x, 2) - param[2])));
			gsl_vector_set(F_dx_dt, 2, (param[9] * (gsl_vector_get(x, 0) - param[0])) + (param[10] * (gsl_vector_get(x, 1) - param[1])) + (param[11] * (gsl_vector_get(x, 2) - param[2])));
			break;
	}
}

/**
* The dF/dx function
* The partial derivative of the jacobian of the DE function with respect to the variable x
* @param param includes at the end the current state estimates in the same order as the states following the model parameters
*/void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){
	switch (regime) {
		case 0:
			gsl_matrix_set(F_dx_dt_dx, 0, 0, param[3]);
			gsl_matrix_set(F_dx_dt_dx, 0, 1, param[4]);
			gsl_matrix_set(F_dx_dt_dx, 0, 2, param[5]);
			gsl_matrix_set(F_dx_dt_dx, 1, 0, param[6]);
			gsl_matrix_set(F_dx_dt_dx, 1, 1, param[7]);
			gsl_matrix_set(F_dx_dt_dx, 1, 2, param[8]);
			gsl_matrix_set(F_dx_dt_dx, 2, 0, param[9]);
			gsl_matrix_set(F_dx_dt_dx, 2, 1, param[10]);
			gsl_matrix_set(F_dx_dt_dx, 2, 2, param[11]);
			break;
	}
}

void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){


	gsl_matrix_set(Ht, 0, 0, 1);
	gsl_matrix_set(Ht, 1, 1, 1);
	gsl_matrix_set(Ht, 2, 2, 1);
 
	gsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);
 
}



void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){


	gsl_matrix_set(eta_noise_cov, 0, 0, param[12]);
	gsl_matrix_set(eta_noise_cov, 1, 0, param[13]);
	gsl_matrix_set(eta_noise_cov, 2, 0, param[14]);
	gsl_matrix_set(eta_noise_cov, 0, 1, param[13]);
	gsl_matrix_set(eta_noise_cov, 1, 1, param[15]);
	gsl_matrix_set(eta_noise_cov, 2, 1, param[16]);
	gsl_matrix_set(eta_noise_cov, 0, 2, param[14]);
	gsl_matrix_set(eta_noise_cov, 1, 2, param[16]);
	gsl_matrix_set(eta_noise_cov, 2, 2, param[17]);

	gsl_matrix_set(y_noise_cov, 0, 0, param[18]);
	gsl_matrix_set(y_noise_cov, 1, 1, param[19]);
	gsl_matrix_set(y_noise_cov, 2, 2, param[20]);
 
}



void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector **pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0, size_t *index_sbj){
	
	gsl_vector *Pvector = gsl_vector_calloc(1);
	gsl_vector *Pintercept = gsl_vector_calloc(1);
	gsl_vector *Padd = gsl_vector_calloc(1);
	gsl_vector *Preset = gsl_vector_calloc(1);
	gsl_vector_set(Pvector, 0, 1);
	gsl_vector_add(Padd, Pvector);
	gsl_vector_add(Padd, Pintercept);
	gsl_vector_add(Preset, Pvector);
	gsl_vector_add(Preset, Pintercept);
	gsl_vector *eta_local = gsl_vector_calloc(3);
	size_t num_regime=pr_0[0]->size;
	size_t dim_latent_var=error_cov_0[0]->size1;
	size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);
	size_t i;
	size_t regime;
	for(regime=0; regime < num_regime; regime++){
		for(i=0; i < num_sbj; i++){
			gsl_vector_set(eta_local, 0, param[21]);
			gsl_vector_set(eta_local, 1, param[22]);
			gsl_vector_set(eta_local, 2, param[23]);
			gsl_vector_set(eta_0[regime], i*dim_latent_var+0, gsl_vector_get(eta_local, 0));
			gsl_vector_set(eta_0[regime], i*dim_latent_var+1, gsl_vector_get(eta_local, 1));
			gsl_vector_set(eta_0[regime], i*dim_latent_var+2, gsl_vector_get(eta_local, 2));
			gsl_vector_set_zero(eta_local);
		}
	gsl_matrix_set((error_cov_0)[regime], 0, 0, param[24]);
	gsl_matrix_set((error_cov_0)[regime], 1, 0, param[25]);
	gsl_matrix_set((error_cov_0)[regime], 2, 0, param[26]);
	gsl_matrix_set((error_cov_0)[regime], 0, 1, param[25]);
	gsl_matrix_set((error_cov_0)[regime], 1, 1, param[27]);
	gsl_matrix_set((error_cov_0)[regime], 2, 1, param[28]);
	gsl_matrix_set((error_cov_0)[regime], 0, 2, param[26]);
	gsl_matrix_set((error_cov_0)[regime], 1, 2, param[28]);
	gsl_matrix_set((error_cov_0)[regime], 2, 2, param[29]);
	}
	for(i=0; i < num_sbj; i++){
		mathfunction_softmax(Padd, pr_0[i]);
	}
	gsl_vector_free(Pvector);
	gsl_vector_free(Pintercept);
	gsl_vector_free(Padd);
	gsl_vector_free(Preset);
	gsl_vector_free(eta_local);
}


/**
 * This function modifies some of the parameters so that it satisfies the model constraint.
 * Do not include parameters in noise_cov matrices 
 */
void function_transform(double *param){
}



void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat){
	gsl_matrix_set_identity(regime_switch_mat);
}



/**
 * The dP/dt function: depend on function_dF_dx, needs to be compiled on the user end
 * but user does not need to modify it or care about it.
 */
void mathfunction_mat_to_vec(const gsl_matrix *mat, gsl_vector *vec){
	size_t i,j;
	size_t nx=mat->size1;
	/*convert matrix to vector*/
	for(i=0; i<nx; i++){
		gsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));
		for (j=i+1;j<nx;j++){
			gsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));
			/*printf("%lu",i+j+nx-1);}*/
		}
	}
}
void mathfunction_vec_to_mat(const gsl_vector *vec, gsl_matrix *mat){
	size_t i,j;
	size_t nx=mat->size1;
	/*convert vector to matrix*/
	for(i=0; i<nx; i++){
		gsl_matrix_set(mat,i,i,gsl_vector_get(vec,i));
		for (j=i+1;j<nx;j++){
			gsl_matrix_set(mat,i,j,gsl_vector_get(vec,i+j+nx-1));
			gsl_matrix_set(mat,j,i,gsl_vector_get(vec,i+j+nx-1));
		}
	}
}
void function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dP_dt){
	
	size_t nx;
	nx = (size_t) floor(sqrt(2*(double) p->size));
	gsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);
	mathfunction_vec_to_mat(p,P_mat);
	gsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);
	function_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);
	gsl_matrix *dFP=gsl_matrix_calloc(nx,nx);
	gsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);
	gsl_matrix_transpose_memcpy(dP_dt, dFP);
	gsl_matrix_add(dP_dt, dFP);
	size_t n_Q_vec=(1+nx)*nx/2;
	gsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);
	size_t i;
	for(i=1;i<=n_Q_vec;i++){
			gsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);
	}
	gsl_matrix *Q_mat=gsl_matrix_calloc(nx,nx);
	mathfunction_vec_to_mat(Q_vec,Q_mat);
	gsl_matrix_add(dP_dt, Q_mat);
	mathfunction_mat_to_vec(dP_dt, F_dP_dt);
	gsl_matrix_free(P_mat);
	gsl_matrix_free(F_dx_dt_dx);
	gsl_matrix_free(dFP);
	gsl_matrix_free(dP_dt);
	gsl_vector_free(Q_vec);
	gsl_matrix_free(Q_mat);
}
