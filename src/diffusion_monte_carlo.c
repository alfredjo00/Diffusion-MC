#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include "diffusion_monte_carlo.h"
#include "dynamic_array.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>

#define SQ(X) ((X) * (X))
#define CB(X) ((X) * (X) * (X))
#define QD(X) ((X) * (X) * (X) * (X))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define PI 3.1415926536
#define N_DIMS 3

const double alpha = 0.15;

const gsl_rng_type *T;
const int seed = 42;

double potenial_fn_1D(double x){
	return 0.5 * SQ(1.0 - exp(-x));
}


double weight_factor_1D(double x, double E_T, double Delta_tau){
	return exp(-(potenial_fn_1D(x) - E_T ) * Delta_tau);
}


double potenial_fn_He(double *R){
	double r1 = sqrt(SQ(R[0]) + SQ(R[1]) + SQ(R[2]));
	double r2 = sqrt(SQ(R[3]) + SQ(R[4]) + SQ(R[5]));
	double r12 = sqrt(SQ(R[0] - R[3]) + SQ(R[1] - R[4]) + SQ(R[2] - R[5]));
	return -2.0/r1 - 2.0/r2 + 1.0/r12;
}


double local_energy(double *R)
{
	double r1[N_DIMS] = {R[0], R[1], R[2]};
	double r2[N_DIMS] = {R[3], R[4], R[5]};

	double r1_unit[N_DIMS], r2_unit[N_DIMS];

	memcpy(r1_unit, r1, sizeof(double[N_DIMS]));
	memcpy(r2_unit, r2, sizeof(double[N_DIMS]));

	normalize_vector(r1_unit, N_DIMS);
	normalize_vector(r2_unit, N_DIMS);

	double r12 = distance_between_vectors(r1, r2, N_DIMS);

	double dot_product = 0.0;

	for(int j = 0; j < N_DIMS; ++j){
		dot_product += (r1_unit[j] - r2_unit[j]) * (r1[j] - r2[j]);
	}

	return -4.0 + dot_product / (r12 * SQ(1.0 + alpha * r12))
					   	- 1.0 / (r12 * CB(1.0 + alpha * r12))
					   	- 1.0 / (4.  * QD(1.0 + alpha * r12))
					   	+ 1.0 / r12;
}


double weight_factor_potential(double *R, double E_T, double Delta_tau){
	return exp(-(potenial_fn_He(R) - E_T ) * Delta_tau);
}


double weight_factor_local_energy(double *R, double E_T, double Delta_tau){
	return exp(-(local_energy(R) - E_T ) * Delta_tau);
}



void drift_v_f(double *R, double *v_f) {
	double r1[N_DIMS] = {R[0], R[1], R[2]};
	double r2[N_DIMS] = {R[3], R[4], R[5]};
	double r12[N_DIMS] = {r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]};

	double r1_unit[N_DIMS], r2_unit[N_DIMS], r12_unit[N_DIMS];

	memcpy(r1_unit, r1, sizeof(double[N_DIMS]));
	memcpy(r2_unit, r2, sizeof(double[N_DIMS]));
	memcpy(r12_unit, r12, sizeof(double[N_DIMS]));

	normalize_vector(r1_unit, N_DIMS);
	normalize_vector(r2_unit, N_DIMS);
	normalize_vector(r12_unit, N_DIMS);

	double r12_norm = vector_norm(r12, N_DIMS);

	v_f[0] = -2.0 * r1_unit[0] - 1/(2. * SQ(1. + alpha * r12_norm)) * r12_unit[0];
	v_f[1] = -2.0 * r1_unit[1] - 1/(2. * SQ(1. + alpha * r12_norm)) * r12_unit[1];
	v_f[2] = -2.0 * r1_unit[2] - 1/(2. * SQ(1. + alpha * r12_norm)) * r12_unit[2];

	v_f[3] = -2.0 * r2_unit[0] - 1/(2. * SQ(1. + alpha * r12_norm)) * r12_unit[0];
	v_f[4] = -2.0 * r2_unit[1] - 1/(2. * SQ(1. + alpha * r12_norm)) * r12_unit[1];
	v_f[5] = -2.0 * r2_unit[2] - 1/(2. * SQ(1. + alpha * r12_norm)) * r12_unit[2];
}


void diffusion_method_1D(Array *save_walkers, double *E_T_array, double *N_walkers_array,
		int N_walkers_0, int N_time, int iterations_save){

	/*
	 * Initialise gsl rng
	 */
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);


	double gamma 		= 0.5;
	double x_limit 		= 5;
	double Delta_tau 	= 0.02;
	double E_T 			= 0.5;
	double weight_x;

	int N_walkers = N_walkers_0;
	int min_walkers = 10;
	int eq_steps = 100;
	int m;

	init_array(save_walkers, (size_t) 2 * iterations_save * N_walkers_0);

	Array array_walkers;
	init_array(&array_walkers, (size_t) N_walkers_0);
	linspace(array_walkers.array, -x_limit, x_limit, N_walkers_0);
	array_walkers.used = N_walkers_0;

	/*
	 * Array to store indices of walkers to be destroyed
	 * at the end of the loop.
	 */
	Array remove_elements;
	size_t rm_elements_size = 1000;
	init_array(&remove_elements, rm_elements_size);
	double E_T_sum = 0;

	for (int i = 0; i < N_time; ++i){
		for (int j = 0; j < N_walkers; ++j){

			array_walkers.array[j] +=  sqrt(Delta_tau) * gsl_ran_gaussian(r, 1.0);
			weight_x = weight_factor_1D(array_walkers.array[j], E_T, Delta_tau);

			m = (int) (weight_x + gsl_rng_uniform(r));

			if (m == 0) {
				insert_array(&remove_elements, j);
			}

			else if (m > 1){
				for (int k = 0; k < m-1; ++k){
					insert_array(&array_walkers, array_walkers.array[j]);
				}
			}

			if (N_time - i <= iterations_save)
				insert_array(save_walkers, array_walkers.array[j]);

			if (j == N_walkers - 1){
				for (int l = 0; l < remove_elements.used; ++l){
					if (array_walkers.used > min_walkers)
						remove_array(&array_walkers, remove_elements.array[remove_elements.used - 1 - l]);
				}
				N_walkers = array_walkers.used;
				reset_array(&remove_elements, rm_elements_size);

				if (i < eq_steps)
					E_T += -gamma * log(((double)N_walkers) / N_walkers_0 );

				else {
					E_T_sum += E_T;
					E_T = E_T_sum / (i - eq_steps + 1) - gamma * log(((double)N_walkers) / N_walkers_0);
				}

				E_T_array[i] 		= E_T;
				N_walkers_array[i] 	= (double) N_walkers;

				break;
			}
		}
	}
	gsl_rng_free(r);
	free_array(&remove_elements);
	free_array(&array_walkers);
}


void initialize_walkers_He(Matrix *m_walkers, int N_walkers_0)
{
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	for (int i = 0; i < N_walkers_0; ++i){
		double r1 		= 0.7 + gsl_rng_uniform(r);
		double theta1 	= acos(2.0 * gsl_rng_uniform(r) - 1.0);
		double phi1 	= 2.0 * PI * gsl_rng_uniform(r);

		double r2 		= 0.7 + gsl_rng_uniform(r);
		double theta2 	= acos(2.0 * gsl_rng_uniform(r) - 1.0);
		double phi2 	= 2.0 * PI * gsl_rng_uniform(r);

		m_walkers->matrix[i][0] = r1 * sin(theta1) * cos(phi1);
		m_walkers->matrix[i][1] = r1 * sin(theta1) * sin(phi1);
		m_walkers->matrix[i][2] = r1 * cos(theta1);

		m_walkers->matrix[i][3] = r2 * sin(theta2) * cos(phi2);
		m_walkers->matrix[i][4] = r2 * sin(theta2) * sin(phi2);
		m_walkers->matrix[i][5] = r2 * cos(theta2);
	}
	gsl_rng_free(r);
}


/*
 * Diffusion method for 2 particles in 3 dims, with importance sampling if drift_order > 0
 */
double diffusion_method_He(Matrix *save_walkers, double *E_T_array,
		double *N_walkers_array, int N_walkers_0, int N_time,
		double Delta_tau, int drift_order, int iterations_save)
{
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	double gamma = 0.5;
	double E_T 	= -2.916;
	double weight_R;
	double v_f[2 * N_DIMS];
	double R1_2[2 * N_DIMS];

	int N_walkers = N_walkers_0;
	int min_walkers = 10;
	int eq_steps = 100;
	int eq_sum = (int) N_time / 10;

	/*
	 * Matrix struct to store walkers.
	 * Dynamically allocates and reallocates memory
	 * based on number of walkers.
	 * Matrix->matrix is a [N][6] dim matrix,
	 * with x1,y1,z1,x2,y2,z2 for all walkers.
	 */
	init_matrix(save_walkers, (size_t) 2 * iterations_save * N_walkers_0);

	Matrix m_walkers;
	init_matrix(&m_walkers, (size_t) (2 * N_walkers_0));
	initialize_walkers_He(&m_walkers, N_walkers_0);
	m_walkers.used = N_walkers_0;

	/*
	 * Array to store indices of the walkers to be destroyed.
	 */
	Array remove_elements;
	size_t rm_elements_size = 10;
	init_array(&remove_elements, rm_elements_size);
	double E_T_mean = 0;
	double E_T_sum = 0;

	/*
	 * Iterates through all walkers for all time steps.
	 *
	 * Does not update number of walkers added or destroyed until the
	 * end of the loop for that time step.
	 *
	 * The mean of the oscillating E_T approaches the ground state energy
	 */

	if (drift_order < 2)
		(void) R1_2;

	for (int i = 0; i < N_time; ++i){
		for (int j = 0; j < N_walkers; ++j){

			if (drift_order >= 1)
				drift_v_f(m_walkers.matrix[j], v_f);

			if (drift_order == 2){
				for (int k = 0; k < 2 * N_DIMS; ++k){
					R1_2[k] =  v_f[k] * Delta_tau * 0.5 + m_walkers.matrix[j][k];
				}
				drift_v_f(R1_2, v_f);
			}

			for (int k = 0; k < 2 * N_DIMS; ++k){
				m_walkers.matrix[j][k] += sqrt(Delta_tau) * gsl_ran_gaussian(r, 1.0);

				if (drift_order >= 1)
					m_walkers.matrix[j][k] += Delta_tau * v_f[k];

			}

			if (drift_order >= 1)
				weight_R = weight_factor_local_energy(m_walkers.matrix[j], E_T, Delta_tau);
			else
				weight_R = weight_factor_potential(m_walkers.matrix[j], E_T, Delta_tau);


			int m = (int) (weight_R + gsl_rng_uniform(r));

			if (m == 0)
				insert_array(&remove_elements, j);

			else if (m > 1){
				for (int k = 0; k < m-1; ++k){
					duplice_row_matrix(&m_walkers, j);
				}
			}

			if (N_time - i <= iterations_save)
				insert_matrix(save_walkers, m_walkers.matrix[j]);

			if (j == N_walkers - 1){
				for (int l = 0; l < remove_elements.used; ++l){
					if (m_walkers.used > min_walkers)
						remove_matrix(&m_walkers, remove_elements.array[remove_elements.used - 1 - l]);
				}

				N_walkers = m_walkers.used;
				reset_array(&remove_elements, rm_elements_size);

				if (i < eq_steps)
					E_T += -gamma * log(((double)N_walkers) / N_walkers_0);

				else {
					E_T_mean = E_T_mean * (i - eq_steps)/(i - eq_steps + 1.0) + E_T / (i - eq_steps + 1.0);

					E_T = E_T_mean - gamma * log(((double)N_walkers) / N_walkers_0);
				}

				E_T_array[i] 		= E_T;
				N_walkers_array[i] 	= (double) N_walkers;
				if (i >= eq_sum){
					E_T_sum += E_T;
				}

				//printf("i %d \t N %d \t E %.12f\n", i, N_walkers, E_T);
				break;
			}
		}
	}
	gsl_rng_free(r);
	free_array(&remove_elements);
	free_matrix(&m_walkers);
	return E_T_sum / (N_time - eq_sum);
}

