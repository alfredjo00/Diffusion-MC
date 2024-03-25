#ifndef INCLUDE_DIFFUSION_MONTE_CARLO_H_
#define INCLUDE_DIFFUSION_MONTE_CARLO_H_

#include "dynamic_array.h"

void diffusion_method_1D(Array *save_walkers, double *E_T_array,
		double *N_walkers_array, int N_walkers_0, int N_time,
		int iterations_save);

double diffusion_method_He(Matrix *save_walkers, double *E_T_array,
		double *N_walkers_array, int N_walkers_0, int N_time,
		double Delta_tau, int drift_order, int iterations_save);

#endif /* INCLUDE_DIFFUSION_MONTE_CARLO_H_ */
