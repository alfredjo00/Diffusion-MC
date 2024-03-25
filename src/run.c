#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "dynamic_array.h"
#include "tools.h"
#include "diffusion_monte_carlo.h"

void run_task_1()
{
	int N_walkers_0 = 200;
	int N_time 		= 10000;
	int iterations_save = 500;

	double *E_T_array 		= (double*) malloc(N_time * sizeof(double));
	double *N_walkers_array = (double*) malloc(N_time * sizeof(double));

	Array walkers;
	diffusion_method_1D(&walkers, E_T_array, N_walkers_array, N_walkers_0,
			N_time, iterations_save);


	double **time_array = malloc(sizeof(double[N_time][2]));
	time_array[0] = E_T_array;
	time_array[1] = N_walkers_array;

	write_to_file("data/task1_E-T_N.txt", time_array, N_time, 2);

	free(E_T_array); 		E_T_array = NULL;
	free(N_walkers_array); 	N_walkers_array = NULL;
	free(time_array);		time_array = NULL;


	double **walker_array = malloc(sizeof(double[walkers.used]));

	walker_array[0] = walkers.array;

	write_to_file("data/task1_walkers.txt", walker_array, walkers.used, 1);
	free(walker_array); walker_array = NULL;

	free_array(&walkers);
}


void run_task_2()
{
	int N_walkers_0 = 1000;
	int N_time 		= 10000;
	int drift_order = 0;
	int iterations_save = 50;
	double Delta_tau = 0.01;

	double *E_T_array 		= (double*) malloc(N_time * sizeof(double));
	double *N_walkers_array = (double*) malloc(N_time * sizeof(double));

	Matrix walkers;
	diffusion_method_He(&walkers, E_T_array, N_walkers_array, N_walkers_0,
			N_time, Delta_tau, drift_order, iterations_save);

	double **time_array = malloc(sizeof(double[N_time][2]));
	time_array[0] = E_T_array;
	time_array[1] = N_walkers_array;

	write_to_file("data/task2_E-T_N.txt", time_array, N_time, 2);

	free(E_T_array); 		E_T_array = NULL;
	free(N_walkers_array); 	N_walkers_array = NULL;
	free(time_array);		time_array = NULL;


	double **walker_array = transpose_array(walkers.used, 6, walkers.matrix);

	write_to_file("data/task2_walkers.txt", walker_array, walkers.used, 6);
	free(walker_array); walker_array = NULL;

	free_matrix(&walkers);
}



void run_task_3(int drift_order)
{
	int N_walkers_0 = 1000;
	int N_time 		= 10000;
	int iterations_save = 50;
	double Delta_tau = 0.2;

	double *E_T_array 		= (double*) malloc(N_time * sizeof(double));
	double *N_walkers_array = (double*) malloc(N_time * sizeof(double));

	Matrix walkers;
	diffusion_method_He(&walkers, E_T_array, N_walkers_array, N_walkers_0,
			N_time, Delta_tau, drift_order, iterations_save);


	double **time_array = malloc(sizeof(double[N_time][2]));
	time_array[0] = E_T_array;
	time_array[1] = N_walkers_array;

	if (drift_order == 1)
		write_to_file("data/task3a_E-T_N.txt", time_array, N_time, 2);
	else
		write_to_file("data/task3b_E-T_N.txt", time_array, N_time, 2);


	free(E_T_array); 		E_T_array = NULL;
	free(N_walkers_array); 	N_walkers_array = NULL;
	free(time_array);		time_array = NULL;


	double **walker_array = transpose_array(walkers.used, 6, walkers.matrix);

	if (drift_order == 1)
		write_to_file("data/task3a_walkers.txt", walker_array, walkers.used, 6);
	else
		write_to_file("data/task3b_walkers.txt", walker_array, walkers.used, 6);

	free(walker_array); walker_array = NULL;

	free_matrix(&walkers);
}

void run_task_4()
{
	int N_walkers_0 = 1000;
	int iterations_save = 0;
	double tau_total = 1000.0;

	int n_points = 40;

	double *E_T_array 		= (double*) malloc(tau_total/0.01 * sizeof(double));
	double *N_walkers_array = (double*) malloc(tau_total/0.01 * sizeof(double));

	Matrix walkers;
	double E_T_1[n_points], E_T_2[n_points];

	double tau[n_points];

	linspace(tau, 0.01, 0.4, n_points);

	for (int i = 0; i < n_points; ++i){
		E_T_1[i] = diffusion_method_He(&walkers, E_T_array, N_walkers_array, N_walkers_0,
				(int) (tau_total / tau[i]), tau[i], 1, iterations_save);
		free_matrix(&walkers);

		E_T_2[i] = diffusion_method_He(&walkers, E_T_array, N_walkers_array, N_walkers_0,
				(int) (tau_total / tau[i]), tau[i], 2, iterations_save);
		free_matrix(&walkers);
		printf("%d \n", i);
	}
	free(E_T_array);
	free(N_walkers_array);

	double **write_array = malloc(sizeof(double[n_points][2]));
	write_array[0] = E_T_1;
	write_array[1] = E_T_2;

	write_to_file("data/task4_energy.txt", write_array, n_points, 2);

	free(write_array); write_array = NULL;
}

int run(int argc, char *argv[])
{
	//run_task_1();
	//run_task_2();
	//run_task_3(1);
	//run_task_3(2);
	run_task_4();
	return 0;
}




