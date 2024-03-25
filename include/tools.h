#ifndef INCLUDE_TOOLS_H_
#define INCLUDE_TOOLS_H_
/* **********************************************
*
* Add v1 and v2 elementwise
* results is stored in res
*
* res should be properly initialised to zero
* for this function to work correctly
*
* **********************************************/
void elementwise_addition(
		double *res,
		double *v1,
		double *v2,
		unsigned int len);

/* **********************************************
*
* Multiply v1 and v2 elementwise
* results is stored in res
*
* res should be properly initialised to zero
* for this function to work correctly
*
* **********************************************/
void elementwise_multiplication(
		double *res,
		double *v1,
		double *v2,
		unsigned int len);

/* **********************************************
*
* Calculate the dot product between
* v1 and v2
*
* the result is returned as a double
*
* **********************************************/
double dot_product(
		double *v1,
		double *v2,
		unsigned int len);

/* **********************************************
*
* Allocate the memory to a 2D array
*
* **********************************************/
void create_2D_array(
		double ***array,
		unsigned int column_size,
		unsigned int row_size);


/* **********************************************
*
* Matrix mult. (not used, not tested)
*
* **********************************************/
void matrix_multiplication(
		double **result,
		double **v1,
		double **v2,
		unsigned int m,
		unsigned int n);

/* **********************************************
*
* Returns vector norm of 1D array
*
* **********************************************/
double vector_norm(
		double *v1,
		unsigned int len);

/* **********************************************
*
* Normalizes vector input.
*
* **********************************************/
void normalize_vector(
		double *v1,
		unsigned int len);

/* **********************************************
*
* Returns arithmetic mean of 1D array.
*
* **********************************************/
double average(
		double *v1,
		unsigned int len);

/* **********************************************
*
* Returns standard deviation of 1D array.
*
* **********************************************/
double standard_deviation(
		double *v1,
		unsigned int len);


/* **********************************************
*
* Returns distance between vectors ==
* norm of the \delta vector
*
* **********************************************/
double distance_between_vectors(
		double *v1,
		double *v2,
		unsigned int len);

/* **********************************************
*
* Creates 1D array with:
* start, start + dt, ..., start + dt * (len_t-1)
*
* **********************************************/
void arange(
		double *array,
		double start,
		int len_t,
		double dt);

/* **********************************************
*
* Writes 2D array to file.
*
* **********************************************/
void write_to_file(
		char *fname,
		double **arr,
		int col_size,
		int row_size);

/* **********************************************
*
* Writes 1D array to a file
*
* **********************************************/
void write_to_file_1D_array(
		char *fname,
		double *arr,
		int N);


/* **********************************************
*
* Reads data from space/tab separated values in
* a file
*
* **********************************************/
void read_data(
		char *fname,
		double **array,
		int cols,
		int rows);

/* **********************************************
*
* 1D array of ones.
*
* **********************************************/
void array_ones(
		double *arr,
		int n);


/* **********************************************
*
* Prints 1D array to console
*
* **********************************************/
void print_array(
		double *arr,
		int n);


/* **********************************************
*
* Trapezoidal integration of function
* starting from t = 0, to (n-1) * dt
*
* **********************************************/
double trapezoidal_int(
		double *f,
		double dt,
		int n);

/* **********************************************
*
* Linear space 1D array
*
* **********************************************/
void linspace(
		double *arr,
		double a,
		double b,
		int n);


/* **********************************************
*
* Total kinetic energy of all particles
*
* **********************************************/
double kinetic_energy(
		double vel[][3],
		int n_atoms,
		double m);


/* **********************************************
*
* Velocity correlation function
* uses matrix of velocities for each particle
*
* **********************************************/
void velocity_correlation(
		double *vfc,
		double vel[][3],
		int n_times,
		int n_atoms);

/* **********************************************
*
* Mean squared displacement
*
* **********************************************/
void mean_squared_displacement(
		double *msd,
		double pos[][3],
		int n_times,
		int n_atoms);

/* **********************************************
*
* Returns a transpose of input 2D array
*
* **********************************************/
double **transpose_array(
		int col_size,
		int row_size,
		double src[][row_size]);

#endif /* INCLUDE_TOOLS_H_ */
