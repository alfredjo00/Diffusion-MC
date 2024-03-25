#ifndef INCLUDE_DYNAMIC_ARRAY_H_
#define INCLUDE_DYNAMIC_ARRAY_H_


#include <stddef.h>


typedef struct Array{
	double *array;
	size_t used;
	size_t size;
} Array;


typedef struct Matrix{
	double (*matrix)[6];
	size_t used;
	size_t size;
} Matrix;


void init_array(Array *a, size_t initial_size);

void init_matrix(Matrix *m, size_t initial_size);

void insert_array(Array *a, double element);

void duplice_row_matrix(Matrix *m, int j);

void insert_matrix(Matrix *m, double a[6]);

void remove_array(Array *a, int index);

void remove_matrix(Matrix *m, int index);

void reset_array(Array *a, size_t reset_size);

void free_array(Array *a);

void free_matrix(Matrix *m);


#endif /* INCLUDE_DYNAMIC_ARRAY_H_ */
