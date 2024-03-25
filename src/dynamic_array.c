#include "dynamic_array.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*
 * structs defined in header file.
 */

/*
 * Buffer of memory to reduce 'realloc' calls.
 */
const int memory_buffer = 1000;

/*
 * Max size of buffer before realloc to reduce buffer, where
 * buffer size = memory_buffer * memory_buffer_ceil_factor
 */
const double memory_buffer_ceil_factor = 5.0;


/*
 * array: dynamic 1D array
 *
 * matrix: dynamic [N] x [6] (2D) array.
 */

void init_array(Array *a, size_t initial_size) {
	a->array = malloc(initial_size * sizeof(double));
	a->used = 0;
	a->size = initial_size;
}

void init_matrix(Matrix *m, size_t initial_size) {
	m->matrix = malloc(sizeof(double[initial_size][6]));
	m->used = 0;
	m->size = initial_size;
}


void insert_array(Array *a, double element) {
	if (a->size <= a->used) {
		a->size += memory_buffer;
		a->array = realloc(a->array, a->size * sizeof(double));
	}
	a->array[a->used++] = element;
}


/*
 * realloc takes a decent amount time so to speed it up
 * allocates a lot more memory then currently necessary.
 *
 */

void duplice_row_matrix(Matrix *m, int j) {
	if (m->size <= m->used) {
		m->size += memory_buffer;
		m->matrix = realloc(m->matrix, sizeof(double[m->size][6]));
	}

	assert(j < m->used);

	for (int i = 0; i < 6; ++i){
		m->matrix[m->used][i] = m->matrix[j][i];
	}
	m->used++;
}

void insert_matrix(Matrix *m, double a[6]) {

	if (m->size <= m->used) {
		m->size += memory_buffer;
		m->matrix = realloc(m->matrix, sizeof(double[m->size][6]));
	}

	for (int i = 0; i < 6; ++i){
		m->matrix[m->used][i] = a[i];
	}
	m->used++;
}


/*
 * only realloc if a lot of memory is wasted.
 */
void remove_array(Array *a, int index) {
	for(int i = index; i < a->used-1; i++){
		a->array[i] = a->array[i + 1];
	}

	a->used--;

	if (a->size >= a->used + memory_buffer_ceil_factor * memory_buffer) {
		a->size = a->used + memory_buffer;
		a->array = realloc(a->array, a->size * sizeof(double));
	}
}

void remove_matrix(Matrix *m, int index) {
	for(int i = index; i < m->used-1; i++){
		for (int j = 0; j < 6; j++){
			m->matrix[i][j] = m->matrix[i + 1][j];
		}
	}

	m->used--;

	if (m->size >= m->used + memory_buffer_ceil_factor * memory_buffer) {
		m->size = m->used + memory_buffer;
		m->matrix = realloc(m->matrix, sizeof(double[m->size][6]));
	}
}


void reset_array(Array *a, size_t reset_size) {
	if (reset_size != a->size)
		a->array = realloc(a->array, reset_size * sizeof(double));

	a->array = memset(a->array, 0, reset_size * sizeof(double));
	a->size = reset_size;
	a->used = 0;
}


void free_array(Array *a) {
	free(a->array);
	a->array = NULL;
	a->used = a->size = 0;
}


void free_matrix(Matrix *m) {
	free(m->matrix);
	m->matrix = NULL;
	m->used = m->size = 0;
}


