#ifndef SYMNMF_H
#define SYMNMF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Constants */
#define EPSILON 1e-4
#define SIGMA 1
#define MAX_ITER 300
#define BETA 0.5

int free_matrix(double ***arr, int rows);
int alloc_matrix(double ***arr, int rows, int cols);

void print_matrix(double **matrix, int rows, int cols);

double **sym(double **data, int rows, int cols);
double **ddg(double **similarity, int size);
double **norm(double **sim_mat, double **diag_mat, int size);

double **symnmf(double **W, double **H, int N, int k);

#endif 
