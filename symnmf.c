#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "symnmf.h"


#define EPSILON 1e-4
#define SIGMA 1
#define MAX_ITER 300
#define BETA 0.5


/**
 * Frees memory allocated for a 2D array of doubles.
 * 
 * @param arr A pointer to the array to be freed.
 * @param rows The number of rows in the array.
 * @return Returns 0 on successful freeing of memory, or another integer if an error occurs.
 */
int free_matrix(double ***arr, int rows) {
    int i;
    if (arr != NULL && *arr != NULL) {
        for (i = 0; i < rows; i++) {
            free((*arr)[i]); 
        }
        free(*arr); 
        *arr = NULL;  
    }
    return 0;
}


/**
 * Copies the content of one matrix to another.
 * Both matrices must already be allocated and have the same dimensions.
 *
 * @param src The source matrix to copy from.
 * @param dest The destination matrix to copy to.
 * @param rows The number of rows in the matrices.
 * @param cols The number of columns in the matrices.
 */
void copy_matrix(double **src, double **dest, int rows, int cols) {
    int i,j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

/**
 * Allocates memory for a 2D array (matrix) of doubles.
 * 
 * @param arr Pointer to the 2D array to be allocated.
 * @param rows The number of rows in the matrix.
 * @param cols The number of columns in the matrix.
 * @return Returns 0 if allocation is successful, 1 if an error occurs.
 */
int alloc_matrix(double ***matrix, int rows, int cols)
{
    int i;
    *matrix = (double **)malloc(rows * sizeof(double *));
    if (*matrix == NULL)
        return 1;  

    for (i = 0; i < rows; i++) { 
        (*matrix)[i] = (double *)calloc(cols, sizeof(double));
        if ((*matrix)[i] == NULL) {
            free_matrix(matrix, i);  
            return 1;  
        }
    }
    return 0;  
}

/**
 * Reads a file that represent the data and initialize the rows and columns.
 * 
 * @param filename The path to the file.
 * @param rows Pointer to store the number of rows.
 * @param cols Pointer to store the number of columns.
 * @return Returns 0 on success, or 1 if an error occurs (e.g., file not found).
 */
int get_matrix_dimensions(const char *filename, int *rows, int *cols) {
    FILE *file;
    char buffer[1024];
    int column_count;
    int row_count;
    char *line, *token;

    column_count = 0;
    row_count = 0;
    file = fopen(filename, "r");
    if (!file) {
        return 1;
    }

    

    if (fgets(buffer, sizeof(buffer), file)) {
        line = buffer;
        token = strtok(line, " ,\t\n");
        while (token != NULL) {
            column_count++;
            token = strtok(NULL, " ,\t\n");
        }
        row_count++;
    }

    while (fgets(buffer, sizeof(buffer), file)) {
        row_count++;
    }

    *rows = row_count;
    *cols = column_count;

    fclose(file);
    return 0;
}



/**
 * Transposes a given matrix, swapping its rows with its columns.
 * This function dynamically allocates memory for the transposed matrix and then fills it by swapping the indices of the original matrix.
 * 
 * @param matrix A pointer to the 2D array of doubles representing the matrix to transpose.
 * @param rows The number of rows in the original matrix.
 * @param cols The number of columns in the original matrix.
 * @return A pointer to the newly allocated transposed 2D array of doubles, or NULL if memory allocation fails.
 */
double** transpose(double **matrix, int rows, int cols) {
    double **transposed_matrix;
    int i,j;

    if (alloc_matrix(&transposed_matrix,cols, rows)==1) return NULL;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            transposed_matrix[j][i] = matrix[i][j];
        }
    }
    return transposed_matrix;
}


/**
 * Multiplies two matrices A and B, and returns the resulting matrix.
 * This function dynamically allocates memory for the result matrix and initializes its elements to zero before performing the multiplication.
 * 
 * @param A A pointer to the 2D array of doubles representing the first matrix (left operand).
 * @param B A pointer to the 2D array of doubles representing the second matrix (right operand).
 * @param A_rows The number of rows in matrix A.
 * @param A_cols The number of columns in matrix A (should match the number of rows in matrix B).
 * @param B_cols The number of columns in matrix B.
 * @return A pointer to the resulting 2D array of doubles after matrix multiplication, or NULL if memory allocation fails.
 */
double** matrix_multiply(double **A, double **B, int A_rows, int A_cols, int B_cols) {
    double ** AB;
    int i,j,k;
    if (alloc_matrix(&AB,A_rows,B_cols) == 1) return NULL;

    for (i = 0; i < A_rows; i++) {
        for (j = 0; j < B_cols; j++) {
            AB[i][j] = 0.0; 
            for (k = 0; k < A_cols; k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }

        }
    }
    return AB;
}

/**
 * Reads a matrix from a file and stores it in a 2D array.
 * 
 * @param filename The path to the file containing the matrix.
 * @param matrix Pointer to the 2D array where the matrix will be stored.
 * @param rows Pointer to the number of rows in the matrix.
 * @param cols Pointer to the number of columns in the matrix.
 * @return Returns 0 on successful reading, or 1 if an error occurs.
 */
int read_matrix_from_file(const char *filename, double ***matrix, int *rows, int *cols) {
    FILE *file;
    int i,j,result;
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file");
        return 1;  
    }

    if (alloc_matrix(matrix, *rows, *cols) != 0) {
        fclose(file);
        return 1;  
    }

    for (i = 0; i < *rows; i++) {
        for (j = 0; j < *cols; j++) {
            result = (j == *cols - 1) ? fscanf(file, "%lf", &(*matrix)[i][j]) 
                                          : fscanf(file, "%lf,", &(*matrix)[i][j]);
            if (result != 1) {
                fprintf(stderr, "Error reading matrix data at row %d, column %d.\n", i, j);
                free_matrix(matrix, i + 1);
                fclose(file);
                return 1;
            }
        }
    }

    fclose(file);
    return 0;
}


/**
 * Prints the contents of a 2D array (matrix) to the standard output, formatting each element to four decimal places.
 * Each row of the matrix is printed on a new line, with elements in the row separated by commas.
 *
 * @param matrix A pointer to the 2D array of doubles representing the matrix to be printed.
 * @param rows The number of rows in the matrix.
 * @param cols The number of columns in the matrix.
 */
void print_matrix(double **matrix, int rows, int cols) {
    int i,j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < cols - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}


/**
 * Computes the similarity matrix for a given dataset using a Gaussian kernel.
 * Allocates memory for the similarity matrix and computes the Gaussian similarity values.
 *
 * @param data Pointer to the 2D array of input data.
 * @param rows The number of rows in the data array.
 * @param cols The number of columns in each row of the data array.
 * @return Pointer to the allocated and computed similarity matrix, or NULL if memory allocation fails.
 */
double** sym(double **data, int rows, int cols) {
    int i,j,k;
    double sigma_sq,diff,sum;
    double **similarity;
    if (alloc_matrix(&similarity, rows, rows) != 0) {
        return NULL; 
    }

    sigma_sq = SIGMA * SIGMA;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < rows; j++) {
            if (i == j) {
                similarity[i][j] = 0.0; 
            } else {
                sum = 0.0;
                for (k = 0; k < cols; k++) {
                    diff = data[i][k] - data[j][k];
                    sum += diff * diff;
                }
                similarity[i][j] = exp(-sum / (2 * sigma_sq)); 
            }
        }
    }
    return similarity; 
}


/**
 * Computes the Diagonal Degree Matrix from a given similarity matrix.
 * Allocates memory for the degree matrix and fills its diagonal with the sum of the corresponding rows in the similarity matrix.
 * 
 * @param similarity A pointer to the 2D array of doubles representing the similarity matrix.
 * @param size The number of rows (and columns, since the matrix is square) of the similarity matrix.
 * @return Pointer to the allocated and filled degree matrix, or NULL if memory allocation fails.
 */
double** ddg(double **similarity, int size) {
    double **degree_matrix;
    int i,j;
    if (alloc_matrix(&degree_matrix, size, size) != 0) {
        return NULL; 
    }

    for (i = 0; i < size; i++) {
        double row_sum = 0.0;
        for (j = 0; j < size; j++) {
            row_sum += similarity[i][j];
        }
        degree_matrix[i][i] = row_sum;
    }
    return degree_matrix; 
}


/**
 * Normalizes the similarity matrix using the diagonal degree matrix and allocates memory for the normalized matrix.
 * 
 * @param sim_mat The similarity matrix.
 * @param diag_mat The diagonal degree matrix.
 * @param size The size of the matrices (number of rows or columns, as the matrices are square).
 * @return Pointer to the allocated and normalized matrix, or NULL if an error occurs (e.g., memory allocation fails).
 */
double** norm(double **sim_mat, double **diag_mat, int size) {
    double **normalized_mat;
    double *D_inv_sqrt;
    int i,j;

    if (alloc_matrix(&normalized_mat, size, size) != 0) {
        return NULL; 
    }

    D_inv_sqrt = (double *)malloc(size * sizeof(double));
    if (D_inv_sqrt == NULL) {
        free_matrix(&normalized_mat, size); 
        return NULL;
    }

    for (i = 0; i < size; i++) {
        if (diag_mat[i][i] > 0)
            D_inv_sqrt[i] = 1.0 / sqrt(diag_mat[i][i]);
        else
            D_inv_sqrt[i] = 0.0;  
    }

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            normalized_mat[i][j] = D_inv_sqrt[i] * sim_mat[i][j] * D_inv_sqrt[j];
        }
    }
    free(D_inv_sqrt);

    return normalized_mat; 
}

/**
 * Calculates the squared Frobenius norm of the difference between two matrices A and B.
 * This function is used to compute the sum of the squares of the element-wise differences between two matrices,
 * which is a measure of the "distance" between them.
 *
 * @param A A pointer to the first matrix (double**), where each element is a double.
 * @param B A pointer to the second matrix (double**), where each element is a double.
 * @param rows The number of rows in both matrices A and B.
 * @param cols The number of columns in both matrices A and B.
 * @return The squared Frobenius norm (double) of the difference between matrices A and B.
 */

double squred_f_norm(double **A, double **B, int rows, int cols) {
    int i,j;
    double sum,diff; 

    sum = 0.0;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            diff = A[i][j] - B[i][j];
            sum += diff * diff;
        }
    }
    return sum;
}


/**
 * Updates the matrix H in-place to minimize the difference between W and the product of H*H^T during Symmetric Non-negative Matrix Factorization (SymNMF).
 * This function performs one iteration of the update process defined in the SymNMF algorithm, which includes matrix multiplications and element-wise updates of H.
 *
 * @param H A pointer to the 2D array of doubles representing the matrix H, which is updated in-place.
 * @param W A pointer to the 2D array of doubles representing the weight matrix W.
 * @param N The number of rows in matrices W and H.
 * @param k The number of columns in matrix H (feature dimension).
 * @return Returns 0 if the update is successful, 1 if an error occurs during matrix operations (e.g., memory allocation failures).
 */

int update_H(double **H, double **W, int N, int k) {
    double **WH ,**Ht,**Ht_H,**H_Ht_H;
    int i,j;
    WH = matrix_multiply(W, H, N, N, k);
    if (WH == NULL) {
        return 1;
    }

    Ht = transpose(H, N, k);
    if (Ht == NULL) {
        free_matrix(&WH, N);
        return 1;
    }

    Ht_H = matrix_multiply(Ht, H, k, N, k);
    if (Ht_H == NULL) {
        free_matrix(&WH, N);
        free_matrix(&Ht, k);
        return 1;
    }

    H_Ht_H = matrix_multiply(H, Ht_H, N, k, k);
    if (H_Ht_H == NULL) {
        free_matrix(&WH, N);
        free_matrix(&Ht, k);
        free_matrix(&Ht_H, k);
        return 1;
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < k; j++) {
            if (H_Ht_H[i][j] != 0)  
                H[i][j] = H[i][j] * (1 - BETA + BETA * (WH[i][j] / H_Ht_H[i][j]));
        }
    }

    free_matrix(&WH, N);
    free_matrix(&Ht, k);
    free_matrix(&Ht_H, k);
    free_matrix(&H_Ht_H, N);
    return 0;
}


/**
 * Performs Symmetric Non-negative Matrix Factorization (SymNMF) on the given weight matrix W using initial matrix H.
 * It iteratively updates matrix H to minimize the difference between W and H*H^T, based on the Frobenius norm.
 *
 * @param W A pointer to the 2D array of doubles representing the weight matrix.
 * @param H A pointer to the 2D array of doubles representing the initial matrix, which is modified in-place.
 * @param N The number of rows (and columns, since W and H are square) in matrices W and H.
 * @param k The number of columns in matrix H (feature dimension).
 * @return A pointer to the 2D array of doubles representing the updated matrix H after convergence, or NULL if an error occurs during the process.
 */
double** symnmf(double **W,double **H,int N,int k) {
    double **prev_H;
    int i;

    if (alloc_matrix(&prev_H,N,k) == 1) return NULL;

    for (i = 0; i < MAX_ITER; i++)
    {
        copy_matrix(H, prev_H, N, k); 
        if(update_H(H,W,N,k) == 1) return NULL;
        if(squred_f_norm(prev_H,H,N,k) < EPSILON) break;

    }
    free_matrix(&prev_H, N);
    return H;
}


int main(int argc, char *argv[]) {
    double **data , **similarity , **deg_matrix , **normalized_matrix;
    int rows, cols;
    char *str_sym ,*str_norm, *str_ddg ,*op;
    data = NULL;
    deg_matrix = NULL;  
    normalized_matrix = NULL;

    str_sym = "sym";
    str_norm = "norm";
    str_ddg = "ddg";
    if (argc != 3) {
        printf("An Error Has Occurred");
        return 1;
    }
    op = argv[1];
    if((strcmp(str_sym , op) != 0 )&&(strcmp(str_ddg , op) != 0 )&&(strcmp(str_norm , op) != 0 )){
        printf("An Error Has Occurred");
        return 1;
    }

    if(get_matrix_dimensions(argv[2], &rows, &cols) != 0){
        printf("An Error Has Occurred");
        return 1;
    }
    if (read_matrix_from_file(argv[2], &data, &rows, &cols) != 0) {
        printf("An Error Has Occurred");
        return 1;
    }

    


    similarity = sym(data, rows, cols);
    if(similarity == NULL){
        printf("An Error Has Occurred");
        free_matrix(&data,rows);
        return 1;
    }

    /* ############### SYM handler ############### */
    if(strcmp(str_sym , op) == 0){
        
        free_matrix(&data,rows);
        print_matrix(similarity,rows,rows);
        free_matrix(&similarity,rows);
        return 0;
    }
    /* ########################################## */

    deg_matrix = ddg(similarity, rows);
    if(deg_matrix == NULL){
        printf("An Error Has Occurred");
        free_matrix(&data,rows);
        free_matrix(&similarity,rows);
        return 1;
    }

    /* ############### DDG handler ############### */
    if (strcmp(str_ddg , op) == 0)
    {
        free_matrix(&data,rows);
        free_matrix(&similarity,rows);
        print_matrix(deg_matrix,rows,rows);
        free_matrix(&deg_matrix,rows);
        return 0;
    }
    /* ########################################## */

    normalized_matrix = norm(similarity,deg_matrix,rows);
    if(normalized_matrix == NULL){
        printf("An Error Has Occurred");
        free_matrix(&data,rows);
        free_matrix(&similarity,rows);
        free_matrix(&deg_matrix,rows);
        return 1;
    }
    /* ############### NORM handler ############### */
    if (strcmp(str_norm , op) == 0)
    {   
        free_matrix(&data,rows);
        free_matrix(&similarity,rows);
        free_matrix(&deg_matrix,rows);
        print_matrix(normalized_matrix,rows,rows);
        free_matrix(&normalized_matrix,rows);
        return 0;
    }
    /* ########################################## */
    return 0;
}
