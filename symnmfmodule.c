#define PY_SSIZE_T_CLEAN

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Python.h"
#include "symnmf.h"


/*
    Parses a Python object - List of Lists to a 2D array in C
    Allocate space for the 2D array and update rows and cols
    on success returns 0.
 */
int py_ListOfList_parser(PyObject **py_LisfOfLists, double ***c_matrix, int *rows, int *cols) {
    int i, j;
    PyObject *py_List, *element;

    *rows = PyObject_Length(*py_LisfOfLists);

    /* Check if the py list is empty or some error equired */
    if (*rows == -1 || *rows == 0) return 1;

    *cols = PyObject_Length(PyList_GetItem(*py_LisfOfLists, 0));

    /* Check for errors in parsing the rows */
    if (*cols == -1) return 1;

    if (alloc_matrix(c_matrix, *rows, *cols) != 0)
        return 1;

    for (i = 0; i < *rows; i++) {
        py_List = PyList_GetItem(*py_LisfOfLists, i);
        /*the passed list from the user should be with even size of rows*/
        if (*cols != PyObject_Length(py_List)) {
            alloc_matrix(c_matrix, *rows,*rows); /*just to free is seccsessfuly after*/
            return 1;
        }
        for (j = 0; j < *cols; j++) {
            element = PyList_GetItem(py_List, j);
            (*c_matrix)[i][j] = PyFloat_AsDouble(element);
        }
    }
    return 0;
}


/*
    Parses a C object -2D array to a List of Lists in Python
    on success returns 0.
 */
int c_Matrix_parser(PyObject **py_LisfOfLists, double ***c_matrix, const int rows, const int cols) {
    PyObject *py_List, *element;
    int i, j;

    *py_LisfOfLists = PyList_New(rows);
    if (*py_LisfOfLists == NULL) return 1;

    for (i = 0; i < rows; i++) {
        py_List = PyList_New(cols);
        if (py_List == NULL)
            return 1;
        for (j = 0; j < cols; j++) {
            element = PyFloat_FromDouble((*c_matrix)[i][j]); /* a double in C is a float in Py*/
            if (element == NULL) return 1;
            PyList_SET_ITEM(py_List, j, element);
        }
        PyList_SET_ITEM(*py_LisfOfLists, i, py_List);
    }
    return 0;
}



/*
    The py_sym function in this Python C extension parses input data from Python to C,
    computes a similarity matrix using the C function sym, and then converts the results back to Python.
    If any error occurs during processing, it cleans up the allocated memory and returns NULL.
*/
static PyObject *py_sym(PyObject *self, PyObject *args) {
    PyObject *py_data, *py_sim_mat;
    double **c_data, **c_sim_mat;
    int rows, cols;

    if (!PyArg_ParseTuple(args, "O", &py_data))
        return NULL;

    if (py_ListOfList_parser(&py_data, &c_data, &rows, &cols) != 0)
        return NULL;
    
    c_sim_mat = sym(c_data,rows,cols);

    if (c_sim_mat == NULL){
        free_matrix(&c_data, rows);
        return NULL;
    }
        

    if (c_Matrix_parser(&py_sim_mat, &c_sim_mat, rows, rows) != 0) {
        free_matrix(&c_data, rows);
        free_matrix(&c_sim_mat, rows);
        return NULL;
    }
    free_matrix(&c_data, rows);
    free_matrix(&c_sim_mat, rows);
    return py_sim_mat;
}


/*
    computes the similarity matrix using sym, calculates the diagonal degree matrix with ddg,
    and then converts the result back to a Python object. If an error occurs during any stage,
    it ensures that all allocated memory is freed before returning NULL.
*/
static PyObject *py_ddg(PyObject *self, PyObject *args) {
    PyObject *py_data, *py_diag_mat;
    double **c_data, **c_sim_mat, **c_diag_mat;
    int rows, cols;
    py_diag_mat = NULL;
    if (!PyArg_ParseTuple(args, "O", &py_data))
        return NULL;

    if (py_ListOfList_parser(&py_data, &c_data, &rows, &cols) != 0)
        return NULL;

    c_sim_mat = sym(c_data,rows,cols);

    if(c_sim_mat == NULL){
        free_matrix(&c_data,rows);
        return NULL;
    }

    c_diag_mat = ddg(c_sim_mat,rows);

    if(c_diag_mat == NULL){
        free_matrix(&c_data,rows);
        free_matrix(&c_sim_mat,rows);
        return NULL;
    }

    if(c_Matrix_parser(&py_diag_mat,&c_diag_mat,rows, rows)!=0){
        free_matrix(&c_data,rows);
        free_matrix(&c_sim_mat,rows);
        free_matrix(&c_diag_mat,rows);
        return NULL;
    }

    free_matrix(&c_data, rows);
    free_matrix(&c_sim_mat, rows);
    free_matrix(&c_diag_mat, rows);

    return py_diag_mat;
}


/*
computes the normalized matrix using the similarity and diagonal degree matrices,
converts the result back to a Python object,and manages memory cleanup in case of errors.
*/
static PyObject *py_norm(PyObject *self, PyObject *args) {
    PyObject *py_data, *py_norm_mat;
    double **c_data, **c_sim_mat, **c_diag_mat, **c_norm_mat;
    int rows, cols;
    py_norm_mat = NULL;

    if (!PyArg_ParseTuple(args, "O", &py_data)) return NULL;

    if (py_ListOfList_parser(&py_data, &c_data, &rows, &cols) != 0) return NULL;

    c_sim_mat = sym(c_data,rows,cols);

    if(c_sim_mat == NULL){
        free_matrix(&c_data,rows);
        return NULL;
    }

    c_diag_mat = ddg(c_sim_mat,rows);

    if(c_diag_mat == NULL){
        free_matrix(&c_data,rows);
        free_matrix(&c_sim_mat,rows);
        return NULL;
    }

    c_norm_mat = norm(c_sim_mat,c_diag_mat,rows);

    if (c_norm_mat == NULL){
        free_matrix(&c_data,rows);
        free_matrix(&c_sim_mat,rows);
        free_matrix(&c_diag_mat,rows);
        return NULL;
    }
        

    if(c_Matrix_parser(&py_norm_mat, &c_norm_mat, rows, rows) != 0){
        free_matrix(&c_data, rows);
        free_matrix(&c_sim_mat,rows);
        free_matrix(&c_diag_mat,rows);
        free_matrix(&c_norm_mat, rows);
        return NULL;
    }

    free_matrix(&c_data, rows);
    free_matrix(&c_sim_mat,rows);
    free_matrix(&c_diag_mat,rows);
    free_matrix(&c_norm_mat, rows);
    return py_norm_mat;
}


/*
computes the symnmf matrix ,it takes two inputs: the initial matrix H and the normalized matrix, both in Python object form.
It parses these into C arrays, computes the SymNMF,and converts the resulting matrix back into a Python object. 
*/
static PyObject *py_symnmf(PyObject *self, PyObject *args) {
    PyObject *py_init_H, *py_final_H, *py_norm_mat;
    double **c_init_H, **c_final_H, **c_norm_mat;
    int N, k;
    py_final_H = NULL;

    if (!PyArg_ParseTuple(args, "OO", &py_init_H, &py_norm_mat))
        return NULL;

    if (py_ListOfList_parser(&py_init_H, &c_init_H, &N, &k) != 0)
        return NULL;
    if (py_ListOfList_parser(&py_norm_mat, &c_norm_mat, &N, &N) != 0)
        return NULL;
    
    c_final_H = symnmf(c_norm_mat,c_init_H,N,k);
    
    if (c_final_H ==NULL) {
        free_matrix(&c_init_H, N);
        free_matrix(&c_norm_mat, N);
        return NULL;
    }

    if (c_Matrix_parser(&py_final_H, &c_final_H, N, k) != 0)
        return NULL;

    free_matrix(&c_init_H, N);
    free_matrix(&c_norm_mat, N);

    return py_final_H;
}


static PyMethodDef symnmfMethods[] = {
        {"symnmf",
                (PyCFunction) py_symnmf,
                     METH_VARARGS,
                PyDoc_STR("Returns the final H matrix")},
        {"sym",
                (PyCFunction) py_sym,
                     METH_VARARGS,
                PyDoc_STR("Returns the similarity matrix")},
        {"ddg",
                (PyCFunction) py_ddg,
                     METH_VARARGS,
                PyDoc_STR("Returns the Diagonal Degree Matrix")},
        {"norm",
                (PyCFunction) py_norm,
                     METH_VARARGS,
                PyDoc_STR("Returns the normalized similarity matrix")},

        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
        PyModuleDef_HEAD_INIT,
        "symnmfmodule",
        NULL,
        -1,
        symnmfMethods
};

PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    return PyModule_Create(&symnmfmodule);
}