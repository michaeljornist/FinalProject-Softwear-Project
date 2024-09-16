# Final Project: Symmetric Non-negative Matrix Factorization (SymNMF)

This project involves implementing a clustering algorithm based on Symmetric Non-negative Matrix Factorization (SymNMF). The project compares the performance of SymNMF to K-means on various datasets and evaluates clustering quality using the silhouette score.

## Files:
- `symnmf.py`: Python interface for the project.
- `symnmf.h`: C header file for the project.
- `symnmf.c`: C implementation of the SymNMF algorithm.
- `symnmfmodule.c`: Python C API wrapper for interfacing the C code with Python.
- `analysis.py`: Script to analyze and compare SymNMF to K-means clustering using the silhouette score.
- `setup.py`: Build file for the Python C extension.
- `Makefile`: Script to build the C code and create the executable.

## Tasks:

### 1. Symmetric Non-negative Matrix Factorization (SymNMF)

The goal of SymNMF is to cluster data by minimizing the Frobenius norm between the normalized similarity matrix and the factorization matrix. The algorithm works in the following steps:

1. **Form the Similarity Matrix**: Create the similarity matrix based on the Euclidean distance between data points.
2. **Diagonal Degree Matrix**: Compute the diagonal degree matrix from the similarity matrix.
3. **Normalized Similarity Matrix**: Calculate the normalized similarity matrix using the diagonal degree matrix.
4. **Optimize H**: Find the matrix `H` that minimizes the objective function using an iterative update rule until convergence.

### 2. C Program

The C implementation calculates various matrices (similarity, diagonal degree, and normalized similarity matrix) and interacts with the Python interface to perform the clustering. The C code must compile cleanly using the provided `Makefile`.

### 3. Python C API

A Python C API is provided to wrap the C code and expose the functions (`symnmf`, `sym`, `ddg`, and `norm`) to Python. This allows seamless integration of the C-based algorithm with Python.

### 4. Analysis Script

The `analysis.py` script compares the performance of SymNMF to K-means clustering on a given dataset. It uses the `silhouette_score` from `sklearn.metrics` to evaluate the clustering quality.

### Compilation and Execution:

#### C Compilation:
To compile the C program and generate the required executables, run:

```bash
make
