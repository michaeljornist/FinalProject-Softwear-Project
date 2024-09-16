import symnmfmodule
import sys
import numpy as np

def printMat(mat):
    """
    prints the matrix in the specified form
    :param mat: matrix of floats to be printed
    :type mat: list of lists
    :return: void
    """
    for row in mat:
        formattedRow = ["%.4f" % num for num in row]
        print(*formattedRow, sep=",")


def read_file(filename):
    try:
        # Load the data without assuming the delimiter, this lets us check the raw content
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Determine if there's more than one line or if a line has multiple comma-separated values
        if len(lines) == 1 and ',' in lines[0]:
            # It's a single row with multiple columns
            np_matrix = np.loadtxt(filename, delimiter=',').reshape(1, -1)
        else:
            # It's either a single column or multiple rows, load normally
            np_matrix = np.loadtxt(filename, delimiter=',')

            # Check if it's genuinely multi-dimensional or just a single column
            if np_matrix.ndim == 1:
                np_matrix = np_matrix.reshape(-1, 1)  # Make it a column vector

        # Return the matrix as a list of lists and its dimensions
        return np_matrix.tolist(), np_matrix.shape[0], np_matrix.shape[1]

    except Exception as e:
        print("An Error Has Occurred")
        sys.exit(1)

def initialize_H(W, k):

    # Calculate the average of all entries in W
    m = np.mean(W)

    # Calculate the upper limit for the uniform distribution
    upper_limit = 2 * np.sqrt(m / k)

    # Initialize H with random values from the interval [0, upper_limit]
    # Assuming W is an n x n matrix, H should be n x k
    n = len(W)
    H = np.random.uniform(0, upper_limit, (n, k))
    if H.ndim == 1:
            H = H.reshape(1, -1)  # Reshape as one row with multiple columns
    return H.tolist()



if __name__ == "__main__":
    np.random.seed(0)

    # test_sym_and_ddg()
     # Check the number of arguments (including the script name)
    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        sys.exit(1)

    # Retrieve the arguments
    try:
        K = int(sys.argv[1])  # Convert K to an integer
        goal = sys.argv[2]
        file_name = sys.argv[3]
    except ValueError:
        print("An Error Has Occurred")
        sys.exit(1)  # Exit the script with an error code


    matrix , rows , cols = read_file(file_name)

    if K >= rows:
        print("An Error Has Occurred")
        sys.exit(1)
    

    ############### SYMNMNF HANDLER ###############
    if goal == "symnmf":
        try:
            W = symnmfmodule.norm(matrix)
            H = initialize_H(W,K)
            symnmf_result = symnmfmodule.symnmf(H,W)
            printMat(symnmf_result)
        except Exception as e:
            print("An Error Has Occurred")
            exit(1)
        
    ############### SYM HANDLER ###############
    elif goal == "sym":
        try:
            sym_result = symnmfmodule.sym(matrix)
            printMat(sym_result)
        except Exception as e:
            print("An Error Has Occurred")
            exit(1)

    ############### DDG HANDLER ###############

    elif goal == "ddg":
        try:
            ddg_result = symnmfmodule.ddg(matrix)
            printMat(ddg_result)
        except Exception as e:
            print("An Error Has Occurred")
            exit(1)
    
    ############### NORM HANDLER ###############
    elif goal == "norm":
        try:
            norm_result = symnmfmodule.norm(matrix)
            printMat(norm_result)
        except Exception as e:
            print("An Error Has Occurred")
            exit(1)
    else:
        print("An Error Has Occurred")
        exit(1)  # Exit the script with an error code


