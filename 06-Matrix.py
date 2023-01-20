from math import *
from random import randint


# Set Matrix Element
def matrix_element(i, j, b, s):
    return sin(pow((i+j),s)/log(i+j+b))


# Create Matrix
def matrix(var,size,b):
    s = (1 + var / 10) / 2
    matrix = [[matrix_element(i, j, b, s) for i in range(1, size + 1)] for j in range(1, size + 1)]
    return matrix


# Create Identity Matrix
def identity_matrix(size):
    id_matrix = [[0 for i in range(size)] for j in range(size)]
    for i in range(size):
        id_matrix[i][i] = 1
    return id_matrix


# Multiplication Of Matrix By Vector
def matrix_multiply_vector(A, b):
    result = []
    for i in range(len(b)):
        temp = 0
        for j in range(len(b)):
            temp += A[i][j] * b[j]
        result.append(temp)
    return result


# Multiplication Of Matrix By Value
def matrix_multiply_value(A, b):
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] *= b
    return A


# Subtraction Of Matrix, A - B
def matrix_subtract(A, B):
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] -= B[i][j]
    return A


#           (A)^N * b
def pow_matrix_N_vector(N, A, b):
    for i in range(N):
        b = matrix_multiply_vector(A, b)
    return b


# Random Normalized Vector
def vector_normalized(V):
    norm = 0
    for i in range(len(V)):
        norm += pow(V[i], 2)
    norm = sqrt(norm)
    for i in range(len(V)):
        V[i] /= norm
    return V


# Absolute Maximum Eigenvalue
def abs_max_eigenvalue(A, achieved_accuracy, prev_val_calc, calc_val):
    N1 = 1
    N2 = N1 * 2
    vector1 = vector2 = vector_normalized([randint(1, 100) for i in range(len(A))])
    while not achieved_accuracy:
        if prev_val_calc:
            eigenvalue1 = calc_val[0]
            del calc_val[0]
            prev_val_calc = False
        else:
            vector1 = vector_normalized(pow_matrix_N_vector(N1, A, vector1))
            i = 0
            while vector1[i] == 0:
                i += 1
            eigenvalue1 = matrix_multiply_vector(A, vector1)[i] / vector1[i]
            prev_val_calc = True

        vector2 = vector_normalized(pow_matrix_N_vector(N2, A, vector2))
        i = 0
        while vector2[i] == 0:
            i += 1
        eigenvalue2 = matrix_multiply_vector(A, vector2)[i] / vector2[i]
        calc_val.append(eigenvalue2)

        if abs(eigenvalue2 - eigenvalue1) > accuracy:
            N1 += 1
            N2 = N1 * 2

        else:
            # print(N2)
            achieved_accuracy = True

    return eigenvalue2


# Minimum Eigenvalue
def min_eigenvalue(A, achieved_accuracy, prev_val_calc, calc_val):
    eigenvalue = abs_max_eigenvalue(A, achieved_accuracy, prev_val_calc, calc_val)
    if eigenvalue > 0:
        temp = eigenvalue
        A = matrix_subtract(A, matrix_multiply_value(identity_matrix(len(A)), temp))
        eigenvalue = abs_max_eigenvalue(A, achieved_accuracy, prev_val_calc, calc_val)
        eigenvalue += temp
    return eigenvalue


# Constants
b = 1
var = 10
matrix_size = 100
accuracy = 1e-5
achieved_accuracy = False
prev_val_calc = False
calc_val = []

# for i in range(1, 19):
#     print(f'Variant {i} = {min_eigenvalue(matrix(i, matrix_size, b), achieved_accuracy, prev_val_calc, calc_val)}')

print(round(min_eigenvalue(matrix(var, matrix_size, b), achieved_accuracy, prev_val_calc, calc_val), 5))

