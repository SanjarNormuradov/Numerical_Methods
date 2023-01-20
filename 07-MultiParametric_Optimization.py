# from math import *
import math
from random import randint

# Create Empty Matrix
def matrix(Nrows, Ncolumns):
    matrix_tmp = [[0 for i in range(Ncolumns)] for j in range(Nrows)]
    return matrix_tmp


# Edge Vector Length Max
def edge_vect_len_max(edg_vec_cmp_mat_tmp):
    max = 0
    for vectorN in range(6):
        temp = edge_vect_len(edg_vec_cmp_mat_tmp, vectorN)
        if temp >= max:
            max = temp
    return max


# Edge Vector Length Min
def edge_vect_len_min(edg_vec_cmp_mat_tmp):
    min = 1e+16
    for vectorN in range(6):
        temp = edge_vect_len(edg_vec_cmp_mat_tmp, vectorN)
        if temp <= min:
            min = temp
    return min


# Calculate Edge Vector Length
def edge_vect_len(edg_vec_cmp_mat_tmp, i):
    edge_vect_len_sum = 0
    for k in range(3):
        edge_vect_len_sum += math.pow(float(edg_vec_cmp_mat_tmp[i][k]), 2)
    edge_vect_len_sum = math.sqrt(edge_vect_len_sum)
    return edge_vect_len_sum


# Calculate Edge Vector's Component Length
def edge_vect_comp_len(edg_vec_cor_mat_tmp, i, j, k):
    return math.fabs(float(edg_vec_cor_mat_tmp[i][k]) - float(edg_vec_cor_mat_tmp[j][k]))


# Calculate Edge Vector's Component Length Matrix 6x3
def edge_vect_comp_mat(edg_vec_cor_mat_tmp):
    edg_vec_cmp_mat_tmp = matrix(6, 3)
    # vector1
    for k in range(3):
        edg_vec_cmp_mat_tmp[0][k] = edge_vect_comp_len(edg_vec_cor_mat_tmp, 0, 1, k)
    # vector2
    for k in range(3):
        edg_vec_cmp_mat_tmp[1][k] = edge_vect_comp_len(edg_vec_cor_mat_tmp, 0, 2, k)
        # print(edg_vec_cmp_mat_tmp[1][k])
    # vector3
    for k in range(3):
        edg_vec_cmp_mat_tmp[2][k] = edge_vect_comp_len(edg_vec_cor_mat_tmp, 0, 3, k)
    # vector4
    for k in range(3):
        edg_vec_cmp_mat_tmp[3][k] = edge_vect_comp_len(edg_vec_cor_mat_tmp, 1, 2, k)
    # vector5
    for k in range(3):
        edg_vec_cmp_mat_tmp[4][k] = edge_vect_comp_len(edg_vec_cor_mat_tmp, 1, 3, k)
    # vector6
    for k in range(3):
        edg_vec_cmp_mat_tmp[5][k] = edge_vect_comp_len(edg_vec_cor_mat_tmp, 2, 3, k)

    return edg_vec_cmp_mat_tmp


def edge_vect_angle(edg_vec_cmp_mat_tmp):
    edg_vec_ang_tmp = []
    # cos(vector1^vector2)
    tmp = 0
    for k in range(3):
        tmp += edg_vec_cmp_mat_tmp[0][k] * edg_vec_cmp_mat_tmp[1][k]
    tmp /= edge_vect_len(edg_vec_cmp_mat_tmp, 0) * edge_vect_len(edg_vec_cmp_mat_tmp, 1)
    edg_vec_ang_tmp.append(tmp)
    # cos(vector1^vector3)
    tmp = 0
    for k in range(3):
        tmp += edg_vec_cmp_mat_tmp[0][k] * edg_vec_cmp_mat_tmp[2][k]
    tmp /= edge_vect_len(edg_vec_cmp_mat_tmp, 0) * edge_vect_len(edg_vec_cmp_mat_tmp, 2)
    edg_vec_ang_tmp.append(tmp)
    # cos(vector2^vector3)
    tmp = 0
    for k in range(3):
        tmp += edg_vec_cmp_mat_tmp[1][k] * edg_vec_cmp_mat_tmp[2][k]
    tmp /= edge_vect_len(edg_vec_cmp_mat_tmp, 1) * edge_vect_len(edg_vec_cmp_mat_tmp, 2)
    edg_vec_ang_tmp.append(tmp)
    # cos(vector4^vector5)
    tmp = 0
    for k in range(3):
        tmp += edg_vec_cmp_mat_tmp[3][k] * edg_vec_cmp_mat_tmp[4][k]
    tmp /= edge_vect_len(edg_vec_cmp_mat_tmp, 3) * edge_vect_len(edg_vec_cmp_mat_tmp, 4)
    edg_vec_ang_tmp.append(tmp)

    return edg_vec_ang_tmp


# Calculate Pyramid Volume
def pyramid_volume(edg_vec_cmp_mat_tmp):
    # print(edg_vec_cmp_mat_tmp)
    diag1_pos = float(edg_vec_cmp_mat_tmp[0][0]) * float(edg_vec_cmp_mat_tmp[1][1]) * float(edg_vec_cmp_mat_tmp[2][2])
    # print(diag1_pos)
    diag2_pos = float(edg_vec_cmp_mat_tmp[0][1]) * float(edg_vec_cmp_mat_tmp[1][2]) * float(edg_vec_cmp_mat_tmp[2][0])
    diag3_pos = float(edg_vec_cmp_mat_tmp[0][2]) * float(edg_vec_cmp_mat_tmp[1][0]) * float(edg_vec_cmp_mat_tmp[2][1])

    diag1_neg = float(edg_vec_cmp_mat_tmp[0][0]) * float(edg_vec_cmp_mat_tmp[1][2]) * float(edg_vec_cmp_mat_tmp[2][1])
    diag2_neg = float(edg_vec_cmp_mat_tmp[0][1]) * float(edg_vec_cmp_mat_tmp[1][0]) * float(edg_vec_cmp_mat_tmp[2][2])
    diag3_neg = float(edg_vec_cmp_mat_tmp[0][2]) * float(edg_vec_cmp_mat_tmp[1][1]) * float(edg_vec_cmp_mat_tmp[2][0])
    pyramid_vol_tmp = (diag1_pos + diag2_pos + diag3_pos - diag1_neg - diag2_neg - diag3_neg) / 6
    return pyramid_vol_tmp


# Calculate One of Pyramid Surface Squares
def pyramid_surface_square(edg_vec_cmp_mat_tmp, edg_vec_ang_tmp, i, j):
    pyr_surf_sqr = edge_vect_len(edg_vec_cmp_mat_tmp, i) * edge_vect_len(edg_vec_cmp_mat_tmp, j)
    if (i + j) == 1:
        pyr_surf_sqr *= math.sqrt(1 - math.pow(float(edg_vec_ang_tmp[0]), 2)) / 2
    elif (i + j) == 2:
        pyr_surf_sqr *= math.sqrt(1 - math.pow(float(edg_vec_ang_tmp[1]), 2)) / 2
    elif (i + j) == 3:
        pyr_surf_sqr *= math.sqrt(1 - math.pow(float(edg_vec_ang_tmp[2]), 2)) / 2
    elif (i + j) == 7:
        pyr_surf_sqr *= math.sqrt(1 - math.pow(float(edg_vec_ang_tmp[3]), 2)) / 2

    return pyr_surf_sqr


# Calculate Pyramid Surface Square Overall
def pyramid_surface_square_total(edg_vec_cmp_mat_tmp):
    edge_vect_angles = edge_vect_angle(edg_vec_cmp_mat_tmp)
    # First Pyramid Surface Square
    pyr_surf_sqr1 = pyramid_surface_square(edg_vec_cmp_mat_tmp, edge_vect_angles, 0, 1)
    # Second Pyramid Surface Square
    pyr_surf_sqr2 = pyramid_surface_square(edg_vec_cmp_mat_tmp, edge_vect_angles, 0, 2)
    # Third Pyramid Surface Square
    pyr_surf_sqr3 = pyramid_surface_square(edg_vec_cmp_mat_tmp, edge_vect_angles, 1, 2)
    # Fourth Pyramid Surface Square
    pyr_surf_sqr4 = pyramid_surface_square(edg_vec_cmp_mat_tmp, edge_vect_angles, 3, 4)
    # Pyramid Surface Square Overall
    pyramid_surf_sqr_tmp = pyr_surf_sqr1 + pyr_surf_sqr2 + pyr_surf_sqr3 + pyr_surf_sqr4
    return pyramid_surf_sqr_tmp


# Calculate Function Value
def func(edg_vec_cor_mat_tmp, var_tmp):
    edg_vec_cmp_mat = edge_vect_comp_mat(edg_vec_cor_mat_tmp)
    # print(edg_vec_cor_mat_tmp)    # ok
    pyramid_vol = pyramid_volume(edg_vec_cmp_mat)
    # print(pyramid_vol)    # bug
    pyramid_surf_sqr = pyramid_surface_square_total(edg_vec_cmp_mat)
    # print(pyramid_surf_sqr)     # ok
    Lmax = edge_vect_len_max(edg_vec_cmp_mat)
    # print(Lmax)   # not ok, same as Lmin
    Lmin = edge_vect_len_min(edg_vec_cmp_mat)
    # print(Lmin)   # not ok, same as Lmax
    func_value = (math.pow(pyramid_vol, 2) / math.pow(pyramid_surf_sqr, 3)) * (Lmax / (Lmax + var_tmp * Lmin))
    # print(func_value)
    return func_value

# Calculate Function Partial Derivative Value
def func_deriv(edg_vec_cor_mat_tmp, var_tmp, i):
    achieved_accuracy_drv = False
    prev_deriv_calc = False
    func_deriv_tmp1 = 0
    func_deriv_tmp2 = 0
    delta_x = randint(1, 10) / 200
    delta_x2 = delta_x / 2
    while not achieved_accuracy_drv:
        if i == 1:
            if not prev_deriv_calc:
                edg_vec_cor_mat_tmp[2][0] += delta_x
                func_deriv_tmp1 += func(edg_vec_cor_mat_tmp, var_tmp)
                edg_vec_cor_mat_tmp[2][0] -= 2 * delta_x
                func_deriv_tmp1 -= func(edg_vec_cor_mat_tmp, var_tmp)
                func_deriv_tmp1 /= 2 * delta_x
                edg_vec_cor_mat_tmp[2][0] += delta_x

            edg_vec_cor_mat_tmp[2][0] += delta_x2
            func_deriv_tmp2 += func(edg_vec_cor_mat_tmp, var_tmp)
            edg_vec_cor_mat_tmp[2][0] -= 2 * delta_x2
            func_deriv_tmp2 -= func(edg_vec_cor_mat_tmp, var_tmp)
            func_deriv_tmp2 /= 2 * delta_x2
            edg_vec_cor_mat_tmp[2][0] += delta_x2
        elif i == 2:
            if not prev_deriv_calc:
                edg_vec_cor_mat_tmp[2][1] += delta_x
                func_deriv_tmp1 += func(edg_vec_cor_mat_tmp, var_tmp)
                edg_vec_cor_mat_tmp[2][1] -= 2 * delta_x
                func_deriv_tmp1 -= func(edg_vec_cor_mat_tmp, var_tmp)
                func_deriv_tmp1 /= 2 * delta_x
                edg_vec_cor_mat_tmp[2][1] += delta_x

            edg_vec_cor_mat_tmp[2][1] += delta_x2
            func_deriv_tmp2 += func(edg_vec_cor_mat_tmp, var_tmp)
            edg_vec_cor_mat_tmp[2][1] -= 2 * delta_x2
            func_deriv_tmp2 -= func(edg_vec_cor_mat_tmp, var_tmp)
            func_deriv_tmp2 /= 2 * delta_x2
            edg_vec_cor_mat_tmp[2][1] += delta_x2
        elif i == 3:
            if not prev_deriv_calc:
                edg_vec_cor_mat_tmp[3][0] += delta_x
                func_deriv_tmp1 += func(edg_vec_cor_mat_tmp, var_tmp)
                edg_vec_cor_mat_tmp[3][0] -= 2 * delta_x
                func_deriv_tmp1 -= func(edg_vec_cor_mat_tmp, var_tmp)
                func_deriv_tmp1 /= 2 * delta_x
                edg_vec_cor_mat_tmp[3][0] += delta_x

            edg_vec_cor_mat_tmp[3][0] += delta_x2
            func_deriv_tmp2 += func(edg_vec_cor_mat_tmp, var_tmp)
            edg_vec_cor_mat_tmp[3][0] -= 2 * delta_x2
            func_deriv_tmp2 -= func(edg_vec_cor_mat_tmp, var_tmp)
            func_deriv_tmp2 /= 2 * delta_x2
            edg_vec_cor_mat_tmp[3][0] += delta_x2
        elif i == 4:
            if not prev_deriv_calc:
                edg_vec_cor_mat_tmp[3][1] += delta_x
                func_deriv_tmp1 += func(edg_vec_cor_mat_tmp, var_tmp)
                edg_vec_cor_mat_tmp[3][1] -= 2 * delta_x
                func_deriv_tmp1 -= func(edg_vec_cor_mat_tmp, var_tmp)
                func_deriv_tmp1 /= 2 * delta_x
                edg_vec_cor_mat_tmp[3][1] += delta_x

            edg_vec_cor_mat_tmp[3][1] += delta_x2
            func_deriv_tmp2 += func(edg_vec_cor_mat_tmp, var_tmp)
            edg_vec_cor_mat_tmp[3][1] -= 2 * delta_x2
            func_deriv_tmp2 -= func(edg_vec_cor_mat_tmp, var_tmp)
            func_deriv_tmp2 /= 2 * delta_x2
            edg_vec_cor_mat_tmp[3][1] += delta_x2
        elif i == 5:
            if not prev_deriv_calc:
                edg_vec_cor_mat_tmp[3][2] += delta_x
                func_deriv_tmp1 += func(edg_vec_cor_mat_tmp, var_tmp)
                edg_vec_cor_mat_tmp[3][2] -= 2 * delta_x
                func_deriv_tmp1 -= func(edg_vec_cor_mat_tmp, var_tmp)
                func_deriv_tmp1 /= 2 * delta_x
                edg_vec_cor_mat_tmp[3][2] += delta_x

            edg_vec_cor_mat_tmp[3][2] += delta_x2
            func_deriv_tmp2 += func(edg_vec_cor_mat_tmp, var_tmp)
            edg_vec_cor_mat_tmp[3][2] -= 2 * delta_x2
            func_deriv_tmp2 -= func(edg_vec_cor_mat_tmp, var_tmp)
            func_deriv_tmp2 /= 2 * delta_x2
            edg_vec_cor_mat_tmp[3][2] += delta_x2
        else: return None
        if math.fabs(func_deriv_tmp2 - func_deriv_tmp1) > accuracy:
            delta_x /= 2
            delta_x2 = delta_x / 2
            func_deriv_tmp1 = func_deriv_tmp2
            prev_deriv_calc = True
            func_deriv_tmp2 = 0
        else:
            achieved_accuracy_drv = True

    return func_deriv_tmp1


# Normalized Function Partial Derivatives
def func_deriv_norm(edg_vec_cor_mat_tmp, var_tmp, j):
    func_deriv_norma = 0
    func_deriv_mat = []
    for k in range(1, 6):
        func_deriv_mat.append(func_deriv(edg_vec_cor_mat_tmp, var_tmp, k))
        func_deriv_norma += math.pow(func_deriv_mat[k-1], 2)
        # print(func_deriv_mat)
        # print(norma)
    func_deriv_norma = math.sqrt(func_deriv_norma)
    if func_deriv_norma == 0:
        func_deriv_normd = 0
    else:
        func_deriv_normd = float(func_deriv_mat[j-1]) / func_deriv_norma
    return func_deriv_normd


# Constants
edg_vec_cor_mat = matrix(4, 3)
var = 10
accuracy = 1e-5
coefficient = 1
global_func_max = 0
Ncount = 0
tmp = 0

while Ncount < 20:
    edg_vec_cor_mat[0][0] = 0  # x1
    edg_vec_cor_mat[0][1] = 0  # y1
    edg_vec_cor_mat[0][2] = 0  # z1
    edg_vec_cor_mat[1][0] = 1  # x2
    edg_vec_cor_mat[1][1] = 0  # y2
    edg_vec_cor_mat[1][2] = 0  # z2
    edg_vec_cor_mat[2][2] = 0  # z3
    n_itter = 0
    achieved_accuracy_func = False
    prev_func_calc = False
    while not achieved_accuracy_func:
        n_itter += 1
        step = max(1e-3, 1 / (n_itter))
        if not prev_func_calc:
            #   Set Initial Random Coordinates (x3, y3, x4, y4, z4)
            for i in range(2, 4):
                for j in range(3):
                    if (i == 2) and (j == 2):
                        continue
                    else:
                        while tmp == 0:
                            tmp = randint(-100, 100) / 10
                        edg_vec_cor_mat[i][j] = tmp
            # edg_vec_cor_mat[2][0] = 0.4998296383806974
            # edg_vec_cor_mat[2][1] = 2.1901062977286836
            # edg_vec_cor_mat[3][0] = 0.4997624922574255
            # edg_vec_cor_mat[3][1] = -0.9468954277365494
            # edg_vec_cor_mat[3][2] = 1.9748454131528452
            # edg_vec_cor_mat[2][0] = 0.4
            # edg_vec_cor_mat[2][1] = 2.2
            # edg_vec_cor_mat[3][0] = 0.5
            # edg_vec_cor_mat[3][1] = -1.4
            # edg_vec_cor_mat[3][2] = 1.2
            #   Increment Coordinates
            for k in range(1, 6):
                if k == 1:
                    edg_vec_cor_mat[2][0] += step * func_deriv_norm(edg_vec_cor_mat, var, k)
                elif k == 2:
                    edg_vec_cor_mat[2][1] += step * func_deriv_norm(edg_vec_cor_mat, var, k)
                elif k == 3:
                    edg_vec_cor_mat[3][0] += step * func_deriv_norm(edg_vec_cor_mat, var, k)
                elif k == 4:
                    edg_vec_cor_mat[3][1] += step * func_deriv_norm(edg_vec_cor_mat, var, k)
                elif k == 5:
                    edg_vec_cor_mat[3][2] += step * func_deriv_norm(edg_vec_cor_mat, var, k)
            #   Calculate Function Value with Increment
            func_val1 = func(edg_vec_cor_mat, var)

        #   Double Increment Coordinates
        for k in range(1, 6):
            if k == 1:
                edg_vec_cor_mat[2][0] += coefficient * step * func_deriv_norm(edg_vec_cor_mat, var, k)
            elif k == 2:
                edg_vec_cor_mat[2][1] += coefficient * step * func_deriv_norm(edg_vec_cor_mat, var, k)
            elif k == 3:
                edg_vec_cor_mat[3][0] += coefficient * step * func_deriv_norm(edg_vec_cor_mat, var, k)
            elif k == 4:
                edg_vec_cor_mat[3][1] += coefficient * step * func_deriv_norm(edg_vec_cor_mat, var, k)
            elif k == 5:
                edg_vec_cor_mat[3][2] += coefficient * step * func_deriv_norm(edg_vec_cor_mat, var, k)
        #   Calculate Function Value with Double Increment
        func_val2 = func(edg_vec_cor_mat, var)

        if math.fabs(func_val1 - func_val2) > accuracy:
            coefficient = 2
            func_val1 = func_val2
            prev_func_calc = True
        else:
            achieved_accuracy_func = True

    if math.fabs(func_val2 - global_func_max) < accuracy:
        Ncount += 1
    else:
        if func_val2 > global_func_max:
            global_func_max = func_val2
            Ncount = 0
        else: continue

print(f"{float(edg_vec_cor_mat[0][0])}  {float(edg_vec_cor_mat[0][1])}  {float(edg_vec_cor_mat[0][2])}")
print(f"{float(edg_vec_cor_mat[1][0])}  {float(edg_vec_cor_mat[1][1])}  {float(edg_vec_cor_mat[1][2])}")
print(f"{round(float(edg_vec_cor_mat[2][0]), 1)}  {round(float(edg_vec_cor_mat[2][1]), 1)}  {float(edg_vec_cor_mat[2][2])}")
print(f"{round(float(edg_vec_cor_mat[3][0]), 1)}  {round(float(edg_vec_cor_mat[3][1]), 1)}  {round(float(edg_vec_cor_mat[3][2]), 1)}")
print(round(global_func_max, 6))
