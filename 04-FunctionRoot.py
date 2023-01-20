import math
from math import *

# 2 * e^(-a * x) = sqrt(x), 0 =< x <= 4

def func_magn(a_parameter, x):
    return 2 * exp(-a_parameter * x) - sqrt(x)

# def aux_func1_magn(x, n_iter):
#     if n_iter == 0:
#         result = x
#     else:
#         result = (x + log(2 / sqrt(aux_func1_magn(x, n_iter - 1))) / 2) / 2
#     return result


def aux_func_magn(x, n_iter):
    if n_iter == 0:
        result = x
    else:
        result = 2626963 * (x + 4 * exp((-2) * a * aux_func_magn(x, n_iter - 1))) / 15714239
    return result


# a = float(input('a = ')
a = 1
seg_begin = 0
seg_end = 4
root = 0
N_itter = 1
calc_magn = []
prev_magn_calc = False
achieved_accuracy = False

# Iterative Method
point = seg_end
while not achieved_accuracy:
    if prev_magn_calc:
        prev_magn = calc_magn[0]
        del calc_magn[0]
        prev_magn_calc = False
    else:
        prev_magn = aux_func_magn(point, N_itter)
        prev_magn_calc = True

    next_magn = aux_func_magn(point, 2*N_itter)
    calc_magn.append(next_magn)
    if (next_magn < seg_begin) or (next_magn > seg_end):
        if point == seg_begin:
            print(f"Limit {round(next_magn, 4)} from starting point {point} is out of segment [{seg_begin}, {seg_end}]")
            print('Iterative Method result in this is Error')
            break
        else:
            point = seg_begin
            calc_magn.clear()
            prev_magn_calc = False
    else:
        if fabs(next_magn - prev_magn) > 1e-8:
            # achieved_accuracy = False
            N_itter += 1
        else:
            achieved_accuracy = True
            print()
            print('****** Iterative Method ******')
            print(f'Root = {round(next_magn, 8)}')
            print(f'N = {N_itter}')
            print(f'Limit starts from x = {point}')

achieved_accuracy = False
N_itter = 0
# Bisection Method
while not achieved_accuracy:
    seg_mid = seg_begin + fabs(seg_begin - seg_end) / 2
    N_itter += 1
    if func_magn(a, x=seg_mid) * func_magn(a, x=seg_begin) > 0:
        seg_begin = seg_mid
    elif func_magn(a, x=seg_mid) * func_magn(a, x=seg_begin) < 0:
        seg_end = seg_mid
    else:
        if func_magn(a, x=seg_mid) == 0:
            root = seg_mid
        else:
            root = seg_begin

    if fabs(seg_begin - seg_end) > 1e-8:
        achieved_accuracy = False
    else:
        achieved_accuracy = True
        root = seg_begin

print()
print('****** Bisection Method ******')
print(f'Root = {round(root, 8)}')
print(f'N = {N_itter}')