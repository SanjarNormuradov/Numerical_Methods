import math
from random import randint

# Calculate Upper Function Value
def upper_func(x, var=10):
    return math.pow(x, 1 - var / 20)

# Calculate Upper Function Value
def lower_func(x, var=10):
    return math.pow(x, 1 + var / 10)

# # Calculate Exact Area
# def exact_area_Sympson(variant=10):
#     N1 = 0
#     segment = 1
#     sum_ends = 1
#     accuracy = 6
#     achieved_accuracy = False
#     calculated_sum_list = []
#     sum_calculated = False
#     while not achieved_accuracy:
#         N1 += 1
#         N2 = N1 * 2
#         if sum_calculated:
#             sum1 = calculated_sum_list[0]
#             del calculated_sum_list[0]
#             sum_calculated = False
#         else:
#             sum1 = sum_ends
#             for i in range(2, N1 + 1):
#                 sum1 += 2 * (upper_func(x=(i-1)*segment/N1, var=variant) - lower_func(x=(i-1)*segment/N1, var=variant))
#             for k in range(1, N1 + 1):
#                 sum1 += 4 * (upper_func(x=(2*k-1)*segment/(2*N1), var=variant) - lower_func(x=(2*k-1)*segment/(2*N1), var=variant))
#             sum1 *= segment / (6 * N1)
#             sum_calculated = True
#
#         sum2 = sum_ends
#         for i in range(2, N2 + 1):
#             sum2 += 2 * (upper_func(x=(i-1)*segment/N2, var=variant) - lower_func(x=(i-1)*segment/N2, var=variant))
#         for k in range(1, N2 + 1):
#             sum2 += 4 * (upper_func(x=(2*k-1)*segment/(2*N2), var=variant) - lower_func(x=(2*k-1)*segment/(2*N2), var=variant))
#         sum2 *= segment / (6 * N2)
#         calculated_sum_list.append(sum2)
#
#         if math.fabs(sum1 - sum2) <= math.pow(10, -accuracy):
#             achieved_accuracy = True
#
#     return round(sum1, accuracy)

# Calculate Exact Area
def exact_area(variant=10):
    return (3 * variant / 20) / ((2-variant/20) * (2+variant/10))

N_accuracy = 3
func_accuracy = 7
var = 10
MC_area_mean_lst = []
MC_area_error_lst = []
N_lst = []
N_cnt = 0
while N_cnt < 7:
    N_cnt += 1
    N1 = 8
    achieved_accuracy = False
    while not achieved_accuracy:
        N1 *= 2
        N2 = N1 * 2
        n1 = round(math.sqrt(N1), 0)
        n2 = round(math.sqrt(N2), 0)
        M1 = 0
        M2 = 0
        for i in range(1, N1 + 1):
            xi = randint(0, n1) / n1
            yi = randint(0, n1) / n1
            if ((yi >= lower_func(xi, var)) and (yi <= upper_func(xi, var))):
                M1 += 1
        MC_area1 = M1 / N1
        for i in range(1, N2 + 1):
            xi = randint(0, n2) / n2
            yi = randint(0, n2) / n2
            if ((yi >= lower_func(xi, var)) and (yi <= upper_func(xi, var))):
                M2 += 1
        MC_area2 = M2 / N2
        if math.fabs(MC_area2 - MC_area1) <= math.pow(10, -N_accuracy):
            achieved_accuracy = True
    N_lst.append(N1)
    t = 3.25
    cnt = 0
    MC_area_mean = 0
    MC_area_lst = []
    while cnt < 10:
        cnt += 1
        M1 = 0
        for i in range(1, N1 + 1):
            xi = randint(0, n1) / n1
            yi = randint(0, n1) / n1
            if ((yi >= lower_func(xi, var)) and (yi <= upper_func(xi, var))):
                M1 += 1
        MC_area1 = M1 / N1
        MC_area_lst.append(round(MC_area1, func_accuracy))
    for i in range(10):
        MC_area_mean += float(MC_area_lst[i])
    MC_area_mean /= 10
    dispersion_sum = 0
    for i in range(10):
        dispersion_sum += math.pow(float(MC_area_lst[i]) - MC_area_mean, 2)
    MC_area_error = round(t * math.sqrt(dispersion_sum / 90), func_accuracy)
    MC_area_mean_lst.append(round(MC_area_mean, func_accuracy))
    MC_area_error_lst.append(MC_area_error)

print(N_lst)
print(MC_area_mean_lst)
print(MC_area_error_lst)
print()
print(f"1) N = {max(N_lst)}  <=  \u0394S = {math.pow(10, -N_accuracy)}")
print(f"2) S = ({max(MC_area_mean_lst)} \u00B1 {min(MC_area_error_lst)}), t = {t}, \u03B1 = 0.99")
print(f"3) S_exact = {round(exact_area(var), func_accuracy)}"
      f"   Difference = {round(math.fabs(round(exact_area(var), func_accuracy) - (float(max(MC_area_mean_lst)) + float(min(MC_area_error_lst)))), func_accuracy)}")