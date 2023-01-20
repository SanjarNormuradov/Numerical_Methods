import math

k = 1
k2 = k * 2
a = float(input('a = '))
sum_0 = 1 / a
sum_1 = math.pow(-1, k) / ((3 * k + 1) * (3 * k + a)) + sum_0
sum_2 = math.pow(-1, k2) / ((3 * k2 + 1) * (3 * k2 + a)) + sum_1
calculated_list = [sum_2]

while math.fabs(sum_1 - sum_2) > 1e-10:
    k += 1
    sum_1 = float(calculated_list[0])# Odd k is not accounted, should be calculated separately
    del calculated_list[0]

    for i in range(k2 + 1, k * 2 + 1):
        sum_2 += math.pow(-1, i) / ((3 * i + 1) * (3 * i + a))
        calculated_list.append(sum_2)

    k2 = k * 2

sum_1 = round(sum_1, 10)
print(f'Sum of the row = {sum_1}')
print(f'The precision was achieved at N = {k}')
