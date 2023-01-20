import math


def func_magnitude(c_parameter, x):
    return math.log1p(c_parameter * x) + math.pow(x, 3) * math.cos(c_parameter * x) / (math.sqrt(x) + math.exp(x))


a = 1
b = 2
c = float(input('c = '))

# Rectangular Method
N = 1
N2 = N * 2
segment = b - a
sum_1 = segment * func_magnitude(c, x=a+segment)
sum_2 = 0
for i in range(1, N2 + 1):
    sum_2 += func_magnitude(c, x=a + (2*i-1) * segment / (2 * N2))
sum_2 *= segment / N2
calculated_list = [sum_2]
sum_counted = True

while math.fabs(sum_1 - sum_2) > 1e-10:
    N += 1
    N2 = N * 2

    if sum_counted:
        sum_1 = float(calculated_list[0])
        del calculated_list[0]
        sum_counted = False
    else:
        sum_1 = 0
        for i in range(1, N + 1):
            sum_1 += func_magnitude(c, x=a + (2*i-1) * segment / (2 * N))
        sum_1 *= segment / N
        sum_counted = True

    sum_2 = 0
    for i in range(1, N2 + 1):
        sum_2 += func_magnitude(c, x=a + (2*i-1) * segment / (2 * N2))
    sum_2 *= segment / N2
    calculated_list.append(sum_2)

sum_1 = round(sum_1, 10)
print('Rectangular Method')
print(f'I = {sum_1}')
print(f'N in rectangular method = {N}')


# Trapezoid Method
N = 1
N2 = N * 2
segment = b - a
sum_whole = func_magnitude(c, x=a) + func_magnitude(c, x=b)
sum_1 = sum_whole
sum_2 = sum_1 + 2 * func_magnitude(c, x=a+segment/2)
sum_1 *= segment / 2
sum_2 *= segment / (2 * N2)
calculated_sum = [sum_2]
sum_counted = True

while math.fabs(sum_1 - sum_2) > 1e-10:
    N += 1
    N2 = N * 2

    if sum_counted:
        sum_1 = calculated_sum[0]
        del calculated_sum[0]
        sum_counted = False
    else:
        sum_1 = sum_whole
        for i in range(2, N + 1):
            sum_1 += 2 * func_magnitude(c, x=a + (i-1) * segment / N)
        sum_1 *= segment / (2 * N)
        sum_counted = True

    sum_2 = sum_whole
    for i in range(2, N2 + 1):
        sum_2 += 2 * func_magnitude(c, x=a + (i-1) * segment / N2)
    sum_2 *= segment / (2 * N2)
    calculated_sum.append(sum_2)

sum_1 = round(sum_1, 10)
print('Trapezoid Method')
print(f'I = {sum_1}')
print(f'N in trapezoid method = {N}')


# Simpson Method

N = 1
N2 = N * 2
segment = b - a
sum_whole = func_magnitude(c, x=a) + func_magnitude(c, x=b)
sum_1 = (sum_whole + 4 * func_magnitude(c, x=a+segment/2)) * segment / 6
sum_2 = sum_whole + 2 * func_magnitude(c, x=a+segment/2)
for k in range(1, N2 + 1):
    sum_2 += 4 * func_magnitude(c, x=a + (2*k-1) * segment / (2 * N2))
sum_2 *= segment / (6 * N2)
calculated_sum = [sum_2]
calculated_item = []
sum_counted = True

while math.fabs(sum_1 - sum_2) > 1e-10:
    N += 1
    N2 = N * 2

    if sum_counted:
        sum_1 = calculated_sum[0]
        del calculated_sum[0]
        sum_counted = False
    else:
        sum_1 = sum_whole
        for i in range(2, N + 1):
            sum_1 += 2 * func_magnitude(c, x=a + (i-1) * segment / N)
        for k in range(1, N + 1):
            sum_1 += 4 * func_magnitude(c, x=a + (2 * k - 1) * segment / (2 * N))
        sum_1 *= segment / (6 * N)
        sum_counted = True

    sum_2 = sum_whole
    for i in range(2, N2 + 1):
        sum_2 += 2 * func_magnitude(c, x=a + (i-1) * segment / N2)
    for k in range(1, N2 + 1):
        sum_2 += 4 * func_magnitude(c, x=a + (2 * k - 1) * segment / (2 * N2))
    sum_2 *= segment / (6 * N2)
    calculated_sum.append(sum_2)

sum_1 = round(sum_1, 10)
print('Sympson Method')
print(f'I = {sum_1}')
print(f'N in Sympson method = {N}')