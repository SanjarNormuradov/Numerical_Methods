import math

# Integral Has Its Limits From 0 To +Inf
# Function Goes Up To +Inf Near Zero


def f_func_magnitude(a_parameter, x):
    return math.log1p(x) * (math.tanh(a_parameter * x) - 1) / math.pow(x, 3/2)


def g_func_magnitude(a_parameter, x):
    return (a_parameter + 1/2) * math.pow(x, 1/2) - math.pow(x, -1/2)


# a = float(input('a = '))
a = 1
K = 1e-2  # K - Upper Of 1st And Lower Limits Of 2nd Sub_Integral, Final Integral More Accurate If K Smaller
L1 = 1e+1  # L1 - Upper Of 2nd And Lower Limits Of 3rd Sub_Integral, Final Integral More Accurate If L1 Larger
L2 = L1 * 2  # L2 - Twice Larger Than Limit L1, To Compare Both Final Integrals For Accuracy At The End
N1 = 1e+5  # N1 - Lower Limit Of 3rd Sub_Integral, Final Integral More Accurate If N1 Larger
N2 = N1 * 2  # N2 - Twice Larger Than Limit N1, To Compare Both Final Integrals For Accuracy At The End
x1 = 1e-4  # x1 - Step In 1st Sub_Integral, Final Integral More Accurate If x Smaller
x2 = x1 / 2  # x2 - Twice Smaller Than Step x1, To Compare Both Final Integrals For Accuracy At The End
X1 = 1e+2  # X1 - Step In 3st Sub_Integral, Final Integral More Accurate If X1 Larger
X2 = X1 * 2  # X2 - Twice Larger Than Step X1, To Compare Both Final Integrals For Accuracy At The End
P = 10000
sub_1_sum_analytic = 2 * ((a + 1/2) * math.pow(K, 3/2) / 3 - math.pow(K, 1/2))
sub_1_sum_1 = sub_1_sum_analytic
sum_1 = sum_2 = sub_2_sum_1 = sub_2_sum_2 = sub_3_sum_1 = sub_3_sum_2 = 0
calculated_sub_1_sum = []
calculated_sub_2_sum = []
calculated_sub_3_sum = []
sub_1_sum_counted = False
sub_2_sum_counted = False
sub_3_sum_counted = False
achieved_accuracy = False

while not achieved_accuracy:
    # 1st Sub_Integrals Calculations For Steps x1 And x2
    sub_1_sum_2 = sub_1_sum_analytic
    if sub_1_sum_counted:
        sub_1_sum_1 = calculated_sub_1_sum[0]
        del calculated_sub_1_sum[0]

        M = int(K / x2)
        step = x2
        sub_1_sum_tiny = (f_func_magnitude(a, x=K-step) + f_func_magnitude(a, x=step)) - (
                g_func_magnitude(a, x=K) + g_func_magnitude(a, x=step))
        for i in range(2, M):
            sub_1_sum_tiny += 2 * (f_func_magnitude(a, x=i * step) - g_func_magnitude(a, x=i * step))
        for j in range(1, M):
            sub_1_sum_tiny += 4 * (f_func_magnitude(a, x=step + (2 * j - 1) * step / 2)
                                   - g_func_magnitude(a, x=step + (2 * j - 1) * step / 2))
        sub_1_sum_tiny *= step / 6
        # 1st Section-Step Of 1st Sub_Integral, From 0 to x1, Rectangular Method
        sub_1_sum_tiny += (f_func_magnitude(a, x=step / 2) - g_func_magnitude(a, x=step / 2)) * step
        sub_1_sum_2 += sub_1_sum_tiny
        calculated_sub_1_sum.append(sub_1_sum_2)
    else:
        for counted_sum in range(1, 3):
            if counted_sum == 1:
                # 1st Sub_Integral With Step x1, 1st Time Counted
                M = int(K / x1)
                step = x1
            else:
                # 1st Sub_Integral With Step x2 (x2 = x1/2), 1st Time Counted
                M = int(K / x2)
                step = x2

            sub_1_sum_tiny = (f_func_magnitude(a, x=K-step) + f_func_magnitude(a, x=step)) - (g_func_magnitude(a, x=K-step) + g_func_magnitude(a, x=step))
            for i in range(2, M):
                sub_1_sum_tiny += 2 * (f_func_magnitude(a, x=i * step) - g_func_magnitude(a, x=i * step))
            for j in range(1, M):
                sub_1_sum_tiny += 4 * (f_func_magnitude(a, x=step + (2 * j - 1) * step / 2) - g_func_magnitude(a, x=step + (2 * j - 1) * step / 2))
            sub_1_sum_tiny *= step / 6
            # 1st Section-Step Of 1st Sub_Integral, From 0 to x1, Rectangular Method
            sub_1_sum_tiny += (f_func_magnitude(a, x=step / 2) - g_func_magnitude(a, x=step / 2)) * step

            if counted_sum == 1:
                sub_1_sum_1 += sub_1_sum_tiny
            else:
                sub_1_sum_2 += sub_1_sum_tiny
                calculated_sub_1_sum.append(sub_1_sum_2)
                sub_1_sum_counted = True

    # 2nd Sub_Integrals Calculations For Limits L1 And L2
    sub_2_sum = 0
    if sub_2_sum_counted:
        sub_2_sum_1 = calculated_sub_2_sum[0]
        del calculated_sub_2_sum[0]

        segment = L2 - K
        M = 2 * P # P is Doubled at the end of every cycle, hence Doubling M is obsolete
        sub_2_sum_2 = f_func_magnitude(a, x=K) + f_func_magnitude(a, x=L2)
        for i in range(2, M):
            sub_2_sum_2 += 2 * f_func_magnitude(a, x=K + (i - 1) * segment / M)
        for j in range(1, M):
            sub_2_sum_2 += 4 * f_func_magnitude(a, x=K + (2 * j - 1) * segment / (2 * M))
        sub_2_sum_2 *= segment / (6 * M)
        calculated_sub_2_sum.append(sub_2_sum_2)
    else:
        for counted_sum in range(1, 3):
            if counted_sum == 1:
                segment = L1 - K
                M = P
                sub_2_sum = f_func_magnitude(a, x=K) + f_func_magnitude(a, x=L1)
            else:
                segment = L2 - K
                M = 2 * P
                sub_2_sum = f_func_magnitude(a, x=K) + f_func_magnitude(a, x=L2)

            for i in range(2, M):
                sub_2_sum += 2 * f_func_magnitude(a, x=K + (i-1) * segment / M)
            for j in range(1, M):
                sub_2_sum += 4 * f_func_magnitude(a, x=K + (2 * j - 1) * segment / (2 * M))
            sub_2_sum *= segment / (6 * M)

            if counted_sum == 1:
                sub_2_sum_1 = sub_2_sum
            else:
                sub_2_sum_2 = sub_2_sum
                calculated_sub_2_sum.append(sub_2_sum_2)
                sub_2_sum_counted = True

    # 3rd Sub_Integrals Calculations For Limits L1:N1 And L2:N2 With Steps X1 And X2
    sub_3_sum = 0
    if sub_3_sum_counted:
        sub_3_sum_1 = calculated_sub_3_sum[0]
        del calculated_sub_3_sum[0]

        segment = N2 - L2
        M = int(segment / X2)
        sub_3_sum_2 = f_func_magnitude(a, x=N2) + f_func_magnitude(a, x=L2+X2)
        for i in range(2, M):
            sub_3_sum_2 += 2 * f_func_magnitude(a, x=L2 + (i - 1) * segment / M)
        for j in range(1, M):
            sub_3_sum_2 += 4 * f_func_magnitude(a, x=L2 + (2 * j - 1) * segment / (2 * M))
        sub_3_sum_2 *= segment / (6 * M)
        calculated_sub_3_sum.append(sub_3_sum_2)
    else:
        for counted_sum in range(1, 3):
            if counted_sum == 1:
                segment = N1 - L1
                L = L1
                M = int(segment/X1)
                sub_3_sum = f_func_magnitude(a, x=N1) + f_func_magnitude(a, x=L1+X1)
            else:
                segment = N2 - L2
                L = L2
                M = int(segment/X2)
                sub_3_sum = f_func_magnitude(a, x=N2) + f_func_magnitude(a, x=L2+X2)

            for i in range(2, M):
                sub_3_sum += 2 * f_func_magnitude(a, x=L + (i - 1) * segment / M)
            for j in range(1, M):
                sub_3_sum += 4 * f_func_magnitude(a, x=L + (2 * j - 1) * segment / (2 * M))
            sub_3_sum *= segment / (6 * M)

            if counted_sum == 1:
                sub_3_sum_1 = sub_3_sum
            else:
                sub_3_sum_2 = sub_3_sum
                calculated_sub_3_sum.append(sub_3_sum_2)
                sub_3_sum_counted = True

    # Parameters Variation
    x2 /= 2
    L2 *= 2
    P *= 2
    N2 *= 2
    X2 *= 2

    # Final Integrals
    sum_1 = sub_1_sum_1 + sub_2_sum_1 + sub_3_sum_1
    sum_2 = sub_1_sum_2 + sub_2_sum_2 + sub_3_sum_2

    print(round(sum_1, 4))
    print(round(sum_2, 4))
    print('')

    if math.fabs(sum_1 - sum_2) > 1e-4:
        achieved_accuracy = False
    else:
        achieved_accuracy = True
print(f'Integral = {round(sum_2, 4)}')