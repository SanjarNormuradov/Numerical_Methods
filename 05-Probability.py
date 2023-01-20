import math
from math import *


def fact_coeff(N1, N2, N3, N4):
    return factorial(N) / (factorial(N1) * factorial(N2) * factorial(N3) * factorial(N4))


def final_probability(P1, P2, P3, P4, N1, N2, N3, N4):
    #   The chance of not getting P1 or getting Q1=(1-P1) is Q1^(N-N2-N3-N4 - N1) or Q1^0=1
    #   Since Situation Q1 Occurs After All Situations (P1, P2, P3, P4) Have Occurred
    return fact_coeff(N1, N2, N3, N4) * (pow(P1, N1) * pow(P2, N2) * pow(P3, N3) * pow(P4, N4))


N = int(input("N = "))
x = int(input("x = "))
y = int(input("y = "))
u = int(input("u = "))
v = int(input("v = "))
a = int(input("a = "))
b = int(input("b = "))

#  No Crossed Segments
if (u > y) | (v < x):
    if N >= (a + b):
        P1 = 0
        N1 = 0
        P2 = (y - x + 1) / 37
        N2 = a
        P3 = (v - u + 1) / 37
        N3 = b
        P4 = (35 - (y - x) - (v - u)) / 37
        N4 = N - N1 - N2 - N3
        P = round(final_probability(P1, P2, P3, P4, N1, N2, N3, N4), 12)
    else:
        P = 0

# Partially Crossed Segments
elif ((u > x) & (y >= u) & (y < v)) | ((y > v) & (v >= x) & (x > u)):
    if ((a + b) <= N):  # (a + b) <= 2 * N
        # 0 =< N1 <= min(a, b)
        P = 0
        for N1 in range(0, min(a, b) + 1):
            N2 = a - N1
            N3 = b - N1
            N4 = N - N1 - N2 - N3
            if ((u > x) & (y >= u) & (y < v)):
                P1 = (y - u + 1) / 37
                P2 = (u - x) / 37
                P3 = (v - y) / 37
                P4 = (36 - (y - u) - (u - x) - (v - y)) / 37
            else:
                P1 = (v - x + 1) / 37
                P2 = (y - v) / 37
                P3 = (x - u) / 37
                P4 = (36 - (y - v) - (x - u) - (v - x)) / 37
            P += final_probability(P1, P2, P3, P4, N1, N2, N3, N4)
        P = round(P, 12)
    else:
        P = 0

# Nested Segments
elif ((u <= x) and (y <= v)) or ((x <= u) and (v <= y)):
    if ((u <= x) and (y <= v)):
        if (a <= b) and ((a + b) <= N):
            N1 = a
            N2 = 0
            N3 = b - N1
            N4 = N - N1 - N2 - N3  # N4 = N - b
            P1 = (y - x + 1) / 37
            P2 = 0
            P3 = ((v - y) + (x - u)) / 37
            P4 = (36 - (y - x) - (v - y) - (x - u)) / 37
            P = round(final_probability(P1, P2, P3, P4, N1, N2, N3, N4), 12)
        else:
            P = 0
    else:
        if (b <= a) and ((a + b) <= N):
            N1 = b
            N2 = a - N1
            N3 = 0
            N4 = N - N1 - N2 - N3  # N4 = N - a
            P1 = (v - u + 1) / 37
            P2 = ((y - v) + (u - x)) / 37
            P3 = 0
            P4 = (36 - (v - u) - (y - v) - (u - x)) / 37
            P = round(final_probability(P1, P2, P3, P4, N1, N2, N3, N4), 12)
        else:
            P = 0

print(f'Probability = {P}')
