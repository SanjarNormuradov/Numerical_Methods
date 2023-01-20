import winsound
import time
import math
from random import randint

# Calculate Absolute Value of Complex Function
def abs_func_value(u, v):
    a1 = 4
    a2 = 9
    b1 = -5
    b2 = 1
    c1 = -5
    c2 = -2
    U = math.pow(u, 3) - 3*u*math.pow(v, 2) + a1*(math.pow(u, 2)-math.pow(v, 2)) - 2*a2*u*v + b1*u - b2*v + c1
    V = 3*math.pow(u, 2)*v - math.pow(v, 3) + 2*a1*v*u + a2*(math.pow(u, 2)-math.pow(v, 2)) + (b1*v+b2*u) + c2
    return math.sqrt(math.pow(U, 2) + math.pow(V, 2))

# Calculate Function Derivative
def func_deriv(u, v):
    func_deriv_comp = []
    a1 = 4
    a2 = 9
    b1 = -5
    b2 = 1
    U = 3*(math.pow(u, 2)-math.pow(v, 2)) + 2*(a1*u-a2*v) + b1
    V = -(6*u*v + 2*(a1*v+a2*u) + b2)
    norm = math.sqrt(math.pow(U, 2) + math.pow(V, 2))
    U /= norm
    V /= norm
    func_deriv_comp.append(U)
    func_deriv_comp.append(V)
    return func_deriv_comp


start_time = time.time()
# Constants
func_accuracy = 3
crd_accuracy = 5
func_values = []
func_root_ReX = []
func_root_ImX = []
func_root_ReX_cmp = []
func_root_ImX_cmp = []
# 6 starting points per quadrant, per each seg_range(= radius of circle)
N_start_pnts = 24
seg_range = 1
N_roots = 2
n_roots = 0

while n_roots < N_roots:
    achieved_func_accuracy = False
    set_start_pnt = False
    set_step = False
    N = 50
    start_pnts = 0
    n_seq = 0
    n_itter = 0
    sgn = 1
    step_psc = 6.00
    # N_psc = 1.1
    while not achieved_func_accuracy:
        if not set_start_pnt:
            # Set Initial Random Coordinates
            ReX = randint(seg_range*(-10), seg_range*(10))/100
            ImX = randint(seg_range*(-10), seg_range*(10))/100
            set_start_pnt = True
            start_pnts += 1

        while n_seq < 2:
            n_seq += 1
            while n_itter < N:
                n_itter += 1
                if not set_step:
                    step = max(1e-5, 1 / (n_itter + 1000))
                func_derivatives = func_deriv(ReX, ImX)
                #   Decrement Coordinates
                ReX += sgn * step * func_derivatives[0]
                ImX += sgn * step * func_derivatives[1]
            if n_seq == 1:
                #   Calculate Function Value After N Decrements
                func_val1 = abs_func_value(ReX, ImX)
                N *= 2
            else:
                #   Calculate Function Value After 2N Decrements
                func_val2 = abs_func_value(ReX, ImX)

        if math.fabs(func_val1 - func_val2) > math.pow(10, -func_accuracy):
            if func_val2 > func_val1:
                sgn *= -1
                step_psc -= 0.01
                step /= step_psc
                set_step = True
                # N_psc += 0.02
            func_val1 = func_val2
            # N *= N_psc
            N *= 1.1
            n_seq = 1
        else:
            if min(round(func_val2, func_accuracy), round(func_val1, func_accuracy)) < 3 * math.pow(10, -func_accuracy):
                achieved_func_accuracy = True
            else:
                set_start_pnt = False
                set_step = False
                n_seq = 0
                n_itter = 0
                sgn = 1
                N = 50
                step_psc = 6.00
                # N_psc = 1.1
                if start_pnts == N_start_pnts:
                    start_pnts = 0
                    seg_range += 1
                    # 6 random starting points per each radius of circle(=seg_range) per quadrant (pi * r^2 / 4)
                    N_start_pnts = (6 * seg_range**2) * 4
    if (n_roots == 0) or ((round(ReX, crd_accuracy-3) not in func_root_ReX_cmp) and (round(ImX, crd_accuracy-3) not in func_root_ImX_cmp)):
        func_values.append(min(round(func_val2, func_accuracy), round(func_val1, func_accuracy)))
        func_root_ReX.append(round(ReX, crd_accuracy))
        func_root_ImX.append(round(ImX, crd_accuracy))
        func_root_ReX_cmp.append(round(ReX, crd_accuracy-3))
        func_root_ImX_cmp.append(round(ImX, crd_accuracy-3))
        n_roots += 1

# print(step_psc)
# print(N_psc)
# print(seg_range)
# print(N_start_pnts)
print(f"Function Values = {func_values}")
print(f"ReX = {func_root_ReX}")
print(f"ImX = {func_root_ImX}")
end_time = time.time()
print(f"Program Runtime = {end_time-start_time}")
winsound.Beep(500, 2000)
