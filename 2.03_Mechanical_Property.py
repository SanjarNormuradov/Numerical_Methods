import math
import time
import winsound
import numpy as np
import matplotlib.pyplot as plt


#   Distance Between Two Atoms
def r_dist(v, c, x_crd, y_crd):
    r = math.sqrt((x_crd[v]-x_crd[c])**2+(y_crd[v]-y_crd[c])**2)
    return r


#   Angle Between 2 Vectors Connecting 3 Atoms
def cos_theta(v, c, w, x_crd, y_crd):
    #   Vector1, Victim and Criminal Atoms
    V1_x = x_crd[c] - x_crd[v]
    V1_y = y_crd[c] - y_crd[v]
    #   Vector2, Victim and Witness Atoms
    V2_x = x_crd[w] - x_crd[v]
    V2_y = y_crd[w] - y_crd[v]
    # Vectors' Module
    V1 = math.sqrt(V1_x**2 + V1_y**2)
    V2 = math.sqrt(V2_x**2 + V2_y**2)
    cos0 = ((V1_x * V2_x) + (V1_y * V2_y))/(V1 * V2)
    return cos0


#   Cut Function
def f_cut(v, c, x_crd, y_crd):
    R = 2.00
    D = 0.15
    r = r_dist(v, c, x_crd, y_crd)
    if r < (R - D):
        func_value = 1
    elif math.fabs(R - r) <= D:
        func_value = (1 - math.sin(math.pi/2 * (r - R)/D))/2
    else:
        func_value = 0
    return func_value


#   Angular Function g(theta)
def g_coef(cos_theta_ijk):
    gamma = 0.11233
    c = 181.910
    d = 6.28433
    h = 0.5556
    g = gamma * (1 + c**2 * (1/d**2 - 1/(d**2 + (h + cos_theta_ijk)**2)))
    return g


#   MultiParticle Element b_ij
def b_coef(v, c, x_crd, y_crd, n):
    khi_vc = 0
    # W - Witness
    for w in range(3 * n - 4):
        if ((w != v) and (w != c)):
                    fcut_vw = f_cut(v, w, x_crd, y_crd)
                    if fcut_vw != 0:
                        cos_theta_vcw = cos_theta(v, c, w, x_crd, y_crd)
                        khi_vc += fcut_vw * g_coef(cos_theta_vcw)
    b = math.pow(1 + khi_vc, -1/2)
    return b


#   Repulsion Energy
def E_rep(v, c, x_crd, y_crd):
    D0 = 6
    r0 = 1.4276
    S = 2.167
    beta = 2.0099
    r = r_dist(v, c, x_crd, y_crd)
    E = D0/(S-1) * math.exp(-beta * math.sqrt(2*S) * (r - r0))
    return E


#   Attraction Energy
def E_att(v, c, x_crd, y_crd):
    D0 = 6
    r0 = 1.4276
    S = 2.167
    beta = 2.0099
    r = r_dist(v, c, x_crd, y_crd)
    E = D0*S/(S-1) * math.exp(-beta * math.sqrt(2/S) * (r - r0))
    return E


start = time.time()
R_OuterCircle = 1.6710941872164469
force_accuracy = 4
result_accuracy = 3
N = 6
#   1. System Building: Set Initial Coordinates
X_crds = []
Y_crds = []
for j in range(3):
    for i in range(N):
        if (j != 0) and ((i == 2) or (i == 3)):
            continue
        else:
            X_crds.append(R_OuterCircle * (math.cos(2 * math.pi * i / N + math.pi / 6) + (2 * j + 1) * math.cos(math.pi / 6)))
            Y_crds.append(R_OuterCircle * (math.sin(2 * math.pi * i / N + math.pi / 6)))

grad_lst = [[0 for k in range(2)] for i in range(3 * N - 4)]
d = math.pow(10, -6)
step_k = 1
step = step_k * math.pow(10, -3)
achieved_accuracy = False

#   2. Structural Relaxation To Find L0 - Initial Optimal Cell Length
while not achieved_accuracy:
    #   Atoms gradF Calculation
    # V - Victim Of Force, C - Criminal Of Force
    for V in range(3 * N - 4):
        for grad_cmp in range(2):
            for div_cmp in range(2):
                if grad_cmp == 0:
                    if div_cmp == 0:
                        #   Increment Coordinate X Of Victim Atom By d/2 To Calculate E1 For gradFx
                        X_crds[V] += d / 2
                    else:
                        #   Decrement Coordinate X Of Victim Atom By (2 * d/2) To Calculate E2 For gradFx
                        X_crds[V] -= d
                else:
                    if div_cmp == 0:
                        #   Increment Coordinate Y Of Victim Atom By d/2 To Calculate E1 For gradFy
                        Y_crds[V] += d / 2
                    else:
                        #   Decrement Coordinate Y Of Victim Atom By (2 * d/2) To Calculate E2 For gradFy
                        Y_crds[V] -= d
                E_VictimAtom = 0
                for C in range(3 * N - 4):
                    if (V != C):
                        fcut_vc = f_cut(V, C, X_crds, Y_crds)
                        if fcut_vc != 0:
                            E_VictimAtom += fcut_vc * (E_rep(V, C, X_crds, Y_crds) - b_coef(V, C, X_crds, Y_crds, N) * E_att(V, C, X_crds, Y_crds))
                if div_cmp == 0:
                    E1 = E_VictimAtom
                else:
                    E2 = E_VictimAtom
            #   Return Coordinates To Initial Values
            if grad_cmp == 0:
                X_crds[V] += d / 2
            else:
                Y_crds[V] += d / 2
            gradF = - (E1 - E2) / d
            grad_lst[V][grad_cmp] = gradF

    #   Coordinates Optimization
    for atm_nmbr in range(3 * N - 4):
        for grd_cmp in range(2):
            if grd_cmp == 0:
                X_crds[atm_nmbr] += step * grad_lst[atm_nmbr][grd_cmp]
            else:
                Y_crds[atm_nmbr] += step * grad_lst[atm_nmbr][grd_cmp]

    #   Checkpoint: |gradF| > accuracy?
    N_grad_accuracy_achieved = 0
    for atm_nmbr in range(3 * N - 4):
        for grd_cmp in range(2):
            if math.fabs(round(grad_lst[atm_nmbr][grd_cmp], force_accuracy)) <= math.pow(10, -force_accuracy):
                N_grad_accuracy_achieved += 1
    if N_grad_accuracy_achieved == 2 * (3 * N - 4):
        achieved_accuracy = True

#   Calculate L0
L0 = round(r_dist(2, 10, X_crds, Y_crds), result_accuracy)

#   3. "Strain Creates Stress" Loop
R = 2.15   #   R - Cutting Radius, Angstrem
X_crds_in = X_crds.copy()
Y_crds_in = Y_crds.copy()
stress_lst = [0]
eps_lst = [0]
for grph_pnt in range(1, 70):
    eps = grph_pnt / 100
    eps_lst.append(eps * 100)
    #   3.1 Put Strain On Atoms: Change X Coordinates
    for i in range(3 * N - 4):
        Y_crds[i] = Y_crds_in[i]
        if (i == 2) or (i == 3):
            continue
        else:
            X_crds[i] = (1 + eps) * X_crds_in[i]
    #   3.2 Structural Relaxation With Fixed Endings
    achieved_accuracy = False
    while not achieved_accuracy:
        #   Atoms gradF Calculation
        # V - Victim Of Force, C - Criminal Of Force
        for V in range(3 * N - 4):
            for grad_cmp in range(2):
                if ((V == 2) or (V == 3) or (V == 10) or (V == 13)) and (grad_cmp == 0):
                    grad_lst[V][grad_cmp] = 0
                else:
                    for div_cmp in range(2):
                        if grad_cmp == 0:
                            if div_cmp == 0:
                                #   Increment Coordinate X Of Victim Atom By d/2 To Calculate E1 For gradFx
                                X_crds[V] += d / 2
                            else:
                                #   Decrement Coordinate X Of Victim Atom By (2 * d/2) To Calculate E2 For gradFx
                                X_crds[V] -= d
                        else:
                            if div_cmp == 0:
                                #   Increment Coordinate Y Of Victim Atom By d/2 To Calculate E1 For gradFy
                                Y_crds[V] += d / 2
                            else:
                                #   Decrement Coordinate Y Of Victim Atom By (2 * d/2) To Calculate E2 For gradFy
                                Y_crds[V] -= d
                        E_VictimAtom = 0
                        for C in range(3 * N - 4):
                            if (V != C):
                                fcut_vc = f_cut(V, C, X_crds, Y_crds)
                                if fcut_vc != 0:
                                    E_VictimAtom += fcut_vc * (E_rep(V, C, X_crds, Y_crds) -
                                                               b_coef(V, C, X_crds, Y_crds, N) * E_att(V, C, X_crds, Y_crds))
                        if div_cmp == 0:
                            E1 = E_VictimAtom
                        else:
                            E2 = E_VictimAtom
                    #   Return Coordinates To Initial Values
                    if grad_cmp == 0:
                        X_crds[V] += d / 2
                    else:
                        Y_crds[V] += d / 2
                    gradF = - (E1 - E2) / d
                    grad_lst[V][grad_cmp] = gradF

        #   Coordinates Optimization
        for atm_nmbr in range(3 * N - 4):
            for grd_cmp in range(2):
                if grd_cmp == 0:
                    X_crds[atm_nmbr] += step * grad_lst[atm_nmbr][grd_cmp]
                else:
                    Y_crds[atm_nmbr] += step * grad_lst[atm_nmbr][grd_cmp]

        #   Checkpoint: |gradF| > accuracy?
        N_grad_accuracy_achieved = 0
        for atm_nmbr in range(3 * N - 4):
            for grd_cmp in range(2):
                if math.fabs(round(grad_lst[atm_nmbr][grd_cmp], force_accuracy)) <= math.pow(10, -force_accuracy):
                    N_grad_accuracy_achieved += 1
        if N_grad_accuracy_achieved == 2 * (3 * N - 4):
            achieved_accuracy = True

    #   3.3 Stress Calculation
    stress = 0
    for frc_cmp in range(4):
        if frc_cmp == 0:
            V = 2
        elif frc_cmp == 1:
            V = 3
        elif frc_cmp == 2:
            V = 10
        else:
            V = 13
        for div_cmp in range(2):
            if div_cmp == 0:
                #   Increment Coordinate X Of Victim Atom By d/2 To Calculate E1 For gradFx
                X_crds[V] += d / 2
            else:
                #   Decrement Coordinate X Of Victim Atom By (2 * d/2) To Calculate E2 For gradFx
                X_crds[V] -= d
            E_VictimAtom = 0
            for C in range(3 * N - 4):
                if (V != C):
                    fcut_vc = f_cut(V, C, X_crds, Y_crds)
                    if fcut_vc != 0:
                        E_VictimAtom += fcut_vc * (E_rep(V, C, X_crds, Y_crds) -
                                                   b_coef(V, C, X_crds, Y_crds, N) * E_att(V, C, X_crds, Y_crds))
            if div_cmp == 0:
                E1 = E_VictimAtom
            else:
                E2 = E_VictimAtom
        #   Return Coordinates To Initial Values
        X_crds[V] += d / 2
        gradF = - (E1 - E2) / d
        stress += math.fabs(gradF)
    #   Stress(GPa) = Force(eV -> J) / Square((h + R) * R), R = 2.15, Ang -> M
    stress /= 4 * (r_dist(1, 4, X_crds, Y_crds) + R) * R
    #   GPa Conversion
    stress *= 160.218
    stress_lst.append(round(stress))

#   4 Young's Modulus Calculations
Eym = 0
for i in range(3, 8):
    Eym += (stress_lst[i] - stress_lst[i - 3]) * 100 / (eps_lst[i] - eps_lst[i - 3])
Eym /= 5

plt.plot(eps_lst, stress_lst, ls='-', c='orange', lw=3, marker="o", markerfacecolor="green", ms=5)
font1 = {'family': 'Times New Roman', 'color': 'blue', 'weight': 'medium', 'size': 12}
font2 = {'family': 'Times New Roman', 'color': 'red', 'weight': 'semibold', 'size': 14}
plt.title("Stress - Strain Graph", loc='center', fontdict=font2)
plt.xlabel("Strain, %", fontdict=font1)
plt.ylabel("Stress, GPa", fontdict=font1)
plt.grid(color='green', linestyle='--', linewidth=1)
print(f"Optimal Length Of L0 = {L0}")
print(f"Young's Modulus = {round(Eym, result_accuracy)} GPa")
end = time.time()
print(f"Program Runtime = {round(end - start, result_accuracy)}")
winsound.Beep(500, 2000)
plt.show()
