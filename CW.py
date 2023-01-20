import time
import winsound
import numpy as np
import matplotlib.pyplot as plt

#   Distance Between Two Atoms
def r_dist(v, c, crds):
    r = np.sqrt((crds[v][0] - crds[c][0])**2 + (crds[v][1] - crds[c][1])**2 + (crds[v][2] - crds[c][2])**2)
    return r


#   Angle Between 2 Vectors Connecting 3 Atoms
def cos_theta(v, c, w, crds):
    #   Vector1, Victim and Criminal Atoms
    V1_x = crds[c][0] - crds[v][0]
    V1_y = crds[c][1] - crds[v][1]
    V1_z = crds[c][2] - crds[v][2]
    #   Vector2, Victim and Witness Atoms
    V2_x = crds[w][0] - crds[v][0]
    V2_y = crds[w][1] - crds[v][1]
    V2_z = crds[w][2] - crds[v][2]
    # Vectors' Module
    V1 = np.sqrt(V1_x**2 + V1_y**2 + V1_z**2)
    V2 = np.sqrt(V2_x**2 + V2_y**2 + V2_z**2)
    cos0 = ((V1_x * V2_x) + (V1_y * V2_y) + (V1_z * V2_z))/(V1 * V2)
    return cos0


#   Cut Function
def f_cut(v, c, crds):
    R = 2.00
    D = 0.15
    r = r_dist(v, c, crds)
    if r < (R - D):
        func_value = 1
    elif np.fabs(R - r) <= D:
        func_value = (1 - np.sin(np.pi/2 * (r - R)/D))/2
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
def b_coef(v, c, crds, n):
    khi_vc = 0
    # W - Witness
    for w in range(n):
        if ((w != v) and (w != c)):
            fcut_vw = f_cut(v, w, crds)
            if fcut_vw != 0:
                cos_theta_vcw = cos_theta(v, c, w, crds)
                khi_vc += fcut_vw * g_coef(cos_theta_vcw)
    b = np.float_power(1 + khi_vc, -1/2)
    return b


#   Repulsion Energy
def E_rep(v, c, crds):
    D0 = 6
    r0 = 1.4276
    S = 2.167
    beta = 2.0099
    r = r_dist(v, c, crds)
    E = D0/(S-1) * np.exp(-beta * np.sqrt(2*S) * (r - r0))
    return E


#   Attraction Energy
def E_att(v, c, crds):
    D0 = 6
    r0 = 1.4276
    S = 2.167
    beta = 2.0099
    r = r_dist(v, c, crds)
    E = D0*S/(S-1) * np.exp(-beta * np.sqrt(2/S) * (r - r0))
    return E


#   System Energy
def E_system(crds, n):
    E_VictimAtom = 0
    #   V - Victim Atom, C - Criminal Atom
    for v in range(n):
        for c in range(n):
            if (v != c):
                fcut_vc = f_cut(v, c, crds)
                if fcut_vc != 0:
                    E_VictimAtom += fcut_vc * (E_rep(v, c, crds) - b_coef(v, c, crds, n) * E_att(v, c, crds))
    E_VictimAtom /= 2
    return E_VictimAtom


def grad_E(crds, n):
    d = np.float_power(10, -6)
    grad_lst = np.zeros((n, 3))
    for v in range(n):
        for grad_cmp in range(3):
            for div_cmp in range(2):
                if div_cmp == 0:
                    #   Increment Coordinate X,Y,Z Of Victim Atom By d/2 To Calculate E1 For gradE
                    crds[v][grad_cmp] += d / 2
                else:
                    #   Decrement Coordinate X,Y,Z Of Victim Atom By (2 * d/2) To Calculate E2 For gradE
                    crds[v][grad_cmp] -= d
                if div_cmp == 0:
                    E1 = E_system(crds, n)
                else:
                    E2 = E_system(crds, n)

            #   Return Coordinates To Initial Values
            crds[v][grad_cmp] += d / 2

            gradE = - (E1 - E2) / d
            grad_lst[v][grad_cmp] = gradE
    return grad_lst

#   Energy Gradient (i.e. Force) Calculation
def struct_relax(crds, n):
    step = np.float_power(10, -3)
    grad_accuracy = 4
    max_grad = 1
    while max_grad > 10**(-grad_accuracy):
        #   1. Atoms gradF Calculation
        #   v - Victim Of Force
        grad_lst = grad_E(crds, n)
        #   2. Coordinates Optimization
        crds += step * grad_lst
        max_grad = np.max(np.abs(grad_lst))
    return


start = time.time()
R = 1.45
Ampl = 0.1
pi = np.pi
#   1.1.  Dimer System Building
Crds_Dmr = np.zeros((2, 3))
for atm_nmbr in range(2):
    for crd_cmp in range(3):
        if crd_cmp == 0:
            Crds_Dmr[atm_nmbr][crd_cmp] = (1 - 2*atm_nmbr) * R / 2
#   1.2. Structural Relaxation
struct_relax(Crds_Dmr, 2)
#   1.3. Set Initial Displacement For Amplitude Of 0.2 A
Crds_Dmr[0][0] += Ampl
Crds_Dmr[1][0] -= Ampl
grad_lst_Dmr = grad_E(Crds_Dmr, 2)

dt = 1
N = 10**3
m = 1243
grad_lst_Dmr_tmp = np.zeros((2, 3))
V_Dmr = np.zeros((2, 3))
V_lst = []
F_lst = []
Crds_lst = []
V_lst.append(0.)
F_lst.append(grad_lst_Dmr[0][0])
Crds_lst.append(round(Crds_Dmr[1][0], 3))
with open("Dimer.xyz", "w") as file:
    for step in range(N):
        file.write("2\n")
        file.write(f"step {step+1}\n")
        #   2.1. Change Coordinates
        Crds_Dmr += V_Dmr * dt + grad_lst_Dmr * dt**2 / (2 * m)
        Crds_lst.append(round(Crds_Dmr[1][0], 3))
        #   2.2. Save Old Forces
        grad_lst_Dmr_tmp = grad_lst_Dmr + 0

        for atm_nmbr in range(2):
            file.write(f"C {Crds_Dmr[atm_nmbr][0]} {Crds_Dmr[atm_nmbr][1]} {Crds_Dmr[atm_nmbr][2]}\n")

        #   2.3. Calculate New Forces
        grad_lst_Dmr = grad_E(Crds_Dmr, 2)
        F_lst.append(grad_lst_Dmr[1][0])

        #   2.4. Change Velocities
        V_Dmr += (grad_lst_Dmr + grad_lst_Dmr_tmp) / 2 * dt / m
        V_lst.append(V_Dmr[1][0])

T = 0
for step in range(N):
    if (V_lst[step] < 0) and T == 0:
        T = 1
    elif (V_lst[step] > 0) and T == 1:
        break

T = dt * step
Vmax = round(max(np.abs(V_lst)) * 10**5)
Vavg = round(sum(np.abs(V_lst[:step])) / T * 10**5)
end = time.time()
winsound.Beep(500, 2000)
print(f"Period = {T} * 10^(-15) sec")
print(f"Max Velocity = {Vmax} m/sec")
print(f"Average Absolute Velocity = {Vavg} m/sec")
print(f"Program Runtime = {round(end - start, 3)} sec")
print(Crds_lst[:step])

