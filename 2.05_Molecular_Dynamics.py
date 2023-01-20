import time
import winsound
import numpy as np

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


#   Energy Gradient (i.e. Force) Calculation
def grad_E(crds, n):
    d = np.float_power(10, -6)
    grad_lst = np.zeros((n, 3))
    for v in range(n):
        for grad_cmp in range(3):
            for div_cmp in range(2):
                if div_cmp == 0:
                    #   Increment Coordinate X,Y,Z Of Victim Atom By d/2 To Calculate E1 For gradE
                    crds[v][grad_cmp] += d / 2
                    E1 = E_system(crds, n)
                else:
                    #   Decrement Coordinate X,Y,Z Of Victim Atom By (2 * d/2) To Calculate E2 For gradE
                    crds[v][grad_cmp] -= d
                    E2 = E_system(crds, n)

            #   Return Coordinates To Initial Values
            crds[v][grad_cmp] += d / 2

            gradE = - (E1 - E2) / d
            grad_lst[v][grad_cmp] = gradE
    return grad_lst


def struct_relax(crds, n):
    step = np.float_power(10, -3)
    grad_accuracy = 3
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
# Roci = 2.1293379427102295
Rcut = 2.15
pi = np.pi
dt = 1
N = 40 * 10**1
m = 1.2 * 10**3
with open("Fullerene_InitCrds.xyz", "r") as file:
    file_txt = file.readlines()

with open("Fullerene.xyz", "w") as file:
    for Vstep in range(1):
        #   1.  Fullerene System Building From Optimized Coordinates
        Crds_Fllrn = np.zeros((40, 3))
        for str_nmbr in range(len(file_txt)):
            if str_nmbr > 1:
                file_str = file_txt[str_nmbr].split()
                for crd_cmp in range(4):
                    if crd_cmp > 0:
                        Crds_Fllrn[str_nmbr - 2][crd_cmp - 1] = file_str[crd_cmp]

        for atm_nmbr in range(20, 40):
            for crd_cmp in range(3):
                if crd_cmp == 0:
                    Crds_Fllrn[atm_nmbr][crd_cmp] = Crds_Fllrn[atm_nmbr - 20][crd_cmp] + Rcut * 3
                else:
                    Crds_Fllrn[atm_nmbr][crd_cmp] = Crds_Fllrn[atm_nmbr - 20][crd_cmp] * 1

        Vinit = 0.117 + Vstep * 0.02
        grad_lst_Fllrn = np.zeros((40, 3))
        grad_lst_Fllrn_tmp = np.zeros((40, 3))
        V_Fllrn = np.zeros((40, 3))
        for atm_nmbr in range(40):
            if atm_nmbr in range(20):
                V_Fllrn[atm_nmbr][0] = Vinit / 2
            else:
                V_Fllrn[atm_nmbr][0] = -Vinit / 2

        for step in range(N):
            file.write("40\n")
            file.write(f"Vinit = {Vinit * 10**5} m/s, step {step+1}\n")
            #   2.1. Change Coordinates
            Crds_Fllrn += V_Fllrn * dt + grad_lst_Fllrn * dt ** 2 / (2 * m)
            #   2.2. Save Old Forces
            grad_lst_Fllrn_tmp = grad_lst_Fllrn * 1.0
            for atm_nmbr in range(40):
                file.write(f"C {Crds_Fllrn[atm_nmbr][0]} {Crds_Fllrn[atm_nmbr][1]} {Crds_Fllrn[atm_nmbr][2]}\n")

            #   2.3. Calculate New Forces
            grad_lst_Fllrn = grad_E(Crds_Fllrn, 40)

            #   2.4. Change Velocities
            V_Fllrn += (grad_lst_Fllrn + grad_lst_Fllrn_tmp) / 2 * dt / m

end = time.time()
winsound.Beep(500, 2000)
print(f"Program Runtime = {round((end - start)/60, 3)} mins")
