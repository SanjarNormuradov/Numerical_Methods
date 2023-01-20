import math
import time
import winsound
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
    for w in range(n):
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
Rnph = 1.6710941872164469
Raz1 = 1.1
Raz2 = Raz1 * math.sin(math.pi / 5) / math.sin(math.pi / 7)
grad_accuracy_Az = 4
grad_accuracy_Nph = 4
energy_accuracy = 2
d = math.pow(10, -6)
# X_crd_Raz1 = []
# Y_crd_Raz1 = []
# X_crd_Raz2 = []
# Y_crd_Raz2 = []
# X_crd_R3 = []
# Y_crd_R3 = []
# for i in range(360):
#     X_crd_Raz1.append(Raz1 * math.cos(2 * math.pi / 360 * i) + Raz1 * (math.cos(math.pi / 5)))
#     Y_crd_Raz1.append(Raz1 * math.sin(2 * math.pi / 360 * i))
#     X_crd_Raz2.append(Raz2 * math.cos(2 * math.pi / 360 * i) - Raz2 * math.cos(math.pi / 7))
#     Y_crd_Raz2.append(Raz2 * math.sin(2 * math.pi / 360 * i))
# X_crd_R3.append(Roc1 * (math.cos(2 * math.pi / 360 * i) + math.cos(math.pi / 6)))
# Y_crd_R3.append(Roc1 * math.sin(2 * math.pi / 360 * i))


#   1.1. Azulene System Building: Set Initial Coordinates
X_crds_Az = []
Y_crds_Az = []
for i in range(5):
    X_crds_Az.append(Raz1*(math.cos(2*math.pi*i/5 + 6*math.pi/5)) + Raz1*math.cos(math.pi/5))
    Y_crds_Az.append(Raz1*(math.sin(2*math.pi*i/5 + 6*math.pi/5)))
for i in range(7):
    if (i == 0) or (i == 1):
        continue
    else:
        X_crds_Az.append(Raz2*(math.cos(2*math.pi*i/7 - math.pi/7)) - Raz2*math.cos(math.pi/7))
        Y_crds_Az.append(Raz2*(math.sin(2*math.pi*i/7 - math.pi/7)))

# plt.grid(color='green', linestyle='--', linewidth=1)
# plt.plot(X_crds_Az, Y_crds_Az, ls='', c='orange', lw=3, marker="o", markerfacecolor="green", ms=5)
# plt.plot(X_crd_Raz1, Y_crd_Raz1, lw=1)
# plt.plot(X_crd_Raz2, Y_crd_Raz2, lw=1)
# plt.show()

#   1.2. Naphthalene System Building: Set Initial Coordinates
X_crds_Nph = []
Y_crds_Nph = []
for j in range(2):
    for i in range(6):
        if (j == 1) and ((i == 0) or (i == 1)):
            continue
        else:
            X_crds_Nph.append(Rnph*(math.cos(2*math.pi*i/6 + (7-8*j)*math.pi/6) - (2*j-1)*math.cos(math.pi/6)))
            Y_crds_Nph.append(Rnph*(math.sin(2*math.pi*i/6 + (7-8*j)*math.pi/6)))

#   2.1. Structural Relaxation To Find Minimal Energy Of Azulene
grad_lst_Az = [[0 for k in range(2)] for i in range(10)]
step = math.pow(10, -3)
achieved_accuracy = False
while not achieved_accuracy:
    #   Atoms gradF Calculation
    # V - Victim Of Force, C - Criminal Of Force
    for V in range(10):
        for grad_cmp in range(2):
            for div_cmp in range(2):
                if grad_cmp == 0:
                    if div_cmp == 0:
                        #   Increment Coordinate X Of Victim Atom By d/2 To Calculate E1 For gradFx
                        X_crds_Az[V] += d / 2
                    else:
                        #   Decrement Coordinate X Of Victim Atom By (2 * d/2) To Calculate E2 For gradFx
                        X_crds_Az[V] -= d
                else:
                    if div_cmp == 0:
                        #   Increment Coordinate Y Of Victim Atom By d/2 To Calculate E1 For gradFy
                        Y_crds_Az[V] += d / 2
                    else:
                        #   Decrement Coordinate Y Of Victim Atom By (2 * d/2) To Calculate E2 For gradFy
                        Y_crds_Az[V] -= d
                E_VictimAtom = 0
                #   Va - Victim Affected
                for Va in range(10):
                    for C in range(10):
                        if (Va != C):
                            fcut_vc = f_cut(Va, C, X_crds_Az, Y_crds_Az)
                            if fcut_vc != 0:
                                E_VictimAtom += fcut_vc * (E_rep(Va, C, X_crds_Az, Y_crds_Az) - b_coef(Va, C, X_crds_Az, Y_crds_Az, 10) * E_att(Va, C, X_crds_Az, Y_crds_Az))
                if div_cmp == 0:
                    E1 = E_VictimAtom / 2
                else:
                    E2 = E_VictimAtom / 2
            #   Return Coordinates To Initial Values
            if grad_cmp == 0:
                X_crds_Az[V] += d / 2
            else:
                Y_crds_Az[V] += d / 2
            gradF = - (E1 - E2) / d
            grad_lst_Az[V][grad_cmp] = gradF

    #   Coordinates Optimization
    for atm_nmbr in range(10):
        for grd_cmp in range(2):
            if grd_cmp == 0:
                X_crds_Az[atm_nmbr] += step * grad_lst_Az[atm_nmbr][grd_cmp]
            else:
                Y_crds_Az[atm_nmbr] += step * grad_lst_Az[atm_nmbr][grd_cmp]

    #   Checkpoint: |gradF| > accuracy?
    N_grad_accuracy_achieved = 0
    for atm_nmbr in range(10):
        for grd_cmp in range(2):
            if math.fabs(round(grad_lst_Az[atm_nmbr][grd_cmp], grad_accuracy_Az)) <= math.pow(10, -grad_accuracy_Az):
                N_grad_accuracy_achieved += 1
    if N_grad_accuracy_achieved == 20:
        achieved_accuracy = True
# plt.grid(color='green', linestyle='--', linewidth=1)
# plt.plot(X_crds_Az, Y_crds_Az, ls='', c='orange', lw=3, marker="o", markerfacecolor="green", ms=5)
# plt.plot(X_crd_Raz1, Y_crd_Raz1, lw=1)
# plt.plot(X_crd_Raz2, Y_crd_Raz2, lw=1)
# plt.show()

#   E_Az - Minimal Energy Of Azulene
E_Az = 0
for V in range(10):
    for C in range(10):
        if V != C:
            fcut_vc = f_cut(V, C, X_crds_Az, Y_crds_Az)
            if fcut_vc != 0:
                E_Az += fcut_vc * (E_rep(V, C, X_crds_Az, Y_crds_Az) - b_coef(V, C, X_crds_Az, Y_crds_Az, 10) * E_att(V, C, X_crds_Az, Y_crds_Az))
E_Az *= 1/2

#   2.2. Structural Relaxation To Find Minimal Energy Of Naphthalene
grad_lst_Nph = [[0 for k in range(2)] for i in range(10)]
step = math.pow(10, -3)
achieved_accuracy = False
while not achieved_accuracy:
    #   Atoms gradF Calculation
    # V - Victim Of Force, C - Criminal Of Force
    for V in range(10):
        for grad_cmp in range(2):
            for div_cmp in range(2):
                if grad_cmp == 0:
                    if div_cmp == 0:
                        #   Increment Coordinate X Of Victim Atom By d/2 To Calculate E1 For gradFx
                        X_crds_Nph[V] += d / 2
                    else:
                        #   Decrement Coordinate X Of Victim Atom By (2 * d/2) To Calculate E2 For gradFx
                        X_crds_Nph[V] -= d
                else:
                    if div_cmp == 0:
                        #   Increment Coordinate Y Of Victim Atom By d/2 To Calculate E1 For gradFy
                        Y_crds_Nph[V] += d / 2
                    else:
                        #   Decrement Coordinate Y Of Victim Atom By (2 * d/2) To Calculate E2 For gradFy
                        Y_crds_Nph[V] -= d
                E_VictimAtom = 0
                #   Va - Victim Affected
                for Va in range(10):
                    for C in range(10):
                        if (Va != C):
                            fcut_vc = f_cut(Va, C, X_crds_Nph, Y_crds_Nph)
                            if fcut_vc != 0:
                                E_VictimAtom += fcut_vc * (E_rep(Va, C, X_crds_Nph, Y_crds_Nph) - b_coef(Va, C, X_crds_Nph, Y_crds_Nph, 10) * E_att(Va, C, X_crds_Nph, Y_crds_Nph))
                if div_cmp == 0:
                    E1 = E_VictimAtom / 2
                else:
                    E2 = E_VictimAtom / 2
            #   Return Coordinates To Initial Values
            if grad_cmp == 0:
                X_crds_Nph[V] += d / 2
            else:
                Y_crds_Nph[V] += d / 2
            gradF = - (E1 - E2) / d
            grad_lst_Nph[V][grad_cmp] = gradF

    #   Coordinates Optimization
    for atm_nmbr in range(10):
        for grd_cmp in range(2):
            if grd_cmp == 0:
                X_crds_Nph[atm_nmbr] += step * grad_lst_Nph[atm_nmbr][grd_cmp]
            else:
                Y_crds_Nph[atm_nmbr] += step * grad_lst_Nph[atm_nmbr][grd_cmp]

    #   Checkpoint: |gradF| > accuracy?
    N_grad_accuracy_achieved = 0
    for atm_nmbr in range(10):
        for grd_cmp in range(2):
            if math.fabs(round(grad_lst_Nph[atm_nmbr][grd_cmp], grad_accuracy_Nph)) <= math.pow(10, -grad_accuracy_Nph):
                N_grad_accuracy_achieved += 1
    if N_grad_accuracy_achieved == 20:
        achieved_accuracy = True
# plt.grid(color='green', linestyle='--', linewidth=1)
# plt.plot(X_crds_Nph, Y_crds_Nph, ls='', c='orange', lw=3, marker="o", markerfacecolor="green", ms=5)
# plt.plot(X_crd_R3, Y_crd_R3, lw=1)
# plt.show()

#   E_Nph - Minimal Energy Of Naphthalene
E_Nph = 0
for V in range(10):
    for C in range(10):
        if V != C:
            fcut_vc = f_cut(V, C, X_crds_Nph, Y_crds_Nph)
            if fcut_vc != 0:
                E_Nph += fcut_vc * (E_rep(V, C, X_crds_Nph, Y_crds_Nph) - b_coef(V, C, X_crds_Nph, Y_crds_Nph, 10) * E_att(V, C, X_crds_Nph, Y_crds_Nph))
E_Nph *= 1/2

#   Etran - Energy Of Transformation
Etran = E_Az - E_Nph

#  3. Direct Road Building, Azulene -> Nephthalene. Draw E(i) Graph, Find Emax To Calculate Eact - Energy Of Activation
N_step = 10**3
energy_lst = []
step_lst = []
X_crds = []
Y_crds = []
for i in range(N_step + 1):
    X_crds.clear()
    Y_crds.clear()
    for atm_nmbr in range(10):
        X_crds.append((i * X_crds_Nph[atm_nmbr] / N_step) + ((N_step - i) * X_crds_Az[atm_nmbr] / N_step))
        Y_crds.append((i * Y_crds_Nph[atm_nmbr] / N_step) + ((N_step - i) * Y_crds_Az[atm_nmbr] / N_step))
    # plt.grid(color='green', linestyle='--', linewidth=1)
    # plt.plot(X_crds_Az, Y_crds_Az, ls='', c='orange', lw=3, marker="o", markerfacecolor="green", ms=6)
    # plt.plot(X_crds_Nph, Y_crds_Nph, ls='', c='orange', lw=3, marker="o", markerfacecolor="yellow", ms=7)
    # plt.plot(X_crds, Y_crds, ls='', c='orange', lw=3, marker="o", markerfacecolor="red", ms=5)
    # plt.show()
    Esys = 0
    for V in range(10):
        for C in range(10):
            if V != C:
                fcut_vc = f_cut(V, C, X_crds, Y_crds)
                if fcut_vc != 0:
                    Esys += fcut_vc * (E_rep(V, C, X_crds, Y_crds) - b_coef(V, C, X_crds, Y_crds, 10) * E_att(V, C, X_crds, Y_crds))
    Esys *= 1/2
    energy_lst.append(Esys)
    step_lst.append(i)
Eact = max(energy_lst) - E_Az

#   4. Transformation Time Estimation When T = 300 K, Temperature Estimation When t = 1 hour
tau = math.pow(10, -12)
k = 1.38 / 1.6 * 10**(-4)
T = 300
t = 1
t1 = tau * math.exp(Eact / (k * T))
t1 /= 3600 * 24 * 365
exp_ten = 0
while (t1 // 10) > 1:
    exp_ten += 1
    t1 /= 10
T1 = Eact / (k * math.log(t / tau))
end = time.time()
print(f"Minimum Energy of Azulene = {round(E_Az, energy_accuracy)} eV")
print(f"Minimum Energy of Naphtalene = {round(E_Nph, energy_accuracy)} eV")
print(f"Energy of Transformation = {round(Etran, energy_accuracy)} eV")
print(f"Energy of Activation = {round(Eact, energy_accuracy)} eV")
print(f"Transformation Time (T = 300 K) = {round(t1)} * 10^{exp_ten} years")
print(f"Transformation Temperature (t = 1 hour) = {round(T1)} K")
print(f"Program Runtime = {round(end - start, energy_accuracy)} sec")
winsound.Beep(500, 2000)
plt.plot(step_lst, energy_lst, ls='-', c='orange', lw=3)
plt.show()
