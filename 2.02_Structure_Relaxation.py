import math
import time
import winsound

#   Distance Between Two Atoms
def r_dist(vi, vj, ci, cj, x_crd, y_crd, z_crd):
    if cj < 3:
        cjm = cj
    elif (cj > 2) and (cj < 6):
        cjm = cj - 3
    else:
        cjm = cj - 6
    r = math.sqrt((x_crd[vj][vi]-x_crd[cjm][ci])**2+(y_crd[vj][vi]-y_crd[cjm][ci])**2+(z_crd[vj][vi]-z_crd[cj][ci])**2)
    return r


#   Angle Between 2 Vectors Connecting 3 Atoms
def cos_theta(vi, vj, ci, cj, wi, wj, x_crd, y_crd, z_crd):
    #   Criminal Coordinates X,Y Limit Check
    if cj < 3:
        cjm = cj
    elif (cj > 2) and (cj < 6):
        cjm = cj - 3
    else:
        cjm = cj - 6
    #   Witness Coordinates X,Y Limit Check
    if wj < 3:
        wjm = wj
    elif (wj > 2) and (wj < 6):
        wjm = wj - 3
    else:
        wjm = wj - 6
    #   Vector1, Victim and Criminal Atoms
    V1_x = x_crd[cjm][ci] - x_crd[vj][vi]
    V1_y = y_crd[cjm][ci] - y_crd[vj][vi]
    V1_z = z_crd[cj][ci] - z_crd[vj][vi]
    #   Vector2, Victim and Witness Atoms
    V2_x = x_crd[wjm][wi] - x_crd[vj][vi]
    V2_y = y_crd[wjm][wi] - y_crd[vj][vi]
    V2_z = z_crd[wj][wi] - z_crd[vj][vi]
    # Vectors' Module
    V1 = math.sqrt(V1_x**2 + V1_y**2 + V1_z**2)
    V2 = math.sqrt(V2_x**2 + V2_y**2 + V2_z**2)
    cos0 = ((V1_x * V2_x) + (V1_y * V2_y) + (V1_z * V2_z))/(V1 * V2)
    return cos0


#   Cut Function
def f_cut(vi, vj, ci, cj, x_crd, y_crd, z_crd):
    R = 2.00
    D = 0.15
    r = r_dist(vi, vj, ci, cj, x_crd, y_crd, z_crd)
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
def b_coef(vi, vj, ci, cj, x_crd, y_crd, z_crd, n):
    khi_ij = 0
    # W - Witness
    for wj in range(9):
        for wi in range(n):
            if ((wi != vi) or (wj != vj)) and ((wi != ci) or (wj != cj)):
                fcut_ik = f_cut(vi, vj, wi, wj, x_crd, y_crd, z_crd)
                if fcut_ik != 0:
                    cos_theta_ijk = cos_theta(vi, vj, ci, cj, wi, wj, x_crd, y_crd, z_crd)
                    khi_ij += fcut_ik * g_coef(cos_theta_ijk)
    b = math.pow(1 + khi_ij, -1/2)
    return b


#   Repulsion Energy
def E_rep(vi, vj, ci, cj, x_crd, y_crd, z_crd):
    D0 = 6
    r0 = 1.4276
    S = 2.167
    beta = 2.0099
    r = r_dist(vi, vj, ci, cj, x_crd, y_crd, z_crd)
    E = D0/(S-1) * math.exp(-beta * math.sqrt(2*S) * (r - r0))
    return E


#   Attraction Energy
def E_att(vi, vj, ci, cj, x_crd, y_crd, z_crd):
    D0 = 6
    r0 = 1.4276
    S = 2.167
    beta = 2.0099
    r = r_dist(vi, vj, ci, cj, x_crd, y_crd, z_crd)
    E = D0*S/(S-1) * math.exp(-beta * math.sqrt(2/S) * (r - r0))
    return E


start = time.time()
result_accuracy = 3
force_accuracy = 4
R_OuterCircle = 1.6710941872164469
l2 = 0.5
var = 5
N = var + 2

# Translational Vectors And Origin
TransVectOr_X = [[0,    0,      0,      10], [0, 0, 0, 0]]
TransVectOr_Y = [[0,    0,      10,     0], [0, 0, 0, 0]]
TransVectOr_Z = [[0,    l2,    0,       0], [0, 0, 0, 0]]

X_crds = [[0 for i in range(N)] for j in range(3)]
Y_crds = [[0 for i in range(N)] for j in range(3)]
Z_crds = [[0 for i in range(N)] for j in range(9)]
#   Set Initial Coordinates Of Mother Cell
for j in range(3):
    for i in range(N):
        X_crds[j][i] = R_OuterCircle * math.cos(2 * math.pi * i / N)
        Y_crds[j][i] = R_OuterCircle * math.sin(2 * math.pi * i / N)
        Z_crds[j][i] = j * l2
        #   Upper Daughter Cell Coordinate Z
        Z_crds[j+3][i] = j * l2 + 3 * float(TransVectOr_Z[0][1])
        #   Lower Daughter Cell Coordinate Z
        Z_crds[j+6][i] = j * l2 - 3 * float(TransVectOr_Z[0][1])

grad_lst = [[[0 for k in range(3)] for i in range(N)] for j in range(3)]
d = math.pow(10, -6)
step_k = 1
step = step_k * math.pow(10, -2)
achieved_accuracy = False
Efin = 0

while not achieved_accuracy:
    #   Atoms gradF Calculation
    # V - Victim Of Force, C - Criminal Of Force
    for Vj in range(3):
        for Vi in range(N):
            for grad_cmp in range(3):
                for div_cmp in range(2):
                    if grad_cmp == 0:
                        if div_cmp == 0:
                            #   Increment Coordinate X Of Victim Atom By d/2 To Calculate E1 For gradFx
                            X_crds[Vj][Vi] += d / 2
                        else:
                            #   Decrement Coordinate X Of Victim Atom By (2 * d/2) To Calculate E2 For gradFx
                            X_crds[Vj][Vi] -= d
                    elif grad_cmp == 1:
                        if div_cmp == 0:
                            #   Increment Coordinate Y Of Victim Atom By d/2 To Calculate E1 For gradFy
                            Y_crds[Vj][Vi] += d / 2
                        else:
                            #   Decrement Coordinate Y Of Victim Atom By (2 * d/2) To Calculate E2 For gradFy
                            Y_crds[Vj][Vi] -= d
                    else:
                        if div_cmp == 0:
                            #   Increment Coordinate Z Of Victim Atom By d/2 To Calculate E1 For gradFz
                            Z_crds[Vj][Vi] += d / 2
                        else:
                            #   Decrement Coordinate Z Of Victim Atom By (2 * d/2) To Calculate E2 For gradFz
                            Z_crds[Vj][Vi] -= d
                    E_VictimAtom = 0
                    for Cj in range(9):
                        for Ci in range(N):
                            if (Vj != Cj) or (Vi != Ci):
                                fcut_ij = f_cut(Vi, Vj, Ci, Cj, X_crds, Y_crds, Z_crds)
                                if fcut_ij != 0:
                                    E_VictimAtom += fcut_ij * (E_rep(Vi, Vj, Ci, Cj, X_crds, Y_crds, Z_crds) - b_coef(Vi, Vj, Ci, Cj, X_crds, Y_crds, Z_crds, N) * E_att(Vi, Vj, Ci, Cj, X_crds, Y_crds, Z_crds))
                    if div_cmp == 0:
                        E1 = E_VictimAtom
                    else:
                        E2 = E_VictimAtom
                        Efin = E2
                #   Return Coordinates To Initial Values
                if grad_cmp == 0:
                    X_crds[Vj][Vi] += d / 2
                elif grad_cmp == 1:
                    Y_crds[Vj][Vi] += d / 2
                else:
                    Z_crds[Vj][Vi] += d / 2
                gradF = - (E1 - E2) / d
                grad_lst[Vj][Vi][grad_cmp] = gradF

    #   Boundary Condition gradFa Calculation
    Orig_j = 0
    Orig_i = 0
    VectAz_j = 0
    VectAz_i = 1
    for div_cmp in range(2):
        if div_cmp == 0:
            #   Increment Vector A Component Az By d/2 To Calculate E1 For gradFa
            TransVectOr_Z[VectAz_j][VectAz_i] += d / 2
        else:
            #   Decrement Vector A Component Az By d To Calculate E2 For gradFa
            TransVectOr_Z[VectAz_j][VectAz_i] -= d
        E_BoundCond = 0
        fcut_ij = f_cut(Orig_i, Orig_j, VectAz_i, VectAz_j, TransVectOr_X, TransVectOr_Y, TransVectOr_Z)
        if fcut_ij != 0:
            E_BoundCond = fcut_ij*(
                    E_rep(Orig_i, Orig_j, VectAz_i, VectAz_j, TransVectOr_X, TransVectOr_Y, TransVectOr_Z) - E_att(Orig_i, Orig_j, VectAz_i, VectAz_j, TransVectOr_X, TransVectOr_Y, TransVectOr_Z))
        if div_cmp == 0:
            E1 = E_BoundCond
        else:
            E2 = E_BoundCond
    #   Return Vector A Component Az To Initial Value
    TransVectOr_Z[VectAz_j][VectAz_i] += d / 2
    gradFa = - (E1 - E2) / d

    #   Coordinates Optimization
    for atms_lyr in range(3):
        for atm_nmbr in range(N):
            for grd_cmp in range(3):
                if grd_cmp == 0:
                    X_crds[atms_lyr][atm_nmbr] += step * grad_lst[atms_lyr][atm_nmbr][grd_cmp]
                elif grd_cmp == 1:
                    Y_crds[atms_lyr][atm_nmbr] += step * grad_lst[atms_lyr][atm_nmbr][grd_cmp]
                else:
                    Z_crds[atms_lyr][atm_nmbr] += step * grad_lst[atms_lyr][atm_nmbr][grd_cmp]
                    Z_crds[atms_lyr + 3][atm_nmbr] += step * grad_lst[atms_lyr][atm_nmbr][grd_cmp]
                    Z_crds[atms_lyr + 6][atm_nmbr] += step * grad_lst[atms_lyr][atm_nmbr][grd_cmp]
    TransVectOr_Z[VectAz_j][VectAz_i] += step * gradFa

    #   Checkpoint: |gradF| > accuracy?
    N_grad_accuracy_achieved = 0
    for atms_lyr in range(3):
        for atm_nmbr in range(N):
            for grd_cmp in range(3):
                if math.fabs(round(grad_lst[atms_lyr][atm_nmbr][grd_cmp], force_accuracy)) <= math.pow(10, -force_accuracy):
                    N_grad_accuracy_achieved += 1
    if math.fabs(round(gradFa, force_accuracy)) <= math.pow(10, -force_accuracy):
        N_grad_accuracy_achieved += 1
    if N_grad_accuracy_achieved == (9 * N + 1):
        achieved_accuracy = True
    # else:
    #     step = (1 - N_grad_accuracy_achieved / (9 * N + 1)) * step_k * math.pow(10, -3)

print(f"Optimal Distance Between Two Neighbor Atoms = {round(r_dist(2, 0, 3, 0, X_crds, Y_crds, Z_crds), result_accuracy)}")
print(f"Optimal Distance Between Two Neighbor Atoms = {round(r_dist(0, 0, 0, 1, X_crds, Y_crds, Z_crds), result_accuracy)}")
print(f"Optimal Length Of Translational Vector A = {round(float(TransVectOr_Z[0][1]), result_accuracy)}")
print(f"Interaction Energy Of An Atom = {round(Efin, result_accuracy)}")
end = time.time()
print(f"Program Runtime = {round(end - start, result_accuracy)}")
winsound.Beep(500, 2000)
