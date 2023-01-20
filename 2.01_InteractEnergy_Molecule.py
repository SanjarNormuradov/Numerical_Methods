import math
import time


#   Distance Between Two Atoms
def r_dist(i, j, x_crd, y_crd, z_crd):
    r = math.sqrt((x_crd[i] - x_crd[j])**2 + (y_crd[i] - y_crd[j])**2 + (z_crd[i] - z_crd[j])**2)
    return r


#   Angle Between 2 Vectors Connecting 3 Atoms
def cos_theta(i, j, k, x_crd, y_crd, z_crd):
    V1_x = x_crd[j] - x_crd[i]
    V1_y = y_crd[j] - y_crd[i]
    V1_z = z_crd[j] - z_crd[i]
    V2_x = x_crd[k] - x_crd[i]
    V2_y = y_crd[k] - y_crd[i]
    V2_z = z_crd[j] - z_crd[i]
    V1 = math.sqrt(V1_x**2 + V1_y**2 + V1_z**2)
    V2 = math.sqrt(V2_x**2 + V2_y**2 + V2_z**2)
    cos0 = ((V1_x * V2_x) + (V1_y * V2_y) + (V1_z * V2_z))/(V1 * V2)
    return cos0


#   Cut Function
def f_cut(i, j, x_crd, y_crd, z_crd):
    R = 2.00
    D = 0.15
    r = r_dist(i, j, x_crd, y_crd, z_crd)
    # if r < R:
    #     func_value = 1
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
def b_coef(i, j, x_crd, y_crd, z_crd, n):
    khi_ij = 0
    for k in range(n):
        if (k != i) and (k != j):
            fcut_ik = f_cut(i, k, x_crd, y_crd, z_crd)
            if fcut_ik != 0:
                cos_theta_ijk = cos_theta(i, j, k, x_crd, y_crd, z_crd)
                khi_ij += fcut_ik * g_coef(cos_theta_ijk)
    b = math.pow(1 + khi_ij, -1/2)
    return b


#   Repulsion Energy
def E_rep(i, j, x_crd, y_crd,  z_crd):
    D0 = 6
    r0 = 1.4276
    S = 2.167
    beta = 2.0099
    r = r_dist(i, j, x_crd, y_crd, z_crd)
    E = D0/(S-1) * math.exp(-beta * math.sqrt(2*S) * (r - r0))
    return E


#   Attraction Energy
def E_att(i, j, x_crd, y_crd, z_crd):
    D0 = 6
    r0 = 1.4276
    S = 2.167
    beta = 2.0099
    r = r_dist(i, j, x_crd, y_crd, z_crd)
    E = D0*S/(S-1) * math.exp(-beta * math.sqrt(2/S) * (r - r0))
    return E


#   System Energy
def System_Energy(Rcircle):
    var = 5
    n = var + 2
    x_crds = []
    y_crds = []
    z_crds = []
    for i in range(n):
        x_crds.append(Rcircle * math.cos(2 * math.pi * i / n))
        y_crds.append(Rcircle * math.sin(2 * math.pi * i / n))
        z_crds.append(0)
    Sys_E = 0
    for i in range(n):
        for j in range(n):
            if j != i:
                fcut_ij = f_cut(i, j, x_crds, y_crds, z_crds)
                if fcut_ij != 0:
                    Sys_E += fcut_ij * (E_rep(i, j, x_crds, y_crds, z_crds) - b_coef(i, j, x_crds, y_crds, z_crds, n) * E_att(i, j, x_crds, y_crds, z_crds))
    Sys_E *= 1/2
    return Sys_E


#   Calculate Interaction Energy From File Data Input
def File_Input_Energy_Output(file_name):
    Ax = 7
    Ay = 0
    Az = 0
    Bx = 0
    By = 7
    Bz = 0
    Cx = 0
    Cy = 0
    Cz = 7
    crd_accuracy = 4
    Energy_accuracy = 5
    with open(file_name + ".txt", 'r') as Molecule:
        N_atoms = int(Molecule.readline())
        Molecule_Name = Molecule.readline().split()
        x_crds = []
        y_crds = []
        z_crds = []
        for i in range(N_atoms):
            atom_info_list = Molecule.readline().split()
            #   Get X coordinate With Space
            x_crds.append(round(float(atom_info_list[1]), crd_accuracy))
            #   Get Y coordinate With Space
            y_crds.append(round(float(atom_info_list[2]), crd_accuracy))
            #   Get Z coordinate With Space
            z_crds.append(round(float(atom_info_list[3]), crd_accuracy))
        Sys_E = 0
        for i in range(N_atoms):
            for j in range(N_atoms):
                if j != i:
                    #   Minimum Distance Search
                    r_min = 10
                    for u in [-1, 0, 1]:
                        for v in [-1, 0, 1]:
                            for w in [-1, 0, 1]:
                                x_crds[j] += u * Ax + v * Bx + w * Cx
                                y_crds[j] += u * Ay + v * By + w * Cy
                                z_crds[j] += u * Az + v * Bz + w * Cz
                                r_tmp = r_dist(i, j, x_crds, y_crds, z_crds)
                                if r_tmp < r_min:
                                    r_min = r_tmp
                                    u_min = u
                                    v_min = v
                                    w_min = w
                                #   Return Coordinates To Initial Values
                                x_crds[j] -= u * Ax + v * Bx + w * Cx
                                y_crds[j] -= u * Ay + v * By + w * Cy
                                z_crds[j] -= u * Az + v * Bz + w * Cz
                    #   Change Coordinates To Found Minimum
                    x_crds[j] += u_min * Ax + v_min * Bx + w_min * Cx
                    y_crds[j] += u_min * Ay + v_min * By + w_min * Cy
                    z_crds[j] += u_min * Az + v_min * Bz + w_min * Cz
                    fcut_ij = f_cut(i, j, x_crds, y_crds, z_crds)
                    if fcut_ij != 0:
                        Sys_E += fcut_ij * (E_rep(i, j, x_crds, y_crds, z_crds) - b_coef(i, j, x_crds, y_crds, z_crds, N_atoms) * E_att(i, j, x_crds, y_crds, z_crds))
                    #   Return Coordinates To Initial Values
                    x_crds[j] -= u_min * Ax + v_min * Bx + w_min * Cx
                    y_crds[j] -= u_min * Ay + v_min * By + w_min * Cy
                    z_crds[j] -= u_min * Az + v_min * Bz + w_min * Cz
        Sys_E *= 1/2
    print(f"Energy of Interaction Between Atoms in {Molecule_Name[0]} = {round(Sys_E, Energy_accuracy)}")
    return None


# def System_Energy_trimer():
#     n = 3
#     x_crds = [0.0, 0.0, 1.4276]
#     y_crds = [1.4276, 0.0, 0.0]
#     Sys_E = 0
#     for i in range(n):
#         for j in range(n):
#             if j != i:
#                 fcut_ij = f_cut(i, j, x_crds, y_crds)
#                 if fcut_ij != 0:
#                     Sys_E += fcut_ij * (E_rep(i,j,x_crds,y_crds) - b_coef(i,j,x_crds,y_crds,n) * E_att(i,j,x_crds,y_crds))
#     Sys_E *= 1/2
#     return Sys_E


start = time.time()
r_accuracy = 3
Energy_accuracy = 5
var = 5
n = var + 2

#   1.  Read .txt file To Calculate Interaction Energy
# File_Input_Energy_Output("TwoAtomMolecule")
# File_Input_Energy_Output("ThreeAtomMolecule")
#   2.  Optimal Radius Of Regular Polygon
achieved_accuracy = False
previous_Energy_calculated = False
R_circle = 2.00
dR = math.pow(10, -r_accuracy)
step = math.pow(10, -r_accuracy)
while not achieved_accuracy:
    if not previous_Energy_calculated:
        System_Energy1 = System_Energy(R_circle)
        previous_Energy_calculated = True
    R_circle -= dR
    System_Energy2 = System_Energy(R_circle)
    divR = (System_Energy2 - System_Energy1) / dR
    R_circle += step * divR
    System_Energy2 = System_Energy(R_circle)
    if math.fabs(System_Energy1 - System_Energy2) > math.pow(10, -Energy_accuracy):
        System_Energy1 = System_Energy2
    else:
        achieved_accuracy = True
Optimal_Length = 2 * R_circle * math.sin(math.pi/n)
print(f"Optimal Length of Regular Polygon = {round(Optimal_Length, r_accuracy)}")
print(R_circle)
end = time.time()
print(f"Program Runtime = {end - start}")
