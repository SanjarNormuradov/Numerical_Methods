import math
import numpy as np
import matplotlib.pyplot as plt

# Function alpha(x)
def alpha(x):
    return math.pow(math.sin(x), 2)

# Function beta(x)
def beta(x):
    return math.cos(4 * x) + 2 * math.sin(x)


func_accuracy = 14
N = 12 * 10**3
dx = 2 * math.pi / (N - 1)
mat_A = [[0 for col in range(N)] for row in range(N)]
mat_A[0][0] = mat_A[0][N-1] = mat_A[1][N-1] = -1
mat_A[0][1] = mat_A[0][N-2] = mat_A[1][0] = 1
for row in range(2, N):
    for col in range(N):
        if (row - 2 == col) or (row == col):
            mat_A[row][col] = 1
        elif row - 1 == col:
            mat_A[row][col] = alpha((row-1)*dx) * dx**2 - 2
col_b = [0 for row in range(N)]
for row in range(2, N):
    col_b[row] = beta((row-1)*dx) * dx**2

inv_mat_A = np.linalg.inv(mat_A)
col_func = np.dot(inv_mat_A, col_b)
found_index = False
low_index = 0
high_index = 0
seg_sum = 0
while not found_index:
    seg_sum += dx
    high_index += 1
    if seg_sum > 1:
        found_index = True
    else:
        low_index += 1
f1 = col_func[high_index] - ((col_func[high_index] - col_func[low_index]) / dx) * (seg_sum - 1)
found_index = False
while not found_index:
    seg_sum += dx
    high_index += 1
    if seg_sum > 2:
        found_index = True
    else:
        low_index += 1
f2 = col_func[high_index] - ((col_func[high_index] - col_func[low_index]) / dx) * (seg_sum - 2)
G = math.pow(col_func[0] - col_func[N-1], 2) + math.pow(((col_func[N-2] - col_func[N-1])/dx) - ((col_func[0] - col_func[1])/dx), 2)
x = [(dx * seg) for seg in range(N)]
plt.plot(x, col_func)
plt.xlabel("x")
plt.ylabel("function")

print(f"A = {round(col_func[0], func_accuracy)}")
print(f"B = {round(col_func[1], func_accuracy)}")
print(f"|f(0) - f(2п)| = {math.fabs(round(col_func[0] - col_func[N-1], func_accuracy))}")
print(f"|f'(0) - f'(2п)| = {math.fabs(round(((col_func[N-2] - col_func[N-1])/dx) - ((col_func[0] - col_func[1])/dx), func_accuracy))}")
print(f"f(1) = {round(f1, func_accuracy)}")
print(f"f(2) = {round(f2, func_accuracy)}")
print(f"G(A,B) = {round(G, func_accuracy)}")
plt.show()

