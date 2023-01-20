import time
from random import randint

# String To Integer Function
def str_to_int(tmp_str):
    return int(tmp_str)

# Create Text File
def create_txt_file(file_name):
    friend_list = []
    friend_sum = 0
    # Erase The Content of Text File Before Appending
    # open(file_name + ".txt", 'w').close()
    with open(file_name + ".txt", 'w') as file:
        for k in range(10**1):
            N_friends = randint(1, 5)
            friend_list.append(N_friends)
            friend_sum += N_friends
        file.write(f"{friend_sum} \n")
        for i in range(10**1):
            file.write('---------------\n')
            found_friends_list = []
            for k in range(friend_list[i]):
                found_friend = False
                while not found_friend:
                    j = randint(1, 10**1)
                    if (j != i+1):
                        if (k == 0) or (j not in found_friends_list):
                            file.write(f"{i+1} {j} 1 \n")
                            found_friends_list.append(j)
                            found_friend = True
    return friend_list

start = time.time()
M1_friend_list = create_txt_file("Matrix1")
M2_friend_list = create_txt_file("Matrix2")
M1 = open("Matrix1.txt", 'r')
M2 = open("Matrix2.txt", 'r')

new_file_list = []
tmp_column_list = []
tmp_row_list = []
N_friends_list = []
new_file_friend_sum = 0
n_row = 0

friend_sum1 = int(M1.readline().strip())
for cnt1 in range(friend_sum1):
    i_row = ""
    get_whole_number = False
    while not get_whole_number:
        str_tmp = M1.readline(1)
        if str_tmp == '-':
            M1.readline()
        elif str_tmp != ' ':
            i_row += str_tmp
        else:
            get_whole_number = True
    k_tmp1 = ""
    get_whole_number = False
    while not get_whole_number:
        str_tmp = M1.readline(1)
        if str_tmp != ' ':
            k_tmp1 += str_tmp
        else:
            get_whole_number = True
    if i_row not in tmp_row_list:
        if len(tmp_row_list) != 0:
            new_file_list.append('------------------\n')
            for j_column in tmp_column_list:
                new_file_list.append(f"{tmp_row_list[n_row-1]} {j_column} {N_friends_list[tmp_column_list.index(j_column)]} \n")
        tmp_row_list.append(i_row)
        n_row += 1
        tmp_column_list.clear()
        N_friends_list.clear()
    if cnt1 != (friend_sum1 - 1):
        # Start Searching From The Next Line For The Next Time
        M1.readline()

    # Skip Number Of Friends in M2
    M2.readline()
    search_start = 0
    cnt3 = 0
    # Sum Of Elements In M2_friend_list Till M2_friend_list[k_tmp1] = Number Of Unnecessary Elements For Comparison
    for cnt3 in range(int(k_tmp1)-1):
        search_start += int(M2_friend_list[cnt3])
    for n in range(search_start + cnt3 + 2):
        M2.readline()
    for cnt2 in range(M2_friend_list[cnt3+1]):
        k_tmp2 = ""
        get_whole_number = False
        while not get_whole_number:
            str_tmp = M2.readline(1)
            if str_tmp == '-':
                M2.readline()
            elif str_tmp != ' ':
                k_tmp2 += str_tmp
            else:
                get_whole_number = True
        if int(k_tmp2) == int(k_tmp1):
            j_column = ""
            get_whole_number = False
            while not get_whole_number:
                str_tmp = M2.readline(1)
                if str_tmp != ' ':
                    j_column += str_tmp
                else:
                    get_whole_number = True
            if j_column not in tmp_column_list:
                tmp_column_list.append(j_column)
                N_friends_list.append(1)
                new_file_friend_sum += 1
            else:
                if cnt2 < (M2_friend_list[cnt3+1] - 1):
                    inc_element_index = tmp_column_list.index(j_column)
                    N_friends_list[inc_element_index] += 1

        if cnt2 == (M2_friend_list[cnt3+1] - 1):
            # Start Searching From The First Line In The Next Cycle
            M2.seek(0)
        else:
            # Start Searching From The Next Line
            M2.readline()
# If M1 Has Rows Of Only One Value
if n_row == 1:
    new_file_list.append('------------------\n')
    tmp_column_list.sort(key=str_to_int)
    for j_column in tmp_column_list:
        new_file_list.append(f"{tmp_row_list[n_row-1]} {j_column} {N_friends_list[tmp_column_list.index(j_column)]} \n")

# M1.seek(0)
# print(M1.read())
# print(M2.read())
M2.close()
M1.close()
with open("Matrix.txt", 'w') as M:
    M.write(f"{new_file_friend_sum} \n")
    M.writelines(new_file_list)
# with open("Matrix.txt", 'r') as M:
#     print(M.read())
end = time.time()
print(f"Program Runtime = {end - start}")
