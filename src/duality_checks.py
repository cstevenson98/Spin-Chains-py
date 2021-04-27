import subprocess
import os
import sys

# This is a script which will check the 6 possible set ups for the
# structured bath xy_model code.

N = int(sys.argv[1])
abs_g = float(sys.argv[2])
abs_omega0 = float(sys.argv[3])
b = float(sys.argv[4])

flag = sys.argv[5]

Lambda = float(sys.argv[6])

I = int(sys.argv[7])

if I == 1:
    #
    my_list = [N, abs_g, abs_omega0, b, flag, Lambda]
    my_list = list(map(str, my_list))
    subprocess.call(['python', 'main.py'] + my_list)
elif I == 2:
    #
    my_list = [N, abs_g, 0.0, b, flag, Lambda]
    my_list = list(map(str, my_list))
    subprocess.call(['python', 'main.py'] + my_list)
elif I == 3:
    #
    my_list = [N, abs_g, -abs_omega0, b, flag, Lambda]
    my_list = list(map(str, my_list))
    subprocess.call(['python', 'main.py'] + my_list)
elif I == 4:
    #
    my_list = [N, -abs_g, abs_omega0, b, flag, Lambda]
    my_list = list(map(str, my_list))
    subprocess.call(['python', 'main.py'] + my_list)
elif I == 5:
    #
    my_list = [N, -abs_g, 0.0, b, flag, Lambda]
    my_list = list(map(str, my_list))
    subprocess.call(['python', 'main.py'] + my_list)
elif I == 6:
    #
    my_list = [N, -abs_g, -abs_omega0, b, flag, Lambda]
    my_list = list(map(str, my_list))
    subprocess.call(['python', 'main.py'] + my_list)
    
