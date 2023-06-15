import numpy as np
file = 'signal_tem.py' 
## this file is our template

import numpy as np

test = 0
for n in np.arange(0,50,1):
    str1 = 'start = %d'%(n*2)
    str2 = 'end = %d'%(n*2+2)
    # Here we employ 50 Python scripts to simultaneously complete 100 noise realizations in parallel.
    str3 = 'whether_print_header = 1'
    if test==0:
        str3 = 'whether_print_header = 0'
    file_data = ""
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if 'start = tem' in line:
                line = line.replace('start = tem',str1)
            if 'end = tem' in line:
                line = line.replace('end = tem',str2)
            if 'whether_print_header = 0' in line:
               line = line.replace('whether_print_header = 0',str3)
            test = 1
            file_data += line
    with open('signal_%d.py'%n ,"w",encoding="utf-8") as f:
        f.write(file_data)
    print(n,str1,str2,str3)
