import numpy as np
from glob import glob
import os
## this file is our template

import numpy as np
file_list = glob('../echo_waves/*.dat')
file_list.sort()
file_list = file_list+file_list[18:24]
amplitude_array = np.concatenate((np.linspace(160,160,6),np.linspace(57,57,6),np.linspace(29,29,6),np.linspace(34,34,6),np.linspace(0,0,6)))

echo_number_list = []
for file in file_list:
    echo_number_list.append(int(file[file.find(".dat")-3:file.find(".dat")]))
echo_number_list = np.array(echo_number_list)
name_list = []
for file in file_list:
    name_list.append(file[file.find('w1'):file.find('.dat')])
name_list = np.array(name_list)
likelihoodlist = ['newlikelihood','oldlikelihood']

file = 'signal_tem.py' 
slurm_file = 'run_tem.slurm'
test = 0
index_file_data = ""
benchmarkindex = 0
run_file_data = ""
python_file_number_list = np.arange(0,50,1)
number_of_noise_realization = 100
sbatch_file_data = ""
clear_file_data =""
for index in np.arange(0,30,1):
    # Here the index is used to label the inject signal file (5 models with 6 durtaions)
    str4 = 'index = %d'%index
    for i in np.arange(0,2):
        # Here the i is used to label the likelihood function (newlikelihood or oldlikelihood)
        str5 =  'likelihood_index = %d'%i
        amplitude_inject = amplitude_array[index]
        if amplitude_inject !=0:
            name = name_list[index]
        else:
            file_name = name_list[index]
            name = 'noise'+file_name[file_name.find('N'):file_name.find('N')+4]
        test = 0
        if not os.path.exists("../python_files_"+name):
            os.mkdir("../python_files_"+name)
            line = "rm -rf python_files_"+name+'\n'
            clear_file_data = clear_file_data + line
        for n in python_file_number_list:
            # Here the n is used to label the number of the python file. It is also the number of CPU cores used in the slurm file.
            # the length of python_file_number_list should be evenly divisible by the number_of_noise_realization
            str1 = 'start = %d'%(n*(number_of_noise_realization//len(python_file_number_list)))
            str2 = 'end = %d'%(n+1)*(number_of_noise_realization//len(python_file_number_list))
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
                    if 'index = tem' in line:
                        line = line.replace('index = tem',str4)
                    if 'likelihood_index=tem' in line:
                        line = line.replace('likelihood_index=tem',str5)
                    test = 1
                    file_data += line
            print(n,str1,str2,str3,str4,str5,name+'_'+likelihoodlist[i]+'_%d.py'%n)
            
            with open("../python_files_"+name+'/'+name+'_'+likelihoodlist[i]+'_%d.py'%n,"w",encoding="utf-8") as f:
                f.write(file_data)
        slurm_file_data = ""
        with open(slurm_file, "r", encoding="utf-8") as f:
            for line in f:
                if '#SBATCH --array=tem' in line:
                    line = line.replace('#SBATCH --array=tem','#SBATCH --array=%d-%d'%(python_file_number_list[0],python_file_number_list[-1]))
                if 'srun python tem${SLURM_ARRAY_TASK_ID}.py' in line:
                    line = line.replace('srun python tem${SLURM_ARRAY_TASK_ID}.py','srun python '+name+'_'+likelihoodlist[i]+'_${SLURM_ARRAY_TASK_ID}.py')
                slurm_file_data = slurm_file_data + line
        with open("../python_files_"+name+'/'+name+'_'+likelihoodlist[i]+'.slurm',"w",encoding="utf-8") as f:
            f.write(slurm_file_data)
        linebefore = 'cd python_files_'+name+'\n'
        lineafter = 'cd .. \n'
        sbatch_file_data = sbatch_file_data + linebefore + 'sbatch '+name+'_'+likelihoodlist[i]+'.slurm\n' + lineafter

with open('../run.sh',"w",encoding="utf-8") as f:
    f.write(sbatch_file_data)    
with open('../clear_python_files.sh',"w",encoding="utf-8") as f:
    f.write(clear_file_data)    
  
