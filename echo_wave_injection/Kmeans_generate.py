import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import bilby
from sklearn.cluster import MiniBatchKMeans
from glob import glob

duration_array_str = glob("echo_waves/*.dat")
duration_array_str.sort()
for i in range(len(duration_array_str)):
    duration_array_str[i] = duration_array_str[i].split("/")[1].split("N")[1].split('.')[0]
duration_array_str = list(dict.fromkeys(duration_array_str))

duration_array = glob("echo_waves/*.dat")
duration_array.sort()
for i in range(len(duration_array)):
    duration_array[i] = int(duration_array[i].split("/")[1].split("N")[1].split('.')[0])
duration_array = list(dict.fromkeys(duration_array))
duration_array.sort()
duration_array = np.array(duration_array)

name_array = glob("echo_waves/*.dat")
name_array.sort()
for i in range(len(name_array)):
    name_array[i] = name_array[i].split("/")[1].split("N")[0]
#delete duplicates in name array
name_array = list(dict.fromkeys(name_array))
name_array.sort()
name_array = name_array[::-1]
name_array.append('noise')

file_array_new = [sorted(glob('%s/%s/*%s*SNR*.json'%(name,string_likelihood,'N'+duration))) for name in name_array for string_likelihood in ['newlikelihood'] for duration in duration_array_str]
file_array_old = [sorted(glob('%s/%s/*%s*SNR*.json'%(name,string_likelihood,'N'+duration))) for name in name_array for string_likelihood in ['oldlikelihood'] for duration in duration_array_str]

names_array = ['B1','B1','B1','B1','B1','B1'
 ,'B2', 'B2', 'B2', 'B2', 'B2', 'B2'
 ,'B3', 'B3', 'B3', 'B3', 'B3', 'B3'
 ,'B4', 'B4', 'B4', 'B4', 'B4', 'B4'
 ,'noise', 'noise', 'noise', 'noise', 'noise', 'noise']


for file_array_new_2,file_array_old_2,name in zip(file_array_new,file_array_old,names_array):
    print(len(file_array_new_2),len(file_array_old_2))
    print(file_array_new_2[0])

    result = bilby.result.read_in_result(file_array_new_2[0])
    priors = bilby.result.read_in_result(file_array_new_2[0]).priors
    posterior_array_new = [bilby.result.read_in_result(result_single).posterior for result_single in file_array_new_2]
    posterior_array_old = [bilby.result.read_in_result(result_single).posterior for result_single in file_array_old_2]
    posterior_array_allplus_new = pd.concat(posterior_array_new,ignore_index=True)
    posterior_array_allplus_old = pd.concat(posterior_array_old,ignore_index=True)
    priors = bilby.result.read_in_result(file_array_new_2[0]).priors

    outdir = file_array_new_2[0].split("/")[0]
    
    priors.to_json(outdir=outdir,label=result.label.replace('i0_SNR','_Kmeans'))
    
    # posterior_new = pd.DataFrame(columns=['width', 'amplitude', 'phase', 'spacing', 'fmin', 'fmax'])
    # posterior_old = pd.DataFrame(columns=['width', 'amplitude', 'phase', 'spacing', 'fmin', 'fmax'])
    
    # for keys in ['width', 'amplitude', 'phase', 'spacing', 'fmin', 'fmax']:
    #     posterior_new[keys] = MiniBatchKMeans(n_clusters=10000).fit(posterior_array_allplus_new[keys].values.reshape(-1, 1)).cluster_centers_.reshape(-1)
    #     posterior_old[keys] = MiniBatchKMeans(n_clusters=10000).fit(posterior_array_allplus_old[keys].values.reshape(-1, 1)).cluster_centers_.reshape(-1)
    # posterior_new.to_json(result.outdir.replace('/newlikelihood','')+'/'+result.label.replace('i0_SNR','_Kmeans.json'))
    # posterior_old.to_json(result.outdir.replace('/newlikelihood','')+'/'+result.label.replace('i0_SNR','_Kmeans.json').replace('newlikelihood','oldlikelihood'))
    
    posterior_new = pd.read_json(outdir+'/'+result.label.replace('i0_SNR','_Kmeans.json'))
    posterior_old = pd.read_json(outdir+'/'+result.label.replace('i0_SNR','_Kmeans.json').replace('newlikelihood','oldlikelihood'))
    
    keys = 'SNR'
    posterior_new[keys] = MiniBatchKMeans(n_clusters=10000).fit(posterior_array_allplus_new[keys].values.reshape(-1, 1)).cluster_centers_.reshape(-1)
    posterior_old[keys] = MiniBatchKMeans(n_clusters=10000).fit(posterior_array_allplus_old[keys].values.reshape(-1, 1)).cluster_centers_.reshape(-1)
    posterior_new.to_json(outdir+'/'+result.label.replace('i0_SNR','_Kmeans.json'))
    posterior_old.to_json(outdir+'/'+result.label.replace('i0_SNR','_Kmeans.json').replace('newlikelihood','oldlikelihood'))
    
