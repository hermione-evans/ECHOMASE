import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv('dataMoreAmplitude/data_all_new.csv')
# data = pd.read_csv('dataMoreDuration/data_all.csv')
# replace the data_all.csv with data_all_new.csv in dataMoreDuration folder to generate the data for duration.


# generate header and csv file.
# Running this part will overwrite the data got before.
headers = [ 'echoamplitude_injected','duration','i','fmin','fmax','SNR',
           'echoamplitude_searched_new','echoamplitude_searched_old',
           'log_likelihood_new','log_likelihood_old',
           'normalized_echoamplitude_injected','normalized_echoamplitude_searched_new','normalized_echoamplitude_searched_old']
with open('dataMoreAmplitude/data_new.csv', 'w') as f:
    f_csv = csv.writer(f)
    f_csv.writerow(headers)
    
    
    
inject_values = np.unique(data.echoamplitude_injected)
duration_array = np.unique(data.duration)
data_all = data


j=0
for echoamplitude in inject_values:
    for duration in duration_array:
        datanew = data_all[(data_all.echoamplitude_injected == echoamplitude)]
        datanew = datanew[(datanew.duration == duration)]
# This 2 lines select data for each "echoamplitude_injected" and "duration"        
        searched_values = np.unique(datanew.echoamplitude_searched.values)
        i_array = np.unique(datanew.i.values)
        for i in i_array:
            datanew2 = datanew[(datanew.i == i)]
            rows = [echoamplitude,duration,i,np.mean(datanew2.fmin),np.mean(datanew2.fmax),np.mean(datanew2.SNR),
                    datanew2.echoamplitude_searched.values[np.argmax(datanew2.log_likelihood_new)],datanew2.echoamplitude_searched.values[np.argmax(datanew2.log_likelihood_old)],
                    np.max(datanew2.log_likelihood_new),np.max(datanew2.log_likelihood_old),
                    np.mean(datanew2.normalized_echoamplitude_injected),datanew2.normalized_echoamplitude_searched.values[np.argmax(datanew2.log_likelihood_new)],datanew2.normalized_echoamplitude_searched.values[np.argmax(datanew2.log_likelihood_old)]]
    # np.argmax(log_likelihood_new) and np.argmax(log_likelihood_old) select the data with maximum likelihood for different "echoamplitude_searched" 
            with open('dataMoreAmplitude/data_new.csv', 'a+') as f:
                f_csv = csv.writer(f)
                f_csv.writerow(rows)
        print(duration,j,len(i_array),'                    ', end="\r",flush=True)
        j=j+1   
# 这部分是对每个噪声提取出likelihood最大的值，之后取likelihood最大值处的searched_value的标准差作为error
# The part is to extract the maximum likelihood value for each noise, and then take the standard deviation of the searched_value at the maximum likelihood value as the error.


data_all_new = pd.read_csv('dataMoreAmplitude/data_new.csv')
inject_values = np.unique(data_all_new.echoamplitude_injected)
duration_array = np.unique(data_all_new.duration)
def error(x):
    return (np.percentile(x,84)-np.percentile(x,16))/2

# generate header and csv file.
# Running this cell will overwrite the data got before.
headers = ['echoamplitude_injected','duration','SNR',
           'echoamplitude_searched_new','echoamplitude_searched_old','echoamplitude_searched_new_error','echoamplitude_searched_old_error',
           'echoamplitude_searched_new_error_up','echoamplitude_searched_new_error_down','echoamplitude_searched_old_error_up','echoamplitude_searched_old_error_down',
           'echoamplitude_searched_new_min','echoamplitude_searched_new_max','echoamplitude_searched_old_min','echoamplitude_searched_old_max',
           'log_likelihood_new','log_likelihood_old','log_likelihood_new_error','log_likelihood_old_error',
           'log_likelihood_new_error_up','log_likelihood_new_error_down','log_likelihood_old_error_up','log_likelihood_old_error_down',
           'normalized_echoamplitude_injected','normalized_echoamplitude_searched_new','normalized_echoamplitude_searched_old',
           'normalized_echoamplitude_searched_new_min','normalized_echoamplitude_searched_new_max','normalized_echoamplitude_searched_old_min','normalized_echoamplitude_searched_old_max',
          ]
with open('dataMoreAmplitude/data_trend_new.csv', 'w') as f:
    f_csv = csv.writer(f)
    f_csv.writerow(headers)
    
data_all = data_all_new
echoamplitude_array = np.unique(data_all.echoamplitude_injected)
for echoamplitude in echoamplitude_array:
    data = data_all[data_all.echoamplitude_injected == echoamplitude]
    duration_array = np.unique(data.duration)
    for duration in duration_array:
        datanew = data[data.duration == duration]
        rows = [echoamplitude,duration,np.mean(datanew.SNR),
                np.mean(datanew.echoamplitude_searched_new),np.mean(datanew.echoamplitude_searched_old),error(datanew.echoamplitude_searched_new),error(datanew.echoamplitude_searched_old),
                np.percentile(datanew.echoamplitude_searched_new,16),np.percentile(datanew.echoamplitude_searched_new,84),np.percentile(datanew.echoamplitude_searched_old,16),np.percentile(datanew.echoamplitude_searched_old,84),
                np.min(datanew.echoamplitude_searched_new),np.max(datanew.echoamplitude_searched_new),np.min(datanew.echoamplitude_searched_old),np.max(datanew.echoamplitude_searched_old),
                np.mean(datanew.log_likelihood_new),np.mean(datanew.log_likelihood_old),error(datanew.log_likelihood_new),error(datanew.log_likelihood_old),
                np.percentile(datanew.log_likelihood_new,16),np.percentile(datanew.log_likelihood_new,84),np.percentile(datanew.log_likelihood_old,16),np.percentile(datanew.log_likelihood_old,84),
                np.mean(datanew.normalized_echoamplitude_injected),np.mean(datanew.normalized_echoamplitude_searched_new),np.mean(datanew.normalized_echoamplitude_searched_old),
                np.min(datanew.normalized_echoamplitude_searched_new),np.max(datanew.normalized_echoamplitude_searched_new),np.min(datanew.normalized_echoamplitude_searched_old),np.max(datanew.normalized_echoamplitude_searched_old)
               ]
        with open('dataMoreAmplitude/data_trend_new.csv', 'a+') as f:
            f_csv = csv.writer(f)
            f_csv.writerow(rows)