import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import scipy.stats

def qnmcomb_cut_1pole(height,loc,width,duration,x):
    ##高度、相位、峰宽度、单个周期长度、频率（输入的x）
    # cut 表示在归一化的峰高小于cut时截断，目前暂时取的是0.2
    '''
    the first parameter is the height of wave.
    the second parameter is phase, should between 0 to spacing.
    When frequency = phase + n * spacing, the comb reaches its peaks.
    As the phase grows larger, the comb goes right.
    the third parameter is width，which describes the width of signal.
    the fourth parameter is spacing, which describes the whole length of a wave.
    the fifth parameter is the time duration of the starin time-domain data.
    the sixth parameter is input frequency
    The unit of frequency is Hz.
    '''
    # omega_p = ((x - phase) / spacing -np.floor((x - phase) / spacing + 1/2) ) *spacing *2 *np.pi
    # here the range of omega_p is [-spacing/2 , spacing/2]*2 *np.pi, the length of omega_p is spacing*2*pi
    # when x = phase ,omega_p = 0
    omega_p = (x - loc) *2 *np.pi
    cut = 0.1
    cutx1 = np.sqrt(1 - cut**2) * width / cut
    cutx2 = 3/duration*2 *np.pi
# here this number*2 means the smallest number for bins used in one cycle.
    # cutx3 = spacing*np.pi
    cutx = np.max([cutx1,cutx2])
# cutx > 2 *pi width/T
    #res = -1j*width /(omega_p-1j*width)
    # res means resonance structure
    #T_dep = 1 - np.exp(-duration*width*(1-1j*omega_p/width))
    # T_dep means the T dependence factor
    y = np.piecewise(omega_p, [omega_p <= - cutx, np.logical_and(-cutx < omega_p, omega_p < cutx),omega_p >= cutx ], 
                   [0, lambda omega_p:width/np.sqrt(width**2+omega_p**2)*np.abs(1 - np.exp(-duration*width*(1-1j*omega_p/width))),0])
    arg = np.piecewise(omega_p, [omega_p <= - cutx, np.logical_and(-cutx < omega_p, omega_p < cutx),omega_p >= cutx ], 
               [np.pi/2, lambda omega_p:np.angle(-1j*width /(omega_p-1j*width)) + np.angle(1 - np.exp(-duration*width*(1-1j*omega_p/width))) ,-np.pi/2])

    return y * height , arg ,cutx

def qnmcombmodel_cut_1pole(frequency, ** params):
    for key in params.keys():
        if isinstance(params[key], np.ndarray):
            params[key] = params[key][:,None]
    '''the params should be a dictionary containing four keys:
    amplitude,phase,width and spacing
    ''' 

    amplitude = params['amplitude']*1.0
    loc = params['loc']*1.0
    width = params['width']*1.0
    duration = params['duration']*1.0
    deltan = params['deltan']*1.0
    tn = params['tn']*1.0
    mu = qnmcomb_cut_1pole(amplitude,loc,width,duration,frequency)[0]
    arg = qnmcomb_cut_1pole(amplitude,loc,width,duration,frequency)[1]+(deltan-2*np.pi*frequency*tn)
    arg = np.mod((arg+np.pi),(2*np.pi))-np.pi
    cutx = qnmcomb_cut_1pole(amplitude,loc,width,duration,frequency)[2]/(2*np.pi)
    return mu,arg,cutx


def QNM_time_duration_arg_range(tau,duration,spacing,numbers,index,random_factor_1,random_factor_2,delta_n0,t_n0,path):
    dt = 1/4096/8
    time_array = np.arange(0,duration,dt)
    time_domain_strain = np.zeros(len(time_array))
    delta_n_array = []
    t_n_array = []
    for f in np.arange(spacing,numbers * spacing + spacing ,spacing):
        #tau = 10 ## the range of tau is between 1 and 40.2
        #time_domain_strain_part = np.real(np.exp(-time_array/tau)*np.exp(1j*2*np.pi*f*time_array))/tau 
        delta_n = delta_n0
        t_n = t_n0
        if random_factor_1 == 1 :
            delta_n = np.random.random()*2*np.pi
        if random_factor_2 == 1 :    
            t_n = np.random.random()/spacing
        #delta_n = 0.0
        #t_n = 0.01
        delta_n_array.append(delta_n)
        t_n_array.append(t_n)
        time_domain_strain_part = np.exp(-(time_array-t_n)/tau)*np.exp(1j*2*np.pi*f*(time_array-t_n))/tau *np.exp(1j*delta_n)*np.heaviside((time_array-t_n),1)
        ## /tau is a normalized factor
        time_domain_strain = time_domain_strain + time_domain_strain_part
        #print(delta_n,t_n)
    if not os.path.exists(path +'arraydata/'):
            os.mkdir(path +'arraydata/')
    np.save(path +'arraydata/'+ 'delta_n_array_index%d'%(index),delta_n_array)
    np.save(path +'arraydata/'+ 't_n_array_index%d'%(index),t_n_array)
    df = 1/duration
    frequency = np.arange(0,1/dt/2,df)
    echorawfre = dt * np.fft.fft(time_domain_strain)[0:len(frequency)]
    #echorawfre = np.loadtxt('test_width/frequency_domain_strain_tau%d.txt'%tau, dtype=np.complex_)

    loc,hei = scipy.signal.find_peaks(np.abs(echorawfre),height=0.5)
    angle_array = np.angle(echorawfre)
    fre_peaks = frequency[loc]
    fre_peaks_gen = np.arange(spacing,numbers * spacing + spacing ,spacing)
    group_goodness = 1
    if len(loc)!=numbers:
        print("error,len loc = %d"%len(loc),"len peaks = %d"%numbers)
        fre_peaks = fre_peaks_gen
        group_goodness = 0
    for loction in fre_peaks:
        params = dict()
        params['amplitude'] = 1
        params['width'] = 1/tau 
        params['duration'] = duration
        params['loc'] = loction
        params['deltan'] = 0
        params['tn'] = 0.0
        angle = qnmcombmodel_cut_1pole(frequency, ** params)[1]
        cutx = qnmcombmodel_cut_1pole(frequency, ** params)[2]
        cutx3 = spacing/2
        cutx = np.min([cutx,cutx3])
        fmine = loction - cutx
        fmaxe = loction + cutx
        angle_array[int(fmine/df):int(fmaxe/df)+2] = angle_array[int(fmine/df):int(fmaxe/df)+2] - angle[int(fmine/df):int(fmaxe/df)+2]
    angle_array = np.mod((angle_array+np.pi),(2*np.pi))-np.pi
    # this method is to drop the phase around 1 pole

    for peak_index in np.arange(0,numbers):
        fcenter = fre_peaks[peak_index]
        if np.abs(fcenter-fre_peaks_gen[peak_index])<=10*df:
            goodness = 1
        else:
            goodness = 0
        fmaxe = fcenter+cutx
        fmine = fcenter-cutx
        x = frequency[int(fmine/df)+1:int(fmaxe/df)+1]
        y = angle_array[int(fmine/df)+1:int(fmaxe/df)+1]
        y = np.unwrap(y)
        # if np.abs(np.argmax(y)-np.argmin(y))<5:
        #     error = np.min(y)-np.max(y)+2*np.pi
        # else:
        #     error = np.max(y)-np.min(y)\
        error = np.max(y)-np.min(y)
        res = scipy.stats.linregress(x,y)

        rows = [index,tau,duration,spacing,numbers,t_n_array[peak_index],delta_n_array[peak_index],fcenter,res.slope,res.rvalue**2,error,goodness,group_goodness]
        with open(path + '.csv', 'a+') as f:
            f_csv = csv.writer(f)
            f_csv.writerow(rows)

path = "data_version2/data_nofix"
headers = ["group_index","tau","duration","spacing","numbers","t_n_injected","delta_n_injected","f_center","t_n_fitted","Rsquared","residual_arg_range",'peaks_goodness','group_goodness']
with open(path + '.csv', 'w') as f:
    f_csv = csv.writer(f)
    f_csv.writerow(headers)

index_i = 0
duration = 50
numbers = 200
for key_factor in [1.25,1.5,2,3,4,5,7.5,10,15,20,40,60,80,100]:
    for spacing in [5,10,20,30,40]:
        tau = key_factor/spacing
        for test_index in np.arange(0,10):
            if spacing * numbers >2048*8:
                continue
            # if key_factor<=40:
            #     index_i = index_i + 1
            #     continue
            QNM_time_duration_arg_range(tau,duration,spacing,numbers,index_i,1,1,0,0,path)
            index_i = index_i + 1
            print(tau,duration,spacing,numbers,test_index)
