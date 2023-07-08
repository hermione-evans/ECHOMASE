 # import packages
from tempfile import tempdir
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
sys.path.append("..")
import os
import scipy.stats
from scipy.interpolate import interp1d
import bilby
from burstevidence_newcomb import burstevidence as burstevidence_old
from burstevidence_newcomb import burstevidence_qnm as burstevidence
from comb_models_version2 import qnmcombmodel_cut as combmodel
from bilby.core.prior import PriorDict
import csv
os.chdir('..')
whether_print_header = 0

begin = tem
end = begin + 2

if(whether_print_header == 0):
    headers = [ 'echoamplitude_injected','echoamplitude_searched','i','duration', 'fmin','fmax',
           'log_likelihood_new','log_likelihood_old','SNR',
           'normalized_echoamplitude_injected','normalized_echoamplitude_searched']
    with open('dataMoreAmplitude/data_all_new.csv', 'w') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(headers)
    whether_print_header = 1

dt = 1/4096
length0 = int(500/dt)

fmin = 110
fmax = 140

duration_array = [4,10,50,100,200,500]

for duration in duration_array:
    df = 1 / duration
    length = int(duration / dt)

    NFFT = int(0.5 / dt) # the segment length is 0.5s
    fs = int(1 / dt)
    psd_window = np.blackman(NFFT)
    NOVL = int(NFFT / 2)

    noiseL0 = np.load('noise100_500second.npy')
    for i in np.arange(begin,end,1):              

        
        param = dict()
        param['width']=0.25
        
        param['phase']= 1.26
        param['spacing']=5
        param['fmin']=fmin
        param['fmax']=fmax
        param['duration']=duration
# above lines are just for parameter setting

        noiseLp = noiseL0[i*length0:(i+1)*length0]
        noiseL = noiseLp[0:length]
        Pxx_strainL, freqs = mlab.psd(noiseLp, Fs = fs, NFFT = NFFT,window = psd_window,noverlap = NOVL)
        psd_strainL = interp1d(freqs, Pxx_strainL)
        frequency = np.arange(0,int(freqs.max()),df)
        noisefreL = dt * np.fft.fft(noiseL)[0:len(frequency)]

        for echoamplitude in [0,0.01,0.03,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
            param['amplitude']= echoamplitude
            
            strainfreL = noisefreL + combmodel(frequency,**param)[0] * np.exp(1j*combmodel(frequency,**param)[1])

            likelihood = burstevidence(x=frequency,y=np.abs(strainfreL),angle = np.angle(strainfreL),sn=psd_strainL(frequency),function = combmodel , df = df )#remember to set df
            likelihood_old = burstevidence_old(x=frequency,y=np.abs(strainfreL),sn=psd_strainL(frequency),function = combmodel , df = df )#remember to set df
            if echoamplitude<0.15:
                searched_array = np.linspace(0, 0.3 , 60)
            else:
                searched_array = np.linspace(0.6*echoamplitude, 1.4*echoamplitude , 60)
            for echoamplitude_searched in searched_array:
                for x in param:
                    likelihood.parameters[x] = param[x]
                    likelihood_old.parameters[x] = param[x]
                likelihood.parameters['amplitude'] = echoamplitude_searched
                likelihood_old.parameters['amplitude'] = echoamplitude_searched
    # this four lines aim to pass parameters into likelihood


                frequencye= frequency[int(fmin/df):int(fmax/df+2)]
                hi = combmodel(frequencye, ** param)[0]
    # hi = combmodel(frequency,**param)[0]*(1-np.exp(-duration/tau)) for time-domain amplitude
                ρoptcomb = np.sqrt(4*df * np.sum(hi**2/psd_strainL(frequencye)))
                normalized_echoamplitude = echoamplitude * np.mean(np.sqrt((4*df)/psd_strainL(frequencye)))
                normalized_echoamplitude_searched = echoamplitude_searched * np.mean(np.sqrt((4*df)/psd_strainL(frequencye)))

                rows = [echoamplitude,echoamplitude_searched,i,duration,fmin,fmax,
                        likelihood.log_likelihood(),likelihood_old.log_likelihood(),ρoptcomb,
                        normalized_echoamplitude,normalized_echoamplitude_searched]

                with open('dataMoreAmplitude/data_all_new.csv', 'a+') as f:
                    f_csv = csv.writer(f)
                    f_csv.writerow(rows)
                print(rows, end="\r")