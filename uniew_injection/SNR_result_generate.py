import pandas as pd
import matplotlib.mlab as mlab
import scipy
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import bilby
import corner
from glob import glob
from scipy.interpolate import interp1d
from comb_models_version2 import qnmcombmodel_cut as combmodel

for duration in [1,2,3,4,8,12,24,52]:
	samplingrate = 4096 
	dt = 1/samplingrate # time resolution in unit of second
	df = 1/duration
	length = int(duration/dt)
	length0 = int(500/dt)

	noise0 = np.load('noise100_500second.npy')

	NFFT = int(0.5 / dt) # the segment length is 0.5s
	fs = int(1/ dt)
	psd_window = np.blackman(NFFT)
	NOVL = int(NFFT/2)

	for i in range(0,100):
		noisep = noise0[i*length0:(i+1)*length0]
		noise = noisep[0:length]
		Pxx_strain, freqs = mlab.psd(noisep, Fs = fs, NFFT = NFFT,window = psd_window,noverlap = NOVL)
		psd_strain = interp1d(freqs, Pxx_strain)
		frequency = np.arange(0,int(freqs.max()),df)
		noisefre = dt * np.fft.fft(noise)[0:len(frequency)]

		#result_name=glob('*duration=%d_*HWP/newlikelihood/*0.077*i%d_result.json'%(duration,i))
		result_name=glob('*duration=%d_*HWP/oldlikelihood/*0.077*i%d_result.json'%(duration,i))
		result = bilby.result.read_in_result(result_name[0])
		posterior = result.posterior
		posterior_SNR = []
		for j in np.arange(len(posterior)):
			inject_params = posterior.iloc[j]
			fmin = inject_params['fmin']
			fmax = inject_params['fmax']
			frequencye = frequency[int(fmin/df):int(fmax/df+2)]
			hi = combmodel(frequencye,**inject_params)[0]
			posterior_SNR.append(np.sqrt(np.sum(hi**2/psd_strain(frequencye))))
		posterior_SNR = np.array(posterior_SNR)	
		posterior['SNR'] = np.sqrt(4*df)*posterior_SNR
		result.posterior = posterior	
		result.label = result.label +'_SNR'
		result.save_to_file()
		#print('new',duration,i)
		print('old',duration,i)