import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import bilby

from bilby.core.prior import PriorDict
import csv
import os

import sys
sys.path.append("..")
from burstevidence_newcomb import burstevidence as burstevidence_old
from burstevidence_newcomb import burstevidence_qnm as burstevidence
from comb_models_version2 import qnmcombmodel_cut as combmodel

start = tem
end = tem 

phase_injected = 1.26
spacing_injected = 5.83
for duration in [1,2,3,4,8,12,24,52]:    
    samplingrate = 4096 
    dt = 1/samplingrate # time resolution in unit of second
    df = 1/duration
    length = int(duration/dt)
    length0 = int(500/dt)

    #     noise = np.random.normal(0, 1,100*length0)
    #     np.save('noise100_500second',noise)
    # These 2 lines aim to generate noise.
    noise0 = np.load('../noise100_500second.npy')


    NFFT = int(0.5 / dt) # the segment length is 0.5s
    fs = int(1/ dt)
    psd_window = np.blackman(NFFT)
    NOVL = int(NFFT/2)

    fmin = 50
    fmax = 339
    # fmax = 339 related to N=49.  While fmax=110 means N=10 (with outdir of "*_n10*")

    inject_params=dict()
    inject_params['width']=1/4

    inject_params['phase'] = phase_injected
    inject_params['spacing'] = spacing_injected
    inject_params['duration'] = duration

    def phase_spacing(parameters):
        
        converted_parameters = parameters.copy()
        converted_parameters['z'] = parameters['spacing'] - parameters['phase']
        return converted_parameters
    priors=PriorDict(conversion_function=phase_spacing)

    priors['width']=bilby.core.prior.LogUniform(1/duration, 2, 'fw')
    # for duration > 4 ,we use priors['width']=bilby.core.prior.LogUniform(1/duration, 0.5, 'fw')
    priors['amplitude']=bilby.core.prior.Uniform(0, 10, 'amplitude')
    priors['phase']=bilby.core.prior.Uniform(0, 10, 'phase')
    priors['spacing']=bilby.core.prior.Uniform(1, 10, 'spacing')
    priors['z'] = bilby.core.prior.Constraint(minimum=0, maximum=10)
    priors['fmin'] = fmin
    priors['fmax'] = fmax
    priors['duration'] = duration

    npoints = 1000
    nact = 10 
    maxmcmc = 10000
    walks = 100

    outdirtemp = '../npoints=1000_duration=%d'%duration+"_4parameter_HWP"
    # or outdirtemp = '../npoints=1000_duration=%d_n10'%duration+"_4parameter_HWP"



    for echoamplitude in [0.077,0]:
        # we tune this echoamplitude to get a suitable SNR
        amplitude_injected = echoamplitude
        inject_params['amplitude'] = amplitude_injected

        outdir = outdirtemp+'/oldlikelihood' # or '/newlikelihood'
        label0 = '4parameter_echoamplitude=%.3f'% echoamplitude+'_oldlikelihood' # or '/newlikelihood'
        whether_print_header = 0
        if(whether_print_header == 0):
            headers = ['i',  'logB', 'maxloglikelihood','SNR_comb_signal','SNR_comb_median', 'SNR_comb_global', 'run_time', 'duration',
                    'width_median', 'width_plus', 'width_minus', 'width_global',
                    'amplitude_median', 'amplitude_plus', 'amplitude_minus', 'amplitude_global', 
                    'phase_median', 'phase_plus', 'phase_minus', 'phase_global',
                    'spacing_median', 'spacing_plus', 'spacing_minus', 'spacing_global',
                    'fmin_median','fmax_median','duration_median']
            if not os.path.exists(outdirtemp):
                os.mkdir(outdirtemp)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            with open(outdir+'/../'+label0+'_all.csv', 'w') as f:
                f_csv = csv.writer(f)
                f_csv.writerow(headers)
            whether_print_header = 1
                
                
        for i in np.arange(start,end,1):

            noisep = noise0[i*length0:(i+1)*length0]
            noise = noisep[0:length]
            label = label0+'i'+str(i)
            Pxx_strain, freqs = mlab.psd(noisep, Fs = fs, NFFT = NFFT,window = psd_window,noverlap = NOVL)
            psd_strain = interp1d(freqs, Pxx_strain)
            
            frequency = np.arange(0,int(freqs.max()),df)

            
            noisefre = dt * np.fft.fft(noise)[0:len(frequency)]
            strainfre = noisefre + combmodel(frequency,**inject_params)[0] * np.exp(1j*combmodel(frequency,**inject_params)[1])
            likelihood=burstevidence_old(x=frequency,y=np.abs(strainfre),
                                        sn=psd_strain(frequency),function = combmodel , df = df )
            # for newlikelihood, we use these code:
            # likelihood=burstevidence(x=frequency,y=np.abs(strainfre),angle = np.angle(strainfre),
            #                             sn=psd_strain(frequency),function = combmodel , df = df )


            result = bilby.run_sampler(
                likelihood=likelihood, priors=priors, sampler='dynesty', npoints=npoints, nact=nact, maxmcmc=maxmcmc, walks=walks, outdir=outdir, label=label,injection_parameters=inject_params)            
            result.plot_corner(filename=outdir+'/'+label+'_corner.png')
            result.save_posterior_samples(filename=outdir+'/'+label+'_samples.dat')

        ###################################################### csv data #####################################################
            logb = result.log_evidence
            maxloglikelihood = result.posterior.max()['log_likelihood']
            #The suffix 'm' means 'median'; the suffix 'g' means 'global'; the suffix 'p' means 'plus'; the suffix 'n' means 'minus'
            #'global' related to the parameter which gives the maximum likelihood
            #'median','plus','minus' related to the median parameter from Bayes posterior
            posterior_params_m=dict()
            for key in result.search_parameter_keys+result.fixed_parameter_keys:
                posterior_params_m[key]=result.get_one_dimensional_median_and_error_bar(key).median
            posterior_params_p=dict()
            for key in result.search_parameter_keys:
                posterior_params_p[key]=result.get_one_dimensional_median_and_error_bar(key).plus
            posterior_params_n=dict()
            for key in result.search_parameter_keys:            
                posterior_params_n[key]=result.get_one_dimensional_median_and_error_bar(key).minus
            posterior_params_g=dict()
            posterior_params_g = result.posterior.loc[[result.posterior.idxmax()['log_likelihood']]].to_dict(orient ='records')[0]

        ############################################################  csv data end  ##############################################
            frequencye= frequency[int(fmin/df):int(fmax/df+2)]

            hi = combmodel(frequencye, ** inject_params)[0]*np.exp(1j*combmodel(frequencye, ** inject_params)[1])
            rho_opt_signal = np.sqrt(4*df * np.sum(np.abs(hi)**2/psd_strain(frequencye)))

            hi_m = combmodel(frequencye, ** posterior_params_m)[0]*np.exp(1j*combmodel(frequencye, ** posterior_params_m)[1])
            ρoptcomb_m = np.sqrt(4*df * np.sum(np.abs(hi_m)**2/psd_strain(frequencye)))
            hi_g = combmodel(frequencye, ** posterior_params_g)[0]*np.exp(1j*combmodel(frequencye, ** posterior_params_g)[1])
            ρoptcomb_g = np.sqrt(4*df * np.sum(np.abs(hi_g)**2/psd_strain(frequencye)))
            
                                
            strainfreN = strainfre/psd_strain(frequency)
            noisefreN = noisefre/psd_strain(frequency)   
            
            axins_xmin1 = 106
            axins_xmax1 = 106.5
            axins_ymax1 = 4


            fig = plt.figure(figsize=(15,10), constrained_layout=True)
            ## Here we add constrained_layout=True to avoid the overlapping of 2 subplots
            ax1=fig.add_subplot(2,1,1)
            ax1.set_title('$SNR_{signal}$='+"%.6g" % rho_opt_signal+'   $SNR_{model}$=' +
                    "%.6g" % ρoptcomb_m+'   $\ln B=$'+"%.6g" % result.log_evidence,fontsize=15)
            ax1.plot(frequency, np.abs(strainfreN)*np.sqrt(psd_strain(frequency)*(4*df)),  c='orange', label='normalized abs(dH/sH+dL/sL)')
            ax1.plot(frequency, np.abs(noisefreN)*np.sqrt(psd_strain(frequency)*(4*df)), c='aquamarine',label='normalized noise data')
            ax1.plot(frequencye,np.abs(hi)*np.sqrt((4*df)/psd_strain(frequencye)) , c='red',alpha=0.8, label='normalized injected signal')
            ax1.plot(frequencye,np.abs(hi_m)*np.sqrt((4*df)/psd_strain(frequencye)), c='blue',alpha=0.8, label='normalized best-fit comb')
            ax1.legend(loc='upper right')
            ax1.axis([fmin, fmax, 0, 7])
            ax1.set_ylabel(''),ax1.set_xlabel('f (Hz)',fontsize=15)
            ax1.tick_params(labelsize =15)

            axins1 = ax1.inset_axes((0.05, 0.7, 0.2, 0.25))
            axins1.plot(frequencye,np.abs(hi)*np.sqrt((4*df)/psd_strain(frequencye)) , c='red', alpha=0.8, label='normalized injected signal')        
            axins1.plot(frequencye,np.abs(hi_m)*np.sqrt((4*df)/psd_strain(frequencye)),c='blue',alpha=0.8, label='normalized best-fit comb')
            axins1.set_xticks(np.linspace(axins_xmin1,axins_xmax1,6))
            axins1.tick_params(labelsize =13)
            axins1.set_xticklabels(['%.1f'%x for x in np.linspace(axins_xmin1,axins_xmax1,6)],fontsize=12)
            axins1.axis([axins_xmin1, axins_xmax1, 0, axins_ymax1])

            ax2=fig.add_subplot(2,1,2)
            ax2.plot(frequency, np.abs(strainfreN)*np.sqrt(psd_strain(frequency)*(4*df)), c='orange', label='normalized abs(dH/sH+dL/sL)')
            ax2.plot(frequency, np.abs(noisefreN)*np.sqrt(psd_strain(frequency)*(4*df)), c='aquamarine', alpha=0.98, label='normalized noise data')
            ax2.plot(frequencye,np.abs(hi)*np.sqrt((4*df)/psd_strain(frequencye)) , c='red', alpha=0.8,label='normalized injected signal')
            ax2.plot(frequencye,np.abs(hi_g)*np.sqrt((4*df)/psd_strain(frequencye)), c='blue',alpha=0.8, label='normalized best-fit comb')
            ax2.set_ylabel(''),ax2.set_xlabel('f (Hz)',fontsize=15)
            ax2.legend(loc='upper right')
            ax2.axis([fmin, fmax, 0, 7])
            ax2.tick_params(labelsize =15)
            ax2.set_title('Global maximum, comb SNR = {0:.2f}'.format(ρoptcomb_g),fontsize=15)

            axins2 = ax2.inset_axes((0.05, 0.7, 0.2, 0.25))
            axins2.plot(frequencye,np.abs(hi)*np.sqrt((4*df)/psd_strain(frequencye)) , c='red', alpha=0.8, label='normalized injected signal')
            axins2.plot(frequencye,np.abs(hi_g)*np.sqrt((4*df)/psd_strain(frequencye)),c='blue',alpha=0.8, label='normalized best-fit comb')
            axins2.set_xticks(np.linspace(axins_xmin1,axins_xmax1,6))
            axins2.tick_params(labelsize =13)
            axins2.set_xticklabels(['%.1f'%x for x in np.linspace(axins_xmin1,axins_xmax1,6)],fontsize=12)
            axins2.axis([axins_xmin1, axins_xmax1, 0, axins_ymax1])

            plt.savefig(outdir+'/'+label+''+'_bestfit_comb.png')
            plt.close()
            
            rows = [i, logb, maxloglikelihood,rho_opt_signal,ρoptcomb_m,ρoptcomb_g, result.sampling_time,duration]
                
            for key in result.search_parameter_keys:
                rows = np.append(rows,[posterior_params_m[key],posterior_params_p[key],posterior_params_n[key],posterior_params_g[key]])
            for key in result.fixed_parameter_keys:
                rows = np.append(rows,[posterior_params_m[key]])
            with open(outdir+'/../'+label0+'_all.csv', 'a+') as f:
                f_csv = csv.writer(f)
                f_csv.writerow(rows)
