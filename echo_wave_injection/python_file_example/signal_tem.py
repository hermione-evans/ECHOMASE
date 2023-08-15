import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import bilby
from glob import glob
from bilby.core.prior import PriorDict
import csv
import os
import sys
sys.path.append("..")
from burstevidence_newcomb import burstevidence as burstevidence_old
from burstevidence_newcomb import burstevidence_qnm as burstevidence
from comb_models_version3 import qnmcombmodel_cut as combmodel


chi = 0.69
f_RD = (1.5251-1.1568*(1-chi)**0.129)/(2*np.pi)
R_bar = 0.00547/(1+1/np.sqrt(1-chi**2))
fmin = 0
fmax = f_RD
# These 5 lines aims to set parameters about echo model.

# tem part
start = tem
end = tem
likelihood_index=tem
index = tem
# start&end tem is the number of noise data we use. In the paper, we use 100 noise data for each duration.
# And we use multiprocessing to run 100 noise data at the same time. (See the file: generate_python_files.py)
# likelihood_index is the index of likelihood we use. 0 means the new likelihood, 1 means the old likelihood.
# index is the index of injected echo model file. We have 30 injected echo model files. (See the folder: echo_waves)

noise0 = np.load('../noise100_500second.npy')
length0 = int(500*4096) # 100*length0 is the length of noise data we generate.
file_list = glob('../echo_waves/*.dat')
file_list.sort()
file_list = file_list+file_list[18:24] # Add noise data to the end of the file_list. We set inject amplitude to 0 to get noise data.
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
likelihood_string = likelihoodlist[likelihood_index]
file = file_list[index]
amplitude_inject = amplitude_array[index]
necho = echo_number_list[index] # The number of echo waves in the injected signal. It is contained in the file name.

if necho >=100:
    npoints = 2000
else:
    npoints = 1000
nact = 10 
maxmcmc = 10000
walks = 100
# These 7 lines aims to set sampling parameters.

if amplitude_inject !=0:
    name = name_list[index]
else:
    name = name_list[index]
    name = 'noise'+ name[name.find('N'):name.find('N')+4]
outdirtemp = "../"+name[0:name.find('N')]
outdir = outdirtemp+'/'+likelihood_string
label0 = name[name.find('N'):name.find('N')+4]+'_npoints=%d_6parameter_echoamplitude=%d'% (npoints,amplitude_inject)+'_'+likelihood_string
whether_print_header = 0
if(whether_print_header == 0):
    headers = ['i',  'logB', 'maxloglikelihood','SNR_inj_signal', 'SNR_comb_global', 'run_time', 'duration',
            'width_median', 'width_plus', 'width_minus', 'width_global',
            'amplitude_median', 'amplitude_plus', 'amplitude_minus', 'amplitude_global', 
            'phase_median', 'phase_plus', 'phase_minus', 'phase_global',
            'spacing_median', 'spacing_plus', 'spacing_minus', 'spacing_global',
            'fmin_median', 'fmin_plus', 'fmin_minus', 'fmin_global',
            'fmax_median', 'fmax_plus', 'fmax_minus', 'fmax_global',
            'duration_median']
    if not os.path.exists(outdirtemp):
        os.mkdir(outdirtemp)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    with open(outdir+'/../'+label0+'_all.csv', 'w') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(headers)
    whether_print_header = 1
# These lines are used to generate the header of the csv file and the folder to save the data.

        
inject_data = np.genfromtxt(file,delimiter="\t",dtype=np.complex128)
omega = np.real(inject_data[0::,0])
frequency = omega/(2*np.pi)
df = np.diff(frequency)[0]
duration = 1/df
dt = duration/len(omega)
frequency = frequency[0:len(frequency)//2]
# These lines are used to read the injected echo model file and get the frequency array.
# The model files contains both the postive and negative frequencies. We only deal with the positive frequencies.

inject_strain_fre = inject_data[0::,1]
inject_strain_fre = (inject_strain_fre + np.conjugate(inject_strain_fre[::-1]))/2
inject_strain_fre = inject_strain_fre[0:len(frequency)]
# Here  we assume the detector response for positive and negative frequencies are the same. See our paper Sec 4 for more details.

length = len(inject_strain_fre)*2


for i in np.arange(start,end):
    label = label0+'i'+str(i)

    noisep = noise0[i*length0:(i+1)*length0]
    noise = noisep[0:length]

    NFFT = 2048
    fs = 1/dt
    psd_window = np.blackman(NFFT)
    NOVL = int(NFFT/2)
    # These lines are used to set the psd parameter of the noise data.
    
    Pxx_strain, freqs = mlab.psd(noisep, Fs = fs, NFFT = NFFT,window = psd_window,noverlap = NOVL)
    psd_strain = interp1d(freqs, Pxx_strain)
    noisefre = dt * np.fft.fft(noise)[0:len(frequency)]
    strainfre = inject_strain_fre * amplitude_inject + noisefre
    # These lines are used to generate the PSD of the noise data and inject signal data.
    
    def frange(parameters):
        converted_parameters = parameters.copy()
        converted_parameters['z'] = parameters['fmax'] - parameters['fmin'] - parameters['spacing']*10
        return converted_parameters
    priors=PriorDict(conversion_function=frange)
    priors['width']=bilby.core.prior.LogUniform(1/duration, R_bar, 'fw')
    priors['amplitude']=bilby.core.prior.Uniform(np.mean(np.abs(noisefre))/100, np.mean(np.abs(noisefre))*10, 'amplitude')
    # In our paper, we use the notation \(\langle \tilde{P} \rangle\) to represent `np.mean(psd_strain(frequency))`, although in the code we actually employ `np.mean(np.abs(noisefre))`. 
    # These two expressions are numerically similar, with `np.mean(np.abs(noisefre))` being approximately 1.2 times the value of `np.mean(psd_strain(frequency))`."

    priors['phase']=bilby.core.prior.Uniform(0, 1, 'phase')
    priors['spacing']=bilby.core.prior.Uniform(R_bar/4, R_bar, 'spacing')
    priors['z'] = bilby.core.prior.Constraint(minimum=0, maximum=f_RD)
    priors['fmin'] = bilby.core.prior.Uniform(0, f_RD, 'fmin')
    priors['fmax'] = bilby.core.prior.Uniform(0, f_RD, 'fmax')
    priors['duration'] = duration
    # These lines are the priors of the parameters. See Table 2 in our paper for more details.
    
    # # dry run test. Just to test the correctness of the code
    # npoints = 10
    # nact = 1
    # maxmcmc = 100
    # walks = 1
    
    if likelihood_index == 0:
        likelihood=burstevidence(x=frequency,y=np.abs(strainfre),angle = -1*np.angle(strainfre),sn=psd_strain(frequency),function = combmodel , df = df )#remember to set df
        # The minus symbal comes from the different conservation of the phase between the Mathematica and the python
        # the echo_waves/*.dat files comes from Mathemtica. The amplitude also have a different conservation, but we have already tuned the amplitude in the amplitude_array to get a suitable SNR_{dom}
    if likelihood_index == 1:
        likelihood=burstevidence_old(x=frequency,y=np.abs(strainfre),sn=psd_strain(frequency),function = combmodel , df = df )#remember to set df
    if os.path.exists(outdir+'/'+label+'_result.json' ):
        result = bilby.result.read_in_result(outdir=outdir, label=label)
    # If the result file exists, we read the result file. Otherwise, we run the sampler.
    else:
        result = bilby.run_sampler(likelihood=likelihood, priors=priors, sampler='dynesty', npoints=npoints, nact=nact, maxmcmc=maxmcmc, walks=walks, outdir=outdir, label=label)            
    # result.plot_corner(filename=outdir+'/'+label+'_corner.png')
    # result.save_posterior_samples(filename=outdir+'/'+label+'_samples.dat')

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


    fmin_m=posterior_params_m['fmin']
    fmax_m=posterior_params_m['fmax']
    fmin_g=posterior_params_g['fmin']
    fmax_g=posterior_params_g['fmax']

    frequencye_m= frequency[int(fmin_m/df):int(fmax_m/df+2)]
    frequencye_g= frequency[int(fmin_g/df):int(fmax_g/df+2)]

    hi_g0 = np.abs(inject_strain_fre * amplitude_inject)[int(fmin_g/df):int(fmax_g/df+2)]#echo signal           
    rho_opt_signal_g = np.sqrt(4*df * np.sum(hi_g0**2/psd_strain(frequencye_g)))
    hi_g = combmodel(frequencye_g, ** posterior_params_g)[0]*np.exp(1j*combmodel(frequencye_g, ** posterior_params_g)[1])
    ρoptcomb_g = np.sqrt(4*df * np.sum(np.abs(hi_g)**2/psd_strain(frequencye_g)))
    # These lines are used to calculate the SNR of the injected signal and the SNR of the best-fit comb.
#######################################################  csv data end  ##############################################                        
    strainfreN = strainfre/psd_strain(frequency)
    noisefreN = noisefre/psd_strain(frequency)   
    
    # We do not draw these figs in slurm run on supercomputer
    # This part is reproduced in SNR_result_generate.py
    # axins_xmin1 = 106
    # axins_xmax1 = 106.5
    # axins_ymax1 = 4

    # axins_xmin2 = 232
    # axins_xmax2 = 233
    # axins_ymax2 = 2


    # fig = plt.figure(figsize=(15,10), constrained_layout=True)
    # ## Here we add constrained_layout=True to avoid the overlapping of 2 subplots
    # ax1=fig.add_subplot(2,1,1)
    # ax1.set_title('$SNR_{signal}$='+"%.6g" % rho_opt_signal_g+'   $SNR_{model}$=' +
    #         "%.6g" % ρoptcomb_g+'   $\ln B=$'+"%.6g" % result.log_evidence,fontsize=15)
    # ax1.plot(frequency, np.abs(strainfreN)*np.sqrt(psd_strain(frequency)*(4*df)),  c='orange', label='normalized abs(dH/sH+dL/sL)')
    # ax1.plot(frequency, np.abs(noisefreN)*np.sqrt(psd_strain(frequency)*(4*df)), c='aquamarine',label='normalized noise data')
    # ax1.plot(frequencye_g,np.abs(hi_g0)*np.sqrt((4*df)/psd_strain(frequencye_g)) , c='red',alpha=0.8, label='normalized injected signal')
    # ax1.plot(frequencye_g,np.abs(hi_g)*np.sqrt((4*df)/psd_strain(frequencye_g)), c='blue',alpha=0.8, label='normalized best-fit comb')
    # ax1.legend(loc='upper right')
    # ax1.axis([fmin_g, fmax_g, 0, 7])
    # ax1.set_ylabel(''),ax1.set_xlabel('Mf ',fontsize=15)
    # ax1.tick_params(labelsize =15)
    # ## Here we change all "plt" into "ax1" for the convenience of add inplot 

    # # axins1 = ax1.inset_axes((0.05, 0.7, 0.2, 0.25))
    # # #axins1.plot(frequency, np.abs(strainfreN)*np.sqrt(psd_strainLH_list*(4*df)), c='orange', label='normalized abs(dH/sH+dL/sL)')
    # # #axins1.plot(frequency, np.abs(noisefreN)*np.sqrt(psd_strainLH_list*(4*df)), c='aquamarine',alpha=0.98,label='normalized noise data')
    # # #axins1.plot(frequencye,hi*np.sqrt((4*df)/psd_strainH(frequencye)+(4*df)/psd_strainL(frequencye)) , c='red', alpha=0.8, label='normalized injected signal')
    # # axins1.plot(frequencye,hi_m*np.sqrt((4*df)/psd_strainH(frequencye)+(4*df)/psd_strainL(frequencye)),c='blue',alpha=0.8, label='normalized best-fit comb')
    # # axins1.set_xticks(np.linspace(axins_xmin1,axins_xmax1,6))
    # # axins1.tick_params(labelsize =13)
    # # axins1.set_xticklabels(['%.1f'%x for x in np.linspace(axins_xmin1,axins_xmax1,6)],fontsize=12)
    # # axins1.axis([axins_xmin1, axins_xmax1, 0, axins_ymax1])

    # # axins2 = ax1.inset_axes((0.45, 0.7, 0.2, 0.25))
    # # #axins2.plot(frequency, np.abs(strainfreN)*np.sqrt(psd_strainLH_list*(4*df)), c='orange', label='normalized abs(dH/sH+dL/sL)')
    # # #axins2.plot(frequency, np.abs(noisefreN)*np.sqrt(psd_strainLH_list*(4*df)), c='aquamarine',alpha=0.98,label='normalized noise data')
    # # axins2.plot(frequencye,hi*np.sqrt((4*df)/psd_strainH(frequencye)+(4*df)/psd_strainL(frequencye)) , c='red', alpha=0.8, label='normalized injected signal')
    # # axins2.plot(frequencye,hi_m*np.sqrt((4*df)/psd_strainH(frequencye)+(4*df)/psd_strainL(frequencye)),c='blue',alpha=0.8, label='normalized best-fit comb')
    # # axins2.set_xticks(np.linspace(axins_xmin2,axins_xmax2,6))
    # # axins2.tick_params(labelsize =13)
    # # axins2.set_xticklabels(['%.1f'%x for x in np.linspace(axins_xmin2,axins_xmax2,6)],fontsize=12)
    # # axins2.axis([axins_xmin2, axins_xmax2, 0, axins_ymax2])
    # # ## All the way to set settings in axins are the same as ax1

    # ax2=fig.add_subplot(2,1,2)
    # ax2.plot(frequency, np.abs(strainfreN)*np.sqrt(psd_strain(frequency)*(4*df)), c='orange', label='normalized abs(dH/sH+dL/sL)')
    # ax2.plot(frequency, np.abs(noisefreN)*np.sqrt(psd_strain(frequency)*(4*df)), c='aquamarine', alpha=0.98, label='normalized noise data')
    # ax2.plot(frequencye_g,np.abs(hi_g0)*np.sqrt((4*df)/psd_strain(frequencye_g)) , c='red', alpha=0.8,label='normalized injected signal')
    # ax2.plot(frequencye_g,np.abs(hi_g)*np.sqrt((4*df)/psd_strain(frequencye_g)), c='blue',alpha=0.8, label='normalized best-fit comb')
    # ax2.set_ylabel(''),ax2.set_xlabel('Mf',fontsize=15)
    # ax2.legend(loc='upper right')
    # ax2.axis([fmin_g, fmax_g, 0, 7])
    # ax2.tick_params(labelsize =15)
    # ax2.set_title('Global maximum, comb SNR = {0:.2f}'.format(ρoptcomb_g),fontsize=15)

    # # axins1 = ax2.inset_axes((0.05, 0.7, 0.2, 0.25))
    # # #axins1.plot(frequency, np.abs(strainfreN)*np.sqrt(psd_strainLH_list*(4*df)), c='orange', label='normalized abs(dH/sH+dL/sL)')
    # # #axins1.plot(frequency, np.abs(noisefreN)*np.sqrt(psd_strainLH_list*(4*df)), c='aquamarine',alpha=0.98,label='normalized noise data')
    # # axins1.plot(frequencye,hi*np.sqrt((4*df)/psd_strainH(frequencye)+(4*df)/psd_strainL(frequencye)) , c='red', alpha=0.8, label='normalized injected signal')
    # # axins1.plot(frequencye,hi_g*np.sqrt((4*df)/psd_strainH(frequencye)+(4*df)/psd_strainL(frequencye)),c='blue',alpha=0.8, label='normalized best-fit comb')
    # # axins1.set_xticks(np.linspace(axins_xmin1,axins_xmax1,6))
    # # axins1.tick_params(labelsize =13)
    # # axins1.set_xticklabels(['%.1f'%x for x in np.linspace(axins_xmin1,axins_xmax1,6)],fontsize=12)
    # # axins1.axis([axins_xmin1, axins_xmax1, 0, axins_ymax1])

    # # axins2 = ax2.inset_axes((0.45, 0.7, 0.2, 0.25))
    # # #axins2.plot(frequency, np.abs(strainfreN)*np.sqrt(psd_strainLH_list*(4*df)), c='orange', label='normalized abs(dH/sH+dL/sL)')
    # # #axins2.plot(frequency, np.abs(noisefreN)*np.sqrt(psd_strainLH_list*(4*df)), c='aquamarine',alpha=0.98,label='normalized noise data')
    # # axins2.plot(frequencye,hi*np.sqrt((4*df)/psd_strainH(frequencye)+(4*df)/psd_strainL(frequencye)) , c='red', alpha=0.8, label='normalized injected signal')
    # # axins2.plot(frequencye,hi_g*np.sqrt((4*df)/psd_strainH(frequencye)+(4*df)/psd_strainL(frequencye)),c='blue',alpha=0.8, label='normalized best-fit comb')
    # # axins2.set_xticks(np.linspace(axins_xmin2,axins_xmax2,6))
    # # axins2.tick_params(labelsize =13)
    # # axins2.set_xticklabels(['%.1f'%x for x in np.linspace(axins_xmin2,axins_xmax2,6)],fontsize=12)
    # # axins2.axis([axins_xmin2, axins_xmax2, 0, axins_ymax2])
    # plt.savefig(outdir+'/'+label+''+'_bestfit_comb.png')
    # plt.close()
    

    rows = [i, logb, maxloglikelihood,rho_opt_signal_g,ρoptcomb_g, result.sampling_time,duration]
        
    for key in result.search_parameter_keys:
        rows = np.append(rows,[posterior_params_m[key],posterior_params_p[key],posterior_params_n[key],posterior_params_g[key]])
    for key in result.fixed_parameter_keys:
        rows = np.append(rows,[posterior_params_m[key]])
    with open(outdir+'/../'+label0+'_all.csv', 'a+') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(rows)
