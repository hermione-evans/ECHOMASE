import bilby
import numpy as np
from scipy.special import ive as Ive
from scipy.special import iv as Iv
import inspect
def logi0(x):
    return np.log(Ive(0,x))+np.abs(x)


class burstevidence(bilby.Likelihood):
    """
    Likelihood for one detector with no detector response
    """

    def __init__(self, x, y, sn, function,df):
        """

        Parameters
        ----------
        x: array-like
            The frequency array
        y: array-like
            The strain data of frequency domain 
        sn: array-like
            noise spectral density. this symbol is coherent with that in (5.2),arXiv:1811.02042v1
        phase: float
            The parametre of comb function, it should be between 0 to 1
        width: float
            describes the width of signal
        spacing: float
            describes the whole length of a wave
        """
        super().__init__(parameters={'fmin': None,'fmax': None})
        
        self.df = df 
        self.x = x 
        self.y = y
        self.sn = sn 
        ## maybe later sn will be replaced by an interpolating function
        self.function = function  
        self.N = len(x)
        
        
        
        parameters = inspect.getargspec(function).args
        parameters.pop(0)
        self.parameters = dict.fromkeys(parameters)
        self.function_keys = self.parameters.keys()
        
    def log_likelihood(self):
        fmin = self.parameters['fmin']
        fmax = self.parameters['fmax']
        xe = self.x[int(fmin/self.df):int(fmax/self.df+2)]
        sne = self.sn[int(fmin/self.df):int(fmax/self.df+2)]
        ye = self.y[int(fmin/self.df):int(fmax/self.df+2)]
        model_parameters = {k: self.parameters[k] for k in self.function_keys}
        hi = self.function(xe,  **model_parameters)[0]
        burst = np.sum(logi0(4*self.df*np.abs(ye)*np.abs(hi)/sne)-2*self.df*np.abs(hi)**2/sne)
        return burst

class burstevidence_old(bilby.Likelihood):
    """
    Likelihood for one detector with no detector response
    """

    def __init__(self, x, y, sn, function,df):
        """

        Parameters
        ----------
        x: array-like
            The frequency array
        y: array-like
            The strain data of frequency domain 
        sn: array-like
            noise spectral density. this symbol is coherent with that in (5.2),arXiv:1811.02042v1
        phase: float
            The parametre of comb function, it should be between 0 to 1
        width: float
            describes the width of signal
        spacing: float
            describes the whole length of a wave
        """
        super().__init__(parameters={'fmin': None,'fmax': None})
        
        self.df = df 
        self.x = x 
        self.y = y
        self.sn = sn 
        ## maybe later sn will be replaced by an interpolating function
        self.function = function  
        self.N = len(x)
        
        
        
        parameters = inspect.getargspec(function).args
        parameters.pop(0)
        self.parameters = dict.fromkeys(parameters)
        self.function_keys = self.parameters.keys()
        
    def log_likelihood(self):
        fmin = self.parameters['fmin']
        fmax = self.parameters['fmax']
        xe = self.x[int(fmin/self.df):int(fmax/self.df+2)]
        sne = self.sn[int(fmin/self.df):int(fmax/self.df+2)]
        ye = self.y[int(fmin/self.df):int(fmax/self.df+2)]
        model_parameters = {k: self.parameters[k] for k in self.function_keys}
        hi = self.function(xe,  **model_parameters)
        burst = np.sum(logi0(4*self.df*np.abs(ye)*np.abs(hi)/sne)-2*self.df*np.abs(hi)**2/sne)
        return burst



class burstevidence_qnm(bilby.Likelihood):
    """
    Likelihood for one detector with no detector response, while keeping the phase associated with the pole structure
    """

    def __init__(self, x, y, angle, sn, function,df):
        """

        Parameters
        ----------
        x: array-like
            The frequency array
        y: array-like
            The abs of strain data of frequency domain 
        angle: array-like
            The angle of strain data of frequency domain 
        sn: array-like
            noise spectral density. this symbol is coherent with that in (5.2),arXiv:1811.02042v1
        phase: float
            The parametre of comb function, it should be between 0 to 1
        width: float
            describes the width of signal
        spacing: float
            describes the whole length of a wave
        """
        super().__init__(parameters={'fmin': None,'fmax': None})

        self.df = df 
        self.x = x 
        self.y = y
        self.angle = angle
        self.sn = sn 
        ## maybe later sn will be replaced by an interpolating function
        self.function = function  
        self.N = len(x)



        parameters = inspect.getfullargspec(function).args
        parameters.pop(0)
        self.parameters = dict.fromkeys(parameters)
        self.function_keys = self.parameters.keys()

    def log_likelihood(self):
        fmin = self.parameters['fmin']
        fmax = self.parameters['fmax']
        xe = self.x[int(fmin/self.df):int(fmax/self.df+2)]
        sne = self.sn[int(fmin/self.df):int(fmax/self.df+2)]
        ye = self.y[int(fmin/self.df):int(fmax/self.df+2)]
        phi = self.angle[int(fmin/self.df):int(fmax/self.df+2)]

        model_parameters = {k: self.parameters[k] for k in self.function_keys}
        hj = self.function(xe,  **model_parameters)[0]
        arg_h = self.function(xe,  **model_parameters)[1]
        qnm_number = self.function(xe,  **model_parameters)[2]
        arg = phi - arg_h

        arg_part = np.split(arg,np.unique(qnm_number, return_index=True)[1][1:])
        hj_part = np.split(hj,np.unique(qnm_number, return_index=True)[1][1:])
        y_part = np.split(ye,np.unique(qnm_number, return_index=True)[1][1:])
        sn_part = np.split(sne,np.unique(qnm_number, return_index=True)[1][1:])
        # here we divide this all array into different parts by qnm_number, which means the array for n-th resonance
        burst_n = []

        for part_index in np.arange(0,len(hj_part)):
            hj_i = hj_part[part_index]
            arg_i = arg_part[part_index]
            y_i = y_part[part_index]
            sn_i = sn_part[part_index]
            burst_n = np.append(burst_n, logi0(np.abs(np.sum(4*self.df*hj_i*y_i*np.exp(1j*arg_i)/sn_i)))-np.sum(2*self.df*np.abs(hj_i)**2/sn_i ))
        burst = np.sum(burst_n)
        return burst