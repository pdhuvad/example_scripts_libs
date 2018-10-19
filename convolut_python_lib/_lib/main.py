# coding: utf-8



__author__ = "Pratik Dhuvad"
__email__ = "tud12582@temple.edu"
__status__ = "Development"
__date__ = "March 20, 2018"



import numpy as np
import os,sys,shutil,subprocess,pickle
from os.path import join

from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen as mg

from scipy import interpolate
from scipy import signal




class mconv:
    '''
    Does convolution with Gaussian or Lorentzian windows
    '''
    
    def __init__(self, datax, datay):
        
        self.datax = datax
        self.datay = datay

        self.NPoints = len(datax)
        self.X = np.array(datax)


            
    def w_gaussian(self, M, wsigma, dx):
        ''' see : https://github.com/scipy/scipy/blob/v0.19.0/scipy/signal/windows.py#L1159-L1219
        M should be odd number'''
        if wsigma <=0 : wsigma = 0.00001
        wsigma = wsigma/dx
        n = np.arange(0, M) - (M - 1.0) / 2.0
        wsigma2 = 2 * wsigma * wsigma
        wg = np.exp(-n ** 2 / wsigma2)    
        return wg

    def w_lorentzian(self, M, wgamma, dx):
        ''' '''
        wgamma = wgamma/dx
        if wgamma <=0 : wgamma = 0.00001
        n = np.arange(0, M) - (M - 1.0) / 2.0
        wl = 1 / ( ((2*n)/wgamma)**2 + 1  )    
        return wl             
           
        
    def Gaussian(self,sigma=None,fwhm=None,saveto=None):
        
        if fwhm:
            if sigma:
                print('ignoring input sigma')
        elif sigma:
            fwhm = sigma * np.sqrt(8 * np.log(2))
        else: 
            raise ValueError('sigma/fwhm was not set....')  
           
        self.sigma = fwhm/np.sqrt(8 * np.log(2))   
        
        M=101
        dx=self.datax[2]-self.datax[1]
        win = self.w_gaussian(M, self.sigma, dx)
        out = signal.convolve(self.datay, win, mode='same') / sum(win)
       
        if saveto:
            if fmt is None:
                fmt="%18.6e %18.6e"                
            of = np.column_stack( (self.X, out) )
            np.savetxt(str(saveto), of, delimiter=" ", fmt=fmt )       
            
        return [self.X, out]
            

    def Lorentzian(self,gamma=None,fwhm=None,saveto=None,M=None):
        # in Lorentzian fwhm is equal to gamma
        
        if fwhm:
            if gamma:
                print('in Lorentzian fwhm is equal to gamma')
        elif gamma:
            fwhm = gamma
        else: 
            raise ValueError('sigma/fwhm was not set....')   
            
        self.gamma = fwhm 
        
        if M is None: M=1001
        dx=self.datax[2]-self.datax[1]    
        win = self.w_lorentzian(M, self.gamma, dx)
        out = signal.convolve(self.datay, win, mode='same') / sum(win)

            
        if saveto:
            if fmt is None:
                fmt="%18.6e %18.6e"                
            of = np.column_stack( (self.X, out) )
            np.savetxt(str(saveto), of, delimiter=" ", fmt=fmt )       
            
        return [self.X, out]    
    
    
    def LorentzianVL(self,saveto=None,fmt=None,M=None,A=None,B=None,offset=None):
        # gamma = A(x-offset) + B   
        
        if A is None: A=0.1
        if B is None: B=0           
        if offset is None: offset=self.datax[0]  
        
        gammas = []    
        for i,d in enumerate(self.datax):
            g = max(0,A*(d-offset)) + B
            gammas.append(g)
  
        if M is None: M=1001   
        dx=self.datax[2]-self.datax[1]
        out = np.zeros(self.NPoints)
        
                
        for i, [gg, x0, y0] in enumerate(zip(gammas,self.datax,self.datay)):
            win = self.w_lorentzian(M, gg, dx)
            c = signal.convolve(self.datay, win, mode='same') / sum(win)
            out[i] = c[i]
                        
        if saveto:
            if fmt is None:
                fmt="%18.6e %18.6e"                
            of = np.column_stack( (self.X, out) )
            np.savetxt(str(saveto), of, delimiter=" ", fmt=fmt )           
            
        #return [self.X, out], [self.datax, gammas]   
        return [self.X, out]

    
    
    def LorentzianVE2(self,saveto=None,fmt=None,M=None,offset=None):
        # gamma = (x-offset)^2 + B   
        
        if A is None: A=0.1
        if B is None: B=0           
        if offset is None: offset=self.datax[0]  
        
        gammas = []    
        for i,d in enumerate(self.datax):
            g = max(0,A*(d-offset)) + B
            gammas.append(g)
  
        if M is None: M=1001   
        dx=self.datax[2]-self.datax[1]
        out = np.zeros(self.NPoints)
        
                
        for i, [gg, x0, y0] in enumerate(zip(gammas,self.datax,self.datay)):
            win = self.w_lorentzian(M, gg, dx)
            c = signal.convolve(self.datay, win, mode='same') / sum(win)
            out[i] = c[i]
                        
        if saveto:
            if fmt is None:
                fmt="%18.6e %18.6e"                
            of = np.column_stack( (self.X, out) )
            np.savetxt(str(saveto), of, delimiter=" ", fmt=fmt )           
            
        #return [self.X, out], [self.datax, gammas]   
        return [self.X, out]



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
class mXANES:
    
    def __init__(self, data=None, data_loadfrom=None, srange=None, structure=None, vcncutoff=5.0, ca=None, Eonset=None, 
                 xanesid=None, source=None, edge=None, multiplicity=1): 

        if data is None:
            if data_loadfrom:
                data = np.loadtxt(data_loadfrom, unpack=True, comments='#', skiprows=0)
                self.E0 = np.array(data[0])
                self.I0 = np.array(data[1])                
            else:
                self.E0 = np.array([])                
                self.I0 = np.array([])                               
        else:
            self.E0 = np.array(data[0])
            self.I0 = np.array(data[1])
            
        if srange:
            sel = (self.E0 >= srange[0]) & (self.E0 <= srange[1])
            self.E0 = self.E0[sel]
            self.I0 = self.I0[sel]            

        # Energy offset
        if Eonset is None:
            self.Eonset = self.E0[0]
        else:
            self.Eonset = Eonset
            
        # XANES id
        if xanesid is None:
            self.xanesid = 'not_set'
        else:
            self.xanesid = xanesid 
            
        # XANES edge
        if edge is None:
            self.edge = 'not_set'
        else:
            self.edge = edge         
            
        # source
        if source is None:
            self.source = 'not_set'
        else:
            self.source = source         
 
        # source
        if ca is None:
            self.ca = 'not_set'
        else:
            self.ca = ca     
            
        # structure
        if structure:
            try:               
                self.ca = structure[0][structure[1]].species_string            
                nnfinder = VoronoiNN(cutoff=vcncutoff,allow_pathological=True)
                vcn = nnfinder.get_cn(structure[0], structure[1], use_weights=True)
                self.vcn = vcn
                self.structure = [structure[0],structure[1],vcn]
            except Exception as exc:
                print(exc)
                self.structure = [structure[0],structure[1],[]]
                self.vcn = []
        else:
            self.structure = [[],[],[]] 
            self.vcn = []
        
        self.peaks = []
        self.multiplicity = multiplicity
        
        self.E = None
        self.I = None
        
                                            
    def Interpolate(self,iterprange,stepsize=0.05):               
        if self.E[0] > iterprange[0]:
            # left padding
            npts = int((self.E[0]-iterprange[0])/stepsize)+1
            x_patch = np.linspace(iterprange[0],self.E[0]-stepsize,npts)
            y_patch=np.empty(len(x_patch)); y_patch.fill(self.I[0])
            self.E = np.concatenate((x_patch,self.E.T), axis=0)
            self.I = np.concatenate((y_patch,self.I.T), axis=0)                       
        if self.E[-1] < iterprange[1]:
            # right padding            
            npts = int((iterprange[1]-self.E[-1])/stepsize)+2
            x_patch = np.linspace(self.E[-1],iterprange[1],npts)
            y_patch=np.empty(len(x_patch)); y_patch.fill(self.I[-1])
            self.E = np.concatenate((self.E.T,x_patch), axis=0)
            self.I = np.concatenate((self.I.T,y_patch), axis=0)             
        f = interpolate.interp1d(self.E,self.I,kind='linear')
        self.E = np.linspace(iterprange[0],iterprange[1], int((iterprange[1]-iterprange[0])/stepsize)+1  )      
        self.I = f(self.E) 
                
    def FindPeaks(self,xin=None,yin=None,srangep=None):        
        if (xin is None) or (yin is None):
            xsearch = self.E0
            ysearch = self.I0
        else:
            xsearch = xin
            ysearch = yin                       
        if srangep:
            sel = (xsearch >= srangep[0]) & (xsearch <= srangep[1])
            xsearch = xsearch[sel]; ysearch = ysearch[sel]        
        ipeaks = signal.argrelextrema(ysearch, np.greater)[0]         
        peaks = []
        for i in ipeaks:
            peaks.append([xsearch[i],ysearch[i]])
        self.peaks = peaks        
        return peaks
  
    def yscale_by(self,yscale):
        self.I = self.I*yscale
        
    def normalize_to(self,nstr):
        if nstr == 'max':            
            self.I = self.I/max(self.I)
        elif nstr == 'tail':
            self.I = self.I/self.I[-1]      
        else:
            self.I = self.I  


            
    def broaden0(self, g_sigma=None,g_fwhm=None, l_gamma=None,l_fwhm=None,
                lvl=None):
        
        E_in, I_in = self.E0, self.I0
  
        if g_sigma:
            E_out, I_out = mconv(E_in, I_in).Gaussian(sigma=g_sigma)
        elif g_fwhm:
            E_out, I_out = mconv(E_in, I_in).Gaussian(fwhm=g_fwhm)
            
        if l_gamma:
            E_out, I_out = mconv(E_in, I_in).Lorentzian(gamma=l_gamma)            
        elif l_fwhm:
            E_out, I_out = mconv(E_in, I_in).Lorentzian(fwhm=l_fwhm)            
              
        if lvl:            
            if len(lvl) == 1:
                E_out, I_out = mconv(E_in, I_in).LorentzianVL(A=lvl[0],B=None,offset=None)
            elif len(lvl) == 2:
                E_out, I_out = mconv(E_in, I_in).LorentzianVL(A=lvl[0],B=lvl[1],offset=None)               
            else:
                E_out, I_out = mconv(E_in, I_in).LorentzianVL(A=lvl[0],B=lvl[1],offset=lvl[2])  
                         
        self.E0, self.I0 = E_out, I_out    



 
    def transform(self,srange=None,irange=None,e0shift=False,
                  y0shift=True,normalize='max',xshift=None):
        if self.E == None:
            self.E = self.E0.copy()
            self.I = self.I0.copy()        
        if e0shift:
            self.E = self.E -self.Eonset                        
        if xshift:
            self.E = self.E + xshift
        if irange:
            self.Interpolate(irange)                  
        if y0shift:
            self.I = self.I -self.I[0]              
        if normalize == 'max':            
            self.I = self.I/max(self.I)
        elif normalize == 'tail':
            self.I = self.I/self.I[-1]
        elif normalize == 'none':
            self.I = self.I         
        else:
            self.I = self.I/max(self.I) 
       
        
            
            







