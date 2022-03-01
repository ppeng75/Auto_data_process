# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 10:17:11 2021

@author: pps
"""
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

class SpecProcess:
    def __init__(self,path,expstyle,filetype='.dat'):
        self.path = path
        self.expstyle = expstyle
        self.filetype = filetype
        self.pathlist = glob.glob(self.path+'/'+'*'+filetype)
        self.filenames = []
        for i in range(0,len(self.pathlist)):
            self.pathlist[i] = self.pathlist[i].replace('\\','/')
            iter1 = re.finditer('/', self.pathlist[i])
            iter2 = re.finditer('_', self.pathlist[i])
            # iter3 = re.finditer('\.dat',name)
            for j in iter1:
                id1 = j
            for k in iter2:
                id2 = k
            # for l in iter3:
                #     id3 = l
            id_backs = id1.end()
            id_unders = id2.start()
            #id_dot = id3.start()
            name = self.pathlist[i][id_backs:id_unders]
            if name in self.filenames:
                continue
            else:
                self.filenames.append(name)              
                
    def loadspecs(self):           
        Xvalue = dict()
        Yvalue = dict()
        Tvalue = dict()
        for i in range(0,len(self.filenames)):
            self.pathname = glob.glob(self.path+'/'+self.filenames[i]+'*'+self.filetype)   
            Nspec = len(self.pathname)
            data0 = np.loadtxt(self.pathname[0])
            Ndata,Ncol = np.shape(data0)
            x_store = np.zeros((Ndata,Nspec))
            y_store = np.zeros((Ndata,Nspec))
            t_store = np.zeros((Ndata,Nspec))
            del Nspec,Ndata
            for j in range(0,len(self.pathname)):
                filenum = str(j+1)
                data = np.loadtxt(self.path+'/'+self.filenames[i]+'_'+filenum+self.filetype)
                x_store[:,j] = data[:,1]
                if Ncol == 3:
                    y_store[:,j] = data[:,2]
                t_store[:,j] = data[:,0]
                del data
            # del self.pathname
            Xvalue[self.filenames[i]] = x_store
            Yvalue[self.filenames[i]] = y_store
            Tvalue[self.filenames[i]] = t_store
        return Xvalue,Yvalue,Tvalue
    
    def avespecs_sam(self,xdata,ydata,tdata):
        xave = dict()
        yave = dict()
        tave = dict()
        for i in range(0,len(self.filenames)):
            xsum = np.sum(xdata[self.filenames[i]],axis=1)
            ysum = np.sum(ydata[self.filenames[i]],axis=1)
            N_scans = np.size(xdata[self.filenames[i]],axis=1)
            xave[self.filenames[i]] = xsum/N_scans
            yave[self.filenames[i]] = ysum/N_scans
            tave[self.filenames[i]] = tdata[self.filenames[i]][:,0]
        return xave,yave,tave
    
    def avespecs_samref(self,xdata,ydata,tdata,Nsam,Nref):
        self.samid = list()
        self.refid = list()
        xave_sam = dict()
        yave_sam = dict()
        tave_sam = dict()
        xave_ref = dict()
        yave_ref = dict()
        tave_ref = dict() 
        for i in range(0,len(self.filenames)):
            N_scans = np.size(xdata[self.filenames[i]],axis=1)
            N_dps = np.size(xdata[self.filenames[i]],axis=0)
            N_cycles = N_scans//(Nsam+Nref)
            xsum_sam = np.zeros((N_dps))
            ysum_sam = np.zeros((N_dps))
            xsum_ref = np.zeros((N_dps))
            ysum_ref = np.zeros((N_dps))
            if N_scans < (Nsam+Nref):
                continue
            else:
                for j in range(0,N_scans):
                    if j % (Nsam+Nref) < Nsam:
                        self.samid.append(j)
                        xsum_sam = xsum_sam+xdata[self.filenames[i]][:,j]
                        ysum_sam = ysum_sam+ydata[self.filenames[i]][:,j]
                    else:
                        self.refid.append(j)
                        xsum_ref = xsum_ref+xdata[self.filenames[i]][:,j]
                        ysum_ref = ysum_ref+ydata[self.filenames[i]][:,j]
                xave_sam[self.filenames[i]] = xsum_sam/N_cycles/Nsam
                yave_sam[self.filenames[i]] = ysum_sam/N_cycles/Nsam
                tave_sam[self.filenames[i]] = tdata[self.filenames[i]][:,0]
                xave_ref[self.filenames[i]] = xsum_ref/N_cycles/Nref
                yave_ref[self.filenames[i]] = ysum_ref/N_cycles/Nref
                tave_ref[self.filenames[i]] = tdata[self.filenames[i]][:,Nsam]
        return xave_sam,yave_sam,tave_sam,xave_ref,yave_ref,tave_ref
    
    def Totalfield(self,xdata,ydata):
        filenames = list(xdata)
        for name in filenames:
            theta = np.arctan(ydata[name]/xdata[name])
            xdata[name] = xdata[name]/np.cos(theta)
        return xdata
                
    
    def plot_dict(self,xvalues,yvalues):
        plt.figure()
        for name in self.filenames:
            plt.plot(xvalues[name],yvalues[name],label=name)
            plt.legend(loc='best',fontsize=5)
            
    def plot_list(self,xvalues,yvalues):
        plt.figure()
        for i in range(len(yvalues)):
            plt.plot(xvalues[i],yvalues[i],label=self.filenames[i])
            plt.legend(loc='best',fontsize=5)
            
    
    
    def polarimetry(self,comp1,comp2,time,idxsam,idxref,idxref_select=[]):
        samx = []
        samsx = []
        samy = []
        samsy = []
        tsam = []
        fsam = []
        refx = []
        refsx = []
        refy = []
        refsy = []
        tref = []
        fref = []
        trans = []
        theta = []
        eta = []
        samnames = []
        refnames = []
        Efield = dict()
        #save sample spectrum data 
        for i in idxsam:
            samx.append(comp1[self.filenames[i]] + comp2[self.filenames[i]])
            samy.append(comp1[self.filenames[i]] - comp2[self.filenames[i]])
            tsam.append(time[self.filenames[i]])
            samnames.append(self.filenames[i])
        #fft sample data and get polarization angle as well as ellipticity  
        for i in range(len(samx)):
            [freq,sx] = self.fftx(samx[i],tsam[i],2)
            [freq,sy] = self.fftx(samy[i],tsam[i],2)
            fsam.append(freq)
            samsx.append(sx)
            samsy.append(sy)
            Ecra = sx + 1j*sy
            Ecri = sx - 1j*sy
            theta.append((np.unwrap(np.angle(Ecra))-np.unwrap(np.angle(Ecri)))/2)
            eta.append((abs(Ecri)-abs(Ecra))/(abs(Ecri)+abs(Ecra)))
            del freq,sx,sy
        #save referene spectrumm data
        if idxref == []:
            for k in idxref:
                k = int(k)
                refx.append(comp1[self.filenames[k]] + comp2[self.filenames[k]])
                refy.append(comp1[self.filenames[k]] - comp2[self.filenames[k]])
                tref.append(time[self.filenames[k]])
                refnames.append(self.filenames[k])
        else:
            for k in idxref_select:
                id1 = int(idxref[k])
                id2 = int(idxref[k+1])
                refx.append(1/2*(comp1[self.filenames[id1]]+comp1[self.filenames[id2]]) + 1/2*(comp2[self.filenames[id1]]+comp2[self.filenames[id2]]))
                refy.append(1/2*(comp1[self.filenames[id1]]+comp1[self.filenames[id2]]) - 1/2*(comp2[self.filenames[id1]]+comp2[self.filenames[id2]]))
                tref.append(time[self.filenames[id1]])
                refnames.append(self.filenames[id1])
        #fft reference data and calculate transmission            
        for i in range(len(tref)):
            [freq,sx] = self.fftx(refx[i],tref[i],2)
            [freq,sy] = self.fftx(refy[i],tref[i],2)
            refsx.append(sx)
            refsy.append(sy)
            fref.append(freq)
            Asam = np.sqrt(abs(samsx[i])**2+abs(samsy[i])**2)
            Aref = np.sqrt(abs(refsx[i])**2+abs(refsy[i])**2)
            trans.append(Asam/Aref)
            del freq,sx,sy
        
        Efield['sam xvalues td'] = samx 
        Efield['sam xvalues fd'] = samsx 
        Efield['sam yvalues td'] = samy
        Efield['sam yvalues fd'] = samsy 
        Efield['sam time axis'] = tsam
        Efield['sam frequency axis'] = fsam
        Efield['ref xvalues td'] = refx 
        Efield['ref xvalues fd'] = refsx 
        Efield['ref yvalues td'] = refy
        Efield['ref yvalues fd'] = refsy 
        Efield['ref time axis'] = tref
        Efield['ref frequency axis'] = fref
        Efield['polarization angles'] = theta
        Efield['Ellipticity angles'] = eta
        Efield['transmission'] = trans
        Efield['sam names'] = samnames
        Efield['ref names'] = refnames
        
        return Efield

    def polarimetry_new(self,comp1,comp2,time,idxsam,idxref,idxref_select=None):
        self.samx = dict()
        self.samsx = dict()
        self.samy = dict()
        self.samsy = dict()
        self.tsam = dict()
        self.fsam = dict()
        self.refx = dict()
        self.refsx = dict()
        self.refy = dict()
        self.refsy = dict()
        self.tref = dict()
        self.fref = dict()
        self.trans = dict()
        self.theta = dict()
        self.eta = dict()
        #save sample spectrum data, do fft and calculate polarization rotation and ellipticity
        for i in idxsam:
            self.samx[self.filenames[i]] = comp1[self.filenames[i]] + comp2[self.filenames[i]]
            self.samy[self.filenames[i]] = comp1[self.filenames[i]] - comp2[self.filenames[i]]
            self.tsam[self.filenames[i]] = time[self.filenames[i]]
            [freq,sx] = self.fftx(self.samx[self.filenames[i]],self.tsam[self.filenames[i]],2)
            [freq,sy] = self.fftx(self.samy[self.filenames[i]],self.tsam[self.filenames[i]],2)
            self.fsam[self.filenames[i]] = freq
            self.samsx[self.filenames[i]] = sx
            self.samsy[self.filenames[i]] = sy
            Ecra = sx + 1j*sy
            Ecri = sx - 1j*sy
            self.theta[self.filenames[i]] = (np.unwrap(np.angle(Ecra))-np.unwrap(np.angle(Ecri))/2)
            self.eta[self.filenames[i]] = (abs(Ecri)-abs(Ecra))/(abs(Ecri)+abs(Ecra))
            del freq,sx,sy
        #save referene spectrumm data, do fft on data
        if idxref_select is None:
            for k in idxref:
                k = int(k)
                self.refx[self.filenames[k]] = comp1[self.filenames[k]] + comp2[self.filenames[k]]
                self.refy[self.filenames[k]] = comp1[self.filenames[k]] - comp2[self.filenames[k]]
                self.tref[self.filenames[k]] = time[self.filenames[k]]
                [freq,sx] = self.fftx(self.refx[self.filenames[i]],self.tref[self.filenames[i]],2)
                [freq,sy] = self.fftx(self.refy[self.filenames[i]],self.tref[self.filenames[i]],2)
                self.fref[self.filenames[i]] = freq
                self.refsx[self.filenames[i]] = sx
                self.refsy[self.filenames[i]] = sy
                del freq,sx,sy
        else:
            for k in idxref_select:
                id1 = int(idxref[k])
                id2 = int(idxref[k+1])
                self.refx.append(1/2*(comp1[self.filenames[id1]]+comp1[self.filenames[id2]]) + 1/2*(comp2[self.filenames[id1]]+comp2[self.filenames[id2]]))
                self.refy[self.filenames[id1]] = 1/2*(comp1[self.filenames[id1]]+comp1[self.filenames[id2]]) - 1/2*(comp2[self.filenames[id1]]+comp2[self.filenames[id2]])
                self.tref[self.filenames[id1]] = time[self.filenames[id1]]
                [freq,sx] = self.fftx(self.refx[self.filenames[i]],self.tref[self.filenames[i]],2)
                [freq,sy] = self.fftx(self.refy[self.filenames[i]],self.tref[self.filenames[i]],2)
                self.fref[self.filenames[i]] = freq
                self.refsx[self.filenames[i]] = sx
                self.refsy[self.filenames[i]] = sy
                del freq,sx,sy
        #calculate transmission
        for i in range(0,len(list(self.samx))):
            samnames = list(self.samx)
            refnames = list(self.refx)
            self.trans[samnames[i]] = np.sqrt(abs(self.samsx[samnames[i]])**2+abs(self.samsy[samnames[i]])**2)/np.sqrt(abs(self.refsx[refnames[i]])**2+abs(self.refsy[refnames[i]])**2)

    def getindex(self,freq,E_sam,E_ref,L,E_echo=None,c=3e8):
        freq = self.formatinput(freq)
        E_sam = self.formatinput(E_sam)
        E_ref = self.formatinput(E_ref)
        fname = list(freq)[0]
        samname = list(E_sam)[0]
        refname = list(E_ref)[0]
        omega = freq[fname]*2*np.pi
        if E_echo is None:
            phi = np.unwrap(np.angle(E_sam[samname]/E_ref[refname]))
            n = dict()
            n[samname] = phi/omega/L*c+1
        else:
            E_echo = self.formatinput(E_echo)
            echoname = list(E_echo)[0]
            l = np.linspace(0.9*L,1.1*L,1000)
            L0 = np.ones(len(l))*L
            phi_sam = np.unwrap(np.angle(E_sam[samname]/E_ref[refname]))
            phi_echo = np.unwrap(np.angle(E_echo[echoname]/E_ref[refname]))
            n_sam = phi_sam[:,np.newaxis]/omega[:,np.newaxis]/l[np.newaxis,:]*c+1
            n_echo = phi_echo[:,np.newaxis]/omega[:,np.newaxis]/l[np.newaxis,:]*c+1                
            dn = n_sam-n_echo
            for i in range(0,len(freq)):
                L0[i] = self.findzeropoints(l,dn[i,:])
            n = dict()
            n[samname] = phi_sam/omega/L0*c+1
        return n
     
    def findzeropoints(self,x,y):
        for i in range(0,len(x)-1):
            x_zeros = []
            if y[i]*y[i+1] <= 0:
               x_zeros.append((x[i]+x[i+1])/2)
        return x_zeros
    
    def specchop(self,x,y,x0,xt):
        idx0 = np.where(x>=x0)[0]
        idxt = np.where(x<=xt)[0]
        x_new = x[idx0:idxt]
        y_new = y[idx0:idxt]
        return x_new,y_new
    
        
    
    def getEfield(self,t0,E0,r_focus,w0,eps0=8.85e-12,c=3e8):
        # r = np.linspace(0,r_focus,100)
        # P = np.zeros((len(t0),len(r)))
        En = E0/max(E0)
        rc = r_focus/1.5
        P = c*eps0*En**2*rc**2/2*(1-np.exp(-r_focus**2/rc**2))
        # for i in range(0,len(r)):
        #     En_r = En*np.exp(-r[i]**2/rc**2)
        #     P[:,i] = c*eps0*En_r**2*np.pi*(r[2]-r[1])*r[i]
        # w_r = np.trapz(P,x=t0,axis=0)
        # w = np.trapz(w_r,x=r)
        w = np.trapz(P,x=t0)
        E_out = np.sqrt(w0/w)/1e5
        return E_out
    
    def formatinput(self,x):
        if type(x) is list:
            x_new = dict()
            for i in range(len(x)):
                name = 'x'+str(i)
                x_new[name] = x[i]
        if type(x) is np.ndarray:
            x_new = dict()
            if np.size(x) == len(x):
                    name = 'x0'
                    x_new[name] = x
            else:
                for i in range(np.size(x,1)):
                    name = 'x'+str(i)
                    x_new[name] = x[:,i]
        if type(x) is dict:
            x_new = x
        return x_new
    
    def transmission(self,sxsam,sxref):
        sxsam = self.formatinput(sxsam)
        sxref = self.formatinput(sxref)
        trans = dict()
        names_sam = list(sxsam)
        names_ref = list(sxref)
        N = len(sxsam)
        for i in range(0,N):
            trans[names_sam[i]] = abs(sxsam[names_sam[i]])/abs(sxref[names_ref[i]])
        return trans
    
    def combinedict(self,xsam,xref,subs='ref'):
        x_new = dict()
        x_new = xsam
        for key,value in xref.items():
            if key in xsam.keys(): 
                key_new = key+subs
                x_new[key_new] = value
            else:
                x_new[key] = value
        return x_new
    
class Basic_functions:
    def __init__(self):
        self.default = 0
        
    def formatinput(self,x):
        if type(x) is list:
            x_new = dict()
            for i in range(len(x)):
                name = 'x'+str(i)
                x_new[name] = x[i]
        if type(x) is np.ndarray:
            x_new = dict()
            if np.size(x) == len(x):
                    name = 'x0'
                    x_new[name] = x
            else:
                for i in range(np.size(x,1)):
                    name = 'x'+str(i)
                    x_new[name] = x[:,i]
        if type(x) is dict:
            x_new = x
        return x_new
    
    def specchop(self,x,y,x0,xt):
        idx0 = np.where(x>=x0)[0][0]
        idxt = np.where(x<=xt)[0][-1]
        x_new = x[idx0:idxt]
        y_new = y[idx0:idxt]
        return x_new,y_new
    
    def findzeropoints(self,x,y):
        for i in range(0,len(x)-1):
            x_zeros = []
            if y[i]*y[i+1] <= 0:
               x_zeros.append((x[i]+x[i+1])/2)
        return x_zeros
    
    def getindex(self,freq,E_sam,E_ref,L,E_echo=None,c=3e8):
        freq = self.formatinput(freq)
        E_sam = self.formatinput(E_sam)
        E_ref = self.formatinput(E_ref)
        fname = list(freq)[0]
        samname = list(E_sam)[0]
        refname = list(E_ref)[0]
        omega = freq[fname]*2*np.pi
        if E_echo is None:
            phi = np.unwrap(np.angle(E_sam[samname]/E_ref[refname]))
            n = dict()
            n[samname] = phi/omega/L*c+1
        else:
            E_echo = self.formatinput(E_echo)
            echoname = list(E_echo)[0]
            l = np.linspace(0.9*L,1.1*L,1000)
            L0 = np.ones(len(l))*L
            phi_sam = np.unwrap(np.angle(E_sam[samname]/E_ref[refname]))
            phi_echo = np.unwrap(np.angle(E_echo[echoname]/E_ref[refname]))
            n_sam = phi_sam[:,np.newaxis]/omega[:,np.newaxis]/l[np.newaxis,:]*c+1
            n_echo = phi_echo[:,np.newaxis]/omega[:,np.newaxis]/l[np.newaxis,:]*c+1                
            dn = n_sam-n_echo
            for i in range(0,len(freq)):
                L0[i] = self.findzeropoints(l,dn[i,:])
            n = dict()
            n[samname] = phi_sam/omega/L0*c+1
        return n
    
    def fftx(self,xvalues,tvalues,pad):
        xvalues = self.formatinput(xvalues)
        tvalues = self.formatinput(tvalues)
        sx = dict()
        freq = dict()
        for name in list(xvalues):
            ts = abs(tvalues[name][1]-tvalues[name][0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            xvalue_reduce = xvalues[name]-np.average(xvalues[name])
            SX = np.fft.fft(xvalue_reduce,n=NFFT)
            sx[name] = SX[0:int(NFFT/2)]
            freq[name] = fs/2*np.linspace(0,1,len(sx[name]))
        return freq,sx     
    
    def array_fftx(self,xvalues,tvalues,pad=2):
        ts = abs(tvalues[1]-tvalues[0])
        NFFT = int(2**(np.ceil(np.log2(len(tvalues)))+pad))
        fs = 1/ts
        xvalue_reduce = xvalues-np.average(xvalues)
        SX = np.fft.fft(xvalue_reduce,n=NFFT)
        sx = SX[0:int(NFFT/2)]
        freq = fs/2*np.linspace(0,1,len(sx))
        return freq,sx 
    
    def ifftx(self,sxvalues,fvalues,pad):
        sxvalues = self.formatinput(sxvalues)
        fvalues = self.formatinput(fvalues)
        fs = abs(fvalues[1]-fvalues[0])
        NFFT = int(2**(np.ceil(np.log2(len(fvalues)))+pad))
        ts = 1/fs
        SX = np.fft.ifft(sxvalues,n=NFFT)
        sx = SX[0:int(NFFT/2)]
        time = ts/2*np.linspace(0,1,len(sx))
        return time,sx
    
    def dict_getabs(self,sxvalues):
        names = list(sxvalues)
        for name in names:
            sxvalues[name] = abs(sxvalues[name])
        return sxvalues
    
    def superosc(self,xvalues):
        xvalues = self.formatinput(xvalues)
        xsum = dict()
        for name in list(xvalues):
            xsum['SO'] = xsum['SO']+xvalues[name]
        return xsum
    
    def array_getderivative(self,xvalues,yvalues):
        xdata = xvalues
        ydata = yvalues
        deri_store = np.zeros(len(xdata))
        for i in range(0,len(xdata)):
            if i == 0:
                deri_store[i] = (ydata[i+1]-ydata[i])/(xdata[i+1]-xdata[i])
            else:
                if i == len(xdata)-1:
                    deri_store[i] = (ydata[i]-ydata[i-1])/(xdata[i]-xdata[i-1])
                else:
                    deri_store[i] = (ydata[i+1]-ydata[i-1])/(xdata[i+1]-xdata[i-1])
        return deri_store
     
    def dict_normalize(self,x):
        x = self.formatinput(x)
        xnew = dict()
        for name in x:
            xnew[name] = x[name]/max(x[name])
        return xnew
                             
         