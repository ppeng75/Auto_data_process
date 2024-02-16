# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 14:04:00 2023

@author: Admin
"""
import numpy as np
import dataprocessor as dp
import matplotlib.pyplot as plt
import plottool_v6 as pt
from numpy import pi
from scipy.signal import stft, spectrogram
import glob
import re
# from numba import njit, jit
from scipy.integrate import simps
from numpy.fft import fft
import tkinter as tk
from tkinter import ttk
import os
from copy import deepcopy
    


def get_fnames(path, filetype='.dat'):
    pathlist = glob.glob(path+'/'+'*'+filetype)
    filenames = []
    for i in range(0,len(pathlist)):
        pathlist[i] = pathlist[i].replace('\\','/')
        iter1 = re.finditer('/', pathlist[i])
        iter2 = re.finditer('_', pathlist[i])
        for j in iter1:  #loop to last / sign
            id1 = j
        for k in iter2: #loop ot last _ sign
            id2 = k
        id_backs = id1.end() #get the index after last /
        id_unders = id2.start() # get the index before last _
        name = pathlist[i][id_backs:id_unders]
        if name in filenames:
            continue
        else:
            filenames.append(name)
    return filenames

def load_single_spec(folder, fname, ftype = '.dat'):
    folder = folder.replace('\\','/')
    flist = glob.glob(folder+'/'+fname+'_*[0-9]'+ftype)
    dim_y = 1
    for i in range(0, len(flist)):
        fnum = str(i+1)
        try:
            data = np.loadtxt(folder+'/'+fname+'_'+fnum+ftype)
        except OSError:
            continue
        except ValueError:
            data = np.loadtxt(folder+'/'+fname+'_'+fnum+ftype, dtype = complex)
        else:
            data = np.loadtxt(folder+'/'+fname+'_'+fnum+ftype)
            if i == 0:
                t_store = data[:, 0]
                x_store = data[:, 1]
                try:
                    y_store = data[:, 2]
                except:
                    dim_y = 0
                    y_store = np.zeros_like(t_store)
                else:
                    y_store = data[:, 2]
            else:
                if dim_y == 0:
                    t_new = data[:, 0]
                    x_new = data[:, 1]
                    y_new = np.zeros_like(t_new)
                    t_store = np.colume_stack((t_store, t_new))
                    x_store = np.column_stack((x_store, x_new))
                    y_store = np.column_stack((y_store, y_new))
                else:
                    t_new = data[:, 0]
                    x_new = data[:, 1]
                    y_new = data[:, 2]
                    t_store = np.column_stack((t_store, t_new))
                    x_store = np.column_stack((x_store, x_new))
                    y_store = np.column_stack((y_store, y_new))
    return t_store, x_store, y_store

def save_data(path, fname, t, x, y = None):
    if os.path.exists(path):
        pass
    else:
        os.makedirs(path)
    if y is None:
        if len(t) == len(x):
            y = np.zeros_like(t)
            spec = np.stack((t,x,y), axis = 1)
            np.savetxt(path+fname, spec)
        else:
            raise Exception('dimension of time and amplitude does not match')
    else:
        if len(t) == len(x) and len(t) == len(y):
            spec = np.stack((t,x,y), axis = 1)
            np.savetxt(path+fname, spec)
        else:
            raise Exception('Dimension of time and amplitude does not match')
            
def str_get_D(string):
    str_seg_all = string.split('_')
    pick_str = ''
    for str_seg in str_seg_all:
        if 'D' in str_seg:
            pick_str = str_seg
            continue
    num = ''
    for char in pick_str:
        if char.isdigit():
            num = num + char
    try:
        num = int(num)
    except ValueError:
        num = 0
    else:
        num = int(num)
    return num




def str_get_Delay(string):
    str_seg_all = string.split('_')
    pick_str = ''
    for str_seg in str_seg_all:
        if 'Delay' in str_seg:
            pick_str = str_seg
            continue
    num = ''
    for char in pick_str:
        if char.isdigit():
            num = num + char
    try:
        num = int(num)
    except ValueError:
        num = 0
    else:
        num = int(num)
    return num

# def str_get_D(string):
#     str_seg_all = string.split('_')
#     pick_str = ''
#     for str_seg in str_seg_all:
#         if 'D' in str_seg:
#             pick_str = str_seg
#             continue
#     num = ''
#     for char in pick_str:
#         if char.isdigit():
#             num = num + char
#     try:
#         num = int(num)
#     except ValueError:
#         num = 0
#     else:
#         num = int(num)
#     return num


# @njit
def interferometric_vis(sx):
    # sx_normal = abs(sx)/max(abs(sx))
    sx_normal = np.absolute(sx)
    i_max = np.max(sx_normal**2)
    i_min = np.min(sx_normal**2)
    ivis = (i_max-i_min)/(i_max+i_min)
    return ivis


def average_spec_noref(t, x):
    x_ave = {}
    t_ave = {}
    for key in list(x):
        try:
            t_ave[key] = t[key][:,0]
        except IndexError:
            t_ave[key] = t[key]
            x_ave[key] = x[key]
        else:
            t_ave[key] = t[key][:,0]
            x_ave[key] = np.average(x[key], axis = 1)
    return t_ave, x_ave

# @njit
def generate_gauss_cos(t, f, sigma = 1):
    return np.cos(2*pi*f*t)*np.exp(-(t)**2/sigma**2)


def dict_ave_samref(x, t, nsam, nref):
    samx = dict()
    samt = dict()
    refx = dict()
    reft = dict()
    comp1 = dict()
    comp2 = dict()
    for key in list(x):
        x_temp = x[key]
        t_temp = t[key]
        samx_sum = 0
        refx_sum = 0
        n_cycle = round(np.size(x_temp, axis = 1)/(nsam+nref))
        for i in range(np.size(x_temp, axis=1)):
            if i % (nsam+nref) < nsam:
                samx_sum = x_temp[:,i] + samx_sum
            else:
                refx_sum = x_temp[:,i] + refx_sum
        samx[key] = samx_sum/n_cycle
        samt[key] = t_temp[:,0]
        reft[key] = t_temp[:,nref]
        refx[key] = refx_sum/n_cycle
        comp1[key] = samx[key]+refx[key]
        comp2[key] = samx[key]-refx[key]
    return samx,samt,refx,reft,comp1,comp2

def array_ave_samref(x, t, nsam, nref):
    samx_sum = 0
    refx_sum = 0
    n_cycle = round(np.size(x, axis = 1)/(nsam+nref))
    for i in range(0, np.size(x, axis=1)):
        if i % (nsam+nref) < nsam:
            samx_sum = x[:,i] + samx_sum
        else:
            refx_sum = x[:,i] + refx_sum
    samx = samx_sum/n_cycle
    samt = t[:,0]
    reft = t[:,nref]
    refx = refx_sum/n_cycle
    comp1 = samx+refx
    comp2 = samx-refx
    return samx,samt,refx,reft,comp1,comp2

def array_cal_J_polar(t, E1x, E1y, E2x, E2y, t_cen = None):
    if len(t) == len(E1x):
        if len(E1x) == len(E2x):
            tw_len = []
            J = []
            if t_cen == None:
                idx0 = int(np.floor(len(E1x)/2))
                idx1 = idx0 + 1
            else:
                idx0 = np.where(t<=t_cen)[0][-1]
                idx1 = idx0 + 1
            while idx0 >= 0 and idx1 <= len(E1x)-1:
                tw_len.append(abs(t[idx1]-t[idx0])/2)
                E1x_w = E1x[idx0:idx1]
                E1y_w = E1y[idx0:idx1]
                E2x_w = E2x[idx0:idx1]
                E2y_w = E2y[idx0:idx1]
                J.append(simps((E1x_w-E2x_w)**2+(E1y_w-E2y_w)**2)/simps(E1x_w**2+E1y_w**2+E2x_w**2+E2y_w**2)*4)
                idx0 = idx0 - 1
                idx1 = idx1 + 1
        else:
            raise Exception('number of data for E1 and E2 does not match')
    else:
        raise Exception('number of data for t and E does not match')
    return np.array(tw_len), np.array(J)
    
        
def FaradayRotate(freq,B):
    '''
    

    Parameters
    ----------
    freq : 1d array
        DESCRIPTION. frequency domain obtained from Fourier transform input waves
    B : TYPE float64
        DESCRIPTION. External magnetic field strength in T

    Returns
    -------
    T_x : 1d array
        DESCRIPTION. complex transmission coefficient in x direction
    T_y : 1d array
        DESCRIPTION. complex transmission coefficient in y direction

    '''
    d=5*1e-4 #sample thickness
    w=2*pi*freq*1e12 #angular frequency 
    #set parameters for calculation
    n_air=1 # index of air
    e=1.6e-19 #electron charge
    c=3e8 # speed of light
    v_e=0.86*1e12 # electron scattering rate
    m_e=0.018*9.1e-31 # mass of electron
    w_t=2*pi*5.90e12 # transverse phonon frequency
    w_l=2*pi*5.54e12 # longitudinal phonon frequency
    epi_b=13.82 # background permittivity
    w_ce=e*B/(m_e) # cyclotron frequency
    w_e=0.28*2*pi*1e12*np.sqrt(epi_b) #plasma frequency is set to constant 0.28
    # w_h=(4*pi*Nh*e^2/m_h*1);
    v_ph=2*pi*1e12 # phonon damping rate
 
   
    epi_ph = epi_b*((w_t**2-w_l**2)/(w_t**2-w**2-1j*v_ph*w))

    ncra=np.sqrt(epi_b-w_e**2/(w*(w+1j*v_e-w_ce))+epi_ph)
    ncri=np.sqrt(epi_b-w_e**2/(w*(w+1j*v_e+w_ce))+epi_ph)

    t_cra_as=2*n_air/(n_air+ncra) #calculate Fresnel transmission coefficient
    t_cra_sa=2*ncra/(n_air+ncra)
    t_cri_as=2*n_air/(n_air+ncri)
    t_cri_sa=2*ncri/(n_air+ncri)
    # calculate complex transmission coefficents
    T_cra =t_cra_as*t_cra_sa*np.exp(1j*(ncra-1)*w*d/c)
    T_cri =(t_cri_as*t_cri_sa*np.exp(1j*(ncri-1)*w*d/c))
    T_x = (T_cra+T_cri)/np.sqrt(2)
    T_y = (T_cra-T_cri)/np.sqrt(2)/1j
        
    return T_x, T_y


def fftx(t, x, pad=2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_reduce = x-np.average(x)
    SX = np.fft.fft(x_reduce,n=NFFT)
    sx = SX[0:int(NFFT/2)]
    freq = fs/2*np.linspace(0,1,len(sx))
    return freq, sx 


def ifftx(t,sx):
    NFFT = 2*len(sx)
    X = np.fft.ifft(sx,n=NFFT)
    x = X[0:len(t)]
    return x


def fftx_blackman(t, x, pad = 2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_reduce = x-np.average(x)
    x_filter = x_reduce*np.blackman(len(x_reduce))
    SX = np.fft.fft(x_filter,n=NFFT)
    sx = SX[0:int(NFFT/2)]
    freq = fs/2*np.linspace(0,1,len(sx))
    return freq, sx


def fftx_stft(t, x, pad = 2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_normal = x-np.average(x)
    f_stft, t_stft, zx = stft(x_normal, fs, nperseg = len(x)/2, noverlap=len(x)/2-1,
                                        window = ('gaussian',14), nfft = NFFT)
    return f_stft, t_stft, zx

def localf_stft(t, x, pad = 2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_normal = x-np.average(x)
    f_stft, t_stft, zxx = stft(x_normal, fs, nperseg = len(x)/2, noverlap=len(x)/2-1,
                                        window = ('gaussian',14), nfft = NFFT)
    nrow, ncol = np.shape(zxx)
    ampl = abs(zxx)
    f_ave = np.zeros(ncol)
    for i in range(0, ncol):
        ampl_sum = sum(ampl[:,i])
        f_ave[i] = (sum(f_stft*ampl[:,i])/ampl_sum) 
    return t_stft, f_ave


def formatinput(x):
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


def spec_shift(t, x, t_shift):
    step = abs(t[1]-t[0])
    move = int(np.floor(t_shift/step))
    if t_shift > 0:
        x_new = np.concatenate((np.zeros(move), x[0:-move]), axis = 0)
    if t_shift < 0:
        x_new = np.concatenate((x[move:], np.zeros(move)))
    return x_new

def field_remove_nan(x):
    x = np.nan_to_num(x)
    for i,data in enumerate(x):
        if data == 0:
            if i+1<=len(x)-1:
                x[i] = x[i+1]
            else:
                x[i] = x[i-1]
    return x

def dict_normalize(x):
    x_in = deepcopy(x)
    x_out = {}
    if type(x) is dict:
        for key in list(x_in):
            x_max = max(x[key])
            x_normal = x[key]/x_max
            x_out[key] = x_normal
    else:
        raise Exception('Input is not a dictionary')
    return x_out
            
    
    


if __name__ == '__main__':
    path = 'C:/Users/Admin/Dropbox/Data/SPEC SO/2023/Jun/06-30-2023'
    all_fnames = get_fnames(path)

    x = {}
    y = {}
    t = {}
    zxx = {}
    f = {}
    t_stft = {}
    fnames_so = []
    f_ave_all = {}
    for file in all_fnames:
        if 'SO' in file:
            fnames_so.append(file)
    fnames_so.sort(key = str_get_D)
    for fname in fnames_so:
        t[fname], x[fname], y[fname] = load_single_spec(path, fname)