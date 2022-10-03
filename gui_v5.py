#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:24:23 2022

@author: ppsapphire
"""

    # -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:33:56 2021

@author: pps

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
import scipy as sp
import pandas as pd
import dataprocessor as dp
import study as st
import plottool_v2 as pt
import curfittool as cft
import get_superoscilation_v3  as gso


class Labtools:        
    def __init__(self):
        self.root = tk.Tk()
        # self.root.columnconfigure(0,weight=1)
        # self.root.rowconfigure(0,weight=1)
        self.nb = ttk.Notebook(self.root)
        self.frame_transmission = ttk.Frame(self.nb)
        self.frame_plotdata = ttk.Frame(self.nb)
        self.frame_index = ttk.Frame(self.nb)
        self.frame_loaddata = ttk.Frame(self.nb)
        self.frame_plotdata_fd = ttk.Frame(self.nb)
        self.frame_plotnormal = ttk.Frame(self.nb)
        self.frame_plotnormal_fd = ttk.Frame(self.nb)
        self.frame_polarimetry = ttk.Frame(self.nb)
        self.frame_getso = ttk.Frame(self.nb)
        self.nb.add(self.frame_loaddata,text='load data')
        self.nb.add(self.frame_transmission,text='transmission')
        self.nb.add(self.frame_index,text='refractive index')
        self.nb.add(self.frame_polarimetry,text='polarimetry')
        self.nb.add(self.frame_plotdata,text='plot data (td)')
        self.nb.add(self.frame_plotdata_fd,text='plot data (fd)')
        self.nb.add(self.frame_plotnormal,text='plot normalized (td)')
        self.nb.add(self.frame_plotnormal_fd,text='plot normalized (fd)')
        self.nb.add(self.frame_getso,text='construct superoscillation')
        self.nb.grid(sticky='nsew')
        # self.nb.columnconfigure(0,weight=1)
        # self.nb.rowconfigure(0,weight=1)
        
        '''
            load data frame
        '''
        self.expstyle_key = tk.StringVar(value='sam+ref')
        mainframe = ttk.Labelframe(self.frame_loaddata,text='load data')
        mainframe.pack(expand=True)
        # mainframe.columnconfigure(0,weight=1)
        # mainframe.columnconfigure(1,weight=1)
        typeentry = ttk.Combobox(mainframe,textvariable=self.expstyle_key,width=20)
        typeentry['values'] = ('sam/ref','sam+ref','polarimetry_so')
        self.nsam = tk.IntVar(value=1)
        self.nref = tk.IntVar(value=1)
        
        button_selectfolder = ttk.Button(mainframe,text='select folder',command=self.selectfolder)
        button1 = ttk.Button(mainframe,text='load data',command=self.load)
        button_avg = ttk.Button(mainframe,text='avg and fft',command=self.basicprocess)
        button_checkdata = ttk.Button(mainframe, text='check data', command=self.checkdata)
        label4 = ttk.Label(mainframe,text='Measurement type: ')
        label1 = ttk.Label(mainframe,text='Sample scan number: ')
        label2 = ttk.Label(mainframe,text='Reference scan number: ')
        self.entry_nsam = ttk.Entry(mainframe,textvariable=self.nsam,state='disabled',width=10)
        self.entry_nref = ttk.Entry(mainframe,textvariable=self.nref,state='disabled',width=10)
        button2 = ttk.Button(mainframe,text='confirm',command=self.state_entry)
   
        
        
        button_selectfolder.grid(column=0,row=0)
        label4.grid(column=0,row=1)
        typeentry.grid(column=1,row=1)
        button2.grid(column=2,row=1)
        self.entry_nsam.grid(column=1,row=2)
        self.entry_nref.grid(column=3,row=2)
        label1.grid(column=0,row=2)
        label2.grid(column=2,row=2)
        button1.grid(column=0,row=3)
        button_checkdata.grid(column=1,row=3)
        button_avg.grid(column=2,row=3)

        
        
        '''status frame'''
        self.actiondisp = tk.StringVar(value='ready')
        window_status = ttk.Labelframe(self.frame_loaddata,text='Status')
        window_status.pack(side=tk.BOTTOM,expand=True, fill=tk.BOTH)
        label_action = ttk.Label(window_status,textvariable=self.actiondisp)
        label_action.grid(column=1,row=6)
        
        
        '''
        study frame
        '''
        frame_study = ttk.Labelframe(self.frame_loaddata,text='studies')
        frame_study.pack(expand=True)
        # frame_study.columnconfigure(0,weight=1)
        button_transmission = ttk.Button(frame_study,text='study transmission',command=self.study_transmission)
        button_index = ttk.Button(frame_study,text='study refractive index',command=self.study_index)
        button_polarimetry = ttk.Button(frame_study,text='study polarimetry',command=self.study_polarimetry)
        button_superosc =  ttk.Button(frame_study,text='get_SO',command=self.SO_spectrum)
        
        
        button_transmission.grid(column=0,row=0,sticky=tk.W)
        button_index.grid(column=0,row=1)
        button_polarimetry.grid(column=0,row=2,sticky=tk.W)
        button_superosc.grid(column=0,row=3)
        
    def selectfolder(self):
        self.folder = tk.filedialog.askdirectory()
    
    def state_entry(self):
        if self.expstyle_key.get() == 'sam+ref' or self.expstyle_key.get() == 'polarimetry' :
            self.entry_nsam.configure(state='normal')
            self.entry_nref.configure(state='normal')
        else:
            self.entry_nsam.configure(state='disabled')
            self.entry_nref.configure(state='disabled')   
        #self.dispActions(self.frame_loaddata, self.expstyle_key.get()+' selected')
        self.actiondisp.set(self.expstyle_key.get()+' selected')
        
    
    def load(self):
        self.path = self.folder
        self.expstyle = self.expstyle_key.get()
        self.specgo = dp.SpecProcess(self.path,self.expstyle)
        self.x,self.y,self.t = self.specgo.loadspecs()
        self.actiondisp.set('file loaded')
    
    def basicprocess(self):
        if self.expstyle == 'sam+ref':
            [xavg_sam,tavg_sam,xavg_ref,tavg_ref,comp1,comp2] = self.specgo.average_samref(self.x, self.t,self.nsam.get(),self.nref.get())
            # self.xsam = self.specgo.Totalfield(xavg_sam, yavg_sam)
            # self.xref = self.specgo.Totalfield(xavg_ref, yavg_ref)
            self.xsam = xavg_sam
            self.xref = xavg_ref
            self.tsam = tavg_sam
            self.tref = tavg_ref
            self.xcomp = comp1
            self.ycomp = comp2
            self.freqsam,self.sxsam = dp.Basic_functions().fftx(self.xsam,self.tsam,2)
            self.freqref,self.sxref = dp.Basic_functions().fftx(self.xref,self.tref,2)
            self.freqxcomp,self.sxcomp = dp.Basic_functions().fftx(self.xcomp,self.tsam,2)
            self.freqycomp,self.sycomp = dp.Basic_functions().fftx(self.ycomp,self.tref,2)
            self.compall = self.specgo.combinedict(self.xcomp, self.ycomp)
            self.freqcompall = self.specgo.combinedict(self.freqxcomp, self.freqycomp)
            self.scompall = self.specgo.combinedict(self.sxcomp, self.sycomp)
            self.xall = self.specgo.combinedict(self.xsam, self.xref)
            self.tall = self.specgo.combinedict(self.tsam, self.tref)
            self.sxall = self.specgo.combinedict(self.sxsam, self.sxref)
            self.sxall_abs = dp.Basic_functions().dict_getabs(self.sxall)
            self.sxall_real = dp.Basic_functions().dict_getreal(self.sxall)
            self.freqall = self.specgo.combinedict(self.freqsam, self.freqref)
            self.xall_normal = dp.Basic_functions().dict_normalize(self.xall)
            self.sxall_abs_normal = dp.Basic_functions().dict_normalize(self.sxall_abs)
            self.plotgo_a = pt.Plottool(self.frame_plotdata,self.tall,self.xall)
            self.plotgo_b = pt.Plottool(self.frame_plotdata_fd, self.freqall, self.sxall)
            self.plotgo_c = pt.Plottool(self.frame_plotnormal, self.tall,self.xall_normal)
            self.plotgo_d = pt.Plottool(self.frame_plotnormal_fd, self.freqall, self.sxall_abs_normal)
            
        
        if self.expstyle == 'sam/ref':
            [xavg,yavg,tavg] = self.specgo.avespecs_sam(self.x,self.y,self.t)
            # self.xall = self.specgo.Totalfield(xavg,yavg)
            self.xall = xavg
            self.xall_smooth = self.specgo.dict_curve_smooth(self.xall)
            self.xall_normal = dp.Basic_functions().dict_normalize(self.xall)
            self.tall = tavg
            self.freqall,self.sxall = dp.Basic_functions().fftx(self.xall, self.tall, 2)
            self.sxall_abs = dp.Basic_functions().dict_getabs(self.sxall)
            self.sxall_intensity = dp.Basic_functions().dict_square(self.sxall_abs)
            self.sxall_abs_normal = dp.Basic_functions().dict_normalize(self.sxall_abs)
            self.sxall_real = dp.Basic_functions().dict_getreal(self.sxall)
            self.plotgo_1 =  pt.Plottool(self.frame_plotdata, self.tall,self.xall)
            self.plotgo_2 = pt.Plottool(self.frame_plotdata_fd, self.freqall, self.sxall_intensity)
            self.plotgo_3 = pt.Plottool(self.frame_plotnormal, self.tall,self.xall_normal)
            self.plotgo_4 = pt.Plottool(self.frame_plotnormal_fd, self.freqall, self.sxall_abs_normal)
            
        if self.expstyle == 'polarimetry_so':
            [self.xsam,self.ysam,self.tsam,self.xref,self.yref,self.tref,self.samc1,self.samc2,self.refc1,self.refc2] = self.specgo.average_polarimetry_so(self.t, self.x)
            self.freqsamx,self.sxsam = dp.Basic_functions().fftx(self.xsam, self.tsam, 2)
            self.freqrefx,self.sxref = dp.Basic_functions().fftx(self.xref, self.tref, 2)
            self.freqsamy,self.sysam = dp.Basic_functions().fftx(self.ysam, self.tsam, 2)
            self.sysam_abs = dp.Basic_functions().dict_getabs(self.sysam)
            self.freqrefy,self.syref = dp.Basic_functions().fftx(self.yref, self.tref, 2)
            self.syref_abs = dp.Basic_functions().dict_getabs(self.syref)
            # self.sEsam = dp.Basic_functions().fftx(self.total_sam, self.tsam, 2)
            # self.sEref = dp.Basic_functions().fftx(self.total_ref,self.tref,2)
            # self.Eall = self.specgo.combinedict(self.ysam, self.yref)
            self.Eall = self.specgo.combinedict(self.samc1, self.samc2)
            self.tall = self.specgo.combinedict(self.tsam, self.tref)
            self.plotgo_1 = pt.Plottool(self.frame_plotdata, self.tall, self.Eall)
            self.plotgo_2 = pt.Plottool(self.frame_plotdata_fd, self.freqsamy, self.sysam_abs)
            
        self.actiondisp.set('average and FFT done, process to plot to check result')    

    
    def checkdata(self):
        pt.Plottool(self.frame_plotdata,self.t,self.x)
        self.actiondisp.set('proceed to plot panel to check loaded data')
        
    '''
    begin transmission panel
    '''
        
    def study_transmission(self):
        if self.expstyle == 'sam+ref':
            self.transgo = st.study_transmission(self.frame_transmission,self.freqsam, self.sxsam, self.sxref)
        if self.expstyle == 'sam/ref':
            self.transgo = st.study_transmission(self.frame_transmission, self.freqsam, self.sxall, self.sxall)
        if self.expstyle == 'polarimetry_so':
            self.transgo = st.study_transmission(self.frame_transmission, self.freqsamy, self.sysam, self.syref)
        
    '''
    begin study index panel
    '''
    def study_index(self):
        self.indexgo = st.study_index(self.frame_index,self.xall,self.tall)
        self.actiondisp.set('Proceed to study index panel')
        
        
    '''
    study polarimetr
    '''    
    def study_polarimetry(self):
        self.polargo = st.study_polarimetry(self.frame_polarimetry, self.freqcompall, self.sxcomp, self.sycomp)
        self.actiondisp.set('Proceed to polarimetry study panel')
        
    def SO_spectrum(self):
        self.build_SO = gso.Display_SO(self.frame_getso)
        
# datago = Labtools()
# datago.root.mainloop()