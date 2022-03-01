# -*- coding: utf-8 -*-
"""
@author: Peisong Peng

This script is a module for dataprocessv2, which does most data process and analysis 
This script creates an interactive plotting UI and can do some basic data analysis such as curve fitting

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
from scipy.optimize import curve_fit
import scipy.signal as ss
              
class Plottool:
    def __init__(self,frame,xvalues,yvalues):
        '''
        

        Parameters
        ----------
        root : tk.Tk()
            DESCRIPTION: tk root to initiate gui base window
        xvalues : dict or list of 1d array or ndarray
            DESCRIPTION: x values to be plotted
        yvalues : TYPE
            DESCRIPTION: y values to be plotted

        Returns
        -------
        plot of desired axis limit, labels and fonts

        '''
        self.x = self.formatinput(xvalues)
        self.y = self.formatinput(yvalues)
        
        '''
        first window for plot editor
        '''
        self.frame = frame
        self.window = ttk.Labelframe(frame,text='plot editor')
        self.window.grid(column=0,row=0)
        self.names_y = list(self.y)
        self.names_x = list(self.x)
        self.listvar_x = tk.StringVar(value=self.names_x)
        self.listvar_y = tk.StringVar(value=self.names_y)
        
        self.labelx = tk.StringVar()
        self.labely = tk.StringVar()
        self.xlim_l = tk.DoubleVar()
        self.xlim_h = tk.DoubleVar()
        self.ylim_l = tk.DoubleVar()
        self.ylim_h = tk.DoubleVar()
        self.fontsize_lg = tk.DoubleVar(value=10)
        self.fontsize_ax = tk.DoubleVar(value=10)
        self.gridswitch = tk.IntVar(value=1)
        self.autoxlim = tk.IntVar(value=1)
        self.autoylim = tk.IntVar(value=1)
        self.legendon = tk.IntVar(value=1)
        self.action = tk.StringVar(value='ready')
        self.ydata_select = tk.Listbox(self.window,listvariable=self.listvar_y,selectmode='extended',width=50)
        self.xdata_select = tk.Listbox(self.window,listvariable=self.listvar_x,selectmode='extended',width=50)
        self.entry_xlabel = ttk.Entry(self.window,textvariable=self.labelx,width=10)
        self.entry_ylabel = ttk.Entry(self.window,textvariable=self.labely,width=10)
        self.entry_fontsize_lg = ttk.Entry(self.window,textvariable=self.fontsize_lg,width=5)
        self.entry_fontsize_ax = ttk.Entry(self.window,textvariable=self.fontsize_ax,width=5)
        entry_xl = ttk.Entry(self.window,textvariable=self.xlim_l,width=5)
        entry_xh = ttk.Entry(self.window,textvariable=self.xlim_h,width=5)
        entry_yl = ttk.Entry(self.window,textvariable=self.ylim_l,width=5)
        entry_yh = ttk.Entry(self.window,textvariable=self.ylim_h,width=5)
        self.checkbutton_grid = ttk.Checkbutton(self.window,variable=self.gridswitch)
        checkbutton_autoxlim = ttk.Checkbutton(self.window,variable=self.autoxlim)
        checkbutton_autoylim = ttk.Checkbutton(self.window,variable=self.autoylim)
        checkbutton_legendon = ttk.Checkbutton(self.window,variable=self.legendon)
        self.label_x = ttk.Label(self.window,text='x-axis name: ')
        self.label_y = ttk.Label(self.window,text='y-axis name: ')
        self.label_font_lg = ttk.Label(self.window,text='legend font size: ')
        self.label_font_ax = ttk.Label(self.window,text='axis font size: ')
        self.label_grid = ttk.Label(self.window,text='grid on: ')
        self.button_plot = ttk.Button(self.window,text='plot',command=self.plotfig)
        button_loadselection_x = ttk.Button(self.window,text='load x selection',command=self.loadselection_x)
        button_loadselection_y = ttk.Button(self.window,text='load y selection',command=self.loadselection_y)
        
        
        label_xlim1 = ttk.Label(self.window,text='x from: ')
        label_xlim2 = ttk.Label(self.window,text='to: ')
        label_ylim1 = ttk.Label(self.window,text='y from: ')
        label_ylim2 = ttk.Label(self.window,text='to: ')
        label_autoxlim = ttk.Label(self.window,text='auto')
        label_autoylim = ttk.Label(self.window,text='auto')
        
        label_legendon = ttk.Label(self.window,text='legend on')
        
        
        
        
        self.xdata_select.grid(column=0,row=0,columnspan=3)
        self.ydata_select.grid(column=3,row=0,columnspan=3)
        button_loadselection_x.grid(column=2,row=1)
        button_loadselection_y.grid(column=5,row=1)
        self.label_x.grid(column=0,row=2)
        self.entry_xlabel.grid(column=1,row=2)
        self.label_y.grid(column=2,row=2)
        self.entry_ylabel.grid(column=3,row=2)
        self.label_font_lg.grid(column=0,row=3)
        self.entry_fontsize_lg.grid(column=1,row=3)
        self.label_font_ax.grid(column=2,row=3)
        self.entry_fontsize_ax.grid(column=3,row=3)
        label_xlim1.grid(column=0,row=4)
        entry_xl.grid(column=1,row=4)  
        label_xlim2.grid(column=2,row=4)
        entry_xh.grid(column=3,row=4) 
        label_autoxlim.grid(column=4,row=4)
        checkbutton_autoxlim.grid(column=5,row=4)
        label_ylim1.grid(column=0,row=5)
        entry_yl.grid(column=1,row=5)
        label_ylim2.grid(column=2,row=5)
        entry_yh.grid(column=3,row=5)
        label_autoylim.grid(column=4,row=5)
        checkbutton_autoylim.grid(column=5,row=5)
        self.checkbutton_grid.grid(column=1,row=6)
        self.button_plot.grid(column=2,row=6)
        self.label_grid.grid(column=0,row=6)
        checkbutton_legendon.grid(column=4,row=6)
        label_legendon.grid(column=3,row=6)
        
        
        
        
        '''
        second window for showing td plots
        '''
        self.window2 = ttk.Labelframe(frame,text='plot')
        self.window2.grid(column=1,row=0,sticky='nsew')
        fig1 = Figure(figsize=(10,5))
        self.canvas0 = FigureCanvasTkAgg(figure=fig1,master=self.window2)
        self.canvas0.get_tk_widget().grid(column=0,row=0,columnspan=10,rowspan=10,sticky='nsew')
        # self.plotall(fig)
        
 
        
        
        '''
        third window for savin plots and corresponding data
        
        '''
        self.window3 = ttk.Labelframe(frame,text='save plot')
        self.window3.grid(column=0,row=2)
        self.legendnames = tk.StringVar()
        self.exportname_fig = tk.StringVar()
        self.exportname_data = tk.StringVar()
        self.userlgname = tk.IntVar(value=0)
        self.userlegendnames = []
        button_exportfig = ttk.Button(self.window3,text='Export Fig',command=self.plot_export)
        self.entry_legendnames = ttk.Entry(self.window3,textvariable=self.legendnames,width=40)
        self.entry_exportname_fig = ttk.Entry(self.window3,textvariable=self.exportname_fig,width=40)
        self.entry_exportname_data = ttk.Entry(self.window3,textvariable=self.exportname_data,width=40)
        label_explain = ttk.Label(self.window3,text='Enter file name including an additional part to original variable name and file extension name')
        label_exportname_fig = ttk.Label(self.window3,text='Enter export figure name: ')
        label_exportname_data = ttk.Label(self.window3,text='Enter export data name: ')
        button_exportdata = ttk.Button(self.window3,text='Export Data',command=self.export_data)
        button_loadlegendnames = ttk.Button(self.window3,text='Load names',command=self.loadnames)
        checkbutton_userlegendname = ttk.Checkbutton(self.window3,variable=self.userlgname)
        label_userlgname = ttk.Label(self.window3,text='user defined legend names')
        checkbutton_userlegendname.grid(column=1,row=0)
        self.entry_legendnames.grid(column=2,row=0,columnspan=2)
        label_userlgname.grid(column=0,row=0)
        button_loadlegendnames.grid(column=4,row=0)
        label_explain.grid(column=0,row=1,columnspan=4)
        label_exportname_fig.grid(column=0,row=2)
        self.entry_exportname_fig.grid(column=1,row=2,columnspan=2)
        button_exportfig.grid(column=3,row=2)
        label_exportname_data.grid(column=0,row=3)
        self.entry_exportname_data.grid(column=1,row=3,columnspan=2)
        button_exportdata.grid(column=3,row=3)
        
        
        
        
        
        '''
        4th window for quantify superoscillation
        '''
        self.window4 = ttk.Labelframe(frame,text='superoscillation quantify')
        self.window4.grid(column=1,row=2)
        self.sonum = tk.DoubleVar(value=0)
        self.holdon = tk.IntVar(value=0)
        self.sqkhb_ax = np.zeros(10)
        label_quantifier = ttk.Label(self.window4,text='Superoscillation Quantifier: ')
        label_sonum = ttk.Label(self.window4,textvariable=self.sonum)
        Button_SO = ttk.Button(self.window4,text='show SO quantify number',command=self.SO_quantify)
        Button_sqk = ttk.Button(self.window4,text='calculate SO quantify k',command=self.SO_quantify_k)
        Button_sqk_hb = ttk.Button(self.window4,text='hb method quantify k',command=self.SO_quaSO_quantify_k_hb)
        Button_sqk_show = ttk.Button(self.window4,text='show SO quantify k plot',command=self.SO_quantify_plot)
        Button_dt = ttk.Button(self.window4,text='show dt',command=self.show_derivative)
        Button_cal = ttk.Button(self.window4,text='calculate dt',command=self.cal_derivative)
        label_quantifier.grid(column=0,row=1,sticky='w')
        label_sonum.grid(column=1,row=1,sticky='w')
        Button_SO.grid(column=0,row=0,sticky='w')
        Button_dt.grid(column=0,row=3,sticky='w')
        Button_cal.grid(column=0,row=2,sticky='w')
        Button_sqk.grid(column=0,row=4,sticky='w')
        Button_sqk_hb.grid(column=1,row=4,sticky='w')
        Button_sqk_show.grid(column=0,row=5,sticky='w')
        
        # button_savespec = ttk.Button(self.window3,text='save this spectrum',command=self.savespec)
        
        '''
        5th window for curve fit
        '''
        
        self.window5 = ttk.LabelFrame(frame,text='curve fit options')
        self.window5.grid(column=1,row=3)
        self.fun_choice = tk.StringVar()
        self.p0 = tk.DoubleVar(value=1)
        self.p1 = tk.DoubleVar(value=1)
        self.p2 = tk.DoubleVar(value=1)
        self.p3 = tk.DoubleVar(value=1)
        self.p4 = tk.DoubleVar(value=1)
        self.p5 = tk.DoubleVar(value=1)
        self.p6 = tk.DoubleVar(value=1)
        self.p7 = tk.DoubleVar(value=1)
        self.y_fit_save = dict()
        self.para_fit_save = dict()
        self.label_fun = tk.StringVar(value='function format display')
        combox_functions = ttk.Combobox(self.window5,textvariable=self.fun_choice,width=8)
        combox_functions['value'] = ['linear','poly2','poly3','poly4','sinu','exponential','gaussian','modcos','cosine']
        button_fit = ttk.Button(self.window5,text='fit',command=self.fitcurve)
        button_savefit = ttk.Button(self.window5,text='save fit curve',command=self.savefit)
        self.entry_p0 = ttk.Entry(self.window5,textvariable=self.p0,state='disabled',width=10)
        self.entry_p1 = ttk.Entry(self.window5,textvariable=self.p1,state='disabled',width=10)
        self.entry_p2 = ttk.Entry(self.window5,textvariable=self.p2,state='disabled',width=10)
        self.entry_p3 = ttk.Entry(self.window5,textvariable=self.p3,state='disabled',width=10)
        self.entry_p4 = ttk.Entry(self.window5,textvariable=self.p4,state='disabled',width=10)
        self.entry_p5 = ttk.Entry(self.window5,textvariable=self.p5,state='disabled',width=10)
        self.entry_p6 = ttk.Entry(self.window5,textvariable=self.p6,state='disabled',width=10)
        self.entry_p7 = ttk.Entry(self.window5,textvariable=self.p7,state='disabled',width=10)
        label_p0 = ttk.Label(self.window5,text='a:')
        label_p1 = ttk.Label(self.window5,text='b:')
        label_p2 = ttk.Label(self.window5,text='c:')
        label_p3 = ttk.Label(self.window5,text='d:')
        label_p4 = ttk.Label(self.window5,text='e:')
        label_p5 = ttk.Label(self.window5,text='f:')
        label_p6 = ttk.Label(self.window5,text='g:')
        label_p7 = ttk.Label(self.window5,text='h:')
        label_parashow = ttk.Label(self.window5,text='fitted paremeters:')
        label_p0_showfit = ttk.Label(self.window5,text='a:')
        label_p1_showfit = ttk.Label(self.window5,text='b:')
        label_p2_showfit = ttk.Label(self.window5,text='c:')
        label_p3_showfit = ttk.Label(self.window5,text='d:')
        label_p4_showfit = ttk.Label(self.window5,text='e:')
        label_p5_showfit = ttk.Label(self.window5,text='f:')
        label_p6_showfit = ttk.Label(self.window5,text='g:')
        label_p7_showfit = ttk.Label(self.window5,text='h:')
        self.label_p0_fitvalue = ttk.Label(self.window5,text='none')
        self.label_p1_fitvalue = ttk.Label(self.window5,text='none')
        self.label_p2_fitvalue = ttk.Label(self.window5,text='none')
        self.label_p3_fitvalue = ttk.Label(self.window5,text='none')
        self.label_p4_fitvalue = ttk.Label(self.window5,text='none')
        self.label_p5_fitvalue = ttk.Label(self.window5,text='none')
        self.label_p6_fitvalue = ttk.Label(self.window5,text='none')
        self.label_p7_fitvalue = ttk.Label(self.window5,text='none')
        label_function = ttk.Label(self.window5,textvariable=self.label_fun)
        combox_functions.grid(column=0,row=0)
        button_fit.grid(column=1,row=0)
        label_p0.grid(column=0,row=1)
        
        self.entry_p0.grid(column=1,row=1)
        label_p1.grid(column=2,row=1)
        self.entry_p1.grid(column=3,row=1)
        label_p2.grid(column=4,row=1)
        self.entry_p2.grid(column=5,row=1)
        label_p3.grid(column=6,row=1)
        self.entry_p3.grid(column=7,row=1)
        label_p4.grid(column=0,row=2)
        self.entry_p4.grid(column=1,row=2)
        label_p5.grid(column=2,row=2)
        self.entry_p5.grid(column=3,row=2)
        label_p6.grid(column=4,row=2)
        self.entry_p6.grid(column=5,row=2)
        label_p7.grid(column=6,row=2)
        self.entry_p7.grid(column=7,row=2)
        label_function.grid(column=0,row=3)
        
        label_parashow.grid(column=0,row=4)
        label_p0_showfit.grid(column=0,row=5)
        self.label_p0_fitvalue.grid(column=1,row=5)
        label_p1_showfit.grid(column=2,row=5)
        self.label_p1_fitvalue.grid(column=3,row=5)
        label_p2_showfit.grid(column=4,row=5)
        self.label_p2_fitvalue.grid(column=5,row=5)
        label_p3_showfit.grid(column=6,row=5)
        self.label_p3_fitvalue.grid(column=7,row=5)
        label_p4_showfit.grid(column=0,row=6)
        self.label_p4_fitvalue.grid(column=1,row=6)
        label_p5_showfit.grid(column=2,row=6)
        self.label_p5_fitvalue.grid(column=3,row=6)
        label_p6_showfit.grid(column=4,row=6)
        self.label_p6_fitvalue.grid(column=5,row=6)
        label_p7_showfit.grid(column=6,row=6)
        self.label_p7_fitvalue.grid(column=7,row=6)
        button_savefit.grid(column=0,row=7)
        
        
        '''
        bottom window for status bar
        '''
        self.window6 = tk.Frame(frame)
        self.window6.grid(column=0,row=3)

        label_status = ttk.Label(self.window6,text='status: ')
        label_showaction = ttk.Label(self.window6,textvariable=self.action)
        
        label_status.grid(column=0,row=0)
        label_showaction.grid(column=1,row=0)
        
        
        
    def fitcurve(self):
        self.canvas0.figure = Figure(figsize=(10,5))
        plot0 = self.canvas0.figure.add_subplot()
        if self.autoxlim.get() == 1:
             self.x_new = self.x[self.names_x[self.idx_selectx[0]]]
             self.y_new = self.y[self.names_y[self.idx_selecty[0]]]
        else:
            idx_start = np.where(self.x[self.names_x[self.idx_selectx[0]]]>self.xlim_l.get())[0][0]
            idx_end = np.where(self.x[self.names_x[self.idx_selectx[0]]]<self.xlim_h.get())[0][-1]
            self.x_new = self.x[self.names_x[self.idx_selectx[0]]][idx_start:idx_end]
            self.y_new = self.y[self.names_y[self.idx_selecty[0]]][idx_start:idx_end]
        plot0.plot(self.x_new,self.y_new,label=self.names_y[self.idx_selecty[0]] ,linestyle='--')
        if self.fun_choice.get() == 'linear':
            self.para_fit = curve_fit(self.linear,self.x_new,self.y_new)[0]
            self.y_fit = self.para_fit[0]*self.x_new+self.para_fit[1]
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
        if self.fun_choice.get() == 'poly2': 
            self.para_fit = curve_fit(self.poly2,self.x_new,self.y_new)[0]
            self.y_fit = self.para_fit[0]*(self.x_new-self.para_fit[1])**2+self.para_fit[2]
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
        if self.fun_choice.get() == 'poly3': 
            self.para_fit = curve_fit(self.poly3,self.x_new,self.y_new)[0]
            self.y_fit = self.para_fit[0]*(self.x_new-self.para_fit[1])**3+self.para_fit[2]
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
        if self.fun_choice.get() == 'poly4': 
            self.entry_p0.configure(state='normal')
            self.entry_p1.configure(state='normal')
            self.entry_p2.configure(state='normal')
            self.entry_p3.configure(state='normal')
            self.entry_p4.configure(state='normal')
            self.label_fun.set('y=a*x**4+b*x**3+c*x**2+d*x**2+e')
            p0 = self.p0.get()
            p1 = self.p1.get()
            p2 = self.p2.get()
            p3 = self.p3.get()
            p4 = self.p4.get()
            p_init = np.array([p0,p1,p2,p3,p4])
            self.para_fit = curve_fit(self.poly4,self.x_new,self.y_new,p0=p_init)[0]
            self.y_fit = self.para_fit[0]*self.x_new**4+self.para_fit[1]*self.x_new**3+self.para_fit[2]*self.x_new**2+self.para_fit[3]*self.x_new+self.para_fit[4]
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
        if self.fun_choice.get() == 'sinu': 
            self.entry_p0.configure(state='normal')
            self.entry_p1.configure(state='normal')
            self.entry_p2.configure(state='normal')
            self.para_fit = curve_fit(self.sinu,self.x_new,self.y_new)[0]
            self.y_fit = self.sinu(self.x_new,self.para_fit[0],self.para_fit[1],self.para_fit[2])
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
        if self.fun_choice.get() == 'exponential': 
            self.entry_p0.configure(state='normal')
            self.entry_p1.configure(state='normal')
            self.entry_p2.configure(state='normal')
            self.label_fun.set('y=a*np.exp(b*x+c)')
            p0 = self.p0.get()
            p1 = self.p1.get()
            p2 = self.p2.get()
            p_init = np.array([p0,p1,p2])
            self.para_fit = curve_fit(self.exponential,self.x_new,self.y_new,p0=p_init)[0]
            self.y_fit = self.exponential(self.x_new,self.para_fit[0],self.para_fit[1],self.para_fit[2])
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
        if self.fun_choice.get() == 'gaussian': 
            self.entry_p0.configure(state='normal')
            self.entry_p1.configure(state='normal')
            p0 = self.p0.get()
            p1 = self.p1.get()
            p_init = np.array([p0,p1])
            self.label_fun.set('y=1/((np.sqrt(np.pi*2)*sigma))*np.exp(-(x-tau)^2/2/sigma^2)')
            self.para_fit = curve_fit(self.gaussian,self.x_new,self.y_new,p0=p_init)[0]
            self.label_p0_fitvalue.config(text=str(self.para_fit[0]))
            self.label_p1_fitvalue.config(text=str(self.para_fit[1]))
            self.y_fit = self.gaussian(self.x_new,self.para_fit[0],self.para_fit[1])
            plot0.plot(self.x_new,self.y_fit,label='fitted data')    
        if self.fun_choice.get() == 'modcos':
            self.entry_p0.configure(state='normal')
            self.entry_p1.configure(state='normal')
            self.entry_p2.configure(state='normal')
            self.entry_p3.configure(state='normal')
            self.entry_p4.configure(state='normal')
            self.entry_p5.configure(state='normal')
            self.entry_p6.configure(state='normal')
            self.label_fun.set('y=cos(2*pi*f*x+b)*(a*x**4+g*x**3+c*x**2+d*x+e)')
            p0 = self.p0.get()
            p1 = self.p1.get()
            p2 = self.p2.get()
            p3 = self.p3.get()
            p4 = self.p4.get()
            p5 = self.p5.get()
            p6 = self.p6.get()
            p_init = np.array([p0,p1,p2,p3,p4,p5,p6])
            self.para_fit = curve_fit(self.modcos,self.x_new,self.y_new,p0=p_init)[0]
            self.label_p0_fitvalue.config(text=str(self.para_fit[0]))
            self.label_p1_fitvalue.config(text=str(self.para_fit[1]))
            self.label_p2_fitvalue.config(text=str(self.para_fit[2]))
            self.label_p3_fitvalue.config(text=str(self.para_fit[3]))
            self.label_p4_fitvalue.config(text=str(self.para_fit[4]))
            self.label_p5_fitvalue.config(text=str(self.para_fit[5]))
            self.label_p6_fitvalue.config(text=str(self.para_fit[6]))
            self.y_fit = self.modcos(self.x_new,self.para_fit[0],self.para_fit[1],self.para_fit[2],self.para_fit[3],self.para_fit[4],self.para_fit[5],self.para_fit[6])
            plot0.plot(self.x_new,self.y_fit,label='fitted data')  
        if self.fun_choice.get() == 'cosine':
            self.entry_p0.configure(state='normal')
            self.entry_p1.configure(state='normal')
            self.entry_p2.configure(state='normal')
            self.entry_p3.configure(state='normal')
            self.label_fun.set('y=a*np.cos(b*x+c)+d')
            p0 = self.p0.get()
            p1 = self.p1.get()
            p2 = self.p2.get()
            p3 = self.p3.get()
            p_init = np.array([p0,p1,p2,p3])
            self.para_fit = curve_fit(self.cosine,self.x_new,self.y_new,p0=p_init)[0]
            self.y_fit = self.cosine(self.x_new,self.para_fit[0],self.para_fit[1],self.para_fit[2],self.para_fit[3])
            plot0.plot(self.x_new,self.y_fit,label='fitted cosine')    
        plot0.grid(self.gridswitch.get())
        plot0.legend(loc='best',fontsize=self.fontsize_lg.get())
        plot0.set_xlabel(self.labelx.get(),fontsize=self.fontsize_ax.get())
        plot0.set_ylabel(self.labely.get(),fontsize=self.fontsize_ax.get())
        self.canvas0.draw()
        
    def savefit(self):
        self.y_fit_save[self.names_x[self.idx_selectx[0]]] = self.y_fit
        self.para_fit_save[self.names_x[self.idx_selectx[0]]] = self.para_fit
        self.action.set('fit saved')

    
    def loadselection_x(self):
        self.idx_selectx = self.xdata_select.curselection()
        self.action.set('x data selected')
        
    def loadselection_y(self):
        self.idx_selecty = self.ydata_select.curselection()
        self.action.set('y data selected')
    
    def loadnames(self):
        self.userlegendnames = self.legendnames.get().split('+')
        self.action.set('legend names loaded')
        
        
    def plotfig(self):
        if self.holdon.get() == 0:
            self.canvas0.figure = Figure(figsize=(10,5))            
        plot0 = self.canvas0.figure.add_subplot()
        if self.userlgname.get() == 1:
            for i in self.idx_selectx:
                for j in self.idx_selecty:
                    if len(self.idx_selectx) == 1:
                        plot0.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.userlegendnames[j])
                    else:
                        if i == j:
                            plot0.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.userlegendnames[j])
        else:
            for i in self.idx_selectx:
                for j in self.idx_selecty:
                    if len(self.idx_selectx) == 1:
                        plot0.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.names_y[j])
                    else:
                        if i == j:
                            plot0.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.names_y[j])
        if self.autoxlim.get() == 0:
            plot0.set_xlim(left=float(self.xlim_l.get()),right=float(self.xlim_h.get()))
        if self.autoylim.get() == 0:
            plot0.set_ylim(bottom=float(self.ylim_l.get()),top=float(self.ylim_h.get()))
        plot0.set_xlabel(self.labelx.get(),fontsize=self.fontsize_ax.get())
        plot0.set_ylabel(self.labely.get(),fontsize=self.fontsize_ax.get())
        if self.legendon.get() == 1:
            plot0.legend(loc='upper right',fontsize=self.fontsize_lg.get())
        plot0.grid(self.gridswitch.get())
        self.canvas0.draw()
        self.action.set('Figure plotted')
      
        
    def plot_export(self):
        path = tk.filedialog.askdirectory()
        plt.figure()
        if self.userlgname.get() == 1:
            for i in self.idx_selectx:
                for j in self.idx_selecty:
                    if len(self.idx_selectx) == 1:
                        plt.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.userlegendnames[j])
                    else:
                        if i == j:
                            plt.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.userlegendnames[j])
        else:
            for i in self.idx_selectx:
                for j in self.idx_selecty:
                    if len(self.idx_selectx) == 1:
                        plt.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.names_y[j])
                    else:
                        if i == j:
                            plt.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.names_y[j])
        if self.legendon.get() == 1:             
            plt.legend(loc='upper right',fontsize=self.fontsize_lg.get())
        plt.xlabel(self.labelx.get(),fontsize=self.fontsize_ax.get())
        plt.ylabel(self.labely.get(),fontsize=self.fontsize_ax.get())
        if self.autoxlim.get() == 0:
            plt.xlim([self.xlim_l.get(),self.xlim_h.get()])
        if self.autoylim.get() == 0:
            plt.ylim([self.ylim_l.get(),self.ylim_h.get()])
        plt.grid(self.gridswitch.get())
        plt.savefig(path+'\\'+self.exportname.get())
        plt.show()
        self.action.set('Figure exported')
    
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
        
        if type(x) is pd.core.frame.DataFrame:
            x_new = x
        
        return x_new
    
    def save_spec(self):
        for i in self.idx_selectx:
            for j in self.idx_selecty:
                if len(self.idx_selectx) != 1:
                    if i==j:
                        newname_x = self.names_x[i] + self.specname.get()
                        newname_y = self.names_y[j] + self.specname.get()
                        self.x[newname_x],self.y[newname_y] = dp.Basic_functions().specchop(self.x[self.names_x[i]],self.y[self.names_y[j]], self.xlim_l.get(), self.xlim_h.get())
                        
                else:
                    newname_x = self.names_x[i] + self.specname.get()
                    newname_y = self.names_y[j] + self.specname.get()
                    self.x[newname_x],self.y[newname_y] = dp.Basic_functions().specchop(self.x[self.names_x[i]],self.y[self.names_y[j]], self.xlim_l.get(), self.xlim_h.get())
        self.action.set('spectrum saved')
        
    def export_data(self):
        path = tk.filedialog.askdirectory()
        for i in self.idx_selectx:
            for j in self.idx_selecty:
                if len(self.idx_selectx) != 1:
                    if i==j:
                        newname_x = self.names_x[i] +'_x'+ self.exportname_data.get()
                        newname_y = self.names_y[j] + '_y'+self.exportname_data.get()
                        if self.autoxlim.get() == 0:
                            self.x[newname_x],self.y[newname_y] = dp.Basic_functions().specchop(self.x[self.names_x[i]],self.y[self.names_y[j]], self.xlim_l.get(), self.xlim_h.get())
                        else:
                            self.x[newname_x] = self.x[self.names_x[i]]
                            self.y[newname_y] = self.y[self.names_y[j]]
                        data_output = np.stack((self.x[newname_x],self.y[newname_y]),axis=1)
                        datasave = open(path+'/'+'export_'+newname_x,'w')
                        np.savetxt(path+'/'+'export_'+newname_x,data_output)
                        datasave.close()
                        del newname_x,newname_y,datasave,data_output
                else:
                    newname_x = self.names_x[i] +'_x' + self.exportname_data.get()
                    newname_y = self.names_y[j] +'_y' + self.exportname_data.get()
                    if self.autoxlim.get() == 0:
                        self.x[newname_x],self.y[newname_y] = dp.Basic_functions().specchop(self.x[self.names_x[i]],self.y[self.names_y[j]], self.xlim_l.get(), self.xlim_h.get())
                    else:
                        self.x[newname_x] = self.x[self.names_x[i]]
                        self.y[newname_y] = self.y[self.names_y[j]]
                    data_output = np.stack((self.x[newname_x],self.y[newname_y]),axis=1)
                    datasave = open(path+'/'+'export_'+newname_x,'w')
                    np.savetxt(path+'/'+'export_'+newname_x,data_output)
                    datasave.close()
                    del newname_x,newname_y,datasave,data_output
        self.action.set('data exported') 
        
        
        
    def plotall(self,figure):
        plot = figure.add_subplot()
        for i in range(0,len(self.names_x)):
            plot.plot(self.x[self.names_x[i]],self.y[self.names_y[i]],label=self.names_x)
        self.canvas0.draw()   
        
    def linear(self,x,a,b):
        return a*x+b
    
    def poly2(self,x,a,b,c):
        return a*(x-b)**2+c
    
    def poly3(self,x,a,b,c):
        return a*(x-b)**3+c
    
    def poly4(self,x,a,b,c,d,e):
        return a*(x)**4+b*x**3+c*x**2+d*x**2+e
    
    def sinu(self,x,a,b,c):
        return a*np.sin(b*x+c)
    
    def exponential(self,x,a,b,c):
        return a*np.exp(b*x+c)
    
    def gaussian(self,x,sigma,tau):
        return 1/((np.sqrt(np.pi*2)*sigma))*np.exp(-(x-tau)**2/2/sigma**2)
    
    def modcos(self,x,a,b,c,d,e,f,g):
        return np.cos(2*np.pi*f*x+b)*(a*x**4+g*x**3+c*x**2+d*x+e)

    def modsin(self,x,a,b,c,d,e,f,g):
        return np.sin(2*np.pi*f*x+b)*(a*x**4+g*x**3+c*x**2+d*x+e)
    
    def cosine(self,x,a,b,c,d):
        return a*np.cos(b*x+c)+d
    
    
    def cal_derivative(self):
        self.deri_value = dict()
        for i in self.idx_selectx:
            for j in self.idx_selecty:
                if self.names_x[i] == self.names_y[j]:
                    self.de_time = self.x[self.names_x[i]]
                    self.de_ampl = self.y[self.names_y[j]]
                    self.deri_value[self.names_x[i]] = dp.Basic_functions().array_getderivative(self.de_time, self.de_ampl)
        self.action.set('derivative calculated')
        return self.deri_value
    
    
    def show_derivative(self):
        self.canvas0.figure = Figure() 
        plot1 = self.canvas0.figure.add_subplot()
        for i in self.idx_selectx:
            for j in self.idx_selecty:
                if self.names_x[i] == self.names_y[j]:
                    plot1.plot(self.x[self.names_x[i]],self.y[self.names_y[j]],label=self.names_x[i],linestyle='--')
                    plot1.plot(self.x[self.names_x[i]],self.deri_value[self.names_x[i]],label = self.names_x[i]+'dt')
        self.canvas0.draw()
        
    
    def get_timewindows(self):
        self.timewindow = dict()
        for name in list(self.deri_value):
            self.timewindow[name] = list()
            self.tw_x = self.deri_value[name]
            self.tw_t = self.x[name]
            for i in range(0,len(self.tw_x)-1):
                tw_x1 = self.tw_x[i]
                tw_x2 = self.tw_x[i+1]
                if tw_x1<0 and tw_x2 > 0:
                    self.timewindow[name].append((tw_x1+tw_x2)/2)
                
                
    
    def SO_quantify(self):
        for i in self.idx_selectx:
            for j in self.idx_selecty:
                if self.names_x[i] == self.names_y[j]:
                    self.y[self.names_y[j]] = self.y[self.names_y[j]]/max(self.y[self.names_y[j]])
                    idx_thresh_high = np.where(self.y[self.names_y[j]]>0.2)[0][-1]
                    idx_thresh_low = np.where(self.y[self.names_y[j]]>0.2)[0][0]
                    f_thresh = np.where(self.x[self.names_x[i]]>=0.8)[0][0]
                    # if self.autoxlim.get() == 0:
                    #     idx_thresh_high = np.where(self.x[self.names_x[i]]<=self.xlim_h.get())[0][-1]
                    # else:
                    #     idx_thresh_high = len(self.x[self.names_x[i]])-1
                    self.field_all = np.trapz(self.y[self.names_y[j]][idx_thresh_low:idx_thresh_high],self.x[self.names_x[i]][idx_thresh_low:idx_thresh_high])
                    self.field_so = np.trapz(self.y[self.names_y[j]][f_thresh:idx_thresh_high],self.x[self.names_x[i]][f_thresh:idx_thresh_high])
                    self.so_quantifier = self.field_so/self.field_all
        self.sonum.set(self.so_quantifier)
    
    def SO_quantify_k_hb(self):
        for i in self.idx_selectx:
            for j in self.idx_selecty:
                if self.names_x[i] == self.names_y[j]:
                    self.sqkhb_ax = (ss.hilbert(self.y[self.names_y[j]])) #hilbert tranform to make the signal complex
                    self.sqkhb_phase = np.unwrap(np.angle(self.sqk_ax),discont=np.pi/2) # calculate phase from the complex signal
                    self.sqkhb_time = self.x[self.names_x[i]]
                    self.sqkhb_k_local = dp.Basic_functions().array_getderivative(self.sqk_time, self.sqk_phase) # get derivative of phase for local k
                    # self.sqk_k_local_normal = self.sqk_k_local/max(self.sqk_k_local)
        self.action.set('hb quantify complete')
        
    def SO_quantify_k(self):
       self.ay = 0+1j*0
       for key in list(self.para_fit_save):
           self.sqk_time = self.x[key] 
           ay_real = self.modcos(self.sqk_time, self.para_fit_save[key][0], self.para_fit_save[key][1], self.para_fit_save[key][2], self.para_fit_save[key][3], self.para_fit_save[key][4],self.para_fit_save[key][5],self.para_fit_save[key][6])
           ay_img = self.modsin(self.sqk_time, self.para_fit_save[key][0], self.para_fit_save[key][1], self.para_fit_save[key][2], self.para_fit_save[key][3], self.para_fit_save[key][4],self.para_fit_save[key][5],self.para_fit_save[key][6])
           self.ay = self.ay + ay_real+1j*ay_img
       self.sqk_phase = np.unwrap(np.angle(self.ay),discont=np.pi/2)
       self.sqk_k_local = dp.Basic_functions().array_getderivative(self.sqk_time, self.sqk_phase) # get derivative of phase for local k
       self.action.set('quantify complete')
                   
                 
                  
    # def SO_quantify_plot(self):
    #     self.canvas0.figure,self.ax = plt.subplots(figsize=(10,5)) 
    #     self.ax2 = self.ax.twinx()
    #     for i in self.idx_selectx:
    #         for j in self.idx_selecty:
    #             if self.names_x[i] == self.names_y[j]:
    #                 if self.x[self.names_x[i]][-1] < 4.8:
    #                     idx_up = len(self.x[self.names_x[i]])-1
    #                 else:
    #                     idx_up = np.where(self.x[self.names_x[i]]>4.8)[0][0]
    #                 if self.x[self.names_x[i]][0] > -4.8:   
    #                     idx_low = 0
    #                 else:
    #                     idx_low = np.where(self.x[self.names_x[i]]<-4.8)[0][-1]
    #                 self.sqp_x = self.x[self.names_x[i]][idx_low:idx_up]
    #                 self.sqp_y1 = self.y[self.names_y[j]][idx_low:idx_up]
    #                 self.sqp_y2 = self.sqk_k_local[idx_low:idx_up]
    #                 self.sqp_y3 = self.sqk_phase[idx_low:idx_up]
    #                 self.ax2.plot(self.sqp_x,self.sqp_y1,label = 'signal',linestyle='--')
    #                 self.ax.plot(self.sqp_x,self.sqp_y2,label = 'k value',color='red')
    #                 self.ax.plot(self.sqp_x,self.sqp_y3, label='phase',linestyle='dashdot',color='green')
    #     self.ax.legend(loc='upper left',fontsize=self.fontsize_lg.get())
    #     self.ax2.legend(loc='upper right',fontsize=self.fontsize_lg.get())
    #     self.ax.grid(self.gridswitch.get())
    #     self.ax.set_xlabel('time (ps)',fontsize=self.fontsize_ax.get())
    #     self.ax.set_ylabel('phase and k value',fontsize=self.fontsize_ax.get())
    #     self.ax2.set_ylabel('signal',fontsize=self.fontsize_ax.get())
    #     self.canvas0.draw()
        
    def SO_quantify_plot(self):
        self.canvas0.figure,self.ax = plt.subplots(figsize=(10,5)) 
        self.ax2 = self.ax.twinx()
        self.sqp_x = self.sqk_time
        self.sqp_y1 = np.real(self.ay)
        self.sqp_y2 = self.sqk_k_local
        self.sqp_y3 = self.sqk_phase
        if sum(self.sqkhb_ax) != 0:
            self.sqp_y4 = self.sqkhb_phase
            self.sqp_y5 = self.sqkhb_k_local
            self.ax2.plot(self.sqp_x,self.sqp_y4,label='hb method phase',linestyle='--',color='yellow')
            self.ax.plot(self.sqp_x,self.sqp_y5,label='hb method local k',color='orange')
        self.ax2.plot(self.sqp_x,self.sqp_y1,label = 'signal',linestyle='--')
        self.ax.plot(self.sqp_x,self.sqp_y2,label = 'local k',color='red')
        self.ax.plot(self.sqp_x,self.sqp_y3, label='phase',linestyle='dashdot',color='green')
        if self.autoxlim.get() == 0:
            self.ax.set_xlim(self.xlim_l.get(),self.xlim_h.get())
        self.ax.legend(loc='upper left',fontsize=self.fontsize_lg.get())
        self.ax2.legend(loc='upper right',fontsize=self.fontsize_lg.get())
        self.ax.grid(self.gridswitch.get())
        self.ax.set_xlabel('time (ps)',fontsize=self.fontsize_ax.get())
        self.ax.set_ylabel('phase and k value',fontsize=self.fontsize_ax.get())
        self.ax2.set_ylabel('signal',fontsize=self.fontsize_ax.get())
        self.canvas0.draw()
        
        
        
# x = np.linspace(1, 10,10)
# y = x**2
# root=tk.Tk()
# window_master=tk.Frame()
# window_master.grid()
# pt=Plottool(window_master, x, y)
# root.mainloop()