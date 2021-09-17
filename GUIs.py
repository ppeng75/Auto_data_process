# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 10:00:49 2021

@author: pps
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
import scipy as sp

class Plot_GUI:  
    def __init__(self,root,x_values,y_values):      
        self.root = root
        self.root.title('Plot')
        self.xvalues = x_values
        self.yvalues = y_values
        self.window0 = ttk.Frame(self.root)
        self.window0.grid()
        self.xlim_l = tk.DoubleVar()
        self.xlim_h = tk.DoubleVar()
        self.ylim_l = tk.DoubleVar()
        self.ylim_h = tk.DoubleVar()
        self.xvars = tk.StringVar()
        self.yvars = tk.StringVar()
        self.setlimit_x = tk.IntVar()
        self.setlimit_y = tk.IntVar()
        self.xidx_list = tk.IntVar(value=999)
        self.yidx_list = tk.IntVar(value=999)
        self.xidx_array = tk.StringVar(value=None)
        self.yidx_array = tk.StringVar(value=None)
        checkbox_x = ttk.Checkbutton(self.window0,text='auto',variable=self.setlimit_x,onvalue=0,offvalue=1)
        checkbox_y = ttk.Checkbutton(self.window0,text='auto',variable=self.setlimit_y,onvalue=0,offvalue=1)
        self.combox_xd = ttk.Combobox(self.window0,textvariable=self.xvars,state='disabled',width=20)
        self.combox_yd = ttk.Combobox(self.window0,textvariable=self.yvars,state='disabled',width=20)
        self.combox_xl = ttk.Combobox(self.window0,textvariable=self.xidx_list,state='disabled',width=10)
        self.combox_yl = ttk.Combobox(self.window0,textvariable=self.yidx_list,state='disabled',width=10)
        self.combox_xa = ttk.Combobox(self.window0,textvariable=self.xidx_array,state='disabled',width=10)
        self.combox_ya = ttk.Combobox(self.window0,textvariable=self.yidx_array,state='disabled',width=10)
        entry_xl = ttk.Entry(self.window0,textvariable=self.xlim_l,width=10)
        entry_xh = ttk.Entry(self.window0,textvariable=self.xlim_h,width=10)
        entry_yl = ttk.Entry(self.window0,textvariable=self.ylim_l,width=10)
        entry_yh = ttk.Entry(self.window0,textvariable=self.ylim_h,width=10)
        labelx1 = ttk.Label(self.window0,text='x from ')
        labelx2 = ttk.Label(self.window0,text='to ')
        labely1 = ttk.Label(self.window0,text='y from ')
        labely2 = ttk.Label(self.window0,text='to ')
        label_cbox_xd = ttk.Label(self.window0,text='x dict: ')
        label_cbox_yd = ttk.Label(self.window0,text='y dict: ')
        label_cbox_xl = ttk.Label(self.window0,text='x list: ')
        label_cbox_yl = ttk.Label(self.window0,text='y list: ')
        label_cbox_xa = ttk.Label(self.window0,text='x array: ')
        label_cbox_ya = ttk.Label(self.window0,text='y array: ')
        button1 = ttk.Button(self.window0,text='Plot',command=self.plotfig)
        self.button_xd = ttk.Button(self.window0,text='get content',command=self.getdictcontent_x,state='disabled')
        self.button_yd = ttk.Button(self.window0,text='get content',command=self.getdictcontent_y,state='disabled')
        self.button_xl = ttk.Button(self.window0,text='get content',command=self.getlistcontent_x,state='disabled')
        self.button_yl = ttk.Button(self.window0,text='get content',command=self.getlistcontent_y,state='disabled')
        self.button_xa = ttk.Button(self.window0,text='get data',command=self.getarraycontent_x,state='disabled')
        self.button_ya = ttk.Button(self.window0,text='get data',command=self.getarraycontent_y,state='disabled')
        #plot controls
        labelx1.grid(column=0,row=0)
        entry_xl.grid(column=1,row=0)
        labelx2.grid(column=2,row=0)
        entry_xh.grid(column=3,row=0)
        checkbox_x.grid(column=4,row=0)
        
        labely1.grid(column=0,row=1)
        entry_yl.grid(column=1,row=1)
        labely2.grid(column=2,row=1)
        entry_yh.grid(column=3,row=1)
        checkbox_y.grid(column=4,row=1)
        
        label_cbox_xd.grid(column=0,row=2)
        self.combox_xd.grid(column=1,row=2)
        self.button_xd.grid(column=2,row=2)
        label_cbox_yd.grid(column=3,row=2)
        self.combox_yd.grid(column=4,row=2)
        self.button_yd.grid(column=5,row=2)
        
        label_cbox_xl.grid(column=0,row=3)
        self.combox_xl.grid(column=1,row=3)
        self.button_xl.grid(column=2,row=3)
        label_cbox_yl.grid(column=3,row=3)
        self.combox_yl.grid(column=4,row=3)
        self.button_yl.grid(column=5,row=3)
        
        label_cbox_xa.grid(column=0,row=4)
        self.combox_xa.grid(column=1,row=4)
        self.button_xa.grid(column=2,row=4)
        label_cbox_ya.grid(column=3,row=4)
        self.combox_ya.grid(column=4,row=4)
        self.button_ya.grid(column=5,row=4)
        
        button1.grid(column=0,row=5)
        
        
        if type(self.xvalues) == dict:
            self.button_xd.configure(state='normal')
            self.combox_xd.configure(state='normal')
            self.combox_xd['values'] = list(self.xvalues)
        if type(self.yvalues) == dict:
            self.button_yd.configure(state='normal')
            self.combox_yd.configure(state='normal')
            self.combox_yd['values'] = list(self.xvalues)
        if type(self.xvalues) == list:
            self.button_xl.configure(state='normal')
            self.combox_xl.configure(state='normal')
            self.combox_xl['values'] = np.arange(0,len(self.xvalues),1)
        if type(self.yvalues) == list:
            self.button_yl.configure(state='normal')
            self.combox_yl.configure(state='normal')
            self.combox_yl['values'] = np.arange(0,len(self.xvalues),1)
        if type(self.xvalues) == np.ndarray:
            self.button_xa.configure(state='normal')
            self.combox_xa.configure(state='normal')
            self.combox_xa['values'] = np.arange(0,len(self.xvalues),1)
        if type(self.yvalues) == np.ndarray:
            self.button_ya.configure(state='normal')
            self.combox_ya.configure(state='normal')
            self.combox_ya['values'] = np.arange(0,len(self.xvalues),1)
            
            
    def getdictcontent_x(self):
        content = self.xvalues[self.xvars.get()]
        if type(content) == list:
            content_idx = np.arange(0,len(content),1)
            self.combox_xl.configure(state='normal')
            self.combox_xl['value'] = list(content_idx)
        if type(content) == np.ndarray:
            if np.size(content) == len(content):
                self.combox_xa.configure(state='normal')
                self.button_xa.configure(state='normal')
                self.combox_xa['value'] = '0'
            else:
                content_idx = np.arange(0,np.size(content,axis=1),1)
                self.combox_xl.configure(state='normal')
                self.button_xl.configure(state='normal')
                self.combox_xl['value'] = list(content_idx)
            
    def getdictcontent_y(self):
        content = self.yvalues[self.yvars.get()]
        if type(content) is list:
            content_idx = np.arange(0,len(content),1)
            self.combox_yl.configure(state='normal')
            self.combox_yl['value'] = list(content_idx)
        if type(content) is np.ndarray:
            if np.size(content) == len(content):
                self.combox_ya.configure(state='normal')
                self.button_ya.configure(state='normal')
                self.combox_ya['value'] = '0'
            else:
                content_idx = np.arange(0,np.size(content,axis=1),1)
                self.combox_yl.configure(state='normal')
                self.button_yl.configure(state='normal')
                self.combox_yl['value'] = list(content_idx)
            
    def getlistcontent_x(self):
        if type(self.xvalues) == dict:
            content = self.xvalues[self.xvars.get()][self.xidx_list.get()]
            if np.size(content) != len(content):
                content_idx = []
                for i in range(0,np.size(content,axis=1)):
                    content_idx.append(str(i))
                self.combox_xa.configure(state='normal')
                self.button_xa.configure(state='normal')
                self.combox_xa['value'] = content_idx
            else:
                self.combox_xa.configure(state='normal')
                self.button_xa.configure(state='normal')
                self.combox_xa['value'] = '0'
        else:
            content = self.xvalues[self.xidx_list.get()]
            if np.size(content) != len(content):
                content_idx = []
                for i in range(0,np.size(content,axis=1)):
                    content_idx.append(str(i))
                self.combox_xa.configure(state='normal')
                self.combox_xa['value'] = content_idx
            else:
                self.combox_xa.configure(state='normal')
                self.combox_xa['value'] = '0'
                
    def getlistcontent_y(self):
        if type(self.yvalues) == dict:
            content = self.yvalues[self.yvars.get()][self.yidx_list.get()]
            if np.size(content) != len(content):
                content_idy = []
                for i in range(0,np.size(content,axis=1)):
                    content_idy.append(str(i))
                self.combox_ya.configure(state='normal')
                self.button_ya.configure(state='normal')
                self.combox_ya['value'] = content_idy
            else:
                self.combox_ya.configure(state='normal')
                self.button_ya.configure(state='normal')
                self.combox_ya['value'] = '0'
        else:
            content = self.yvalues[self.yidx_list.get()]
            if np.size(content) != len(content):
                content_idy = []
                for i in range(0,np.size(content,axis=1)):
                    content_idy.append(str(i))
                self.combox_ya.configure(state='normal')
                self.button_ya.configure(state='normal')
                self.combox_ya['value'] = content_idy
            else:
                self.combox_ya.configure(state='normal')
                self.button_ya.configure(state='normal')
                self.combox_ya['value'] = '0'
            
    def getarraycontent_x(self):
        if type(self.xvalues) is dict:
            if self.xidx_list.get() != 999:
                if self.xidx_array.get() != '0':
                    self.xdata = self.xvalues[self.xvars.get()][self.xidx_list.get()][:,int(self.xidx_array.get())]
                else:
                    self.xdata = self.xvalues[self.xvars.get()][self.xidx_list.get()][:]
            else:
                if self.xidx_array.get() != '0':
                    self.xdata = self.xvalues[self.xvars.get()][:,int(self.xidx_array.get())]
                else:
                    self.xdata = self.xvalues[self.xvars.get()][:]
        else:
            if self.xidx_list.get() != 999:
                if self.xidx_array.get() != '0':
                    self.xdata = self.xvalues[self.xidx_list.get()][:,int(self.xidx_array.get())]
                else:
                    self.xdata = self.xvalues[self.xidx_list.get()][:]
            else:
                if self.xidx_array.get() != '0':
                    self.xdata = self.xvalues[:,int(self.xidx_array.get())]
                else:
                    self.xdata = self.xvalues[:]
            
      
    def getarraycontent_y(self):
        if type(self.yvalues) is dict:
            if self.yidx_list.get() != 999:
                if self.yidx_array.get() != '0':
                    self.ydata = self.yvalues[self.yvars.get()][self.yidx_list.get()][:,int(self.xidx_array.get())]
                else:
                    self.ydata = self.yvalues[self.yvars.get()][self.yidx_list.get()][:]
            else:
                if self.yidx_array.get() != '0':
                    self.ydata = self.yvalues[self.yvars.get()][:,int(self.yidx_array.get())]
                else:
                    self.ydata = self.yvalues[self.yvars.get()][:]
        else:
            if self.yidx_list.get() != 999:
                if self.yidx_array.get() != '0':
                    self.ydata = self.yvalues[self.yidx_list.get()][:,int(self.yidx_array.get())]
                else:
                    self.ydata = self.yvalues[self.yidx_list.get()][:]
            else:
                if self.yidx_array.get() != '0':
                    self.ydata = self.yvalues[:,int(self.yidx_array.get())]
                else:
                    self.ydata = self.yvalues[:]
            
    def plotfig(self):
        self.fig = Figure()
        plot0 = self.fig.add_subplot()
        plot0.plot(self.xdata,self.ydata)
        if self.setlimit_x.get() == 1:
            plot0.set_xlim(left=float(self.xlim_l.get()),right=float(self.xlim_h.get()))
        if self.setlimit_y.get() == 1:
            plot0.set_ylim(bottom=float(self.ylim_l.get()),top=float(self.ylim_h.get()))
        canvas0 = FigureCanvasTkAgg(figure=self.fig,master=self.window0)
        canvas0.draw()
        canvas0.get_tk_widget().grid(column=0,row=6,columnspan=5)       
            
            
class Polarimetry_GUI:
    def __init__(self,names):
        self.root = tk.Tk()
        self.filenames = names
        self.mainwindow = ttk.Frame(self.root)
        self.mainwindow.grid()
        self.namedisplay = str()
        self.idxchar_sam=tk.StringVar()
        self.idxchar_ref=tk.StringVar()
        for i in range(len(self.filenames)):
            self.namedisplay = self.namedisplay+str(i)+':'+self.filenames[i]+'\n '
        canvas1 = tk.Canvas(self.mainwindow)    
        canvas1.create_text(10,10,text=self.namedisplay,anchor=tk.NW)
        label_sam = ttk.Label(self.mainwindow,text='sample scans index')
        label_ref = ttk.Label(self.mainwindow,text='reference scans index')
        entry1 = ttk.Entry(self.mainwindow,textvariable=self.idxchar_sam,width=15)
        entry2 = ttk.Entry(self.mainwindow,textvariable=self.idxchar_ref,width=15)
        button1 = ttk.Button(self.mainwindow,text='next',command=self.preparedata)
        canvas1.grid(column=0,row=0,columnspan=4)
        label_sam.grid(column=0,row=1,columnspan=2)
        label_ref.grid(column=2,row=1,columnspan=2)
        entry1.grid(column=0,row=2,columnspan=2)
        entry2.grid(column=2,row=2,columnspan=2)
        button1.grid(column=3,row=3)
        
    def preparedata(self):
        self.idxint_sam = []
        self.idxint_ref = []
        for idx in self.idxchar_sam.get().split(','):
            self.idxint_sam.append(int(idx))
        for idx in self.idxchar_ref.get().split(','):
            self.idxint_ref.append(int(idx))
        listvar = tk.StringVar(value=self.idxint_ref)
        label1 = ttk.Label(self.mainwindow,text='Select the index of reference scans \n that will be kept:')
        self.listbox1 = tk.Listbox(self.mainwindow,listvariable=listvar,selectmode='extended')
        button1 = ttk.Button(self.mainwindow,text='save selection',command=self.saveselection)
        button2 = ttk.Button(self.mainwindow,text='next',command=self.root.destroy)
        self.listbox1.grid(column=0,row=5,columnspan=2)
        label1.grid(column=0,row=4,columnspan=2)
        button1.grid(column=2,row=5)
        button2.grid(column=4,row=6)
        
    def saveselection(self):
        self.idx_select = self.listbox1.curselection()

class Loaddata_GUI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title('Data Process')
        mainframe = ttk.Frame(self.root)
        # mainframe['padding'] = '5,5,5,5'
        mainframe['borderwidth'] = 2
        mainframe.grid()
        self.selectfolder = tk.filedialog.askdirectory()
        self.expstyle_key = tk.StringVar()
        typeentry = ttk.Combobox(mainframe,textvariable=self.expstyle_key,width=20)
        typeentry['values'] = ('sam/ref','sam+ref','polarimetry')
        self.nsam = tk.IntVar()
        self.nref = tk.IntVar()
        button1 = ttk.Button(mainframe,text='ok',command=self.root.destroy)
        label4 = ttk.Label(mainframe,text='Measurement type: ')
        label1 = ttk.Label(mainframe,text='Sample scan number: ')
        label2 = ttk.Label(mainframe,text='Reference scan number: ')
        self.entry_nsam = ttk.Entry(mainframe,textvariable=self.nsam,state='disabled',width=10)
        self.entry_nref = ttk.Entry(mainframe,textvariable=self.nref,state='disabled',width=10)
        button2 = ttk.Button(mainframe,text='confirm',command=self.state_entry)
        label4.grid(column=0,row=0)
        typeentry.grid(column=1,row=0)
        button2.grid(column=2,row=0)
        self.entry_nsam.grid(column=1,row=1)
        self.entry_nref.grid(column=3,row=1)
        label1.grid(column=0,row=1)
        label2.grid(column=2,row=1)
        button1.grid(column=0,row=2)
        
    def state_entry(self):
        if self.expstyle_key.get() == 'sam+ref' or self.expstyle_key.get() == 'polarimetry' :
            self.entry_nsam.configure(state='normal')
            self.entry_nref.configure(state='normal')
        else:
            self.entry_nsam.configure(state='disabled')
            self.entry_nref.configure(state='disabled')

class Editplot:
    def __init__(self,xvalues,yvalues):
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
        self.root = tk.Tk()
        self.window = ttk.Frame(self.root)
        self.window.grid()
        self.names = list(self.x)
        self.listvar = tk.StringVar(value=self.names)
        self.labelx = tk.StringVar()
        self.labely = tk.StringVar()
        self.xlim_l = tk.DoubleVar()
        self.xlim_h = tk.DoubleVar()
        self.ylim_l = tk.DoubleVar()
        self.ylim_h = tk.DoubleVar()
        self.fontsize_lg = tk.DoubleVar()
        self.fontsize_ax = tk.DoubleVar()
        self.gridswitch = tk.IntVar(value=0)
        self.autoxlim = tk.IntVar(value=1)
        self.autoylim = tk.IntVar(value=1)
        self.xdata_select = tk.Listbox(self.window,listvariable=self.listvar,selectmode='extended',width=100)
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
        self.label_x = ttk.Label(self.window,text='x-axis name: ')
        self.label_y = ttk.Label(self.window,text='y-axis name: ')
        self.label_font_lg = ttk.Label(self.window,text='legend font size: ')
        self.label_font_ax = ttk.Label(self.window,text='axis font size: ')
        self.label_grid = ttk.Label(self.window,text='grid on: ')
        self.button_plot = ttk.Button(self.window,text='plot',command=self.plotfig)
        button_loadselection = ttk.Button(self.window,text='load selection',command=self.loadselection)
        button_exportfig = ttk.Button(self.window,text='Export',command=self.plot_export)
        label_xlim1 = ttk.Label(self.window,text='x from: ')
        label_xlim2 = ttk.Label(self.window,text='to: ')
        label_ylim1 = ttk.Label(self.window,text='y from: ')
        label_ylim2 = ttk.Label(self.window,text='to: ')
        label_autoxlim = ttk.Label(self.window,text='auto')
        label_autoylim = ttk.Label(self.window,text='auto')
        self.xdata_select.grid(column=0,row=0,columnspan=5)
        button_loadselection.grid(column=5,row=0)
        self.label_x.grid(column=0,row=1)
        self.entry_xlabel.grid(column=1,row=1)
        self.label_y.grid(column=2,row=1)
        self.entry_ylabel.grid(column=3,row=1)
        self.label_font_lg.grid(column=0,row=2)
        self.entry_fontsize_lg.grid(column=1,row=2)
        self.label_font_ax.grid(column=2,row=2)
        self.entry_fontsize_ax.grid(column=3,row=2)
        label_xlim1.grid(column=0,row=3)
        entry_xl.grid(column=1,row=3)  
        label_xlim2.grid(column=2,row=3)
        entry_xh.grid(column=3,row=3) 
        label_autoxlim.grid(column=4,row=3)
        checkbutton_autoxlim.grid(column=5,row=3)
        label_ylim1.grid(column=0,row=4)
        entry_yl.grid(column=1,row=4)
        label_ylim2.grid(column=2,row=4)
        entry_yh.grid(column=3,row=4)
        label_autoylim.grid(column=4,row=4)
        checkbutton_autoylim.grid(column=5,row=4)
        self.checkbutton_grid.grid(column=1,row=5)
        self.button_plot.grid(column=2,row=5)
        self.label_grid.grid(column=0,row=5)
        button_exportfig.grid(column=0,row=7)
        
    def loadselection(self):
        self.idx_select = self.xdata_select.curselection()
        
    def plotfig(self):
        self.fig = Figure()
        plot0 = self.fig.add_subplot()
        for i in self.idx_select:
            plot0.plot(self.x[self.names[i]],self.y[self.names[i]],label=self.names[i])
        if self.autoxlim.get() == 0:
            plot0.set_xlim(left=float(self.xlim_l.get()),right=float(self.xlim_h.get()))
        if self.autoylim.get() == 0:
            plot0.set_ylim(bottom=float(self.ylim_l.get()),top=float(self.ylim_h.get()))
        plot0.set_xlabel(self.labelx.get(),fontsize=self.fontsize_ax.get())
        plot0.set_ylabel(self.labely.get(),fontsize=self.fontsize_ax.get())
        plot0.legend(loc='best',fontsize=self.fontsize_lg.get())
        plot0.grid(self.gridswitch)
        canvas0 = FigureCanvasTkAgg(figure=self.fig,master=self.window)
        canvas0.draw()
        canvas0.get_tk_widget().grid(column=0,row=6,columnspan=6)  
        
    def plot_export(self):
        plt.figure()
        for i in self.idx_select:
            plt.plot(self.x[self.names[i]],self.y[self.names[i]],label=self.names[i])
        plt.legend(loc='best',fontsize=self.fontsize_lg.get())
        plt.xlabel(self.labelx.get(),fontsize=self.fontsize_ax.get())
        plt.ylabel(self.labely.get(),fontsize=self.fontsize_ax.get())
        if self.autoxlim.get() == 0:
            plt.xlim([self.xlim_l.get(),self.xlim_h.get()])
        if self.autoylim.get() == 0:
            plt.ylim([self.ylim_l.get(),self.ylim_h.get()])
        plt.grid(self.gridswitch.get())
        plt.show()
        self.root.destroy()
    
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
    

# class GetIndex:
#     def __init__(self,x,y):
#         self.x = x
#         self.y = y
#         self.root = tk.Tk()
#         self.window = ttk.Frame(self.root)
#         self.window.grid()
#         self.names = list(x)
#         self.listnames = tk.StringVar(value=self.names)
        
#         listbox_sam = tk.Listbox(self.window,listvariable=self.listnames,selectmode='extended',width=100)
#         listbox_ref = tk.Listbox(self.window,listvariable=self.listnames,selectmode='extended',width=100)
        