import tkinter as tk
from tkinter import ttk
import dataprocessor as dp
import numpy as np
import plottool_v3 as pt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import copy
import matplotlib.pyplot as plt
import configparser as cfp


class study_transmission:
    def __init__(self,frame,freq,sxsam,sxref):
        self.frame_tr = frame
        self.freq = freq
        self.sxsam = sxsam
        self.sxref = sxref
        self.window1 = ttk.Labelframe(self.frame_tr,text='main')
        self.window1.grid(column=0,row=0)
        self.name_sam = list(sxsam)
        self.name_ref = list(sxref)
        self.list_sam = tk.StringVar(value=self.name_sam)
        self.list_ref = tk.StringVar(value=self.name_ref)
        self.action = tk.StringVar(value='ready')
        self.listbox_sam= tk.Listbox(self.window1,listvariable=self.list_sam,selectmode='extended')
        self.listbox_ref= tk.Listbox(self.window1,listvariable=self.list_ref,selectmode='extended')
        button_select_sam = ttk.Button(self.window1,text='use selection as sample',command=self.select_sam)
        button_select_ref = ttk.Button(self.window1,text='use selection as reference',command=self.select_ref)
        button_calculate = ttk.Button(self.window1,text='calculate',command=self.calculate)
        button_select_combine_ref = ttk.Button(self.window1,text='select and combine',command=self.combine_select_ref)
        button_checktrans = ttk.Button(self.window1,text='check result',command=self.checktrans)
        label_status = ttk.Label(self.window1,text='status: ')
        label_action = ttk.Label(self.window1,textvariable=self.action)
        
        
        self.listbox_sam.grid(column=0,row=0,columnspan=2)
        self.listbox_ref.grid(column=2,row=0,columnspan=2) 
        button_select_sam.grid(column=0,row=1)
        button_select_ref.grid(column=2,row=1)
        button_select_combine_ref.grid(column=3,row=1)
        button_calculate.grid(column=4,row=2)
        button_checktrans.grid(column=0,row=3)
        label_status.grid(column=0,row=4)
        label_action.grid(column=1,row=4)
        
        
        self.window2 = ttk.Labelframe(self.frame_tr,text='plot')
        self.window2.grid(column=0,row=1,sticky='nsew')
        
        
    def select_sam(self):
        self.idx_sam = self.listbox_sam.curselection()
        self.action.set('sample scans selected')
        
    def select_ref(self):
        self.idx_ref = self.listbox_ref.curselection()
        self.action.set('reference scans selected')
    
    def calculate(self):
        self.transmission = dict()
        if len(self.idx_sam) == len(self.idx_ref) and len(self.idx_ref) != 1:
            for i in range(0,len(self.idx_sam)):
               self.transmission[self.name_sam[self.idx_sam[i]]] = abs(self.sxsam[self.name_sam[self.idx_sam[i]]])/abs(self.sxref[self.name_ref[self.idx_ref[i]]])
        
        if len(self.idx_ref) == 1 and len(self.idx_sam) == 1:
            self.transmission[self.name_sam[self.idx_sam[0]]] = abs(self.sxsam[self.name_sam[self.idx_sam[0]]])/abs(self.sxref[self.name_ref[self.idx_ref[0]]])
        self.action.set('transmission calculated')
                
    def combine_select_ref(self):
        self.idx_ref = self.listbox_datas.curselection()
        for i in self.idx_ref:
            self.sxref[self.name_ref[i]] = (self.sxref[self.name_ref[i]] + self.sxref[self.name_ref[i+1]])/2
        self.action.set('reference combined')
            
    def checktrans(self):
        self.st_plotgo = pt.Plottool(self.window2, self.freq, self.transmission, xlimits = (0.2, 2), ylimits = (0,1))


class study_index:
    def __init__(self,frame,x,t):
        self.x = copy.deepcopy(x)
        self.t = copy.deepcopy(t)
        self.x_backup = copy.deepcopy(x)
        self.t_backup = copy.deepcopy(t)
        # self.x = dp.Basic_functions().normalize_signal_samref(self.x)
        self.frame_si = frame
        self.window_1 = ttk.Labelframe(self.frame_si,text='main')
        self.window_2 = ttk.Labelframe(self.frame_si,text='display')
        self.window_3 = ttk.LabelFrame(self.frame_si,text='status')
        self.window_4 = ttk.Labelframe(self.frame_si, text = 'save results')
        self.window_1.grid(column=0,row=0)
        self.window_2.grid(column=1,row=0)
        self.window_3.grid(column=0,row=2)
        self.window_4.grid(column = 0, row = 1)
        self.name = list(x)
        self.list = tk.StringVar(value=self.name)
        # self.si_action = tk.StringVar(value='ready')
        self.si_listbox= tk.Listbox(self.window_1,listvariable=self.list,selectmode='extended')
        self.si_idx_echo = None
        try:
            open('study.ini')
        except:
            self.si_action = tk.StringVar(value = 'Ready with no default values loaded') 
            self.choppoints_sam = tk.StringVar()
            self.choppoints_ref = tk.StringVar()
            self.choppoints_echo = tk.StringVar()
            self.L = tk.DoubleVar(value = 0.00108)
        else: 
            config = cfp.ConfigParser()
            config.read('study.ini')
            self.si_action = tk.StringVar(value = 'Ready with default values loaded') 
            self.choppoints_sam = tk.StringVar(value = config['si_default_values']['Choppoints_sam'])
            self.choppoints_ref = tk.StringVar(value = config['si_default_values']['Choppoints_ref'])
            self.choppoints_echo = tk.StringVar(value = config['si_default_values']['Choppoints_echo'])
            self.L = tk.DoubleVar(value = float(config['si_default_values']['Sample_thickness']))
            
        
        self.f0 = tk.DoubleVar(value = 0.25)
        self.f1 = tk.DoubleVar(value = 1.6)
        # self.index_save = {}
        # self.absorption_save = {}
        self.optical_constants = {}
        self.freq_save = {}
        
        
        '''window content for window_1: operate spectrums'''
        
        button_select_sam = ttk.Button(self.window_1,text='use selection as sample',command=self.si_select_sam)
        button_select_ref = ttk.Button(self.window_1,text='use selection as reference',command=self.si_select_ref)
        button_select_echo = ttk.Button(self.window_1,text='use selection as echo',command=self.si_select_echo)
        button_calculate = ttk.Button(self.window_1,text='calculate',command=self.calculate_index)
        button_show = ttk.Button(self.window_1,text='display result',command=self.display)
        
        label_L = ttk.Label(self.window_1,text='Enter estimated thickness: ')
        entry_L = ttk.Entry(self.window_1,textvariable=self.L,width=10)
        label_choppoints_sam = ttk.Label(self.window_1, text='Enter times to chop sample (start, stop): ')
        entry_choppoints_sam = ttk.Entry(self.window_1, textvariable = self.choppoints_sam, width=10)
        label_choppoints_ref = ttk.Label(self.window_1, text='Enter times to chop ref (start, stop): ')
        entry_choppoints_ref = ttk.Entry(self.window_1, textvariable = self.choppoints_ref, width=10)
        label_choppoints_echo = ttk.Label(self.window_1, text='Enter times to build echo from sample (start, stop): ')
        entry_choppoints_echo = ttk.Entry(self.window_1, textvariable = self.choppoints_echo, width=10)
        button_chopspec_sam = ttk.Button(self.window_1,text='chop sample spec',command=self.chop_sam_spec)
        button_chopspec_ref = ttk.Button(self.window_1,text='chop reference spec',command=self.chop_ref_spec)
        button_buildpec_echo = ttk.Button(self.window_1,text='chop echo spec',command=self.build_echo_spec)
        button_padspec = ttk.Button(self.window_1,text='pad specs', command=self.pad_spec)
        button_save_default = ttk.Button(self.window_1, text = 'save as default', command = self.build_config_file)
        button_show_selected = ttk.Button(self.window_1, text = 'show selected', command = self.plot_selected)
        button_restore = ttk.Button(self.window_1, text = 'restore spectrums', command = self.restore_backups)
        
        self.si_listbox.grid(column=0,row=0,columnspan=3)
        button_select_sam.grid(column=0,row=1)
        button_select_ref.grid(column=1,row=1)
        button_select_echo.grid(column=2,row=1)
        button_show_selected.grid(column = 3, row = 1)
        label_L.grid(column=0,row=2)
        entry_L.grid(column=1,row=2)
        button_calculate.grid(column=2,row=3)
        button_show.grid(column=3,row=3)
        label_choppoints_sam.grid(column = 0,row = 4)
        entry_choppoints_sam.grid(column = 1,row = 4)
        button_chopspec_sam.grid(column = 2, row = 4)
        button_padspec.grid(column = 3,row = 4)
        label_choppoints_ref.grid(column = 0,row = 5)
        entry_choppoints_ref.grid(column = 1,row = 5)
        button_chopspec_ref.grid(column = 2, row = 5)
        label_choppoints_echo.grid(column = 0,row = 6)
        entry_choppoints_echo.grid(column = 1,row = 6)
        button_buildpec_echo.grid(column = 2, row = 6)
        button_save_default.grid(column = 0, row = 7)
        button_restore.grid(column = 0, row = 8)
        
        
        
        '''window content for save results window'''
        
        # button_save = ttk.Button(self.window_4, text = 'save results', command = self.save)
        button_plotsave = ttk.Button(self.window_4, text = 'plot saved results', command = self.Edit_plot_save)
        # entry_f0 = ttk.Entry(self.window_4, textvariable = self.f0, width = 10)
        # entry_f1 = ttk.Entry(self.window_4, textvariable = self.f1, width = 10)
        # label_f0 = ttk.Label(self.window_4, text = 'Freq start from: ')
        # label_f1 = ttk.Label(self.window_4, text = 'Freq stop at: ')
        # label_f0.grid(column = 0, row = 0)
        # entry_f0.grid(column = 1, row = 0)
        # label_f1.grid(column = 2, row = 0)
        # entry_f1.grid(column = 3, row = 0)
        # button_save.grid(column = 0, row = 1)
        button_plotsave.grid(column = 0, row = 0)
        # button_save.grid(column = 0, row = 5)
        self.fig, self.ax = plt.subplots(1,1, figsize = (16,9))
        self.canvas0 = FigureCanvasTkAgg(figure = self.fig, master = self.window_2)
        self.canvas0.get_tk_widget().grid(column=0,row=0,columnspan=6)
        button_examplot = ttk.Button(self.window_2,text='exam plot',command=self.exam_plot)
        button_examplot.grid(column=0,row=1)
        
        '''window for status bar'''
        
        
        label_action = ttk.Label(self.window_3, textvariable = self.si_action)
        label_action.grid(column=0,row=0)
        
        self.bf = dp.Basic_functions()
    
    def si_select_sam(self):
        self.si_idx_sam = self.si_listbox.curselection()[0]
        self.si_action.set('sample spectrum selected')
        
    def si_select_ref(self):
        self.si_idx_ref = self.si_listbox.curselection()[0]
        self.si_action.set('reference spectrum selected')
        
    def si_select_echo(self):
        self.si_idx_echo = self.si_listbox.curselection()[0]
        self.si_action.set('echo spectrum selected')
    
    def restore_backups(self):
        self.x = copy.deepcopy(self.x_backup)
        self.t = copy.deepcopy(self.t_backup)
    
    def build_config_file(self):
        si_config = cfp.ConfigParser()
        si_config['si_default_values'] = {'Sample_thickness': str(self.L.get()),
                                           'Choppoints_sam': self.choppoints_sam.get(),
                                           'Choppoints_ref': self.choppoints_ref.get(),
                                           'Choppoints_echo': self.choppoints_echo.get()}
        with open('study.ini', 'w') as si_default_values:
            si_config.write(si_default_values)
        self.si_action.set('Current entry values are saved as default values')
            
    def plot_selected(self):
        key_sam = self.name[self.si_idx_sam]
        key_ref = self.name[self.si_idx_ref]
        self.ax.plot(self.t[key_sam], self.x[key_sam], label = 'sample')
        self.ax.plot(self.t[key_ref], self.x[key_ref], label = 'reference')
        self.ax.legend()
        self.canvas0.draw()
        
    def calculate_index(self):  
        # f_bot = self.f0.get()
        
        freq,sx = dp.Basic_functions().fftx(self.x, self.t, 2)
        freq_sam = freq[self.name[self.si_idx_sam]]
        # freq_low = freq_sam[np.where(freq_sam <= 2)]
        # self.freq_sam = freq_sam
        sx_sam = sx[self.name[self.si_idx_sam]]
        # sx_sam = sx_sam[np.where(freq_sam <= 2)]
        # freq_ref = freq[self.name[self.si_idx_ref]]
        sx_ref = sx[self.name[self.si_idx_ref]]
        # sx_ref = sx_ref[np.where(freq_sam <= 2)]
        if self.si_idx_echo == None:
            freq_new, sx_sam_new, sx_ref_new, self.index, self.absorption, self.phi_sam, self.phi_ref = self.bf.getindex_array(freq_sam, sx_sam,
                                                                                                             sx_ref,self.L.get())
        else:
            # freq_echo = freq[self.name[self.si_idx_echo]]
            sx_echo = sx[self.name[self.si_idx_echo]]
            sx_echo = sx_echo[np.where(freq_sam <= 2)]
            freq_new, sx_sam_new, sx_ref_new, self.index, self.absorption, self.phi_sam, self.phi_ref = self.bf.getindex_array(freq_sam, sx_sam, 
                                                                                                             sx_ref, self.L.get(),
                                                                                                             E_echo=sx_echo)
        self.freq = freq_new
        self.canvas0.figure = Figure(figsize=(15,8))
        plot0 = self.canvas0.figure.add_subplot()
        # plot0.plot(freq_new,abs(sx_sam))
        plot0.plot(freq_new, abs(sx_sam_new), freq_new, abs(sx_ref_new))
        # plot0.plot(freq_new, self.phi_sam, freq_new, self.phi_ref)
        if self.si_idx_echo != None:
            plot0.plot(freq_new,abs(sx_echo))
        # plot0.set_xlim(left = freq_new[0],right = 2)
        # plot0.set_ylim(bottom = 0, top = 1.5)
        plot0.set_xlabel('frequency (THz)',fontsize=15)
        plot0.set_ylabel('amplitude (V)',fontsize=15)
        plot0.grid(1)
        self.canvas0.draw()
        
        self.freq_save[self.name[self.si_idx_sam]] = self.freq
        self.optical_constants[self.name[self.si_idx_sam]+'_index'] = self.index
        self.optical_constants[self.name[self.si_idx_sam]+'_absorb'] = self.absorption
        # self.index_save[self.name[self.si_idx_sam]] = self.index
        # self.absorption_save[self.name[self.si_idx_sam]] = self.absorption
        
        self.si_action.set('computation complete')
    
    
    def display(self):
        self.canvas0.figure,self.ax = plt.subplots(figsize=(15,8)) 
        self.ax2 = self.ax.twinx()
        self.ax.plot(self.freq,self.index, color = 'red', label = self.name[self.si_idx_sam])
        self.ax.set_xlim(left=0,right=1.65)
        self.ax.set_xlabel('frequency (THz)',fontsize=15)
        self.ax.set_ylabel('refractive index',fontsize=15, color  = 'red')
        self.ax.legend(loc = 'best', fontsize = 15)
        self.ax.set_ylim(bottom = 0, top = 2)
        self.ax2.plot(self.freq,self.absorption, color = 'blue')
        self.ax2.set_xlim(left=0,right=1.65)
        self.ax2.set_xlabel('frequency (THz)',fontsize=15)
        self.ax2.set_ylabel('absorption',fontsize=15, color = 'blue')
        self.canvas0.draw()
        
        
    def save(self):
        self.freq_save[self.name[self.si_idx_sam]] = self.freq
        self.optical_constants[self.name[self.si_idx_sam]+'_index'] = self.index
        self.optical_constants[self.name[self.si_idx_sam]+'_absorb'] = self.absorption
        # self.index_save[self.name[self.si_idx_sam]] = self.index
        # self.absorption_save[self.name[self.si_idx_sam]] = self.absorption
        self.si_action.set('result saved')
        
    def plot_save(self):
        self.canvas0.figure,self.ax = plt.subplots(2,1,figsize=(15,8)) 

        for key in list(self.index_save):
            self.ax[0].plot(self.freq_save[key],self.index_save[key], label = key)
            self.ax[1].plot(self.freq_save[key],self.absorption_save[key], label = key)
        # self.ax[0].set_xlim(left=0,right=2)
        self.ax[0].set_xlabel('frequency (THz)',fontsize=15)
        self.ax[0].set_ylabel('refractive index',fontsize=15)
        self.ax[0].legend(loc = 'best', fontsize = 15)
        # self.ax[0].set_ylim(bottom = 0, top = 2)
        # self.ax[1].set_xlim(left=0,right=2)
        self.ax[1].set_xlabel('frequency (THz)',fontsize=15)
        self.ax[1].set_ylabel('absorption',fontsize=15)
        self.ax[1].legend(loc = 'best', fontsize = 15)
        self.canvas0.draw()
        self.si_action.set('Figure plotted')
        
    def Edit_plot_save(self):
        # self.canvas0.figure,self.ax = plt.subplots(2,1,figsize=(15,8)) 
        window_temp1 = tk.Toplevel()
        plotgo = pt.Plottool(window_temp1, self.freq_save, self.optical_constants, nrow = 2, ncol = 1)
        
        
    def exam_plot(self):
        root = tk.Tk()
        frame = ttk.Frame(root)
        frame.pack(fill=tk.BOTH,expand=True)
        self.canvas1 = FigureCanvasTkAgg(self.canvas0.figure,master=frame)
        self.canvas1.figure = self.canvas0.figure
        self.canvas1.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.canvas1.draw()
        self.toolbar = NavigationToolbar2Tk(self.canvas1,frame)
        
        
        
        
    def build_echo_spec(self):
        key_sam = self.name[self.si_idx_sam]
        t_points = self.choppoints_sam.get().split(',')
        t0 = float(t_points[0])
        t1 = float(t_points[1])
        idx0 = np.where(self.t[key_sam] >= t0)[0][0]
        idx1 = np.where(self.t[key_sam] <= t1)[0][-1]
        key_echo = key_sam + '_echo'
        self.x[key_echo] = self.x[key_sam][idx0:idx1]
        self.t[key_echo] = self.t[key_sam][idx0:idx1]
        self.name = list(self.x)
        self.list.set(self.name)
        self.si_action.set('Echo spectrum built')
        
        
    def chop_ref_spec(self):
        key = self.name[self.si_idx_ref]
        t_points = self.choppoints_ref.get().split(',')
        t0 = float(t_points[0])
        t1 = float(t_points[1])
        idx0 = np.where(self.t[key] >= t0)[0][0]
        idx1 = np.where(self.t[key] <= t1)[0][-1]
        self.x[key] = self.x[key][idx0:idx1]
        self.t[key] = self.t[key][idx0:idx1]
        self.si_action.set('reference spectrum chopped')
        
    def chop_sam_spec(self):
        key = self.name[self.si_idx_sam]
        t_points = self.choppoints_sam.get().split(',')
        t0 = float(t_points[0])
        t1 = float(t_points[1])
        idx0 = np.where(self.t[key] >= t0)[0][0]
        idx1 = np.where(self.t[key] <= t1)[0][-1]
        self.x[key] = self.x[key][idx0:idx1]
        self.t[key] = self.t[key][idx0:idx1]
        self.si_action.set('sample spectrum chopped')
        
    def pad_spec(self):
        key_sam = self.name[self.si_idx_sam]
        key_ref = self.name[self.si_idx_ref]
        t_sam = np.round(self.t[key_sam], 2)
        t_ref = np.round(self.t[key_ref], 2)
        # self.t_ref = t_ref
        # x_ref = self.t[key_ref]
        t_max = max(t_sam[-1], t_ref[-1])
        # self.t_max = t_max
        t_min = min(t_sam[0], t_ref[0])
        if self.si_idx_echo != None:
            key_echo = self.name[self.si_idx_echo]
            t_echo = np.round(self.t[key_echo], 2)
            t_max = max(t_max, t_echo[-1])
            t_min = min(t_min, t_echo[0])
        
        timestep = t_sam[1] - t_sam[0]
        t_all = np.round(np.arange(t_min - 5, t_max + timestep/2, step = timestep), 2)
        self.t[key_sam] = t_all
        self.t[key_ref] = t_all
        
        self.t_sam_0 = t_all[np.where(t_all < t_sam[0])]
        self.t_sam_1 = t_all[np.where(t_all > t_sam[-1])]
        t_sam_0 = t_all[np.where(t_all < t_sam[0])]
        t_sam_1 = t_all[np.where(t_all > t_sam[-1])]
        pad_sam_0 = np.zeros_like(t_sam_0)
        pad_sam_1 = np.zeros_like(t_sam_1)
        self.x[key_sam] = np.concatenate((pad_sam_0, self.x[key_sam], pad_sam_1), axis = 0)
        
        t_ref_0 = t_all[np.where(t_all < t_ref[0])]
        t_ref_1 = t_all[np.where(t_all > t_ref[-1])]
        pad_ref_0 = np.zeros_like(t_ref_0)
        pad_ref_1 = np.zeros_like(t_ref_1)
        self.x[key_ref] = np.concatenate((pad_ref_0, self.x[key_ref], pad_ref_1), axis = 0)
        
        if self.si_idx_echo != None:
            t_echo_0 = t_all[np.where(t_all < t_echo[0])]
            t_echo_1 = t_all[np.where(t_all > t_echo[-1])]
            pad_echo_0 = np.zeros_like(t_echo_0)
            pad_echo_1 = np.zeros_like(t_echo_1)
            self.x[key_echo] = np.concatenate((pad_echo_0, self.x[key_echo], pad_echo_1), axis = 0)
            self.t[key_echo] = t_all
        
        
            
        self.canvas0.figure = Figure(figsize=(15,8))
        
        
        
        plot1 = self.canvas0.figure.add_subplot(311)
        plot2 = self.canvas0.figure.add_subplot(312)
        plot3 = self.canvas0.figure.add_subplot(313)
        plot1.plot(self.t[key_sam] , self.x[key_sam],label = 'sample')
        plot2.plot(self.t[key_ref] , self.x[key_ref],label = 'reference')
        if self.si_idx_echo != None:
            plot3.plot(self.t[key_echo], self.x[key_echo],label = 'echo')
            plot3.set_xlabel('time (ps)')
            plot3.set_ylabel('amplitude (V)')
            plot3.legend(loc='best')
        plot1.set_xlabel('time (ps)')
        plot1.set_ylabel('amplitude (V)')
        plot1.legend(loc='best')
        plot2.set_xlabel('time (ps)')
        plot2.set_ylabel('amplitude (V)')
        plot2.legend(loc='best')
       
        self.canvas0.draw()
        
        
        
class study_polarimetry:
    def __init__(self,frame,freq,sx,sy):
        self.frame_po = frame
        self.freq = freq
        self.sx = sx
        self.sy = sy
        self.sx_abs = dp.Basic_functions().dict_getabs(self.sx)
        self.sy_abs = dp.Basic_functions().dict_getabs(self.sy)
        self.sxy = dp.Basic_functions().dict_total_field(self.sx_abs, self.sy_abs)
        self.window1 = ttk.Labelframe(self.frame_po,text='main')
        self.window1.grid(column=0,row=0)
        self.name = list(sx)
        self.list_name = tk.StringVar(value=self.name)
        self.selection_sam = tk.StringVar(value='Selected sample spectrums: \n')
        self.selection_ref = tk.StringVar(value='Selected reference spectrums: \n')
        self.action = tk.StringVar(value='ready')
        self.listbox= tk.Listbox(self.window1,listvariable=self.list_name,selectmode='extended')
        label_select_sam = ttk.Label(self.window1,textvariable=self.selection_sam)
        label_select_ref = ttk.Label(self.window1,textvariable=self.selection_ref)
        button_select_sam = ttk.Button(self.window1,text='use selection as sample',command=self.select_sam)
        button_select_ref = ttk.Button(self.window1,text='use selection as reference',command=self.select_ref)
        button_calculate = ttk.Button(self.window1,text='calculate',command=self.calculate_transmission)
        button_calculate_polarimetry = ttk.Button(self.window1,text='calculate polarimetry',command=self.calculate_polarimetry)
        button_select_combine_ref = ttk.Button(self.window1,text='select and combine',command=self.combine_select_ref)
        button_checktrans = ttk.Button(self.window1,text='check result',command=self.checktrans)
        button_show_polarangle = ttk.Button(self.window1,text='show polarangle',command=self.show_polarangle)
        button_show_ellipticity = ttk.Button(self.window1,text='show ellipticity',command=self.show_ellipticity)
        label_status = ttk.Label(self.window1,text='status: ')
        label_action = ttk.Label(self.window1,textvariable=self.action)
        
        
        self.listbox.grid(column=0,row=0,columnspan=2)
        label_select_sam.grid(column=2,row=0,columnspan=2)
        label_select_ref.grid(column=4,row=0,columnspan=2)
        button_select_sam.grid(column=0,row=1)
        button_select_ref.grid(column=2,row=1)
        button_select_combine_ref.grid(column=3,row=1)
        button_calculate.grid(column=4,row=2)
        button_calculate_polarimetry.grid(column=2,row=2)
        button_show_polarangle.grid(column=1,row=3)
        button_checktrans.grid(column=0,row=3)
        button_show_ellipticity.grid(column=2,row=3)
        label_status.grid(column=0,row=4)
        label_action.grid(column=1,row=4)
        
        
        self.window2 = ttk.Labelframe(self.frame_po,text='plot')
        self.window2.grid(column=0,row=1,sticky='nsew')
        
        self.window3 = ttk.LabelFrame(self.frame_po,text='plot')
        self.window3.grid(column=1,row=1,sticky='nsew')
        
        
    def select_sam(self):
        self.idx_sam = self.listbox.curselection()
        self.selection_sam.set('Selected sample spectrums: \n')
        for i in self.idx_sam:
            self.selection_sam.set(self.selection_sam.get()+self.name[i]+'\n')
        self.action.set('sample scans selected')
        
    def select_ref(self):
        self.idx_ref = self.listbox.curselection()
        self.selection_ref.set('Selected reference spectrums: \n')
        for i in self.idx_ref:
            self.selection_ref.set(self.selection_ref.get()+self.name[i]+'\n')
        self.action.set('reference scans selected')
    
    def calculate_transmission(self):
        self.transmission_x = dict()
        self.transmission_y = dict()
        if len(self.idx_sam) == len(self.idx_ref) and len(self.idx_ref) != 1:
            for i in self.idx_sam:
                for j in self.idx_ref:
                    if self.name[i] in list(self.sy):
                        self.transmission_x[self.name[i]] = abs(self.sx[self.name[i]])/abs(self.sx[self.name[j]])
                        self.transmission_y[self.name[i]] = abs(self.sy[self.name[i]])/abs(self.sy[self.name[j]])
        if len(self.idx_ref) == 1:
            for i in self.idx_sam:
                self.transmission_x[self.name[i]] = abs(self.sx[self.name[i]])/abs(self.sx[self.name[self.idx_ref[0]]])
                self.transmission_y[self.name[i]] = abs(self.sy[self.name[i]])/abs(self.sy[self.name[self.idx_ref[0]]])
        self.action.set('transmission calculated')
        
    def calculate_polarimetry(self):
        self.polarangle = dict()
        self.ellipticity = dict()
        for i in self.idx_sam:
            Ecra = self.sx[self.name[i]]+1j*self.sy[self.name[i]]
            Ecri = self.sx[self.name[i]]-1j*self.sy[self.name[i]]
            self.ellipticity[self.name[i]] = (abs(Ecri)-abs(Ecra))/(abs(Ecri)+abs(Ecra))
            self.polarangle[self.name[i]] = (np.unwrap(np.angle(Ecra))-np.unwrap(np.angle(Ecri)))/2
        
                
    def combine_select_ref(self):
        self.idx_ref = self.listbox_datas.curselection()
        for i in self.idx_ref:
            self.sxref[self.name_ref[i]] = (self.sxref[self.name_ref[i]] + self.sxref[self.name_ref[i+1]])/2
        self.action.set('reference combined')
            
    def checktrans(self):
        pt.Plottool(self.window2, self.freq, self.transmission)
        
    def show_polarangle(self):
        pt.Plottool(self.window2, self.freq, self.polarangle)
    def show_ellipticity(self):
        pt.Plottool(self.window3, self.freq, self.ellipticity)