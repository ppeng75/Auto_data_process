# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 11:14:37 2021

@author: pps
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import tkinter as tk
from tkinter import ttk
from numpy import pi
from scipy.optimize import curve_fit
from inspect import signature


# class startup:
#     def __init__(self):
#         a = 0


class Curvefit:
    def __init__(self, root, xvalues,yvalues):
        self.root = root
        self.x = xvalues
        self.y = yvalues
        self.window0 = tk.Frame(self.root)
        self.window0.grid(column=0,row=0,sticky='nw')
        self.xliml = tk.DoubleVar(value=min(self.x))
        self.xlimr = tk.DoubleVar(value=max(self.x))
        self.fun_choice = tk.StringVar()
        self.autoxlim = tk.IntVar(value=1)
        self.gridswitch = tk.IntVar(value=1)
        
        
        label_xliml = ttk.Label(self.window0,text='x from: ')
        label_xlimr = ttk.Label(self.window0,text='to: ')
        label_funchoice = ttk.Label(self.window0,text='fit type: ')
        label_autoxlim = ttk.Label(self.window0,text='auto')
        entry_xliml = ttk.Entry(self.window0,textvariable = self.xliml, width=10)
        entry_xlimr = ttk.Entry(self.window0,textvariable = self.xlimr, width=10)
        combox_functions = ttk.Combobox(self.window0,textvariable=self.fun_choice,width=10)
        combox_functions['value'] = ['linear','poly2','poly3','poly4','sinu','exponential','InSb trans']
        checkbutton_autoxlim = ttk.Checkbutton(self.window0,variable=self.autoxlim)
        label_grid = ttk.Label(self.window0, text = 'grid on')
        checkbutton_grid = ttk.Checkbutton(self.window0, variable = self.gridswitch)
        button_fit = ttk.Button(self.window0,text='fit',command=self.fitcurve)
        
        
        
    
        
        
        # widges for window1 
        
        self.window1 = ttk.Labelframe(self.root,text='plot')
        self.window1.grid(column=1,row=0,columnspan = 3, rowspan = 3, sticky='nw')
        fig1 = Figure(figsize=(12,6))
        self.canvas0 = FigureCanvasTkAgg(figure=fig1,master=self.window1)
        self.canvas0.get_tk_widget().grid(column=0,row=0,columnspan=10,rowspan=10,sticky='nsew')
        button_operate_plot = ttk.Button(self.window1,text='examine',command=self.exam_plot)
        
        plot = self.canvas0.figure.add_subplot()
        plot.plot(self.x,self.y)
        

        button_operate_plot.grid(column=0,row=12)
        
        label_xliml.grid(column=0,row=0)
        entry_xliml.grid(column=1,row=0)
        label_xlimr.grid(column=2,row=0)
        entry_xlimr.grid(column=3,row=0)
        label_autoxlim.grid(column=4,row=0)
        checkbutton_autoxlim.grid(column=5,row=0)
        label_funchoice.grid(column=0,row=1)
        combox_functions.grid(column=1,row=1,columnspan=2)
        button_fit.grid(column=0, row=2)
        label_grid.grid(column=1, row=2)
        checkbutton_grid.grid(column=2, row=2)
        
        
        
        # widges for window2
        
        self.window2 = ttk.Labelframe(self.root, text = 'parameter adjust ')
        self.window2.grid(column=0,row=1)
        
        self.a = tk.DoubleVar(value = 1)
        self.b = tk.DoubleVar(value = 1)
        self.c = tk.DoubleVar(value = 1)
        self.d = tk.DoubleVar(value = 1)
        self.e = tk.DoubleVar(value = 1)
        self.f = tk.DoubleVar(value = 1)
        self.g = tk.DoubleVar(value = 1)
     
 
        
        label_p1 = ttk.Label(self.window2, text = 'a: ')
        label_p2 = ttk.Label(self.window2, text = 'b: ')
        label_p3 = ttk.Label(self.window2, text = 'c: ')
        label_p4 = ttk.Label(self.window2, text = 'd: ')
        label_p5 = ttk.Label(self.window2, text = 'e: ')
        label_p6 = ttk.Label(self.window2, text = 'f: ')
        label_p7 = ttk.Label(self.window2, text = 'g: ')
        
        self.entry_a = ttk.Entry(self.window2, textvariable = self.a, state = 'disabled')
        self.entry_b = ttk.Entry(self.window2, textvariable = self.b, state = 'disabled')
        self.entry_c = ttk.Entry(self.window2, textvariable = self.c, state = 'disabled')
        self.entry_d = ttk.Entry(self.window2, textvariable = self.d, state = 'disabled')
        self.entry_e = ttk.Entry(self.window2, textvariable = self.e, state = 'disabled')
        self.entry_f = ttk.Entry(self.window2, textvariable = self.f, state = 'disabled')
        self.entry_g = ttk.Entry(self.window2, textvariable = self.g, state = 'disabled')
        
        label_p1.grid(column=0,row=0)
        self.entry_a.grid(column=1,row=0)
        label_p2.grid(column=2,row=0)
        self.entry_b.grid(column=3,row=0)
        label_p3.grid(column=0,row=1)
        self.entry_c.grid(column=1,row=1)
        label_p4.grid(column=2,row=1)
        self.entry_d.grid(column=3,row=1)
        label_p5.grid(column=0,row=2)
        self.entry_e.grid(column=1,row=2)
        label_p6.grid(column=2,row=2)
        self.entry_f.grid(column=3,row=2)
        label_p7.grid(column=0, row=3)
        self.entry_g.grid(column=1, row=3)
        
        # widgese for window3: show fitted parameters
        self.window3 = ttk.Labelframe(self.root, text = 'fitting results')     
        self.window3.grid(column = 0, row = 2)
        
        self.function_expression = tk.StringVar()
        label_fun_expression = ttk.Label(self.window3, textvariable = self.function_expression)
        label_fun_expression.grid(column = 0, row = 4)
        
        self.a_fitvalue = tk.StringVar(value='nan  ')
        self.b_fitvalue = tk.StringVar(value='nan  ')
        self.c_fitvalue = tk.StringVar(value='nan  ')
        self.d_fitvalue = tk.StringVar(value='nan  ')
        self.e_fitvalue = tk.StringVar(value='nan  ')
        self.f_fitvalue = tk.StringVar(value='nan  ')
        self.g_fitvalue = tk.StringVar(value='nan  ')
        
        self.label_a_value = ttk.Label(self.window3, textvariable = self.a_fitvalue)
        self.label_b_value = ttk.Label(self.window3, textvariable = self.b_fitvalue)
        self.label_c_value = ttk.Label(self.window3, textvariable = self.c_fitvalue)
        self.label_d_value = ttk.Label(self.window3, textvariable = self.d_fitvalue)
        self.label_e_value = ttk.Label(self.window3, textvariable = self.e_fitvalue)
        self.label_f_value = ttk.Label(self.window3, textvariable = self.f_fitvalue)
        self.label_g_value = ttk.Label(self.window3, textvariable = self.g_fitvalue)
        
        label_p1_fit = ttk.Label(self.window3, text = 'a: ')
        label_p2_fit = ttk.Label(self.window3, text = 'b: ')
        label_p3_fit = ttk.Label(self.window3, text = 'c: ')
        label_p4_fit = ttk.Label(self.window3, text = 'd: ')
        label_p5_fit = ttk.Label(self.window3, text = 'e: ')
        label_p6_fit = ttk.Label(self.window3, text = 'f: ')
        label_p7_fit = ttk.Label(self.window3, text = 'g: ')
        
        label_p1_fit.grid(column = 0, row = 0)
        self.label_a_value.grid(column = 1, row = 0)
        label_p2_fit.grid(column = 2, row = 0)
        self.label_b_value.grid(column = 3, row = 0)
        label_p3_fit.grid(column = 0, row = 1)
        self.label_c_value.grid(column = 1, row = 1)
        label_p4_fit.grid(column = 2, row = 1)
        self.label_d_value.grid(column = 3, row = 1)
        label_p5_fit.grid(column = 0, row = 2)
        self.label_e_value.grid(column = 1, row = 2)
        label_p6_fit.grid(column = 2, row = 2)
        self.label_f_value.grid(column = 3, row = 2)
        label_p7_fit.grid(column = 0, row = 3)
        self.label_g_value.grid(column = 1, row = 3)
        
        
        
        
        
        
    def savefit(self):
        self.y_fit_save[self.names_x[self.idx_selectx[0]]] = self.y_fit
        self.para_fit_save[self.names_x[self.idx_selectx[0]]] = self.para_fit
        self.action.set('fit saved')    
        
    def fitcurve(self):
        self.canvas0.figure = Figure(figsize=(12,6))
        plot0 = self.canvas0.figure.add_subplot()
        if self.autoxlim.get() == 1:
             self.x_new = self.x
             self.y_new = self.y
        else:
            idx0 = np.where(self.x > self.xliml.get())[0][0]
            idx1 = np.where(self.x < self.xlimr.get())[0][-1]
            self.x_new = self.x[idx0:idx1]
            self.y_new = self.y[idx0:idx1]
        plot0.plot(self.x_new,self.y_new,label='original' ,linestyle='--')
        
        if self.fun_choice.get() == 'linear':
            self.entry_a.configure(state = 'normal')
            self.entry_b.configure(state = 'normal')
            self.p_init = np.array([self.a.get(), self.b.get()])
            self.para_fit = curve_fit(self.linear,self.x_new,self.y_new, p0 = self.p_init)[0]
            self.a_fitvalue.set(str(self.para_fit[0])+'  ')
            self.b_fitvalue.set(str(self.para_fit[1])+'  ')
            self.y_fit = self.para_fit[0]*self.x_new+self.para_fit[1]
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
            self.function_expression.set('y = a*x + b')
            
        if self.fun_choice.get() == 'poly2': 
            self.entry_a.configure(state='normal')
            self.entry_b.configure(state = 'normal')
            self.entry_c.configure(state = 'normal')
            self.p_init = np.array([self.a.get(), self.b.get(), self.c.get()])
            self.para_fit = curve_fit(self.poly2,self.x_new,self.y_new, p0 = self.p_init)[0]
            self.a_fitvalue.set(str(self.para_fit[0])+'  ')
            self.b_fitvalue.set(str(self.para_fit[1])+'  ')
            self.c_fitvalue.set(str(self.para_fit[2])+'  ')
            self.y_fit = self.para_fit[0]*(self.x_new-self.para_fit[1])**2+self.para_fit[2]
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
            self.function_expression.set('y = a*(x-b)**2 + c')
            
        if self.fun_choice.get() == 'poly3': 
            self.entry_a.configure(state='normal')
            self.entry_b.configure(state = 'normal')
            self.entry_c.configure(state = 'normal')
            self.entry_d.configure(state = 'normal')
            self.p_init = np.array([self.a.get(), self.b.get(), self.c.get(), self.d.get()])
            self.para_fit = curve_fit(self.poly3,self.x_new,self.y_new, p0 = self.p_init)[0]
            self.a_fitvalue.set(str(self.para_fit[0])+'  ')
            self.b_fitvalue.set(str(self.para_fit[1])+'  ')
            self.c_fitvalue.set(str(self.para_fit[2])+'  ')
            self.d_fitvalue.set(str(self.para_fit[3])+'  ')
            self.y_fit = self.para_fit[0]*self.x**3+self.para_fit[1]*self.x**2+self.para_fit[2]*self.x+self.para_fit[3]*self.x
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
            self.function_expression.set('y = ax^3+bx^2+cx+d')
            
            
            
        if self.fun_choice.get() == 'poly4': 
            self.entry_a.configure(state='normal')
            self.entry_b.configure(state='normal')
            self.entry_c.configure(state='normal')
            self.entry_d.configure(state='normal')
            self.entry_e.configure(state='normal')
            self.function_expression.set('y=a*x^4+b*x^3+c*x^2+d*x^2+e')
            self.p_init = np.array([self.a.get(), self.b.get(), self.c.get(), self.d.get(), self.e.get()])
            self.para_fit = curve_fit(self.poly4,self.x_new,self.y_new,p0 = self.p_init)[0]
            self.y_fit = self.para_fit[0]*self.x_new**4+self.para_fit[1]*self.x_new**3+self.para_fit[2]*self.x_new**2+self.para_fit[3]*self.x_new+self.para_fit[4]
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
            self.a_fitvalue.set(str(self.para_fit[0])+'  ')
            self.b_fitvalue.set(str(self.para_fit[1])+'  ')
            self.c_fitvalue.set(str(self.para_fit[2])+'  ')
            self.d_fitvalue.set(str(self.para_fit[3])+'  ')
            self.e_fitvalue.set(str(self.para_fit[4])+'  ')
            
        if self.fun_choice.get() == 'sinu': 
            self.entry_a.configure(state='normal')
            self.entry_b.configure(state='normal')
            self.entry_c.configure(state='normal')
            self.p_init = np.array([self.a.get(), self.b.get(), self.c.get(), self.d.get()])
            self.para_fit = curve_fit(self.sinu,self.x_new,self.y_new, p0 = self.p_init)[0]
            self.y_fit = self.sinu(self.x_new,self.para_fit[0],self.para_fit[1],self.para_fit[2],self.para_fit[3])
            self.a_fitvalue.set(str(self.para_fit[0])+'  ')
            self.b_fitvalue.set(str(self.para_fit[1])+'  ')
            self.c_fitvalue.set(str(self.para_fit[2])+'  ')
            self.d_fitvalue.set(str(self.para_fit[3])+'  ')
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
            self.function_expression.set('y = a*sin(2*pi*b*x-c)+d')
            
            
            
        if self.fun_choice.get() == 'exponential': 
            self.entry_a.configure(state='normal')
            self.entry_b.configure(state='normal')
            self.entry_c.configure(state='normal')
            self.function_expression.set('y=a*np.exp(b*x+c)')
            self.p_init = np.array([self.a.get(), self.b.get(), self.c.get()])
            self.para_fit = curve_fit(self.exponential,self.x_new,self.y_new,p0=self.p_init)[0]
            self.a_fitvalue.set(str(self.para_fit[0])+'  ')
            self.b_fitvalue.set(str(self.para_fit[1])+'  ')
            self.c_fitvalue.set(str(self.para_fit[2])+'  ')
            self.y_fit = self.exponential(self.x_new,self.para_fit[0],self.para_fit[1],self.para_fit[2])
            plot0.plot(self.x_new,self.y_fit,label='fitted data')
            
            
        if self.fun_choice.get() == 'gaussian': 
            self.entry_a.configure(state='normal')
            self.entry_b.configure(state='normal')
            self.p_init = np.array([self.a.get(), self.b.get(), self.c.get()])
            self.function_expression.set('y=1/((sqrt(pi*2)*c))*a*np.exp(-(x-b)^2/(2*b^2))')
            self.para_fit = curve_fit(self.gaussian,self.x_new,self.y_new,p0=self.p_init)[0]
            self.a_fitvalue.set(str(self.para_fit[0])+'  ')
            self.b_fitvalue.set(str(self.para_fit[1])+'  ')
            self.c_fitvalue.set(str(self.para_fit[2])+'  ')
            self.y_fit = self.gaussian(self.x_new,self.para_fit[0],self.para_fit[1],self.para_fit[2])
            plot0.plot(self.x_new,self.y_fit,label='fitted data')  
            
            
        if self.fun_choice.get() == 'modcos':
            self.entry_a.configure(state='normal')
            self.entry_b.configure(state='normal')
            self.entry_c.configure(state='normal')
            self.entry_d.configure(state='normal')
            self.entry_e.configure(state='normal')
            self.entry_f.configure(state='normal')
            self.entry_g.configure(state='normal')
            self.label_fun.set('y=cos(2*pi*f*x+b)*(a*x**4+g*x**3+c*x**2+d*x+e)')
            self.p_init = np.array([self.a.get(), self.b.get(), self.c.get(), self.d.get(), self.e.get(), self.f.get(), self.g.get()])
            self.para_fit = curve_fit(self.modcos,self.x_new,self.y_new,p0=self.p_init)[0]
            self.label_p0_fitvalue.config(text=str(self.para_fit[0]))
            self.label_p1_fitvalue.config(text=str(self.para_fit[1]))
            self.label_p2_fitvalue.config(text=str(self.para_fit[2]))
            self.label_p3_fitvalue.config(text=str(self.para_fit[3]))
            self.label_p4_fitvalue.config(text=str(self.para_fit[4]))
            self.label_p5_fitvalue.config(text=str(self.para_fit[5]))
            self.label_p6_fitvalue.config(text=str(self.para_fit[6]))
            self.label_p6_fitvalue.config(text=str(self.para_fit[7]))
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
        plot0.legend(loc='best',fontsize=15)
        self.canvas0.draw()
        
    def specchop(self,x0,xt):
        self.idx0 = np.where(self.x>=x0)
        self.idxt = np.where(self.x<=xt)
        x_new = self.x[self.idx0[0][0]:self.idxt[0][-1]]
        y_new = self.y[self.idx0[0][0]:self.idxt[0][-1]]
        return x_new,y_new
    
    def exam_plot(self):
        root = tk.Tk()
        frame = ttk.Frame(root)
        frame.pack(fill=tk.BOTH,expand=True)
        self.canvas1 = FigureCanvasTkAgg(self.canvas0.figure,master=frame)
        self.canvas1.figure = self.canvas0.figure
        self.canvas1.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.canvas1.draw()
        self.toolbar = NavigationToolbar2Tk(self.canvas1,frame)
    
    def linear(self,x,a,b):
        return a*x+b
    
    def poly2(self,x,a,b,c):
        return a*(x-b)**2+c
    
    def poly3(self,x,a,b,c,d):
        return a*x**3+b*x**2+c*x+d
    
    def poly4(self,x,a,b,c,d,e):
        return a*x**4 + b*x**3 + c*x**2 + d*x + e
    
    def sinu(self,x,a,b,c,d):
        return a*np.sin(2*pi*b*x-c)+d
    
    def exponential(self,x,a,b,c):
        return a*np.exp(b*x+c)
    
    def gaussian(self, x, a, b, c):
        return 1/((np.sqrt(np.pi*2)*c))*a*np.exp(-(x-b)**2/2/c**2)
    
    def InSb_trans(self,freq,ve,epi_b):
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
        B = 0
        d=5*1e-4 #sample thickness
        w=2*pi*freq*1e12 #angular frequency 
        #set parameters for calculation
        n_air=1 # index of air
        e=1.6e-19 #electron charge
        c=3e8 # speed of light
        v_e=ve*1e12 # electron scattering rate
        m_e=0.018*9.1e-31 # mass of electron
        w_t=2*pi*5.9*1e12 # transverse phonon frequency
        w_l=2*pi*5.54*1e12 # longitudinal phonon frequency
        # epi_b=16 # background permittivity
        w_ce=e*B/(m_e) # cyclotron frequency
        w_e=0.28*2*pi*1e12*np.sqrt(epi_b) #plasma frequency is set to constant 0.28
        # w_e=0.28*2*pi*1e12
        # w_h=(4*pi*Nh*e^2/m_h*1);
        v_ph=1*2*pi*1e12 # phonon damping rate
     
       
        epi_ph = epi_b*((w_t**2-w_l**2)/(w_t**2-w**2-1j*v_ph*w))

        ncra=np.sqrt(epi_b-w_e**2/(w*(w+1j*v_e-w_ce))+epi_ph)
        ncri=np.sqrt(epi_b-w_e**2/(w*(w+1j*v_e+w_ce))+epi_ph)

        t_cra_as=2*n_air/(n_air+ncra) #calculate Fresnel transmission coefficient
        t_cra_sa=2*ncra/(n_air+ncra)
        t_cri_as=2*n_air/(n_air+ncri)
        t_cri_sa=2*ncri/(n_air+ncri)
        # calculate complex transmission coefficents
        T_cra =t_cra_as*t_cra_sa*np.exp(1j*(ncra-1)*w*d/c)
            
        return abs(T_cra)
# x = np.linspace(0, 10,100)
# y = x**2
# fitgo = Curvefit(x, y)
# fitgo.root.mainloop()