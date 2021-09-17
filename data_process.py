# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:47:00 2021

@author: pps
"""


import GUIs as gui
import dataprocessor as dp
import numpy as np


loadgo = gui.Loaddata_GUI() 
loadgo.root.mainloop()    
path0 = loadgo.selectfolder
expstyle0 = loadgo.expstyle_key.get()
spec2 = dp.SpecProcess(path0,expstyle0) 
[x,y,t] = spec2.loadspecs()

if expstyle0 == 'sam+ref':
    [xavg_sam,yavg_sam,tavg_sam,xavg_ref,yavg_ref,tavg_ref] = spec2.avespecs_samref(x, y, t, loadgo.nsam.get(),loadgo.nref.get())

if expstyle0 == 'sam/ref':
    [xavg,yavg,tavg] = spec2.avespecs_sam(x,y,t)
    x_out = spec2.Totalfield(xavg,yavg)
    t_out = tavg
    # plotgo = gui.Editplot(t_out,x_out)
    # plotgo.root.mainloop()
    Efield = spec2.getEfield(t_out['BCSO3_ref_0deg_n45deg_0T_14K_1']*1e-12,x_out['BCSO3_ref_0deg_n45deg_0T_14K_1'],1.3e-3/2,1.5e-4)

if expstyle0 == 'polarimetry':
    [xavg_c1,yavg_c1,tavg_c1,xavg_c2,yavg_c2,tavg_c2] = spec2.avespecs_samref(x, y, t, loadgo.nsam.get(),loadgo.nref.get())
    x_c1 = spec2.Totalfield(xavg_c1,yavg_c1)
    x_c2 = spec2.Totalfield(xavg_c2,yavg_c2)
    t_c1c2 = tavg_c1
    polargo = gui.Polarimetry_GUI(spec2.filenames)
    polargo.root.mainloop()
    spec2.polarimetry_new(x_c1, x_c2, t_c1c2, polargo.idxint_sam, polargo.idxint_ref,polargo.idx_select)
