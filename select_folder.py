# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:39:17 2024

@author: Admin
"""

import numpy as np

import Fast_data_process as fdp

import tkinter as tk
from tkinter import ttk

class select_folder:
    def __init__(self):
        root = tk.Tk()
        window = ttk.Frame(master = root)
        window.grid(column=0, row=0)
        label = ttk.Label(master = window, text = 'Please Select Folder: ')
        button_selectfolder = ttk.Button(master = window, text = 'select', command = self.select_folder)
        button_confirm = ttk.Button(master = window, text = 'confirm', command = root.destroy)
        
        label.grid(column=0,row=0)
        button_selectfolder.grid(column=1,row=0)
        button_confirm.grid(column=2,row=0)
        
        root.mainloop()

    def select_folder(self):
        self.folder = tk.filedialog.askdirectory()
        

if __name__ == '__main__':
    sf = select_folder()
    folder = sf.folder
    