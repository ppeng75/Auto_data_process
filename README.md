# Lab-dataprocess
# dataprocessv2 starts the data process UI, the other py files are modules to realize data process functions
# This set of codes runs under python environment and does basic data process such as FFT and index calculation for a Terahertz time domain spectroscopy system. Input data should be .dat/.txt format. 
Attention: 
1. Raw data needs to be well organized and stored in a good matrix format, where the first column is x-variable inputs and following columns are different y variables.
2. Currently the program can load data that has multiple y variable inputs, but only the first y variable is used in further processing.
3. all data file name should be in the format of '*_n.dat' or where n is an integer number.
4. sample .dat files are also included
# It has a built-in UI for plot operation and curve fitting
