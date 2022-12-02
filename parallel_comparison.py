from re import T
import sys
import os
from datetime import datetime
from tkinter import TRUE
import h5py
import csv
import re
import numpy as np
from swmmtoolbox import swmmtoolbox
from numpy import inf
from tabulate import tabulate

def get_array_from_dset(file_name,dset_name):
    #returns the data from H5 file for data_set_name
    file = h5py.File(file_name,'r')
    all_dset_names=file.keys()
    if dset_name not in all_dset_names:
        print("---------------------- ERROR IN get_array_from_dset ----------------------")
        print("passed dataset name: "+dset_name+" is not in output.h5, only the following datasets are in output.h5")
        print("------------------------------------------------------------------")
        print(all_dset_names)   
        exit()
        return 0
    
    else:
        dset = file[dset_name]
        return dset

#-----------------------------------------------------------------------------------
# USER SETTING CONTROL
Qtolerance = 25.0            # percentage flow tolerance for comparing between the norms
Ytolerance = 25.0            # percentage depth tolerance
Q_balance_error_tol = 10.0   # flow balance tolerance error for nodes
Htolerance = 1.0             # absolute tolerance for head (m)
recompile_swmmC  = False    # logical for recompiling swmmC
print_timeseries = True     # logical to print individual swmm5 vs swmm5+ link and node results
unit             = 'CFS'     # unit of original swmm-c input file. options 'CFS' or 'SI'
#-----------------------------------------------------------------------------------

# Getting current working directory and time when the script is ran so that we can create a new folder
cwd = os.getcwd()
time_now = str(datetime.now())
time_now = time_now.replace(' ', '_')
num_processors = [1,2,3,4]
settings_path  = ""
h5_file_lists = []

volume = []
#checking if a input file is given
if(len(sys.argv) < 2):
    print('no local path to input file provided')
    exit(1)

if(sys.argv[1] == '-h'):
    print("--------------USEFUL INFO FOR RUNNING SCRIPT---------------")
    print("The format for running this script is python comparison_script.py -i *local path to input file* -s *local path to json setting file* -n *num of processors* ")
    print("If no setting file given it will use the default settings of swmm5_plus ")
    print("If no processor amount given, default is 1 processor")
    exit(1)

has_output_path = False

if((len(sys.argv)%2) != 0):

    for arg_id in range(1,len(sys.argv),2):
        arg = sys.argv[arg_id+1] 
        if(sys.argv[arg_id] == "-s"):
            settings_path  = arg
        if(sys.argv[arg_id] == "-i"):
            if(str.rfind(sys.argv[arg_id+1],'/') == -1):
                inp_name = sys.argv[arg_id+1][::len(sys.argv[1])-4]
            else:
                index = str.rfind(sys.argv[arg_id+1],'/')
                inp_name = sys.argv[arg_id+1][index+1:len(sys.argv[arg_id+1])-4]
            inp_path = cwd + '/' + arg
        if(sys.argv[arg_id] == "-o"):
            output_path = arg
            has_output_path = True
else:
    print("incorrect amount of arguments given")
    print("run python comparison_script.py -h for info on using the script")
    exit(1)

# allow for different find of output paths
if has_output_path:
    output_dir = output_path+inp_name+"_parallel_comparison"
    output_dir_timestamped = output_dir+'/'+time_now+'/'
else:
    #setting the output directory
    output_dir = inp_name+"_parallel_comparison"
    output_dir_timestamped = output_dir+'/'+time_now+'/'

#creates new output_dir for the test_case and inside of it a timestamped version 
os.system('mkdir '+ output_dir)
os.system('cd ' + output_dir)
os.system('cd ' + output_dir+ '\n  mkdir '+time_now)

#setting the input, output and report paths needed for running SWMM5_C 
#inp_path = cwd + '/' + sys.argv[1][::len(sys.argv)-1]
out_path = output_dir_timestamped + inp_name +'.out'
rpt_path = output_dir_timestamped + inp_name +'.rpt'


for ii in range(len(num_processors)):    

    #build and run swmm5_plus
    print(' ')
    print('**************************************************************')
    print('Running SWMM5+ with,',num_processors[ii],'Processors')
    print('**************************************************************')
    print(' ')
    os.system('cd build \n make \n mv SWMM ..')
    # set the number of processors 
    os.environ["FOR_COARRAY_NUM_IMAGES"] = str(num_processors[ii])

    if(settings_path==""):
        os.system('./SWMM -i ' + inp_path + ' -o ' + output_dir_timestamped)
    else:
        os.system('./SWMM -i ' + inp_path + ' -s ' + settings_path + ' -o '+ output_dir_timestamped)

    # locating the swmm5_plus output directory inside of the timestamped folder
    # We have to loop because when swmm5_plus runs it also names the output with a timestamped folder so we don't know it before runtime
    for x in os.listdir(output_dir_timestamped):
        if(str.rfind(x,'.') == -1):
            swmm5_plus_dir = output_dir_timestamped+'' + x

    x = os.listdir(swmm5_plus_dir)[0]
    swmm5_plus_dir = swmm5_plus_dir + '/' + x
        
    # now we have the location of the h5 file, and the list of all the datasets in the h5 file
    h5_file = h5py.File(swmm5_plus_dir+'/output.h5','r')
    all_dset_names=h5_file.keys()

    h5_file_name = 'processor_'+str(num_processors[ii])
    h5_file_lists.append(swmm5_plus_dir+'/output.h5')
    
print(h5_file_lists)


   