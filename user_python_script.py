import pandas as pd
import os
import sys
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.transforms as transforms
import matplotlib.patches as patches
from swmmtoolbox import swmmtoolbox

matplotlib.use("webagg")


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

def convert_dset_to_csv(file_name,dset_name):
    #converts hd5f data set to csv
    file = h5py.File(file_name,'r')
    all_dset_names=file.keys()
    if dset_name not in all_dset_names:
        print("---------------------- ERROR IN convert_dset_to_csv ----------------------")
        print("passed dataset name: "+dset_name+" is not in output.h5, only the following datasets are in output.h5")
        print("------------------------------------------------------------------")
        print(all_dset_names)
        exit()
        return 0

    else:
        dset = file[dset_name]
        csv_name = dset_name+'.csv'
        with open(csv_name, 'w',newline='') as csv_file:
            csvwriter = csv.writer(csv_file)
            csvwriter.writerows((dset.attrs['header_data'][:,:].astype('U13')))
            csvwriter.writerows(np.round(dset[:,:],6))

def get_index_from_data_array(array,data_name):
    # returns the index of specific data in the link/node dataset in h5
    # this function doesnot work for linkFV/nodeFV dataset

    # find the header in the data array and convert them to string
    header = array.attrs['header_data'][1,:].astype('U101')
    # find the colum index in the dataset matching the attribute name
    idx = [i for i, s in enumerate(header) if data_name == s]
    if idx == []:
        print("---------------------- ERROR IN get_index_from_data_array ----------------------")
        print("passed dataset name: "+data_name+" is not in output.h5")
        print("please turn on "+data_name+" output from SWMM5+ settings file")
        print("-------------------------------------------------------------------------------")
        print(header)   
        exit()
        return 0
    else:
        return idx[0]

def get_unit_from_data_array(array,data_name):
    # returns units of a specific data in the link/node dataset in the h5 file
    
    
    header = array.attrs['header_data'][1,:].astype('U101')
    units  = array.attrs['header_data'][2,:].astype('U101')
    # find the colum index in the dataset matching the attribute name
    idx = [i for i, s in enumerate(header) if data_name == s]
    if idx == []:
        print("---------------------- ERROR IN get_index_from_data_array ----------------------")
        print("passed dataset name: "+data_name+" is not in output.h5")
        print("please turn on "+data_name+" output from SWMM5+ settings file")
        print("-------------------------------------------------------------------------------")
        print(header)   
        exit()
        return 0
    else:
        return units[idx[0]]
    

#Used for plotting time series data of various varibles.
#This Function takes the path including the .h5 file as the file_name, 
#dset_name is the name of the data set you are acessing, data_name is the actual varible you will be acessing
#For example if you wanted to plot the depth for a link named "Link_X" the call would be the following
#emsg = time_series_plot_variable("<path_to_hdf5_file.h5>", "Link_X","Depth")
def time_series_plot_variable(file_name,output_plot_path,dset_name,data_name):

    #stores the data set for the name given
    dset = get_array_from_dset(file_name,dset_name)
    
    #stores the index in the data_set where the given data type is stored
    dset_index = get_index_from_data_array(dset,data_name)
    
    #stores the units for the data type given 
    datatype_units = get_unit_from_data_array(dset,data_name)
    
    #stores the time units for this time series 
    time_units = get_unit_from_data_array(dset,"Time")

    #storing index of the time column to be used for the plot, which is always no matter the dataset
    time_index = 0 
    
    #Formats the X axis to be longer than the Y, for readabilty 
    plt.figure(figsize=(10, 4))
    
    #Plots the data against the time 
    plt.plot(dset[:,time_index],dset[:,dset_index], linestyle="-")
    
    #Adds a grid the the plot
    plt.grid()

    #Labeling X,Y and Title for the plot
    plt.xlabel("Time" + " ("+time_units+")")
    plt.ylabel(data_name +" ("+datatype_units+")")
    plt.title(dset_name+" "+data_name+" Time Series")

    #Saves plot to plot_path given    
    plt.savefig(output_plot_path+dset_name+"_"+data_name+'_time_series.png',dpi=500)
    
    #Open the plot in a broswer for the user to see 
    plt.show()
    plt.clf()
    plt.close()

    #Don't return anything instead use pass
    pass 


#=================================================================================================================================
#-------------------------------FUNCTION DEFINITIONS FINISHED : START OF ACTUAL SCRIPT--------------------------------------------
#=================================================================================================================================

#Used for checking if paths to the output.h5 file and user defined plot path are given
has_output_path = False
has_plot_path   = False 

# This reads in the command line output which is used specify which folder to read the output.h5 file and what path to output plots to
# -o is used to flag which folder stores the output.h5 file, and -p is the path to fold which is used to store plots when created
if((len(sys.argv)%2) != 0):
    for arg_id in range(1,len(sys.argv),2):
        arg = sys.argv[arg_id+1]
        if(sys.argv[arg_id] == "-o"):
            output_path = arg
            has_output_path = True
        if(sys.argv[arg_id] == "-p"):
            plot_path = arg 
            has_plot_path = True
            
# if no arguments are put on the command line then the script will return this and quit
else:
    print("---------------------------USEFUL INFO FOR RUNNING THE User SCRIPT------------------------------")
    print("The format for running this script is python user_python_script.py -o <local path to SWMM5+ output folder>")
    print("SEE EXAMPLE RUN BELOW")
    print("python profile_animate.py -o SN_rectangular_variable_inflow_output/test_20230417_1749 -p Output_Plots/ ")
    print("")
    quit()

# if no plot path is specified then it will write to the current directory
if(has_plot_path == False):
    plot_path = ""

# retrieve test case name from the directory
test_case = output_path.split("/")[-2][:-7]

# sets the output h5 file
output_file = output_path+'/output.h5'
h5_file = h5py.File(output_file,'r')


#creates a time series plot the flowrate in the node "node_elem_15009" and will output to given plot_path
time_series_plot_variable(output_file, plot_path,"node_elem_15009","Flowrate")