from re import T
import sys
import os
from datetime import datetime
from tkinter import TRUE
import h5py
import csv
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

#-----------------------------------------------------------------------------------
# USER SETTING CONTROL
tol = 5.0                   # tolerance for comparing between the norms
recompile_swmmC  = False    # logical for recompiling swmmC
print_timeseries = True     # logical to print individual swmm5 vs swmm5+ link and node results
#-----------------------------------------------------------------------------------

# Getting current working directory and time when the script is ran so that we can create a new folder
cwd = os.getcwd()
time_now = str(datetime.now())
time_now = time_now.replace(' ', '_')
num_processors = 1
settings_path  = ""

# removes the SWMM5_C code from last comparison run and rebuilds it if needed 
if  recompile_swmmC or os.system('find swmm5_C'):
    os.system('rm -rf swmm5_C')
    os.system('cd interface \n cp Makefile_swmm5 src/')
    os.system('cd interface/src \n make -f Makefile_swmm5 ')

#checking if a input file is given

if(len(sys.argv) < 2):
    print('no local path to input file provided')
    exit()

if(sys.argv[1] == '-h'):
    print("--------------USEFUL INFO FOR RUNNING SCRIPT---------------")
    print("The format for running this script is python comparison_script.py -i *local path to input file* -s *local path to json setting file* -n *num of processors* ")
    print("If no setting file given it will use the default settings of swmm5_plus ")
    print("If no processor amount given, default is 1 processor")
    exit()

has_output_path = False

if((len(sys.argv)%2) != 0):

    for arg_id in range(1,len(sys.argv),2):
        arg = sys.argv[arg_id+1] 
        if(sys.argv[arg_id] == "-s"):
            settings_path  = arg
        if(sys.argv[arg_id] == "-n"):
            num_processors = arg
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
    exit()

# allow for different find of output paths
if has_output_path:
    output_dir = output_path+'/'+inp_name+"_comparison"
    output_dir_timestamped = output_dir+'/'+time_now+'/'
else:
    #setting the output directory
    output_dir = inp_name+"_comparison"
    output_dir_timestamped = output_dir+'/'+time_now+'/'

#creates new output_dir for the test_case and inside of it a timestamped version 
os.system('mkdir '+ output_dir)
os.system('cd ' + output_dir)
os.system('cd ' + output_dir+ '\n  mkdir '+time_now)

#setting the input, output and report paths needed for running SWMM5_C 
#inp_path = cwd + '/' + sys.argv[1][::len(sys.argv)-1]
out_path = output_dir_timestamped + inp_name +'.out'
rpt_path = output_dir_timestamped + inp_name +'.rpt'


#run the swmm5_C
os.system('./swmm5_C '+inp_path+' '+rpt_path+' '+out_path)

#build and run swmm5_plus
os.system("export FOR_COARRAY_NUM_IMAGES="+str(num_processors))
os.system('cd build \n make \n mv SWMM ..')

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

# this will be used to keep a running list of which links and nodes are now within given tolerances
list_of_errors =[]

# Loop through all of the data set names 
for x in all_dset_names:
    
    # Check if the data set is a link
    if(x[0:5]=='link_'):

        # ... store link name 
        link_name = x[5::]
        # ... extract SWMM5-C data
        # extract the flowrates from the swmm5_C .out file and convert to numpy array to store
        swmmC_link_Q = swmmtoolbox.extract(out_path,"link,"+link_name+',Flow_rate').to_numpy().ravel()
        # extract the flowrates from the swmm5_C .out file and convert to numpy array to store
        swmmC_link_Y = swmmtoolbox.extract(out_path,"link,"+link_name+',Flow_depth').to_numpy().ravel()

        # ... extract SWMM5+ data
        z = get_array_from_dset(swmm5_plus_dir+'/output.h5',x)
        # extract the flowrates from the swmm5_plus .h5 file
        swmmF_link_Q = z[1:,3]
        # extract the depths from the swmm5_plus .h5 file
        swmmF_link_Y = z[1:,2]
        # extract the timestamp
        time = z[1:,0]

        # print link flowrate and depth data
        rsme_link_Q = np.linalg.norm(swmmC_link_Q - swmmF_link_Q) / np.sqrt(len(swmmC_link_Q))
        rsme_link_Y = np.linalg.norm(swmmC_link_Y - swmmF_link_Y) / np.sqrt(len(swmmC_link_Y))
        print(' ')
        print('-------------------------------------------------------------------------------')
        print('*** SWMM5-C to SWMM5+ link :', link_name,' result comparison ***')
        if print_timeseries:
            link_col_headers = ["Time (hrs.)","SWMM-C Q (cms)", "SWMM5+ Q (cms)", "SWMM-C Y (m)", "SWMM5+ Y (m)"]
            link_merged_array = np.array([time,swmmC_link_Q, swmmF_link_Q, swmmC_link_Y, swmmF_link_Y]).T
            link_table = tabulate(link_merged_array , link_col_headers,floatfmt = ".3f")
            print(' ') 
            print(link_table)
            print(' ')
        print('Flowrate   RMSE SWMM-C vs SWMM5+ :',"%.3f" %rsme_link_Q)
        print('Flow depth RMSE SWMM-C vs SWMM5+ :',"%.3f" %rsme_link_Y)
        print('-------------------------------------------------------------------------------')

        # ... Calculate the norms
        # calculate L1,L2,Linf norms for the swmm_c link flowrates
        swmmC_link_Q_l1 = np.linalg.norm(swmmC_link_Q,1)
        swmmC_link_Q_l2 = np.linalg.norm(swmmC_link_Q)
        swmmC_link_Q_linf = np.linalg.norm(swmmC_link_Q,inf)
        # calculate L1,L2,Linf norms for the swmm_plus link flowrates
        swmmF_link_Q_l1 = np.linalg.norm(swmmF_link_Q,1)
        swmmF_link_Q_l2 = np.linalg.norm(swmmF_link_Q)
        swmmF_link_Q_linf = np.linalg.norm(swmmF_link_Q,inf)
        # calculate L1,L2,Linf norms for the swmm_c link depths
        swmmC_link_Y_l1 = np.linalg.norm(swmmC_link_Y,1)
        swmmC_link_Y_l2 = np.linalg.norm(swmmC_link_Y)
        swmmC_link_Y_linf = np.linalg.norm(swmmC_link_Y,inf)
        # calculate L1,L2,Linf norms for the swmm_plus link flowrate
        swmmF_link_Y_l1 = np.linalg.norm(swmmF_link_Y,1)
        swmmF_link_Y_l2 = np.linalg.norm(swmmF_link_Y)
        swmmF_link_Y_linf = np.linalg.norm(swmmF_link_Y,inf)

        #check if the L1, L2, Linf norms are within a given range and if not append to list of errors
        if(abs(swmmC_link_Q_l1 - swmmF_link_Q_l1) > tol):
            list_of_errors.append('link: '+link_name+" flowrates are not within given L1 range")
        if(abs(swmmC_link_Q_l2 - swmmF_link_Q_l2) > tol):
            list_of_errors.append('link: '+link_name+" flowrates are not within given L2 range")
        if(abs(swmmC_link_Q_linf - swmmF_link_Q_linf) > tol):
            list_of_errors.append('link: '+link_name+" flowrates are not within given Linf range")
        if(abs(swmmC_link_Y_l1 - swmmF_link_Y_l1) > tol):
            list_of_errors.append('link: '+link_name+" depths are not within given L1 range")
        if(abs(swmmC_link_Y_l2 - swmmF_link_Y_l2) > tol):
            list_of_errors.append('link: '+link_name+" depths are not within given L2 range")
        if(abs(swmmC_link_Y_linf - swmmF_link_Y_linf) > tol):
            list_of_errors.append('link: '+link_name+" depths are not within given Linf range")

    # Check if the data set is a node
    if(x[0:10]=='node_face_'):

        # ... store node name
        node_name = x[10::]

        # ... extract swmmC node data
        # extract the depths from the swmm5_C .out file and convert to numpy array to store
        swmmC_node_H = swmmtoolbox.extract(out_path,"node,"+node_name+',Hydraulic_head').to_numpy().ravel()

        # ... extract swmm5plus node data
        # extract the flowrates from the swmm5_plus .h5 file
        z = get_array_from_dset(swmm5_plus_dir+'/output.h5',x)
        swmmF_node_H = (z[1:,5] + z[1:,6])/2. # averaging the u/s and d/s peizometric heads
        # extract the timestamp
        time = z[1:,0]

        rsme_node_Y = np.linalg.norm(swmmC_node_H - swmmF_node_H) / np.sqrt(len(swmmC_node_H))
        print(' ')
        print('-------------------------------------------------------------------------------')
        print('*** SWMM5-C to SWMM5+ node :', node_name,' result comparison ***')
        # print node depth data
        if print_timeseries:
            node_col_headers = ["Time (hrs.)", "SWMMC H (m)", "SWMM5+ H (m)"]
            node_merged_array = np.array([time,swmmC_node_H, swmmF_node_H]).T
            node_table = tabulate(node_merged_array , node_col_headers,floatfmt = ".3f")
            print(' ')
            print(node_table)
            print(' ')
        print('Head RMSE SWMM-C vs SWMM5+ :',"%.3f" %rsme_node_Y)
        print('-------------------------------------------------------------------------------')

        # calculate L1,L2,Linf norms for the swmm_c output
        swmmC_node_Y_l1 = np.linalg.norm(swmmC_node_H,1)
        swmmC_node_Y_l2 = np.linalg.norm(swmmC_node_H)
        swmmC_node_Y_linf = np.linalg.norm(swmmC_node_H,inf)

        # calculate L1,L2,Linf norms for the swmm_plus output
        swmmF_node_Y_l1 = np.linalg.norm(swmmF_node_H,1)
        swmmF_node_Y_l2 = np.linalg.norm(swmmF_node_H)
        swmmF_node_Y_linf = np.linalg.norm(swmmF_node_H,inf)

        # check if the L1, L2, Linf norms are within a given range and if not append to list of errors
        if(abs(swmmC_node_Y_l1 - swmmF_node_Y_l1) > tol):
            list_of_errors.append('node: '+node_name+" depths are not within given L1 range")
        if(abs(swmmC_node_Y_l2 - swmmF_node_Y_l2) > tol):
            list_of_errors.append('node: '+node_name+" depths are not within given L2 range")
        if(abs(swmmC_node_Y_linf - swmmF_node_Y_linf) > tol):
            list_of_errors.append('node: '+node_name+" depths are not within given Linf range")
        

print(' ')
if(len(list_of_errors) == 0):
    print("no links or nodes are out of the L0, L1, and Linf norm range for given tolerance", tol)
else:
    print("Issues: some links or nodes are out of the L0, L1, and Linf norm range for given tolerance", tol)
    print(list_of_errors)

print(' ')
print("-------------------------------End of comparison-------------------------------")
print(' ')