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
import matplotlib.pyplot as plt
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
Qtolerance = 25.0            # percentage flow tolerance for comparing between the norms
Ytolerance = 25.0            # percentage depth tolerance
Q_balance_error_tol = 10.0   # flow balance tolerance error for nodes
Htolerance = 1.0             # absolute tolerance for head (m)
recompile_swmmC  = False    # logical for recompiling swmmC
print_timeseries = True     # logical to print individual swmm5 vs swmm5+ link and node results
print_plots      = True     # logical to print plots at the comparison results folder
unit             = 'CFS'     # unit of original swmm-c input file. options 'CFS' or 'SI'
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
    exit(1)

# allow for different find of output paths
if has_output_path:
    output_dir = output_path+inp_name+"_comparison"
    output_dir_timestamped = output_dir+'/'+time_now+'/'
    plot_dir = output_dir_timestamped+'/'+inp_name+'_plots'
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

# read the report file (rpt) to get the units in swmm5c
file = open(rpt_path, 'r')
# set the line counter to zero
ii = 0
for line in file.readlines():
    # counter for line number
    ii = ii+1
    if re.search('Flow Units', line, re.I):
        unit = line[29:32]
    elif re.search('Node Inflow Summary', line, re.I):
        # node summary data starts from the 9th line from 'Node Inflow Summary'
        node_summary_start = ii + 9
    elif re.search('Node Surcharge Summary', line, re.I):
        # node summary data ends at the 4th line above 'Node Surcharge Summary'
        node_summary_end = ii - 4
        # node surcharge data starts from the 9th line from 'Node Surcharge Summary'
        node_surcharge_start = ii + 9
    elif re.search('Node Flooding Summary', line, re.I):
        # node summary data ends at the 4th line above 'Node Flooding Summary'
        node_surcharge_end = ii - 4
    elif re.search('Conduit Surcharge Summary', line, re.I):
        # conduit surcharge data starts from the 8th line from 'Conduit Surcharge Summary'
        link_surcharge_start = ii + 8
    elif re.search('Analysis begun on', line, re.I):
        # conduit surcharge data ends at the 3th line above 'Analysis begun on'
        link_surcharge_end = ii - 3

# Go back to position 0
# Or beginning of file
file.seek(0, 0)
# go through line by line again
jj = 0
# empty list for the problematic nodes
problem_nodes = []
surcharged_links = []
surcharged_nodes = []
flow_balance_error_percents = []
for line in file.readlines():
    # counter for line number
    jj = jj+1
    # remove the whitespace from the start and end of the line
    line = line.strip()

    # appending nodes with mass balance error
    if jj >= node_summary_start and jj <= node_summary_end:
        # convert last 6 strings to float
        flow_balance_error = float(line[-6:])
        if abs(flow_balance_error) >= Q_balance_error_tol:
            problem_nodes.append(line.split()[0])
            flow_balance_error_percents.append(flow_balance_error)
    # appending surcharged nodes
    elif jj >= node_surcharge_start and jj <= node_surcharge_end:
        surcharged_nodes.append(line.split()[0])
    # appending surcharged links
    elif jj >= link_surcharge_start and jj <= link_surcharge_end:
        surcharged_links.append(line.split()[0])

# this will be used to keep a running list of which links and nodes are now within given tolerances
list_of_errors =[]

if unit == 'CFS':
    Yf = 3.28084
    Qf = 35.3147
    Yunit = '(ft)'
    Qunit = '(cfs)'
    Htolerance = Htolerance * Yf
elif unit == 'GPM':
    Yf = 3.28084
    Qf = 15850.37
    Yunit = '(ft)'
    Qunit = '(gpm)'
    Htolerance = Htolerance * Yf
elif unit == 'MGD': 
    Yf = 3.28084
    Qf = 22.8245
    Yunit = '(ft)'
    Qunit = '(mgd)'
    Htolerance = Htolerance * Yf  
elif unit == 'CMS':
    Yf = 1.
    Qf = 1.
    Yunit = '(m)'
    Qunit = '(cms)'
elif unit == 'LPS':
    Yf = 1.
    Qf = 1000.000
    Yunit = '(m)'
    Qunit = '(lps)'
elif unit == 'MLD':
    Yf = 1.
    Qf = 84.6000
    Yunit = '(m)'
    Qunit = '(mld)'
else:
    print('Worng unit type seletced')
    exit(1)
# make the plot directory if required
if print_plots:
    os.system('mkdir '+ plot_dir)  

# Loop through all of the data set names 
for x in all_dset_names:
    
    # Check if the data set is a link
    if(x[0:5]=='link_'):

        # ... store link name 
        link_name = x[5::]
        # find if any of the nodes exceeds the mass balance tolerance
        if any([link in link_name for link in surcharged_links]):
            idx = surcharged_links.index(link_name)
            print(' ')
            print('-------------------------------------------------------------------------------')
            print(link_name, 'is surcharged and has been excluded from comparison')
            print('-------------------------------------------------------------------------------')
            print(' ')
        else:
            # ... extract SWMM5-C data
            # extract the flowrates from the swmm5_C .out file and convert to numpy array to store
            swmmC_link_Q = swmmtoolbox.extract(out_path,"link,"+link_name+',Flow_rate').to_numpy().ravel()
            # extract the flowrates from the swmm5_C .out file and convert to numpy array to store
            swmmC_link_Y = swmmtoolbox.extract(out_path,"link,"+link_name+',Flow_depth').to_numpy().ravel()

            # ... extract SWMM5+ data
            z = get_array_from_dset(swmm5_plus_dir+'/output.h5',x)
            # extract the flowrates from the swmm5_plus .h5 file
            swmmF_link_Q = z[:,3] * Qf
            # extract the depths from the swmm5_plus .h5 file
            swmmF_link_Y = z[:,2] * Yf
            # extract the timestamp
            time = z[:,0]
            array_len_Q = len(swmmC_link_Q)
            array_len_Y = len(swmmC_link_Y)
            # print link flowrate and depth data

            # RMSE
            rmse_link_Q = np.linalg.norm(swmmC_link_Q - swmmF_link_Q[:array_len_Q]) / np.sqrt(len(swmmC_link_Q))
            rmse_link_Y = np.linalg.norm(swmmC_link_Y - swmmF_link_Y[:array_len_Y]) / np.sqrt(len(swmmC_link_Y))
            # ... RMSE % error normalized by either the maximum flowrate, depth or 0.1
            norm_rmse_link_Q = 100*rmse_link_Q  / np.maximum(0.1, np.maximum( np.amax(swmmF_link_Q), np.amax(swmmC_link_Q)))
            norm_rmse_link_Y = 100*rmse_link_Y  / np.maximum(0.1, np.maximum( np.amax(swmmF_link_Y), np.amax(swmmC_link_Y)))
            print(' ')
            print('-------------------------------------------------------------------------------')
            print('*** SWMM5-C to SWMM5+ link :', link_name,' result comparison ***')
            if print_timeseries:
                link_col_headers = ["Time (hrs.)","SWMM-C Q "+Qunit, "SWMM5+ Q "+Qunit, "SWMM-C Y "+Yunit, "SWMM5+ Y "+Yunit]
                link_merged_array = np.array([time[:array_len_Q],swmmC_link_Q, swmmF_link_Q[:array_len_Q], swmmC_link_Y, swmmF_link_Y[:array_len_Y]]).T
                link_table = tabulate(link_merged_array , link_col_headers,floatfmt = ".3f")
                print(' ') 
                print(link_table)
                print(' ')
            print('Flowrate   RMSE:',"%.3f" %rmse_link_Q,Qunit,'; or ',"%.2f%%" %norm_rmse_link_Q, ' normalized by ',"%.3f" %np.maximum(0.1, np.amax(swmmC_link_Q)),Qunit)
            print('Flow depth RMSE:',"%.3f" %rmse_link_Y,Yunit,' ; or ',"%.2f%%" %norm_rmse_link_Y, ' normalized by ',"%.3f" %np.maximum(0.1, np.amax(swmmC_link_Y)),Yunit)
            print('-------------------------------------------------------------------------------')
            print(' ')
            # ... Calculate the norms
            # calculate L1,L2,Linf norms for the swmm_c link flowrates
            swmmC_link_Q_l1   = np.linalg.norm(swmmC_link_Q,1)
            swmmC_link_Q_l2   = np.linalg.norm(swmmC_link_Q)
            swmmC_link_Q_linf = np.linalg.norm(swmmC_link_Q,inf)
            # calculate L1,L2,Linf norms for the swmm_plus link flowrates
            swmmF_link_Q_l1   = np.linalg.norm(swmmF_link_Q[:array_len_Q],1)
            swmmF_link_Q_l2   = np.linalg.norm(swmmF_link_Q[:array_len_Q])
            swmmF_link_Q_linf = np.linalg.norm(swmmF_link_Q[:array_len_Q],inf)
            # calculate L1,L2,Linf norms for the swmm_c link depths
            swmmC_link_Y_l1   = np.linalg.norm(swmmC_link_Y,1)
            swmmC_link_Y_l2   = np.linalg.norm(swmmC_link_Y)
            swmmC_link_Y_linf = np.linalg.norm(swmmC_link_Y,inf)
            # calculate L1,L2,Linf norms for the swmm_plus link flowrate
            swmmF_link_Y_l1 = np.linalg.norm(swmmF_link_Y[:array_len_Y],1)
            swmmF_link_Y_l2 = np.linalg.norm(swmmF_link_Y[:array_len_Y])
            swmmF_link_Y_linf = np.linalg.norm(swmmF_link_Y[:array_len_Y],inf)

            if print_plots:
                Q_plot_name = plot_dir+'/'+'link_'+link_name+'_Q.png'
                Y_plot_name = plot_dir+'/'+'link_'+link_name+'_Y.png'
                # link flowrate plot
                plt.figure(1)
                plt.plot(figsize=(8,6))
                plt.plot(time[:array_len_Q],swmmC_link_Q[:array_len_Q],'o-',color='xkcd:red',markersize=4,label='SWMM-C')
                plt.plot(time[:array_len_Q],swmmF_link_Q[:array_len_Q],'o-',color='xkcd:blue',markersize=4,label='SWMM5+')
                plt.xlabel('Time (hrs.)')
                plt.ylabel('Flowrate '+Qunit)
                plt.legend(loc='best',facecolor='white',framealpha=1.0)
                plt.savefig(Q_plot_name,bbox_inches = 'tight',pad_inches=0, format='png')
                plt.close()
                # link depth plot
                plt.figure(2)
                plt.plot(figsize=(8,6))
                plt.plot(time[:array_len_Y],swmmC_link_Y[:array_len_Y],'o-',color='xkcd:red',markersize=4,label='SWMM-C')
                plt.plot(time[:array_len_Y],swmmF_link_Y[:array_len_Y],'o-',color='xkcd:blue',markersize=4,label='SWMM5+')
                plt.xlabel('Time (hrs.)')
                plt.ylabel('Depth '+Yunit)
                plt.legend(loc='best',facecolor='white',framealpha=1.0)
                plt.savefig(Y_plot_name,bbox_inches = 'tight',pad_inches=0, format='png')
                plt.close()

            # ... Find the normalized errors for fail detection MOVED UP BY BRH
            #norm_rmse_link_Q = np.linalg.norm(1 - abs(swmmF_link_Q[:array_len_Q]/swmmC_link_Q)) / np.sqrt(len(swmmC_link_Q))

            # Fail check
            if(norm_rmse_link_Q > Qtolerance):
                print('Failed link',link_name,'. Normalized rmse Q = ',norm_rmse_link_Q,'%')
                list_of_errors.append('link: '+link_name+" flowrates are not within given error tolerance range")

            if(norm_rmse_link_Y > Ytolerance):
                if (swmmC_link_Q_l2 < 0.001):
                    print('near-zero flow, so depth in links is not valid in SWMM-C')
                else:
                    print('Failed link',link_name,'. Normalized rmse Y = ',norm_rmse_link_Y,'%')
                    list_of_errors.append('link: '+link_name+" depths are not within given error tolerance range")    

    # Check if the data set is a node
    if((x[0:10]=='node_face_') or (x[0:10]=='node_elem_')):
        
        if (x[0:10]=='node_face_'):
            is_nJ2 = True
        else:
            is_nJ2 = False
        # ... store node name
        node_name = x[10::]

        # find if any of the nodes exceeds the mass balance tolerance
        if any([node in node_name for node in problem_nodes]):
            idx = problem_nodes.index(node_name)
            print(' ')
            print('-------------------------------------------------------------------------------')
            print(node_name, 'has a flow balance percent error of', flow_balance_error_percents[idx], 'which is greater than the given tolerance',Q_balance_error_tol)
            print('Thus node, ', node_name, 'has been excluded from the compariosn')
            print('-------------------------------------------------------------------------------')
            print(' ')
        elif any([node in node_name for node in surcharged_nodes]):
            print(' ')
            print('-------------------------------------------------------------------------------')
            print(node_name, 'is surcharged and has been excluded from comparison')
            print('-------------------------------------------------------------------------------')
            print(' ')
        else:
            # ... extract swmmC node data
            # extract the depths from the swmm5_C .out file and convert to numpy array to store
            swmmC_node_H = swmmtoolbox.extract(out_path,"node,"+node_name+',Hydraulic_head').to_numpy().ravel()

            # ... extract swmm5plus node data
            # extract the flowrates from the swmm5_plus .h5 file
            z = get_array_from_dset(swmm5_plus_dir+'/output.h5',x)
            if is_nJ2:  
                swmmF_node_H = ((z[:,5] + z[:,6])/2.) * Yf # averaging the u/s and d/s peizometric heads
            else:
                swmmF_node_H = (z[:,5]) * Yf # take the JM peizometric head
            # extract the timestamp
            time = z[:,0]

            array_len_H = len(swmmC_node_H)

            # --- RMSE of head
            rmse_node_H = np.linalg.norm(swmmC_node_H - swmmF_node_H[:array_len_H]) / np.sqrt(len(swmmC_node_H))
            print(' ')
            print('-------------------------------------------------------------------------------')
            print('*** SWMM5-C to SWMM5+ node :', node_name,' result comparison ***')
            # print node depth data
            if print_timeseries:
                node_col_headers = ["Time (hrs.)", "SWMMC H "+Yunit, "SWMM5+ H "+Yunit]
                node_merged_array = np.array([time[:array_len_H],swmmC_node_H, swmmF_node_H[:array_len_H]]).T
                node_table = tabulate(node_merged_array , node_col_headers,floatfmt = ".3f")
                print(' ')
                print(node_table)
                print(' ')
            print('Head RMSE:',"%.3f" %rmse_node_H,Yunit)
            print('-------------------------------------------------------------------------------')

            # calculate L1,L2,Linf norms for the swmm_c output
            swmmC_node_Y_l1 = np.linalg.norm(swmmC_node_H,1)
            swmmC_node_Y_l2 = np.linalg.norm(swmmC_node_H)
            swmmC_node_Y_linf = np.linalg.norm(swmmC_node_H,inf)

            # calculate L1,L2,Linf norms for the swmm_plus output
            swmmF_node_Y_l1 = np.linalg.norm(swmmF_node_H[:array_len_H],1)
            swmmF_node_Y_l2 = np.linalg.norm(swmmF_node_H[:array_len_H])
            swmmF_node_Y_linf = np.linalg.norm(swmmF_node_H[:array_len_H],inf)

            if print_plots:
                H_plot_name = plot_dir+'/'+'node_'+node_name+'_H.png'
                # node head plot
                plt.figure(3)
                plt.plot(figsize=(8,6))
                plt.plot(time[:array_len_H],swmmC_node_H[:array_len_H],'o-',color='xkcd:red',markersize=4,label='SWMM-C')
                plt.plot(time[:array_len_H],swmmF_node_H[:array_len_H],'o-',color='xkcd:blue',markersize=4,label='SWMM5+')
                plt.xlabel('Time (hrs.)')
                plt.ylabel('Head '+Yunit)
                plt.legend(loc='best',facecolor='white',framealpha=1.0)
                plt.savefig(H_plot_name,bbox_inches = 'tight',pad_inches=0, format='png')
                plt.close()

            # Fail check (uses absolute error for head)
            if(rmse_node_H > Htolerance):
                print('Failed node',node_name,'. RSME Head = ',rmse_node_H, Yunit)
                list_of_errors.append('node: '+node_name+" head errors are not within given tolerance range")  

print(' ')
if(len(list_of_errors) == 0):
    print("no links or nodes are out of rangee given % tolerance", Qtolerance, Ytolerance, Htolerance)
else:
    print("Issues: some links or nodes are out of the L0, L1, and Linf norm range for givem tolerance", Qtolerance, Ytolerance, Htolerance)
    print(list_of_errors)

print(' ')
print("-------------------------------End of comparison-------------------------------")
print(' ')
print(' ')

if(len(list_of_errors) == 0):
    print("no links or nodes are out of range for given % tolerance", Qtolerance, Ytolerance)
    sys.stdout.write("no links or nodes are out of the L0, L1, and Linf norm range for given tolerance {Qtolerance,Ytolerance}")
    exit(0)
else:
    print("Issues: some links or nodes are out of the L0, L1, and Linf norm range for given tolerance {Qtol,Ytol}")
    print(list_of_errors)
    sys.stderr.write(list_of_errors)
    exit(1)