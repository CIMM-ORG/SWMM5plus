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
unit             = 'CMS'     # unit of original swmm-c input file. options 'CFS' or 'SI'
Verbose_output   = False
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
else:
    #setting the output directory
    output_dir = inp_name+"_parallel_comparison"
#creates new output_dir for the test_case and inside of it a timestamped version 
os.system('mkdir '+ output_dir)
os.system('cd ' + output_dir)
os.system('cd ' + output_dir+ '\n  mkdir '+time_now)

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

    # allow for different find of output paths
    if has_output_path:
        output_dir_timestamped = output_dir+'/'+time_now+'_processor_'+str(num_processors[ii])+'/'
    else:
        #setting the output directory
        output_dir_timestamped = output_dir+'/'+time_now+'_processor_'+str(num_processors[ii])+'/'
    
    #setting the input, output and report paths needed for running SWMM5_C 
    #inp_path = cwd + '/' + sys.argv[1][::len(sys.argv)-1]
    out_path = output_dir_timestamped + inp_name +'.out'
    rpt_path = output_dir_timestamped + inp_name +'.rpt'

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
    
#print(h5_file_lists)

#print(h5_file_lists[0])

h5_file_num_image_1 = h5py.File(h5_file_lists[0],'r')

#Creating dictionary with opened hdf5 for each number of processor 
# "num_images" : file
dict_of_h5_files = dict([])

#Creating list of hdf5 for each processor and opening them
ii = 0
for ii in range(0,len(h5_file_lists)):
    dict_of_h5_files[num_processors[ii]] = h5py.File(h5_file_lists[ii],'r')

#get all of the data_set_names
all_dset_names=dict_of_h5_files[1].keys()
# this will be used to keep a running list of which links and nodes are now within given tolerances
list_of_errors =[]

#set units for comparison 
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


Q_row_to_append = []
Y_row_to_append = []
Q_List_Final   = []
Y_List_Final    = []
H_row_to_append = []
H_List_Final   = []

for x in all_dset_names:

    if(x[0:5]=='link_'):

        link_name = x[5::]

        # Q_row_to_append.append("Q RMSE")
        Q_row_to_append.append(link_name)
        Q_row_to_append.append("0")

        # Y_row_to_append.append("Y RMSE")
        Y_row_to_append.append(link_name)
        Y_row_to_append.append("0")

        #Store flow and depth data from the serial run of the code
        z_serial = get_array_from_dset(h5_file_lists[0],x)
        link_Q_serial = z_serial[:,3] 
        link_Y_serial = z_serial[:,2]
        time = z_serial[:,0]
       
        #loop through all of the parellel output files and perform link comparison 
        for parallel_file in h5_file_lists[1:]:
            z_parallel = get_array_from_dset(parallel_file,x)
            link_Q_parallel = z_parallel[:,3] 
            link_Y_parallel = z_parallel[:,2] 


            rmse_link_Q = np.linalg.norm(link_Q_parallel - link_Q_serial) / np.sqrt(len(link_Q_parallel))
            rmse_link_Y = np.linalg.norm(link_Y_parallel - link_Y_serial) / np.sqrt(len(link_Y_parallel))

            norm_rmse_link_Q = 100*rmse_link_Q  / np.maximum(0.1, np.maximum( np.amax(link_Q_serial), np.amax(link_Q_parallel)))
            norm_rmse_link_Y = 100*rmse_link_Y  / np.maximum(0.1, np.maximum( np.amax(link_Y_serial), np.amax(link_Y_parallel)))

            Q_row_to_append.append(norm_rmse_link_Q)
            Y_row_to_append.append(norm_rmse_link_Y)

            if(Verbose_output): 
                print(' ')
                print('-------------------------------------------------------------------------------')
                print('*** Parallel SWMM5+(SWMM5+P) to Serial SWMM5+ NUM_Images :', num_processors[h5_file_lists.index(parallel_file)], 'link :', link_name,' result comparison ***')
                if print_timeseries:
                    link_col_headers = ["Time (hrs.)","SWMM5+P Q "+Qunit, "SWMM5+S Q "+Qunit, "SWMM5+P Y "+Yunit, "SWMM5+S Y "+Yunit]
                    link_merged_array = np.array([time,link_Q_parallel, link_Q_serial, link_Y_parallel, link_Y_serial]).T
                    link_table = tabulate(link_merged_array , link_col_headers,floatfmt = ".3f")
                    print(' ') 
                    print(link_table)
                    print(' ')
                print('Flowrate   RMSE:',"%.3f" %rmse_link_Q,Qunit,'; or ',"%.2f%%" %norm_rmse_link_Q, ' normalized by ',"%.3f" %np.maximum(0.1, np.amax(link_Q_parallel)),Qunit)
                print('Flow depth RMSE:',"%.3f" %rmse_link_Y,Yunit,' ; or ',"%.2f%%" %norm_rmse_link_Y, ' normalized by ',"%.3f" %np.maximum(0.1, np.amax(link_Y_parallel)),Yunit)
                print('-------------------------------------------------------------------------------')
                print(' ')

                # ... Calculate the norms
                # calculate L1,L2,Linf norms for the swmm5+P link flowrates
                link_Q_l1_parallel   = np.linalg.norm(link_Q_parallel,1)
                link_Q_l2_parallel   = np.linalg.norm(link_Q_parallel)
                link_Q_linf_parallel = np.linalg.norm(link_Q_parallel,inf)
                # calculate L1,L2,Linf norms for the swmm5+S link flowrates
                link_Q_l1_serial   = np.linalg.norm(link_Q_serial,1)
                link_Q_l2_serial   = np.linalg.norm(link_Q_serial)
                link_Q_linf_serial = np.linalg.norm(link_Q_serial,inf)
                # calculate L1,L2,Linf norms for the swmm5+P link depths
                link_Y_l1_parallel   = np.linalg.norm(link_Y_parallel,1)
                link_Y_l2_parallel   = np.linalg.norm(link_Y_parallel)
                link_Y_linf_parallel = np.linalg.norm(link_Y_parallel,inf)
                # calculate L1,L2,Linf norms for the swmm5+S link flowrate
                link_Y_l1_serial = np.linalg.norm(link_Y_serial,1)
                link_Y_l2_serial = np.linalg.norm(link_Y_serial)
                link_Y_linf_serial = np.linalg.norm(link_Y_serial,inf)
        
                # ... Find the normalized errors for fail detection MOVED UP BY BRH
                #norm_rmse_link_Q = np.linalg.norm(1 - abs(swmmF_link_Q[:array_len_Q]/swmmC_link_Q)) / np.sqrt(len(swmmC_link_Q))

                # Fail check
                if(norm_rmse_link_Q > Qtolerance):
                    print('Failed link',link_name,'. Normalized rmse Q = ',norm_rmse_link_Q,'%')
                    list_of_errors.append('link: '+link_name+" flowrates are not within given error tolerance range")

                if(norm_rmse_link_Y > Ytolerance):
                    if (link_Q_l2_parallel < 0.001):
                        print('near-zero flow, so depth in links is not valid in SWMM-C')
                    else:
                        print('Failed link',link_name,'. Normalized rmse Y = ',norm_rmse_link_Y,'%')
                        list_of_errors.append('link: '+link_name+" depths are not within given error tolerance range")   
     
        Q_List_Final.append(Q_row_to_append[:])
        # Q_List_Final.append(Y_row_to_append[:])
        Y_List_Final.append(Y_row_to_append[:])
        del Q_row_to_append[:]
        del Y_row_to_append[:]        

    # Check if the data set is a node
    if((x[0:10]=='node_face_') or (x[0:10]=='node_elem_')):
        
        if (x[0:10]=='node_face_'):
            is_nJ2 = True
        else:
            is_nJ2 = False
        # ... store node name
        node_name = x[10::]

        z_serial = get_array_from_dset(h5_file_lists[0],x)
        if is_nJ2:  
                node_H_serial = ((z_serial[:,5] + z_serial[:,6])/2.) * Yf # averaging the u/s and d/s peizometric heads
        else:
            node_H_serial = (z_serial[:,5]) * Yf # take the JM peizometric head
        
        # extract the timestamp
        time = z_serial[:,0]
        # H_row_to_append.append("H RMSE")
        H_row_to_append.append(node_name)
        H_row_to_append.append("0")

        #loop through all of the parellel output files and perform link comparison 
        for parallel_file in h5_file_lists[1:]:

            z_parallel = get_array_from_dset(parallel_file,x)
            if is_nJ2:  
                node_H_parallel = ((z_parallel[:,5] + z_parallel[:,6])/2.) * Yf # averaging the u/s and d/s peizometric heads
            else:
                node_H_parallel = (z_parallel[:,5]) * Yf # take the JM peizometric head

            # --- RMSE of head
            rmse_node_H = np.linalg.norm(node_H_parallel - node_H_serial) / np.sqrt(len(node_H_parallel))
            H_row_to_append.append(rmse_node_H)

            if(Verbose_output): 
                print(' ')
                print('-------------------------------------------------------------------------------')
                print('*** SWMM5+P to SWMM5+S NUM_Images :', num_processors[h5_file_lists.index(parallel_file)], 'node :', node_name,' result comparison ***')
                # print node depth data
                if print_timeseries:
                    node_col_headers = ["Time (hrs.)", "SWMM5+P H "+Yunit, "SWMM5+S H "+Yunit]
                    node_merged_array = np.array([time, node_H_parallel, node_H_serial]).T
                    node_table = tabulate(node_merged_array , node_col_headers,floatfmt = ".3f")
                    print(' ')
                    print(node_table)
                    print(' ')
                print('Head RMSE:',"%.3f" %rmse_node_H,Yunit)
                print('-------------------------------------------------------------------------------')
                
                # calculate L1,L2,Linf norms for the swmm_c output ** should this be H instead of Y?
                node_Y_l1_parallel = np.linalg.norm(node_H_parallel,1)
                node_Y_l2_parallel = np.linalg.norm(node_H_parallel)
                node_Y_linf_parallel = np.linalg.norm(node_H_parallel,inf)

                # calculate L1,L2,Linf norms for the swmm_plus output
                node_Y_l1_serial   = np.linalg.norm(node_H_serial,1)
                node_Y_l2_serial   = np.linalg.norm(node_H_serial)
                node_Y_linf_serial = np.linalg.norm(node_H_serial,inf)

                # Fail check (uses absolute error for head)
                if(rmse_node_H > Htolerance):
                    print('Failed node',node_name,'. RSME Head = ',rmse_node_H, Yunit)
                    list_of_errors.append('node: '+node_name+" head errors are not within given tolerance range")  
        H_List_Final.append(H_row_to_append[:])
        del H_row_to_append[:]

if(Verbose_output):
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
else: 
    #np.set_printoptions(threshold=sys.maxsize)
    link_headers = ["Link Name"]
    node_headers = ["Node Name"]
    for processor in num_processors:
        link_headers.append("Processor "+str(processor) + " (%)")
        node_headers.append("Processor "+str(processor) + " (%)")

    # print(np.array(Q_List_Final))
    Q_table_final = tabulate(Q_List_Final,link_headers,floatfmt = ".3f")
    Y_table_final  = tabulate(Y_List_Final,link_headers,floatfmt = ".3f")
    H_table_final  = tabulate(H_List_Final,node_headers,floatfmt = ".3f")
    print(' ')
    print("--------------------------------------------------------------------")
    print('Link FLows -- Normalized Root Mean Square Error (%)')
    print("--------------------------------------------------------------------")
    print(Q_table_final)
    print("--------------------------------------------------------------------")
    print(' ')
    print("--------------------------------------------------------------------")
    print('Link Depths -- Normalized Root Mean Square Error (%)')
    print("--------------------------------------------------------------------")
    print(Y_table_final)
    print("--------------------------------------------------------------------")
    print(' ')
    print("--------------------------------------------------------------------")
    print('Node Heads -- Normalized Root Mean Square Error (%)')
    print("--------------------------------------------------------------------")
    print(H_table_final)
    print("--------------------------------------------------------------------")
    exit(0)