import sys
import os
from datetime import datetime
import h5py
import csv
import numpy as np
from swmmtoolbox import swmmtoolbox
from numpy import inf

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



#Getting current working directory and time when the script is ran so that we can create a new folder
cwd = os.getcwd()
time_now = str(datetime.now())
time_now = time_now.replace(' ', '_')
num_processors = 1
settings_path  = ""
#REMOVE COMMENTS TO RUN FULL CODE AGAIN !!!!!!!!!!!!!!!!!!!!
#removes the SWMM5_C code from last comparison run and rebuilds it
os.system('rm -rf swmm5_C')
os.system('cd interface/src \n make -f Makefile_swmm5 ')

#checking if a input file is given
if(len(sys.argv) < 2):
    print('no local path to input file provided')
    exit()

#assuming the inp file is in the same directory as script 
#triming the input provided to store the input file name
#if(str.rfind(sys.argv[1],'/') == -1):
#    inp_name = sys.argv[1][::len(sys.argv[1])-4]

#else:    
#    index = str.rfind(sys.argv[1],'/')
#    inp_name = sys.argv[1][index+1:len(sys.argv[1])-4]

if(sys.argv[1] == '-h'):
    print("--------------USEFUL INFO FOR RUNNING SCRIPT---------------")
    print("The format for running this script is python comparison_script.py -i *local path to input file* -s *local path to json setting file* -n *num of processors* ")
    print("If no setting file given it will use the default settings of swmm5_plus ")
    print("If no processor amount given, default is 1 processor")
    exit()

if((len(sys.argv)%2) != 0):

    for arg_id in range(1,len(sys.argv),2):
        
        #if(len(sys.argv)-1 == arg_id+1):
        #    print("inside of last arg")
        #    arg = sys.argv[arg_id+1][:len(sys.argv[arg_id+1])-1]
        #    print(arg)
        #else:
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
else:
    print("incorrect amount of arguments given")
    print("run python comparison_script.py -h for info on using the script")
    exit()

print(inp_name)
print(inp_path)
#print(num_processors)
#print(settings_path)


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
print(output_dir_timestamped)

if(settings_path==""):
    os.system('./SWMM -i ' + inp_path + ' -o '+cwd+'/'+output_dir_timestamped)
else:
    os.system('./SWMM -i ' + inp_path + ' -s ' + settings_path + ' -o '+cwd+'/'+output_dir_timestamped)

#locating the swmm5_plus output directory inside of the timestamped folder
#We have to loop because when swmm5_plus runs it also names the output with a timestamped folder so we don't know it before runtime
for x in os.listdir(output_dir_timestamped):
    print(x)

    if(str.rfind(x,'.') == -1):
        swmm5_plus_dir = cwd+'/'+output_dir_timestamped+'' + x
        print(swmm5_plus_dir)

x = os.listdir(swmm5_plus_dir)[0]
swmm5_plus_dir = swmm5_plus_dir + '/' + x
    
     
#now we have the location of the h5 file, and the list of all the datasets in the h5 file
h5_file = h5py.File(swmm5_plus_dir+'/output.h5','r')
all_dset_names=h5_file.keys()

#this will be used to keep a running list of which links and nodes are now within given tolerances
list_of_errors =[]

#Loop through all of the data set names 
for x in all_dset_names:
    
    #Check if the data set is a link
    if(x[0:5]=='link_'):

        #store link name 
        link_name = x[5::]
        
        #extract the flowrates from the swmm5_C .out file
        y = swmmtoolbox.extract(out_path,"link,"+link_name+',Flow_rate')
        
        #convert to numpy array and store
        swmm_c_flowrates = y.to_numpy()

        #calculate L1,L2,Linf norms for the swmm_c output
        swmm_c_l1 = np.linalg.norm(swmm_c_flowrates,1)
        swmm_c_l2 = np.linalg.norm(swmm_c_flowrates)
        swmm_c_linf = np.linalg.norm(swmm_c_flowrates,inf)
        #print((swmm_c_l1))
        #print(swmm_c_l2)
        #print(swmm_c_linf)
        
        #extract the flowrates from the swmm5_plus .h5 file
        z = get_array_from_dset(swmm5_plus_dir+'/output.h5',x)
        swmm_plus_flowrates = z[1:,2]

        #calculate L1,L2,Linf norms for the swmm_plus output
        swmm_plus_l1 = np.linalg.norm(swmm_plus_flowrates,1)
        swmm_plus_l2 = np.linalg.norm(swmm_plus_flowrates)
        swmm_plus_linf = np.linalg.norm(swmm_plus_flowrates,inf)
        #print((swmm_plus_l1))
        #print(swmm_plus_l2)
        #print(swmm_plus_linf)

        #check if the L1, L2, Linf norms are within a given range and if not append to list of errors
        if(abs(swmm_c_l1 - swmm_plus_l1) > .01):
            list_of_errors.append('link: '+link_name+" not with in give L1 range")
        if(abs(swmm_c_l2 - swmm_plus_l2) > .01):
            list_of_errors.append('link: '+link_name+" not with in give L2 range")
        if(abs(swmm_c_linf - swmm_plus_linf) > .01):
            list_of_errors.append('link: '+link_name+" not with in give Linf range")

        

    if(x[0:10]=='node_face_'):

        #store node name
        node_name = x[10::]

        #extract the depths from the swmm5_C .out file
        y = swmmtoolbox.extract(out_path,"node,"+node_name+',Depth_above_invert')

        #convert to numpy array
        swmm_c_depths = y.to_numpy()
        
         #extract the flowrates from the swmm5_plus .h5 file
        z = get_array_from_dset(swmm5_plus_dir+'/output.h5',x)
        swmm_plus_depths = z[1:,3]

         #calculate L1,L2,Linf norms for the swmm_c output
        swmm_c_l1 = np.linalg.norm(swmm_c_depths,1)
        swmm_c_l2 = np.linalg.norm(swmm_c_depths)
        swmm_c_linf = np.linalg.norm(swmm_c_depths,inf)
        #print((swmm_c_l1))
        #print(swmm_c_l2)
        #print(swmm_c_linf)

        #calculate L1,L2,Linf norms for the swmm_plus output
        swmm_plus_l1 = np.linalg.norm(swmm_plus_depths,1)
        swmm_plus_l2 = np.linalg.norm(swmm_plus_depths)
        swmm_plus_linf = np.linalg.norm(swmm_plus_depths,inf)
        #print((swmm_plus_l1))
        #print(swmm_plus_l2)
        #print(swmm_plus_linf)

         #check if the L1, L2, Linf norms are within a given range and if not append to list of errors
        if(abs(swmm_c_l1 - swmm_plus_l1) > .001):
            list_of_errors.append('node: '+node_name+" not with in give L1 range")
        if(abs(swmm_c_l2 - swmm_plus_l2) > .001):
            list_of_errors.append('node: '+node_name+" not with in give L2 range")
        if(abs(swmm_c_linf - swmm_plus_linf) > .001):
            list_of_errors.append('node: '+node_name+" not with in give Linf range")


print("-----------------------------------------------------")
print("------------------End of comparison------------------")
print("-----------------------------------------------------")
if(len(list_of_errors) == 0):
    print("no links or nodes are out of the given L0, L1, L2, and Linf range")
else:
    print(list_of_errors)