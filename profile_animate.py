# dependencies and packages
import pandas as pd
import os
import sys
from datetime import datetime
from tkinter import TRUE
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# important functions
def get_index_from_data_array(array,array_name,data_name):
    # returns the index of specific data in the link/node dataset in h5
    # this function doesnot work for linkFV/nodeFV dataset

    # find the header in the data array and convert them to string
    header = array.attrs['header_data'][1,:].astype('U101')
    # find the colum index in the dataset matching the attribute name
    idx = [i for i, s in enumerate(header) if data_name == s]
    if idx == []:
        print("---------------------- ERROR IN get_index_from_data_array ----------------------")
        print("passed dataset name: "+data_name+" is not in "+array_name)
        print("please turn on "+data_name+" output from the SWMM5+ settings .json file")
        print("Tips:")
        print("these output data are required to be turned on for the animation script to work,")
        print("\"isHeadOut\"               : true,")
        print("\"isElemZBottomOut\"        : true,")
        print("\"isElemZCrownOut\"         : true,")
        print("\"isElemLengthOut\"         : true,")
        print("\"isHeadOut\"               : true ")
        print("-------------------------------------------------------------------------------") 
        exit()
        exit(1)
    else:
        return idx[0]

# keys for element types
CC   = 1
nJm  = 2
diag = 3

#-----------------------------------------------------------------------------------
# USER SETTING CONTROL
# select what unit to produce the animation to
unit = 'CFS'
# save animation
save_animation = False
# saved animation dpi
_dpi = 100
#-----------------------------------------------------------------------------------

# Getting current working directory and time when the script is ran so that we can create a new folder
cwd = os.getcwd()

# check all path to the SWMM5+ output path
if((len(sys.argv)%2) != 0):
    for arg_id in range(1,len(sys.argv),2):
        arg = sys.argv[arg_id+1]
        if(sys.argv[arg_id] == "-o"):
            output_path = arg
            has_output_path = True
else:
    print("---------------------------USEFUL INFO FOR RUNNING THE ANIMATION SCRIPT------------------------------")
    print("The format for running this script is python profile_animate.py -o *local path to SWMM5+ output folder")
    exit(1)

# retrieve test case name from the directory
test_case = output_path.split("/")[-2][:-7]
# animation file name
animation_name = test_case+'.gif'
# output h5 file
output_file = output_path+'/output.h5'

if os.path.isfile(output_file):
    # now we have the location of the h5 file, and the list of all the datasets in the h5 file
    h5_file = h5py.File(output_file,'r')
else:
    print("--------------ERROR: Could not find the h5 file in the SWMM5+ output folder---------------")
    print("           Make sure to turn on the hdf5 output option from SWMM5+ settings file          ")
    exit(1)
    
# get to unit conversion factors to change the units of SWMM5+ output to match SWMM5_c output
if unit == 'CFS':
    Yf = 3.28084
    Qf = 35.3147
    Yunit = '(ft)'
    Qunit = '(cfs)'
elif unit == 'GPM':
    Yf = 3.28084
    Qf = 15850.37
    Yunit = '(ft)'
    Qunit = '(gpm)'
elif unit == 'MGD': 
    Yf = 3.28084
    Qf = 22.8245
    Yunit = '(ft)'
    Qunit = '(mgd)' 
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

# find the name of all the dataset
all_dset_keys=h5_file.keys()

# Declaring empty lists for link and node names
link_name = []
node_name = []
node_type = []

# Loop through all of the data set 
# to find the link and node names
for keys in all_dset_keys:
    # Check if the data set is a link
    if(keys[0:5]=='link_'):
        # ... store link name 
        link_name.append(keys[5::])
    # Check if the data set is a node
    if((keys[0:10]=='node_face_') or (keys[0:10]=='node_elem_')):
        # ... store node name
        node_name.append(keys[10::])
        if (keys[0:10]=='node_face_'):
            node_type.append('nJ2')
        else:
            node_type.append('nJm')

# get the profile from the h5 file
all_attribute_names=h5_file.attrs.keys()
if 'profiles' not in all_attribute_names:
    print("----------------------------- ERROR IN GETTING PROFILES --------------------------------")
    print("profiles is not in the output.h5; please add profiles at the end of SWMM5 input file as:")
    print("                                                                                        ")
    print("[PROFILES]                                                                              ")
    print(";;Name           Links                                                                  ")
    print(";;-------------- ----------                                                             ")
    print("\"profile_name    \" link_1 link_2 link_3 ... ...                                       ")
    print("                                                                                        ")
    print("----------------------------------------------------------------------------------------") 
    exit()

else:
    # say we have a profile like below from the hdf file
    profile = h5_file.attrs['profiles'][0,:].astype('U13')


# say we have a profile like below from the hdf file
profile  = h5_file.attrs['profiles'][0,:].astype('U13')

element_head     = []
element_length   = []
element_zbottom  = []
element_zcrown   = []
element_type     = []
element_name     = []
feature_type     = []
feature_name     = []
feature_length   = []
feature_zcrown   = []
feature_zbottom  = []

nFeatures = 0

for feature_no,feature in enumerate(profile):
    
    if(feature_no % 2 == 0):
        # find the name of the Piezometric head dataset in hdf5
        # depending on the node type
        
        # find the index of the feature in the nodelist
        # to find the node type
        nidx = node_name.index(feature)
        
        if node_type[nidx] == 'nJm':
            # find the name of the Piezometric head dataset in the hdf5
            # using the node ID
            node_data_set = 'nodeFV_'+feature+'_PiezometricHead'
            node_elem_heads = h5_file[node_data_set][:,1]
            # append the link heads to a compined list
            element_head.append(node_elem_heads)
            # now retrieve static dataset for geometry
            # find the name of the static file in the hdf5 dataset
            static_data_set = 'nodeFV_static_'+feature
            static_data     = h5_file[static_data_set]
            # now find the index of zBottom in the link static dataset 
            index = get_index_from_data_array(static_data,static_data_set,'Element Z Bottom')
            elem_zbottom = static_data[0,index]
            element_zbottom.append(elem_zbottom)
            # now find the index of zCrown in the link static dataset 
            index = get_index_from_data_array(static_data,static_data_set,'Element Z Crown')
            elem_zcrown = static_data[0,index]
            element_zcrown.append(elem_zcrown)
            # now get the length of the node
            index = get_index_from_data_array(static_data,static_data_set,'Element Length')
            elem_lengths = static_data[0,index]
            element_length.append(elem_lengths)
            # append the node length
            feature_length.append(elem_lengths)
            # append element type
            element_type.append(nJm)
            # append the name of the feature
            feature_name.append(feature)
            # append feature type
            feature_type.append(nJm)
            # get node zcrown and bottom
            feature_zcrown.append(elem_zcrown)
            feature_zbottom.append(elem_zbottom)
            # append the name of the features the element comes from
            element_name.append(feature)
            
    else:
        # find the name of the Piezometric head dataset in the hdf5
        # using the link ID
        head_data_set   = 'linkFV_'+feature+'_PiezometricHead'
        # retrieve the actual piezometric head data from hdf5
        link_elem_heads = h5_file[head_data_set][:,1:]
        # append the link heads to a compined list
        element_head.append(link_elem_heads)
        # now retrieve static dataset for geometry
        # find the name of the static file in the hdf5 dataset
        static_data_set = 'linkFV_static_'+feature
        static_data     = h5_file[static_data_set]
        # now find the index of zBottom in the link static dataset 
        index = get_index_from_data_array(static_data,static_data_set,'Element Z Bottom')
        elem_zbottom = static_data[:,index]
        element_zbottom.append(elem_zbottom)
        # now find the index of zCrown in the link static dataset 
        index = get_index_from_data_array(static_data,static_data_set,'Element Z Crown')
        elem_zcrown = static_data[:,index]
        element_zcrown.append(elem_zcrown)
        # now find the index of element length in the link static dataset 
        index = get_index_from_data_array(static_data,static_data_set,'Element Length')
        elem_lengths = static_data[:,index]
        element_length.append(elem_lengths)
        # find the lngth of the whole link
        feature_length.append(sum(elem_lengths))
        # append element type
        element_type.append(CC * np.ones(elem_lengths.shape[0]))
        # append the name of the feature
        feature_name.append(feature)
        # append feature type
        feature_type.append(CC)
        # append link zbottom and crown
        feature_zcrown.append(0.5*(min(elem_zcrown)+max(elem_zcrown)))
        feature_zbottom.append(0.5*(min(elem_zbottom)+max(elem_zbottom)))
        # append the name of the features the element comes from
        element_name.append(np.repeat(feature,elem_lengths.shape[0]))
    
# find the simulation time
sim_time = h5_file[head_data_set][:,0]
nTimeSteps = h5_file[head_data_set].shape[0]

# combine and stack all the lists in to numpy arrays
element_head    = np.column_stack(element_head) * Yf
element_zbottom = np.hstack(element_zbottom) * Yf
element_zcrown  = np.hstack(element_zcrown) * Yf 
element_length  = np.hstack(element_length) * Yf
element_type    = np.hstack(element_type)
element_name    = np.hstack(element_name)
feature_length  = np.hstack(feature_length) * Yf
feature_zcrown  = np.hstack(feature_zcrown) * Yf
feature_zbottom = np.hstack(feature_zbottom) * Yf

 
# find nelem
nelem = element_length.shape[0]
# find the x val to plot result
xval = np.zeros(nelem, dtype=np.float64) 
length = 0


for index,item in enumerate(element_length):
    xval[index] = length + item * 0.5
    length      = length + item
            
# set the location of link node labels
xval_feature = np.zeros(feature_length.shape[0], dtype=np.float64) 
length = 0


for index,item in enumerate(feature_length):
    xval_feature[index] = length + item * 0.5
    length              = length + item
    
# animation plot
fig, ax = plt.subplots()
plt.rcParams['figure.figsize'] = [10, 4]
plt.rcParams.update({'font.size': 11})
fig.tight_layout(pad=2)
x = 0
# take 10% of the gradient difference as buffer for plotting
max_zrown = max(element_zcrown)
max_head  = np.amax(element_head)
min_zbottom = min(element_zbottom)
buffer = 0.1 * (abs(max(max_zrown,max_head))-abs(min_zbottom))

# geometry profile in the animation
for indx, feature in enumerate(feature_name):
    
    # plot the conduits first
    if feature_type[indx] == CC:
        # plot the z bottom
        ax.plot(xval[element_name==feature], element_zbottom[element_name==feature], '-k') 
        # plot the crown elevation
        ax.plot(xval[element_name==feature], element_zcrown[element_name==feature], '-k') 
        # plot vlines to seperate the links
        ax.axvline(x=np.cumsum(feature_length)[indx],color='b',linestyle='-',linewidth=0.15)
        # link labels
        ax.text(xval_feature[indx],feature_zcrown[indx]+buffer,feature_name[indx],
                rotation=90,ha='center',va='bottom', fontsize='small')
        
    # plot junction geometry
    elif feature_type[indx] == nJm:
        # JM is not at the start or end of the network
        if indx > 0 or indx < feature_type.shape[0]:
            # find the index of the JM
            J_idx = np.where(element_name==feature)[0]
            # store the index of the elem u/s, JM and elem d/s 
            J_plot_idx = np.array([J_idx-1,J_idx,J_idx+1])
            # plot the zbottom of the junction by combining elem u/s, JM and elem d/s zbottom
            ax.plot(xval[J_plot_idx], element_zbottom[J_plot_idx], '-k')
            # plot a vertical line to show the z crown at the u/s of the junction
            plt.vlines(x = xval[J_plot_idx[0]], ymin = element_zcrown[J_plot_idx[0]],
                       ymax = element_zcrown[J_plot_idx[1]], color='black')
            # plot a vertical line to show the z crown at the d/s of the junction
            plt.vlines(x = xval[J_plot_idx[2]], ymin = element_zcrown[J_plot_idx[2]],
                       ymax = element_zcrown[J_plot_idx[1]], color='black')
            
            # node labels
            ax.text(xval_feature[indx],feature_zbottom[indx]-buffer,feature_name[indx],
                rotation=90,ha='center',va='top', fontsize='small')
        
          
line,  = ax.plot(xval, element_head[x,:], '-o', markersize=2.0, color='xkcd:cerulean')
time_text = ax.text(0.98, 0.95,'',ha='right',va='top',transform=plt.gca().transAxes,fontsize='small')


plt.xlabel('Length along the profile '+Yunit)
plt.ylabel('Piezometric Head '+Yunit)
plt.xlim(min(xval),max(xval))
plt.ylim(min_zbottom-buffer,max(max_zrown,max_head)+buffer)


def animate(ii):
    line.set_ydata(element_head[x+ii,:])  # update the data.
    time_text.set_text('Time = %.1f hr.' %(sim_time[ii]))
    return line,


ani = animation.FuncAnimation(fig, animate, frames = nTimeSteps, interval=50, blit=False)
plt.show()

if save_animation:
    ani.save(animation_name, dpi = _dpi)
