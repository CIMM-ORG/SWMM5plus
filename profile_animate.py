# dependencies and packages
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
# matplotlib.use("webagg")
# matplotlib.rcParams['webagg.address'] = "0.0.0.0"

# important functions
def get_index_from_data_array(array,array_name,data_name):
    # returns the index of specific data in the link/node dataset in h5
    # this function doesnot work for linkFV/nodeFV dataset

    # find the header in the data array and convert them to string
    header = array.attrs['header_data'][1,:].astype('U')
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
        print("-------------------------------------------------------------------------------") 
        quit()
        
    else:
        return idx[0]

# keys for element type
# all the links (including diags are considered CC in this animation code)
CC   = 1
nJ2  = 2
nJm  = 3

#-----------------------------------------------------------------------------------
# USER SETTING CONTROL
# select what unit to produce the animation to
unit = 'CMS'
# save animation
save_animation = False
#-----------------------------------------------------------------------------------
# profile name given by arg
arg_profile_name = ""
has_arg_profile = False

# Getting current working directory and time when the script is ran so that we can create a new folder
cwd = os.getcwd()

# check all path to the SWMM5+ output path
# check for single profile to output
if((len(sys.argv)%2) != 0):
    for arg_id in range(1,len(sys.argv),2):
        arg = sys.argv[arg_id+1]
        if(sys.argv[arg_id] == "-o"):
            output_path = arg
            has_output_path = True
        if(sys.argv[arg_id] == "-s"):
            arg_profile_name = arg
            has_arg_profile = True
else:
    print("---------------------------USEFUL INFO FOR RUNNING THE ANIMATION SCRIPT------------------------------")
    print("The format for running this script is python profile_animate.py -o <local path to SWMM5+ output folder> -s <name of specified profile>")
    print("python profile_animate.py -o SN_rectangular_variable_inflow_output/test_20230417_1749 -s test_2")
    print("")
    quit()

# retrieve test case name from the directory
if os.name == "posix":  # POSIX systems (Linux, macOS, etc.)
    test_case = output_path.split("/")[-2][:-7]
elif os.name == "nt":  # Windows
    test_case = output_path.split("\\")[-2][:-7]   
else:
    print("Unsupported operating system.")
    quit()

# output h5 file
output_file = output_path+'/output.h5'



if os.path.isfile(output_file):
    # now we have the location of the h5 file, and the list of all the datasets in the h5 file
    h5_file = h5py.File(output_file,'r')
else:
    print("--------------ERROR: Could not find the h5 file in the SWMM5+ output folder---------------")
    print("           Make sure to turn on the hdf5 output option from SWMM5+ settings file          ")
    quit()
    
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
    print('Wrong unit type selected')
    quit()

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
    if(keys[0:5]=='link_' and  keys[0:11] != 'link_static'):
        # ... store link name 
        link_name.append(keys[5::])
    # Check if the data set is a node
    if((keys[0:10]=='node_face_') or (keys[0:10]=='node_elem_')):
        # ... store node name
        node_name.append(keys[10::])
        if keys[0:10]=='node_face_':
            node_type.append('nJ2')
        elif keys[0:10]=='node_elem_':
            node_type.append('nJm')

# get the profile from the h5 file
all_attribute_names=h5_file.attrs.keys()
if not all_attribute_names:
    print("------------------------ MISSING ATTRIBUTE SET IN .H5 FILE -------------------------")
    print("missing attributes in the output.h5; add profiles at the end of SWMM5 input file as:")
    print("                                                                                    ")
    print("[PROFILES]                                                                          ")
    print(";;Name           Links                                                              ")
    print(";;-------------- ----------                                                         ")
    print("\"profile_name    \" link_1 link_2 link_3 ... ...                                   ")
    print("                                                                                    ")
    print("------------------------------------------------------------------------------------") 
    quit()

for profile_name_test in all_attribute_names:
    if not all_attribute_names:
        print("----------------------------- ERROR IN GETTING PROFILES --------------------------------")
        print("profiles is not in the output.h5; please add profiles at the end of SWMM5 input file as:")
        print("                                                                                        ")
        print("[PROFILES]                                                                              ")
        print(";;Name           Links                                                                  ")
        print(";;-------------- ----------                                                             ")
        print("\"profile_name    \" link_1 link_2 link_3 ... ...                                       ")
        print("                                                                                        ")
        print("----------------------------------------------------------------------------------------") 
        quit()

    else:
        # say we have a profile like below from the hdf file
        if(has_arg_profile):
            if(arg_profile_name.strip() not in all_attribute_names):
                print("----------------------------- ERROR IN GETTING SPECIFIED PROFILE --------------------------------")
                print("--------------------------- BELOW IS LIST OF PROFILES IN .H5 FILE -----------------==------------")
                print(all_attribute_names)
                quit()

            else:
                profile_name = arg_profile_name
                profile = h5_file.attrs[arg_profile_name][0,:].astype('U')
        else:
            profile_name = profile_name_test
            profile = h5_file.attrs[profile_name_test][0,:].astype('U')

    # say we have a profile like below from the hdf file
    element_head     = []
    element_length   = []
    element_zbottom  = []
    element_zcrown   = []
    element_type     = []
    element_name     = []
    element_rank     = []
    feature_type     = []
    feature_length   = []

    nFeatures = 0

    for feature_no,feature in enumerate(profile):

        if(feature_no % 2 == 0):
            # find the name of the Piezometric head dataset in hdf5
            # depending on the node type

            # find the index of the feature in the nodelist
            # to find the node type
            nidx = node_name.index(feature)

            if node_type[nidx] == 'nJ2':
                # find the head of the nJ2 face
                node_data_set = 'node_face_'+feature
                node_data = h5_file[node_data_set]
                index_1 = get_index_from_data_array(node_data,node_data_set,'PiezometricHeadUpstream')
                index_2 = get_index_from_data_array(node_data,node_data_set,'PiezometricHeadDownstream')
                node_face_heads = 0.5 * (node_data[:,index_1] + node_data[:,index_2])
                # append the nJ2 heads to a compined list
                element_head.append(node_face_heads)
                # append rest of the geometry as zero (updated later)
                element_zbottom.append(0)
                element_zcrown.append(0)
                element_length.append(0)
                feature_length.append(0)
                # append element type
                element_type.append(nJ2)
                # append feature type
                feature_type.append(nJ2)
                # append the name of the features the element comes from
                element_name.append(feature)
                element_rank.append(feature_no)


            elif node_type[nidx] == 'nJm':
                # find the name of the Piezometric head dataset in the hdf5
                # using the node ID
                node_data_set = 'nodeFV_'+feature+'_PiezometricHead'
                node_elem_heads = h5_file[node_data_set][:,1]
                # append the nJm heads to a compined list
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
                # append feature type
                feature_type.append(nJm)
                # append the name of the features the element comes from
                element_name.append(feature)
                element_rank.append(feature_no)

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
            # append feature type
            feature_type.append(CC)
            # append the name of the features the element comes from
            element_name.append(np.repeat(feature,elem_lengths.shape[0]))
            element_rank.append(feature_no * np.ones(elem_lengths.shape[0]))

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
    element_rank    = np.hstack(element_rank)
    feature_length  = np.hstack(feature_length) * Yf
    feature_type    = np.hstack(feature_type)
    
    # find nelem
    nelem = element_length.shape[0]

    # find the x val to plot heads
    xval = np.zeros(nelem, dtype=np.float64)
    length_counter = 0

    for idx,item in enumerate(element_length):

        if element_type[idx] == nJ2:
            # xval is for the head plot 
            xval[idx] = length_counter

        elif element_type[idx] == nJm:
            # xval is for the head plot
            xval[idx] = length_counter + item * 0.5
            # advance for the next element        
            length_counter = length_counter + item 

        elif element_type[idx] == CC:
            # xval is for the head plot
            xval[idx] = length_counter + item * 0.5
            # advance for the next element        
            length_counter = length_counter + item 

    # set the location of link node labels
    xval_feature = np.zeros(feature_length.shape[0], dtype=np.float64) 
    length = 0


    for index,item in enumerate(feature_length):
        xval_feature[index] = length + item * 0.5
        length              = length + item


    # animation plot
    plt.rcParams['figure.figsize'] = [10, 4]
    plt.rcParams.update({'font.size': 11})
    fig, ax = plt.subplots()
    fig.tight_layout(pad=2)
    x = 0
    # take 10% of the gradient difference as buffer for plotting
    max_zrown = max(element_zcrown)
    max_head  = np.amax(element_head)
    if np.all(element_zbottom == 0):
        min_zbottom = 0
    else:
        min_zbottom = np.min(element_zbottom[np.nonzero(element_zbottom)])
    # buffer for plotting and labels
    buffer = 0.1 * (abs(max_head)-abs(min_zbottom))


    # transpormation for plotting
    trans = transforms.blended_transform_factory(
        ax.transData, ax.transAxes)
    # geometry profile in the animation

    for indx, feature in enumerate(profile):

        # plot the conduits first
        if feature_type[indx] == CC:
            # plot the z bottom
            ax.plot(xval[element_rank==indx], element_zbottom[element_rank==indx], '-k') 
            # plot the crown elevation
            ax.plot(xval[element_rank==indx], element_zcrown[element_rank==indx], '-k') 

            # calculations fro link labels
            feature_zbottom = 0.5 * (min(element_zbottom[element_rank==indx])
                                    +max(element_zbottom[element_rank==indx]))
            # # link labels
            # ax.text(xval_feature[indx],feature_zbottom-0.5*buffer,feature,
            #         rotation=90,ha='center',va='top', fontsize='small')

        elif feature_type[indx] == nJ2:
            # find the index of the J2 face
            J_idx = np.where(element_rank==indx)[0]
            if indx == 0:
                # store the index of the elem d/s 
                J_plot_idx = np.array(J_idx+1)
                vline_xval = xval[J_plot_idx]
            elif indx == profile.shape[0]-1:
                # store the index of the elem u/s 
                J_plot_idx = np.array(J_idx-1)
                vline_xval = xval[J_plot_idx]
            else:
                # store the index of the elem u/s, elem d/s 
                J_plot_idx = np.array([J_idx-1,J_idx+1])
                vline_xval = sum(xval[J_plot_idx]) / 2.

            # plot the zbottom of the junction by combining elem u/s, JM and elem d/s zbottom
            ax.plot(xval[J_plot_idx], element_zbottom[J_plot_idx],'-k')
            # plot the zcrown  of the junction by combining elem u/s, JM and elem d/s zbottom
            ax.plot(xval[J_plot_idx], element_zcrown[J_plot_idx],'-k')
            # plot vlines to distinguish nodes
            ax.axvline(x=vline_xval,color='b',linestyle='-',linewidth=0.15)

        # plot junction geometry
        elif feature_type[indx] == nJm:
            # find the index of the JM
            J_idx = np.where(element_rank==indx)[0]

            # JM at the start of the profile
            if indx == 0:
                # plot horizontoal lines to represent zbottom   
                plt.hlines(y=element_zbottom[J_idx],xmin=xval[J_idx],xmax=xval[J_idx+1], color='black')
                # plot vertical line to show the offsets
                plt.vlines(x = xval[J_idx+1], ymin = element_zbottom[J_idx],
                           ymax = element_zbottom[J_idx+1], color='black')

                # plot a vertical line to show the z crown at the d/s of the junction
                plt.vlines(x = xval[J_idx+1], ymin = element_zcrown[J_idx+1],
                           ymax = element_zcrown[J_idx], color='black')

                # plot another horizontal line to close up the JM
                plt.hlines(y=element_zcrown[J_idx],xmin=xval[J_idx],xmax=xval[J_idx+1], color='black')

                # plot vlines to distinguish nodes
                ax.axvline(x=xval[J_idx+1],color='b',linestyle='-',linewidth=0.15)

            # JM at the end of the profile
            elif indx == profile.shape[0]-1:
                # plot horizontoal lines to represent zbottom   
                plt.hlines(y=element_zbottom[J_idx],xmin=xval[J_idx-1],xmax=xval[J_idx], color='black')
                # plot vertical line to show the offsets
                plt.vlines(x = xval[J_idx-1], ymin = element_zbottom[J_idx],
                           ymax = element_zbottom[J_idx-1], color='black')

                # plot a vertical line to show the z crown at the u/s of the junction
                plt.vlines(x = xval[J_idx-1], ymin = element_zcrown[J_idx-1],
                           ymax = element_zcrown[J_idx], color='black')

                # plot another horizontal line to close up the JM
                plt.hlines(y=element_zcrown[J_idx],xmin=xval[J_idx-1],xmax=xval[J_idx], color='black')

                # plot vlines to distinguish nodes
                ax.axvline(x=xval[J_idx-1],color='b',linestyle='-',linewidth=0.15)

            # JM at the middle of the profile
            else:
                # plot horizontoal lines to represent zbottom   
                plt.hlines(y=element_zbottom[J_idx],xmin=xval[J_idx-1],xmax=xval[J_idx+1], color='black')
                # plot vertical line to show the offsets
                plt.vlines(x = xval[J_idx-1], ymin = element_zbottom[J_idx],
                           ymax = element_zbottom[J_idx-1], color='black')
                plt.vlines(x = xval[J_idx+1], ymin = element_zbottom[J_idx],
                           ymax = element_zbottom[J_idx+1], color='black')

                # plot a vertical line to show the z crown at the u/s of the junction
                plt.vlines(x = xval[J_idx-1], ymin = element_zcrown[J_idx-1],
                           ymax = element_zcrown[J_idx], color='black')
                # plot a vertical line to show the z crown at the d/s of the junction
                plt.vlines(x = xval[J_idx+1], ymin = element_zcrown[J_idx+1],
                           ymax = element_zcrown[J_idx], color='black')

                # plot another horizontal line to close up the JM
                plt.hlines(y=element_zcrown[J_idx],xmin=xval[J_idx-1],xmax=xval[J_idx+1], color='black')

                # plot vlines to distinguish nodes
                ax.axvline(x=xval[J_idx],color='b',linestyle='-',linewidth=0.15)

                # node labels
    #             ax.text(xval_feature[indx],feature_zbottom[indx]-buffer,feature_name[indx],
    #                 rotation=90,ha='center',va='top', fontsize='small')

    # animation line          
    line,  = ax.plot(xval, element_head[x,:], '-o', markersize=2.0, color='xkcd:cerulean')
    time_text = ax.text(0.98, 0.95,'',ha='right',va='top',transform=plt.gca().transAxes,fontsize='small')

    #labeling the plots, using profile name and units for length and head
    plt.title(profile_name)
    plt.xlabel('Length along the profile '+Yunit)
    plt.ylabel('Piezometric Head '+Yunit)
    plt.xlim(min(xval),max(xval))
    plt.ylim(min_zbottom-buffer,max_head+buffer)
    # plt.ylim(0,30)
    
    #this automatically helps make sure that the labels aren't cutoff and that the layout is correctly formated 
    plt.tight_layout()

    def animate(ii):
        line.set_ydata(element_head[x+ii,:])  # update the data.
        time_text.set_text('Time = %.1f hr.' %(sim_time[ii]))
        return line,


    ani = animation.FuncAnimation(fig, animate, frames = nTimeSteps, interval=10, blit=False)

    #saving the animation before showing it
    #MIGHT NEED TO BE CHANGED TO THE OUTPUT FOLDER RATHER THAN DUMBING TO CURRENT DIRECTORY
    if save_animation:
        # animation file name
        animation_name = output_path+""+test_case+"_"+ profile_name+'.gif' 
        writergif = animation.PillowWriter(fps=30)
        ani.save(animation_name,writer=writergif)


    # show the animation plot
   
    plt.show()
    plt.clf()
    plt.close()
    # break if a certian profile is animated
    if (has_arg_profile):
        break
else:
    "NO MORE PROFILES TO OUTPUT"