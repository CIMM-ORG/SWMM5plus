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

time_now = str(datetime.now())
time_now = time_now.replace(' ', '_')
batch_testing_csv = 'test.csv'
output_dir = ""
has_output_path = False

#reading in to check if an output path for the batch output is supplied 
if((len(sys.argv)%2) != 0):

    for arg_id in range(1,len(sys.argv),2):
        arg = sys.argv[arg_id+1] 

        if(sys.argv[arg_id] == "-o"):
                print("inside output")
                output_path = arg
                has_output_path = True


#create the path to the output depending if an output was given 
if has_output_path:
    output_dir_timestamped = output_path+'/'+time_now+'/'
    os.system('mkdir '+output_path)
else:
    #setting the output directory
    output_path = "default_batch_comparison"
    output_dir_timestamped = output_path+'/'+time_now

#create the timestamped directory for the start of the batch comparison
os.system('mkdir '+ output_dir_timestamped)

#read through the csv of input files 
#each row is a different test case with the format
#<local_path_to_inp_file>,<local_path_to_json_setting_file>,<name_of_test>
#the last column of the csv is for the names of the test allow for you to run the same inp_file and name it something different 
with open(batch_testing_csv,'r' ,newline = '') as csv_file:
    csv_row_reader = csv.reader(csv_file)
    for row in csv_row_reader:

        inp_file  = row[0].strip()
        json_file = row[1].strip()
        test_case_name = row[2].strip()
        os.system('mkdir '+output_dir_timestamped+"/"+test_case_name.strip())
        #print("python comparison_script.py -i "+inp_file+" -s "+json_file)
        os.system("python compare.py -i "+inp_file+" -s "+json_file+ " -o "+ output_dir_timestamped+"/"+test_case_name.strip() +" -nc false " +" -b True")

         
