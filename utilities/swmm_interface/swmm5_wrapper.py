import math, re # Used to create .rpt and .out paths
import six
import platform

from ctypes import c_double, CDLL, c_float, pointer, c_char_p # Required to handle with DLL variable
from os import sep
from time import time # Required to get computational times.

# ------------------------ CONSTANTS ---------------------------

# Types of objects
JUNCTION = 0
SUBCATCH = 1
NODE = 2
LINK = 3
STORAGE = 4
ORIFICE = 414
OUTFALL = 417

# Unit system
US = 0
SI = 1
DIMENSIONLESS = 0

# Attributes
DEPTH = 200			# [LINK] && [NODE]
VOLUME = 201		# [LINK] && [NODE]
FLOW = 202			# [LINK]
SETTING = 203		# [LINK]
FROUDE = 204		# [LINK]
INFLOW = 205		# [NODE]
FLOODING = 206		# [NODE]
PRECIPITATION = 207	# [SUBCATCHMENT]
RUNOFF = 208		# [SUBCATCHMENT]
LINK_AREA = 209		# [LINK]

# Start options
NO_REPORT = 0
WRITE_REPORT = 1

# Input file constants
INVERT = 400		# [NODE]
DEPTH_SIZE = 401	# [LINK] [NODE]
STORAGE_A = 402		# [NODE]
STORAGE_B = 403		# [NODE]
STORAGE_C = 404		# [NODE]
LENGTH = 405		# [LINK]
ROUGHNESS = 406		# [LINK]
IN_OFFSET = 407		# [LINK]
OUT_OFFSET = 408	# [LINK]
AREA = 409			# [SUBCATCHMENTS]
IMPERV = 410		# [SUBCATCHMENTS]
WIDTH = 411			# [SUBCATCHMENTS]
SLOPE = 412			# [SUBCATCHMENTS]
OUTLET = 413		# [SUBCATCHMENTS]
FROM_NODE = 415		# [LINK]
TO_NODE = 416		# [LINK]

# ------------------- GLOBAL PRIVATE CONSTANTS -----------------

_ERROR_PATH = -300
_ERROR_ATR = -299
_ERROR_TYPE = -298
_ERROR_NFOUND = -297
_ERROR_INCOHERENT = -296
_ERROR_IS_NUMERIC = -295
_ERROR_INVALID_TOPOLOGY = -400
_ERROR_SYS = -500

_ERROR_MSG_NFOUND = AttributeError("Error: Object not found")
_ERROR_MSG_TYPE = AttributeError("Error: Type of object not compatible")
_ERROR_MSG_ATR = AttributeError("Error: Attribute not compatible")
_ERROR_MSG_PATH = AttributeError("Error: Incorrect file path")
_ERROR_MSG_INCOHERENT = TypeError("Error: Incoherent parameter")
_ERROR_MSG_IS_NUMERIC = TypeError("Error: This function just handle numerical attributes")
_ERROR_MSG_INVALID_TOPOLOGY = AttributeError("Error: Invalid network, a node can have 3 upstream and 3 downstream links at most")
_ERROR_MSG_SYS = SystemError("Error: The system crashed, check the API")

# ------------------- GLOBAL PRIVATE VARIABLES -----------------

_SWMM5 = None # C library
_elapsedTime = c_double(0.000001) # Elapsed time in decimal days
_ptrTime = pointer( _elapsedTime ) # Pointer to elapsed time
_start_time = time() # Simulation start time
_end_time = time() # Simulation end time
_type_constants = (JUNCTION, SUBCATCH, NODE, LINK, STORAGE, ORIFICE, OUTFALL,)
_unit_constants = (US, SI, DIMENSIONLESS, )
_attribute_constants = (DEPTH, VOLUME, FLOW, SETTING, FROUDE, INFLOW, FLOODING, PRECIPITATION, RUNOFF, LINK_AREA,)
_report_constants = (NO_REPORT, WRITE_REPORT,)
_input_file_constants = (INVERT, DEPTH_SIZE, STORAGE_A, STORAGE_B, STORAGE_C, LENGTH, ROUGHNESS, IN_OFFSET, OUT_OFFSET,
						AREA, IMPERV, WIDTH, SLOPE, OUTLET, FROM_NODE, TO_NODE,)

DEBUG = True

#############################################################################################
# SWMM lib extra functionality
#############################################################################################

def print_info(inp, units):

	'''
	Inputs:  inp (str) = Path to the input file .inp
			 units (int) = unit system (US, SI)
			 libpath (str) = path to 
	Outputs: None
	Purpose: creates CSV files with information for SWMM6 FORTRAN engine

	Creates two files:

		nodes_info.csv
			(int) n_left: number of nodes left in the list 
				(the first row tells the total number of nodes)
			(string) node_id: id of the node
			(int) ni_idx: index of the node in SWMM's data structure for nodes
			(int) ni_node_type: code for node type according to SWMM
				(0) JUNCTION
				(1) OUTFALL
				(2) STORAGE
				(3) DIVIDER
			(int) ni_N_link_u: number of links connected upstream to the node
			(int) ni_N_link_d: number of links connected downstream to the node
			(int) ni_Mlink_u1: id of first link connected upstream to the node
			(int) ni_Mlink_u2: id of second link connected upstream to the node
			(int) ni_Mlink_u3: id of thrid link connected upstream to the node
			(int) ni_Mlink_d1: id of first link connected downstream to the node
			(int) ni_Mlink_d2: id of second link connected downstream to the node
			(int) ni_Mlink_d3: id of third link connected downstream to the node

		links_info.csv
			(int) l_left: number of links left in the list 
			(string) link_id
			(int) li_idx
			(int) li_link_type
			(int) li_geometry
			(int) li_Mnode_u
			(int) li_Mnode_d
			(float) lr_Length: length of the link (0 if type != CONDUIT)
			(float) lr_Slope: average slope of the link, estimated with extreme points
			(float) lr_Roughness: Manning coefficient of the link (0 if type != CONDUIT)
			(float) lr_InitialFlowrate: initial flow rate
			(float) lr_InitialUpstreamDepth: initial upstream depth
			(float) lr_InitialDnstreamDepth: initial downstream depth
	'''

	global _SWMM5
	

	if platform.system() == 'Windows':
		libpath =  'src' + sep + 'VS2017-DLL' + sep + 'x64' + sep + 'Release' + sep + 'libswmm5.dll'
	else:
		libpath =  'src' + sep + 'libswmm5.so'

	# C library is loaded
	_SWMM5 = CDLL(libpath)
	open_file(inp)  # Step 1
	start(WRITE_REPORT)  # Step 2
	error = _SWMM5.swmm_printInfo(units)
	if (error == _ERROR_INVALID_TOPOLOGY):
		raise(_ERROR_MSG_INVALID_TOPOLOGY)
	elif (error == _ERROR_INCOHERENT):
		raise(_ERROR_MSG_INCOHERENT)
	elif (error == _ERROR_SYS):
		raise(_ERROR_MSG_SYS)
	if DEBUG:
		print ("\nPrinting FORTRAN file -  OK")



#############################################################################################
# SWMM lib default functionality
#############################################################################################

def open_file(inp, msg=False):

	'''
	Inputs:  inp (str) -> Path to the input file .inp
			 msg (Bool)-> Display message in the terminal if True.
	Outputs: None
	Purpose: opens the files required to run a SWMM simulation
	'''

	# Creates paths for the report and the output files
	rpt = inp.replace('.inp', '.rpt')
	out = inp.replace('.inp', '.out')

	error = _SWMM5.swmm_open(c_char_p(six.b(inp)), c_char_p(six.b(rpt)), c_char_p(six.b(out)))
	if (error != 0):
		raise _ERROR_MSG_PATH
	if DEBUG:
		print ("\nOpenning SWMM  -  OK")


def start(write_report, msg=False):

	'''
	Inputs:  write_report (int) -> swmm.py constant related to the write report
			 file option.
			 msg (Bool)-> Display message in the terminal if True.
	Outputs: None
	Purpose: starts a SWMM simulation. Raise Exception if there is an error.
	'''

	# Parameter Error
	if write_report not in _report_constants:
		raise _ERROR_MSG_INCOHERENT

	_start_time = time()
	error = _SWMM5.swmm_start(write_report)
	if (error != 0):
		raise SystemError ("Error %d occured during the initialization of the simulation" % error)
	if DEBUG:
		print ("\nInitializing SWMM  -  OK")


def run_step():

	'''
	Inputs:  None
	Outputs: None
	Purpose: advances the simulation by one routing time step. Raise Exception
			 if there is an error.
	'''

	error = _SWMM5.swmm_step(_ptrTime)
	if (error != 0):
		raise SystemError ("Error %d ocurred at time %.2f" % (error, _elapsedTime.value))


def end(msg=False):

	'''
	Inputs:  msg (Bool)-> Display message in the terminal if True.
	Outputs: None
	Purpose: ends a SWMM simulation. Raise Exception if there is an error.
	'''

	error = _SWMM5.swmm_end()

	if (error != 0):
		raise SystemError ("Error %d: The simulation can not be ended" % error)
	_end_time = time()
	if DEBUG:
		print( ("Correctly Ended in %.2f seconds!" % (_end_time - _start_time)))


def save_report(msg=False):

	'''
	Inputs:  msg (Bool)-> Display message in the terminal if True.
	Outputs: None
	Purpose: writes simulation results to report file. Raise Exception if there is an error.
	'''

	error = _SWMM5.swmm_report()

	if (error != 0):
		raise SystemError ("Error %d: The report file could not be written correctly" % error)
	if DEBUG:
		print( ("Report file correctly written!"))


def close(msg=False):

	'''
	Inputs:  msg (Bool)-> Display message in the terminal if True.
	Outputs: None
	Purpose: closes a SWMM project. Raise Exception if there is an error.
	'''

	error = _SWMM5.swmm_close()
	if (error != 0):
		raise SystemError ("Error %d: The file can not be closed correctly" % error)
	if DEBUG:
		print( ("Correctly Closed!"))


def get_mass_bal_error():

	'''
	Inputs: None.
	Outputs: _ (tuple) -> Values of the errors related to mass balance.
			 			  [0] -> Runoff error
			 			  [1] -> Flow error
			 			  [2] -> Quality error
	Purpose: gets the mass balance errors of the simulation.
	'''

	runOffErr = c_float(0)
	flowErr = c_float(0)
	qualErr = c_float(0)
	ptrRunoff = pointer(runOffErr)
	ptrFlow = pointer(flowErr)
	ptrQual = pointer(qualErr)

	error = _SWMM5.swmm_getMassBalErr(ptrRunoff, ptrFlow, ptrQual)
	if (error != 0):
		raise SystemError ("Error %d: The errors can not be retrieved" % error)

	return (runOffErr.value, flowErr.value, qualErr.value)
