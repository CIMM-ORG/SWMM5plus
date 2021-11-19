#ifndef API_H
#define API_H

// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
  #define WINDOWS
#endif
#ifdef __WIN32__
  #define WINDOWS
#endif

// --- define DLLEXPORT

#ifdef WINDOWS
    #define DLLEXPORT __declspec(dllexport) __stdcall
#else
    #define DLLEXPORT
#endif

#define UPSTREAM 0
#define DOWNSTREAM 1
#define SSIGN(X) (X > 0) - (X < 0)
#define CFTOCM(cf) cf*0.0283168466 // Cubic feet to cubic meters
#define FT2TOM2(sft) sft*0.09290304 // Square feet to square meters
#define FTTOM(ft) ft*0.3048 // Feet to meters
#define API_NULL_VALUE_I -998877
#define NUM_API_INT_VARS 0
#define NUM_API_DOUBLE_VARS 2
#define NUM_API_TABLES 1

// Enums written in caps are extracted from native
// EPA-SWMM, whereas lower case vars are added to EPA-SWMM

enum api_output_node_attribute {
  output_node_depth = 0,
  output_node_volume,
  output_node_latflow,
  output_node_inflow,
  MAX_API_OUTPUT_NODE_ATTR
};

enum api_output_link_attribute {
  output_link_depth = 0,
  output_link_flow,
  output_link_volume,
  output_link_direction,
  MAX_API_OUTPUT_LINK_ATTR
};

enum api_node_attributes {
  node_ID = 1,
  node_type,
  node_outfall_type,
  node_invertElev,
  node_initDepth,
  node_StorageConstant,
  node_StorageCoeff,
  node_StorageExponent,
  node_StorageCurveID,
  node_extInflow_tSeries,
  node_extInflow_tSeries_x1,
  node_extInflow_tSeries_x2,
  node_extInflow_basePat,
  node_extInflow_basePat_type,
  node_extInflow_baseline,
  node_extInflow_sFactor,
  node_has_extInflow,
  node_dwfInflow_monthly_pattern,
  node_dwfInflow_daily_pattern,
  node_dwfInflow_hourly_pattern,
  node_dwfInflow_weekend_pattern,
  node_dwfInflow_avgvalue,
  node_has_dwfInflow,
  node_depth,
  node_fullDepth,
  node_inflow,
  node_volume,
  node_overflow
};

enum api_link_attributes {
  link_ID = 1,
  link_subIndex,
  link_node1,
  link_node2,
  link_offset1,
  link_offset2,
  link_q0,
  link_flow,
  link_depth,
  link_volume,
  link_froude,
  link_setting,
  link_left_slope,
  link_right_slope,
  weir_end_contractions,
  discharge_coeff1,
  discharge_coeff2,
  conduit_roughness,
  conduit_length,
  // --- special elements attributes
  link_type,
  weir_type,
  orifice_type,
  pump_type,
  // --- xsect attributes
  link_xsect_type,
  link_geometry,
  link_xsect_wMax,
  link_xsect_yBot,
  link_xsect_yFull,
};

// API vars are those necessary for external applications
//   but have not been stored in the original SWMM data structures
//   These variables are found in the input file but are either
//   discarded or summarized

// Number of objects computed for
// interface purposes (starts in 1000)
// # define API_START_INDEX 1000
// enum api_num_objects {
// = API_START_INDEX;
// };

// Temporal variables for interface
// purposes
// enum api_int_vars {

// };

enum api_table_attributes {
  table_ID = 1,
  table_type,
  table_refers_to,
};

enum api_double_vars {
  api_left_slope,
  api_right_slope,
};

enum api_tables {
  api_time_series
};

typedef struct {
  int IsInitialized;
  double elapsedTime;
  // int* int_vars[NUM_API_INT_VARS];
  double* double_vars[NUM_API_DOUBLE_VARS];
} Interface;

// --- use "C" linkage for C++ programs

#ifdef __cplusplus
extern "C" {
#endif

// --- Simulation

int DLLEXPORT api_initialize(char* f1, char* f2, char* f3, int run_routing);
void DLLEXPORT api_finalize();
double DLLEXPORT api_run_step();

// --- Property-extraction

// * During Simulation

int DLLEXPORT api_get_node_results(char* node_name, float* inflow, float* overflow, float* depth, float* volume);
int DLLEXPORT api_get_link_results(char* link_name, float* flow, float* depth, float* volume);

// * After Initialization
double DLLEXPORT api_get_start_datetime();
double DLLEXPORT api_get_end_datetime();
int DLLEXPORT api_get_flowBC(int node_idx, double current_datetime, double* flowBC);
int DLLEXPORT api_get_headBC(int node_idx, double current_datetime, double* headBC);
int DLLEXPORT api_get_report_times(double * report_start_datetime, int * report_step, int * hydrology_step);
int DLLEXPORT api_get_node_attribute(int node_idx, int attr, double* value);
int DLLEXPORT api_get_link_attribute(int link_idx, int attr, double* value);
int DLLEXPORT api_get_num_objects(int object_type);
int DLLEXPORT api_get_object_name(int object_idx, char* object_name, int object_type);
int DLLEXPORT api_get_object_name_len(int object_idx, int object_type, int* len);
int DLLEXPORT api_get_num_table_entries(int table_idx, int table_type, int* num_entries);
int DLLEXPORT api_get_table_attribute(int table_idx, int attr, double* value);
int DLLEXPORT api_get_first_entry_table(int table_idx, int table_type, double* x, double* y);
int DLLEXPORT api_get_next_entry_table(int table_idx, int table_type, double* x, double* y);
int DLLEXPORT api_get_next_entry_tseries(int tseries_idx);

// Output fcns
int DLLEXPORT api_write_output_line(double t);
int DLLEXPORT api_update_nodeResult(int node_idx, int resultType, double newNodeResult);
int DLLEXPORT api_update_linkResult(int link_idx, int resultType, double newLinkResult);

// --- Print-out
int add_link(int li_idx, int ni_idx, int direction, int* ni_N_link_u, int* ni_Mlink_u1, int* ni_Mlink_u2, int* ni_Mlink_u3, int* ni_N_link_d, int* ni_Mlink_d1, int* ni_Mlink_d2, int* ni_Mlink_d3);
int DLLEXPORT api_export_linknode_properties(int units);
int DLLEXPORT api_export_link_results(int link_idx);
int DLLEXPORT api_export_node_results(int node_idx);

// --- Utils
int DLLEXPORT api_find_object(int object_type, char *id);
int check_api_is_initialized();
int api_load_vars();
int getTokens(char *s);

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif