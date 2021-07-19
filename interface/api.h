#define UPSTREAM 0
#define DOWNSTREAM 1
#define SSIGN(X) (X > 0) - (X < 0)
#define CFTOCM(cf) cf*0.0283168466 // Cubic feet to cubic meters
#define FT2TOM2(sft) sft*0.09290304 // Square feet to square meters
#define FTTOM(ft) ft*0.3048 // Feet to meters
#define nullvalueI -998877

#ifndef INTERFACE_H
#define INTERFACE_H

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


#define nullvalue -998877022E8
#define NUM_API_INT_VARS 0
#define NUM_API_DOUBLE_VARS 2
#define NUM_API_TABLES 1

// Interface error codes:
#define ERROR_FEATURE_NOT_COMPATIBLE 100001

enum api_node_attributes {
  node_ID = 1,
  node_type,
  node_subtype,
  node_invertElev,
  node_initDepth,
  node_extInflow_tSeries,
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
  node_inflow,
  node_volume,
  node_overflow
};

enum api_link_attributes {
  link_ID = 1,
  link_subIndex,
  link_node1,
  link_node2,
  link_q0,
  link_flow,
  link_depth,
  link_volume,
  link_froude,
  link_setting,
  link_left_slope,
  link_right_slope,
  conduit_roughness,
  conduit_length,
  // --- xsect attributes
  link_type,
  link_xsect_type,
  link_geometry,
  link_xsect_wMax,
  link_xsect_yBot,
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

void* DLLEXPORT api_initialize(char* f1, char* f2, char* f3);
void DLLEXPORT api_finalize(void* f_api);

// --- Property-extraction

// * During Simulation

int DLLEXPORT api_get_node_results(void* f_api, char* node_name, float* inflow, float* overflow, float* depth, float* volume);
int DLLEXPORT api_get_link_results(void* f_api, char* link_name, float* flow, float* depth, float* volume);

// * After Initialization

int DLLEXPORT api_get_node_attribute(void* f_api, int k, int attr, double* value);
int DLLEXPORT api_get_link_attribute(void* f_api, int k, int attr, double* value);
int DLLEXPORT api_get_num_objects(void* f_api, int object_type);
int DLLEXPORT api_get_object_name(void* f_api, int k, char* object_name, int object_type);
double DLLEXPORT api_get_flowBC(void* f_api, int node_idx, double current_datetime);
double DLLEXPORT api_get_start_datetime();
double DLLEXPORT api_get_end_datetime();
int DLLEXPORT api_get_object_name_len(void* f_api, int k, int object_type);
int DLLEXPORT api_get_object_name(void* f_api, int k, char* object_name, int object_type);

// --- Print-out
void DLLEXPORT api_print_object_name(int k, int object_type);
int add_link(int li_idx, int ni_idx, int direction, int* ni_N_link_u, int* ni_Mlink_u1, int* ni_Mlink_u2, int* ni_Mlink_u3, int* ni_N_link_d, int* ni_Mlink_d1, int* ni_Mlink_d2, int* ni_Mlink_d3);
int DLLEXPORT interface_export_linknode_properties(void* f_api, int units);
int DLLEXPORT interface_export_link_results(void* f_api, char* link_name);
int DLLEXPORT interface_export_node_results(void* f_api, char* node_name);

// --- Utils
int check_api_is_initialized(Interface* api);
int api_load_vars(void * f_api);
int getTokens(char *s);
int api_findObject(int type, char *id);

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif