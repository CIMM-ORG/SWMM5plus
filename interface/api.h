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

// these "nodef" are identical to the fortran api_nodef_... values
enum api_nodef_attributes {
  nodef_ID = 1,
  nodef_type,            // 2
  nodef_outfall_type,    // 3
  nodef_invertElev,      // 4
  nodef_initDepth,       // 5
  nodef_StorageConstant,    // 6
  nodef_StorageCoeff,       // 7
  nodef_StorageExponent,    // 8
  nodef_StorageCurveID,     // 9
  nodef_extInflow_tSeries,   // 10
  nodef_extInflow_tSeries_x1,  // 11
  nodef_extInflow_tSeries_x2,   // 12
  nodef_extInflow_basePat,      // 13
  nodef_extInflow_basePat_type,  // 14
  nodef_extInflow_baseline,      // 15
  nodef_extInflow_sFactor,       // 16
  nodef_has_extInflow,             // 17
  nodef_dwfInflow_monthly_pattern,  // 18
  nodef_dwfInflow_daily_pattern,    // 19
  nodef_dwfInflow_hourly_pattern,   // 20
  nodef_dwfInflow_weekend_pattern,  // 21
  nodef_dwfInflow_avgvalue,         // 22
  nodef_has_dwfInflow,              // 23
  // brh20211207s
  //node_depth,                      // xx
  nodef_newDepth,                   // 24
  // brh20211207e
  nodef_fullDepth,                  // 25
  nodef_inflow,                     // 26
  nodef_volume,                     // 27
  // brh20211207s
  //node_overflow                    // 28
  nodef_overflow,                    // 28
  nodef_rptFlag                      // 29
};

// these "linkf" are identical to the fortran api_linkf_... values
enum api_linkf_attributes {
  linkf_ID = 1,
  linkf_subIndex,  // 2 *
  linkf_node1,     // 3 *
  linkf_node2,     // 4 *
  linkf_offset1,   // 5 *
  linkf_offset2,   // 6 *
  linkf_q0,        // 7 *
  linkf_flow,      // 8 *
  linkf_depth,     // 9 *
  linkf_volume,    // 10 *
  linkf_froude,    // 11 *
  linkf_setting,   // 12 *
  linkf_left_slope,        // 13 *
  linkf_right_slope,       // 14 *
  linkf_weir_end_contractions,  // 15 *
  linkf_weir_side_slope,        // 16 *
  linkf_curveid,           // 17 *
  linkf_discharge_coeff1,       // 18 *
  linkf_discharge_coeff2,       // 19 *
  linkf_conduit_roughness,      // 20 *
  linkf_conduit_length,         // 21 *
  // brh 20211207s
  linkf_rptFlag,           // 22 new in api.c
  // brh 20211207s
  // --- special elements attributes
  linkf_type,              // 23 *
  linkf_weir_type,              // 24 *
  linkf_orifice_type,           // 25 *
  linkf_outlet_type,            // 26 *
  linkf_pump_type,              // 27 *
  // --- xsect attributes
  linkf_xsect_type,        // 28 *
  linkf_geometry,          // 29 missing in api.c
  linkf_xsect_wMax,        // 30 *
  linkf_xsect_yBot,        // 31 *
  linkf_xsect_yFull        // 32 *
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
//brh20211208s
//int DLLEXPORT api_get_report_times(double * report_start_datetime, int * report_step, int * hydrology_step);
int DLLEXPORT api_get_report_times(
  double * report_start_datetime, 
  int * report_step, 
  int * hydrology_step,
  int * hydrology_dry_step,
  double * hydraulic_step);
  double DLLEXPORT api_get_NewRunoffTime();
//brh20211208e

int DLLEXPORT api_get_nodef_attribute(int node_idx, int attr, double* value);
int DLLEXPORT api_get_linkf_attribute(int link_idx, int attr, double* value);
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

// --- Hydrology 
int DLLEXPORT api_call_runoff_execute();
int DLLEXPORT api_get_subcatch_runoff(int id, double *runoff);

// --- Utils
int DLLEXPORT api_find_object(int object_type, char *id);
int check_api_is_initialized();
int api_load_vars();
int getTokens(char *s);



#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif