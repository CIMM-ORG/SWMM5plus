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
#define MTOFT(m) m/0.3048  // meters to feet
#define CMTOCFT(cm) cm / 0.0283168466  // cubic meters to cubic feet
#define API_NULL_VALUE_I -998877
#define NUM_API_INT_VARS 0
#define NUM_API_DOUBLE_VARS 2
#define NUM_API_TABLES 1

//-----------------------------------------------------------------------------
//  SI Unit conversion factors
//-----------------------------------------------------------------------------
const double UCF_SI[10][2] = 
      {//  US      SI
      {141732.3,  3600000.0 },         // RAINFALL (in/hr, mm/hr --> m/sec)
      {39.3701,   1000.0    },         // RAINDEPTH (in, mm --> m)
      {3402000.0, 86400000.0},         // EVAPRATE (in/day, mm/day --> m/sec)
      {3.28084,   1.0       },         // LENGTH (ft, m --> m)
      {0.0002471, 1.0e-4    },         // LANDAREA (ac, ha --> m2)
      {35.3147,   1.0       },         // VOLUME (ft3, m3 --> m3)
      {0.62,      1.0       },         // WINDSPEED (mph, km/hr --> km/hr)
      {1.0,       1.8       },         // TEMPERATURE (deg F, deg C --> deg F)
      {2.203e-6,  1.0e-6    },         // MASS (lb, kg --> mg)
      {43560.0,   3048.0    }          // GWFLOW (cfs/ac, cms/ha --> ft/sec)
      };
const double QCF_SI[6] =                  // Flow Conversion Factors:
    {35.3147, 15850.37, 22.8245,          // cfs, gpm, mgd --> cms
     1.00000, 1000.000, 84.6000};         // cms, lps, mld --> cms

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
// in define_api_keys.f90
enum api_nodef_attributes {
  nodef_ID = 2,
  nodef_type,                       // 3
  nodef_outfall_type,               // 4
  nodef_outfall_idx,                // 5
  nodef_invertElev,                 // 6
  nodef_initDepth,                  // 7
  nodef_StorageConstant,            // 8
  nodef_StorageCoeff,               // 9
  nodef_StorageExponent,            // 10
  nodef_StorageCurveID,             // 11
  nodef_StorageFevap,               // 12
  nodef_extInflow_tSeries,          // 13
  nodef_extInflow_tSeries_x1,       // 14
  nodef_extInflow_tSeries_x2,       // 15
  nodef_extInflow_basePat_idx,      // 16
  nodef_extInflow_basePat_type,     // 17
  nodef_extInflow_baseline,         // 18
  nodef_extInflow_sFactor,          // 19
  nodef_has_extInflow,              // 20
  nodef_dwfInflow_monthly_pattern,  // 21
  nodef_dwfInflow_daily_pattern,    // 22
  nodef_dwfInflow_hourly_pattern,   // 23
  nodef_dwfInflow_weekend_pattern,  // 24
  nodef_dwfInflow_avgvalue,         // 25
  nodef_has_dwfInflow,              // 26
  nodef_newDepth,                   // 27
  nodef_fullDepth,                  // 28
  nodef_surDepth,                   // 29
  nodef_inflow,                     // 30
  nodef_volume,                     // 31
  nodef_overflow,                   // 32
  nodef_pondedarea,                 // 33
  nodef_rptFlag,                    // 34
  nodef_hasFlapGate,                // 35
  nodef_RouteTo,                    // 36
  nodef_head_tSeries,               // 37
  nodef_head_tSeries_x1,            // 38
  nodef_head_tSeries_x2,            // 39
  nodef_has_extHead                 // 40
};
// skip 2 numbers for index end and start flags in 
//                               // 41
//                               // 42
// these "linkf" are identical to the fortran api_linkf_... values in define_api_keys.f90
enum api_linkf_attributes {
  linkf_ID = 43,                // 43  <======= COMMON ERROR POINT (mismatch between = and #)
  linkf_subIndex,               // 44 *
  linkf_direction,              // 45 *
  linkf_node1,                  // 46 *
  linkf_node2,                  // 47 *
  linkf_offset1,                // 48 *
  linkf_offset2,                // 49 *
  linkf_q0,                     // 50 *
  linkf_qlimit,                 // 51
  linkf_flow,                   // 52 *
  linkf_depth,                  // 53 *
  linkf_volume,                 // 54 *
  linkf_froude,                 // 55 *
  linkf_setting,                // 56 
  linkf_targetsetting,          // 57
  linkf_timelastset,            // 58 *
  linkf_left_slope,             // 59 *
  linkf_right_slope,            // 60 *
  linkf_weir_end_contractions,  // 61 *
  linkf_weir_side_slope,        // 62 *
  linkf_weir_road_width,        // 63 *
  linkf_weir_road_surface,      // 64 *
  linkf_curveid,                // 65 *
  linkf_discharge_coeff1,       // 66 *
  linkf_discharge_coeff2,       // 67 *
  linkf_initSetting,            // 68 *
  linkf_yOn,                    // 69 *
  linkf_yOff,                   // 70 *
  linkf_conduit_roughness,      // 71 *
  linkf_conduit_length,         // 72 *
  linkf_conduit_barrels,        // 73
  linkf_culvertCode,            // 74
  linkf_rptFlag,                // 75
  linkf_hasFlapGate,            // 76
  linkf_cLossInlet,             // 77
  linkf_cLossOutlet,            // 78
  linkf_cLossAvg,               // 79
  linkf_seepRate,               // 80
  linkf_commonBreak,            // 81
  // --- special elements attribute
  linkf_type,                   // 82 *
  linkf_sub_type,               // 83 *
  linkf_typeBreak,              // 84
  // --- xsect attributes
  linkf_xsect_type,         // 85 *
  linkf_geometry,           // 86 
  linkf_xsect_wMax,         // 87 *
  linkf_xsect_yBot,         // 88 *
  linkf_xsect_yFull,        // 89 *
  linkf_xsect_aFull,        // 90 *
  linkf_xsect_rFull,        // 91 *
  linkf_xsect_rBot,         // 92
  linkf_transectid,         // 93
  linkf_forcemain_coef     // 94
};
// skip 2 numbers for index start and end flags
// end flag                  // 95
// start flag                // 96
// these are identical to transect values in define_api_keys.f90
enum api_transectf_attributes {
  transectf_ID = 97,       // 97 <======= COMMON ERROR POINT (mismatch between = and #)
  transectf_yFull,         // 98
  transectf_aFull,         // 99
  transectf_rFull,         // 100
  transectf_wMax,          // 101
  transectf_ywMax,         // 102
  transectf_sMax,          // 103
  transectf_aMax,          // 104
  transectf_lengthFactor,  // 105
  transectf_roughness      // 106
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

// --- Controls

int DLLEXPORT api_teststuff();
int DLLEXPORT api_controls_count(int* nRules, int* nPremise, int* nThenAction, int* nElseAction);

int DLLEXPORT api_controls_get_premise_data(
    int* locationL,        int* locationR,
    int* linknodesimTypeL, int* linknodesimTypeR,
    int* attributeL,       int* attributeR, 
    int* thisPremiseLevel, int rIdx);

int DLLEXPORT api_controls_get_action_data(
    int* location,     
    int* attribute,
    int* thisActionLevel, int rIdx, int isThen);    
 
int DLLEXPORT api_controls_transfer_monitor_data(     
    double Depth, double Volume, double Inflow, double Flow, 
    double StatusSetting, double TimeLastSet, int LinkNodeIdx, int linknodesimType);

int DLLEXPORT api_controls_execute(
    double currentTimeEpoch, double ElapsedDays, double dtDays);

// --- Simulation

int DLLEXPORT api_initialize(char* f1, char* f2, char* f3, int run_routing);
int DLLEXPORT api_finalize();
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
int DLLEXPORT api_get_SWMM_times(
  double * starttime_epoch,
  double * endtime_epoch,
  double * report_start_datetime, 
  int    * report_step, 
  int    * hydrology_wet_step,
  int    * hydrology_dry_step,
  int    * sweep_start_dayofyear,
  int    * sweep_end_dayofyear,
  int    * dry_days,
  double * hydraulic_step,
  double * total_duration);

  double DLLEXPORT api_get_NewRunoffTime();
//brh20211208e

int DLLEXPORT api_get_SWMM_setup(
    int*    flow_units,
    int*    infiltration_model_type,
    int*    flow_routing_model_type,
    int*    link_offset_type,
    int*    force_main_equation,
    int*    ignore_rainfall,
    int*    ignore_snowmelt,
    int*    ignore_groundwater,
    int*    ignore_rdii,
    int*    ignore_routing,
    int*    ignore_quality,
    int*    allow_ponding,
    int*    steadystate_skip,
    double* steadystate_system_flow_tolerance,
    double* steadystate_lateral_flow_tolerance,
    double* routing_step_lengthening_time,
    double* routing_step_Courant_factor,
    double* routing_step_minimum, 
    int*    inertial_damping_type,
    int*    normal_flow_limiter_type,
    double* minimum_surface_area,
    double* mininum_conduit_slope,
    int*    maximum_number_of_trials,
    double* head_convergence_tolerance,
    int*    number_parallel_threads,
    int*    tempdir_provided,
    int*    control_rule_step,
    int*    surcharge_method
   );

int DLLEXPORT api_get_nodef_attribute(int node_idx, int attr, double* value);
int DLLEXPORT api_get_linkf_attribute(int link_idx, int attr, double* value);

int DLLEXPORT api_get_transectf_attribute(int transect_idx, int attr, double* value);
int DLLEXPORT api_get_N_TRANSECT_TBL();
int DLLEXPORT api_get_transect_table(
  int transect_idx, int table_len,
  double* tarea, double* twidth, double* thydradius);

int DLLEXPORT api_get_num_objects(int object_type);
int DLLEXPORT api_get_object_name(int object_idx, char* object_name, int object_type);
int DLLEXPORT api_get_object_name_len(int object_idx, int object_type, int* len);
int DLLEXPORT api_get_num_table_entries(int table_idx, int table_type, int* num_entries);
int DLLEXPORT api_get_table_attribute(int table_idx, int attr, double* value);
int DLLEXPORT api_get_first_entry_table(int table_idx, int table_type, double* x, double* y);
int DLLEXPORT api_get_next_entry_table(int table_idx, int table_type, double* x, double* y);
int DLLEXPORT api_get_next_entry_tseries(int tseries_idx, double timemax);

int DLLEXPORT api_reset_timeseries_to_start(int tseries_idx);

// Output fcns
int DLLEXPORT api_write_output_line(double t);
int DLLEXPORT api_update_nodeResult(int node_idx, int resultType, double newNodeResult);
int DLLEXPORT api_update_linkResult(int link_idx, int resultType, double newLinkResult);

// --- Printout
int add_link(int li_idx, 
              int ni_idx, 
              int direction, 
              int* ni_N_link_u, 
              int* ni_Mlink_u1, 
              int* ni_Mlink_u2, 
              int* ni_Mlink_u3, 
              int* ni_N_link_d, 
              int* ni_Mlink_d1, 
              int* ni_Mlink_d2, 
              int* ni_Mlink_d3);

int DLLEXPORT api_export_linknode_properties(int units);
int DLLEXPORT api_export_link_results(int link_idx);
int DLLEXPORT api_export_node_results(int node_idx);

// --- Hydrology 
int DLLEXPORT api_export_runon_volume(int kk, double flowrate);
int DLLEXPORT api_call_runoff_execute();
int DLLEXPORT api_get_subcatch_runoff(int id, double *runoff);
int DLLEXPORT api_getNumRdiiFlows(double thisDateTime, int *nRDII);
int DLLEXPORT api_getRdiiFlows(int rdiiIdx, int *nodeIdx, double *flowrate);

// --- Climate
int DLLEXPORT api_call_climate_setState(double thisDate);
int DLLEXPORT api_get_evaporation_rate(double *evapRate);

// --- Utils
int DLLEXPORT api_find_object(int object_type, char *id);
int check_api_is_initialized();
int api_load_vars();
int getTokens(char *s);

//=============================================================================
//   General purpose functions
//=============================================================================

double SI_Unit_Conversion(int u)
//
//  Input:   u = integer code of quantity being converted
//  Output:  returns a units conversion factor
//  Purpose: computes a conversion factor from SWMM's internal
//           units to user's units
//
{
    if ( u < FLOW ) return UCF_SI[u][UnitSystem];
    else            return QCF_SI[FlowUnits];
}

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif