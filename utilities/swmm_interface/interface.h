#ifndef INTERFACE_H_
#define INTERFACE_H_

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
#define NUM_API_VARS 2

enum api_node_attributes {
  node_ID = 1,
  node_type,
  node_invertElev,
  node_initDepth,
  node_extInflow_tSeries,
  node_extInflow_basePat,
  node_extInflow_baseline,
  node_depth,
  node_inflow,
  node_volume,
  node_overflow
};

enum api_link_attributes {
  link_ID = 1,
  link_subIndex,
  link_type,
  link_node1,
  link_node2,
  link_xsect_type,
  link_xsect_wMax,
  link_xsect_yBot,
  link_q0,
  link_geometry,
  conduit_roughness,
  conduit_length,
  link_flow,
  link_depth,
  link_volume,
  link_froude,
  link_setting,
  link_left_slope,
  link_right_slope
};

// API vars are those necessary for external applications
//   but have not been stored in the original SWMM data structures
//   These variables are found in the input file but are either
//   discarded or summarized
enum api_vars {
  api_left_slope,
  api_right_slope
};

typedef struct {
  int IsInitialized;
  double* vars[NUM_API_VARS];
} Interface;

// --- use "C" linkage for C++ programs

#ifdef __cplusplus
extern "C" {
#endif

void* DLLEXPORT api_initialize(char* f1, char* f2, char* f3);
void DLLEXPORT api_finalize(void* f_api);
double DLLEXPORT api_get_node_attribute (void* f_api, int k, int attr);
double DLLEXPORT api_get_link_attribute (void* f_api, int k, int attr);
int DLLEXPORT api_num_links ();
int DLLEXPORT api_num_nodes ();
int api_load_vars (void* f_api);
int  getTokens(char *s);

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif