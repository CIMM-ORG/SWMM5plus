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

enum api_node_attributes {
  node_ID = 1,
  node_type,
  node_invertElev,
  node_initDepth,
  node_extInflow_tSeries,
  node_extInflow_basePat,
  node_extInflow_baseline
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
  conduit_length
};

typedef struct {
  int IsInitialized;
} Interface;

#define nullvalue -998877

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

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif