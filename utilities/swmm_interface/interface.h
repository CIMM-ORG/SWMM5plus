#include <stdlib.h>
#include <stdio.h>

#define UPSTREAM 0
#define DOWNSTREAM 1
#define SSIGN(X) (X > 0) - (X < 0) // SIGN function of X
#define CFTOCM(m) m*0.0283168466 // Cubic feet to cubic meters
#define FT2TOM2(m) m*0.09290304 // Square feet to square meters
#define FTTOM(m) m*0.3048 // Feet to meters
#define nullvalueI -998877

// FORTRAN engine constants
// data types for linkI(:,li_link_type)
#define L_CHANNEL 1
#define L_PIPE 2

// data types for linkI(:,li_geometry)
#define L_NOT_COMPATIBLE -1
#define L_RECTANGULAR 1
#define L_PARABOLIC 2
#define L_TRAPEZOIDAL 3
#define L_WIDTH_DEPTH 4

// data types for elemI(:,ei_roughness_type)
#define E_MANNINGS_N 1
#define E_CD 2

// data types for nodeI(:,ni_node_type)
#define NJ2 1
#define NJM 2
#define N_BC_DN 3
#define N_BC_UP 4

// Fortran starting index
// Fortran indexes start in FIDX, therefore
// k \in [FIDX, Nobjects[.]]
#define FIDX 1

// Row in Nodes (Fortran)
enum node_attributes {
  ni_node_type,
  ni_N_link_u,
  ni_N_link_d,
  ni_Mlink_u1,
  ni_Mlink_u2,
  ni_Mlink_u3,
  ni_Mlink_d1,
  ni_Mlink_d2,
  ni_Mlink_d3,
  nr_Zbottom,
  num_node_attributes // this is always the last
};

// Row in Nodes (Fortran)
enum link_attributes {
  li_link_type,
  li_roughness_type,
  li_geometry,
  li_Mnode_u,
  li_Mnode_d,
  li_InitialDepthType,
  lr_Length,
  lr_BreadthScale,
  lr_Slope,
  lr_Roughness,
  lr_InitialFlowrate,
  lr_InitialDepth,
  lr_InitialUpstreamDepth,
  lr_InitialDnstreamDepth,
  num_link_attributes // this is always the last
};

typedef struct
{
  int unit_system;
  int length_units;
  int flow_units;
  int manning_units;
} InterfaceUnits;

typedef struct
{
  float ** node_attributes;
  int num_nodes;
  int num_links;
  InterfaceUnits units;
} Interface;

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

// --- use "C" linkage for C++ programs

#ifdef __cplusplus
extern "C" {
#endif

int add_link (Interface* api, int li_idx, int ni_idx, int direction);
float get_node_type (int k, float total_n_links);
float get_link_xsect_attrs (int k, int attr, float length_units);
void* DLLEXPORT api_initialize (char* f1, char* f2, char* f3, int unit_system);
int DLLEXPORT api_finalize (void* fapi);
float DLLEXPORT api_get_node_attribute (void* fapi, int k, int attr);
float DLLEXPORT api_get_link_attribute (void* fapi, int k, int attr);
void DLLEXPORT api_print_info (void* fapi);
int DLLEXPORT api_num_links (void* fapi);
int DLLEXPORT api_num_nodes (void* fapi);
#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif