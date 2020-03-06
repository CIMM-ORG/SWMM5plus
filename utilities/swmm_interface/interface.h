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
#define FIDX 1

// Row in Nodes (Fortran)
enum node_attributes {
  e_ni_node_type,
  e_ni_N_link_u,
  e_ni_N_link_d,
  e_ni_Mlink_u1,
  e_ni_Mlink_u2,
  e_ni_Mlink_u3,
  e_ni_Mlink_d1,
  e_ni_Mlink_d2,
  e_ni_Mlink_d3,
  e_nr_Zbottom,
  num_node_attributes // this is always the last
};

// Row in Nodes (Fortran)
enum link_attributes {
  e_li_link_type,
  e_li_roughness_type,
  e_li_geometry,
  e_li_Mnode_u,
  e_li_Mnode_d,
  e_li_InitialDepthType,
  e_lr_Length,
  e_lr_BreadthScale,
  e_lr_Slope,
  e_lr_Roughness,
  e_lr_InitialFlowrate,
  e_lr_InitialDepth,
  e_lr_InitialUpstreamDepth,
  e_lr_InitialDnstreamDepth,
  num_link_attributes // this is always the last
};

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

int add_link (int li_idx, int ni_idx, int direction, float * ni_N_link_u, float * ni_Mlink_u1, float * ni_Mlink_u2, float * ni_Mlink_u3, float * ni_N_link_d, float * ni_Mlink_d1, float * ni_Mlink_d2, float * ni_Mlink_d3);
float get_node_type (int k, float total_n_links);
float get_link_xsect_attrs (int k, int attr, float length_units);
int init_node_tmp_table (float * ni_N_link_u, float * ni_Mlink_u1, float * ni_Mlink_u2, float * ni_Mlink_u3, float * ni_N_link_d, float * ni_Mlink_d1, float * ni_Mlink_d2, float * ni_Mlink_d3);
int get_node_attrs (int k, int length_units, float * attrs, float * ni_N_link_u, float * ni_Mlink_u1, float * ni_Mlink_u2, float * ni_Mlink_u3, float * ni_N_link_d, float * ni_Mlink_d1, float * ni_Mlink_d2, float * ni_Mlink_d3);
float get_link_attribute (int k, int attr, int units);
int DLLEXPORT num_links();
int DLLEXPORT num_nodes();
void DLLEXPORT record_nodes_data (double * nodeMatrix, int units);
void DLLEXPORT record_links_data (double * linkMatrix, int units);
void DLLEXPORT print_info (int);

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif