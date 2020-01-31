#include <stdlib.h>
#include <stdio.h>

// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
  #define WINDOWS
#endif
#ifdef __WIN32__
  #define WINDOWS
#endif

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
#define N_NOT_COMPATIBLE -1

int interface_print_info(int);