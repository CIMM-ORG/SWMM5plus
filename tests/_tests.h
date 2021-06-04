
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

// --- use "C" linkage for C++ programs

#ifdef __cplusplus
extern "C" {
#endif

void DLLEXPORT interface_print_inflow(char* node_name);
void DLLEXPORT interface_print_tseries_props(int k);
void DLLEXPORT print_tseries_name(int k);

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif