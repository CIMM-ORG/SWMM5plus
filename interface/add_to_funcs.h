// FILE add_to_functions.c
//==========================================================================
// SWMM5+ release, version 1.0.0
// 20230608
// Hydraulics engine that links with EPA SWMM-C
// June 8, 2023
//
// Description:
// FUNCTIONS FOR SWMM5+ ACCESS TO EPA-SWMM CONTROLS
// The following should be appended to the bottom of the functions.h file
// prior to compiling with SWMM5+
//
// Methods
// Developed from controls.c algorithms in EPA-SWMM to return data needed
// for the SWMM5+ API to access control information. A controlled element
// in SWMM5+ on one processor might have to access a control location on
// another processor.
//=============================================================================
int     controls_display(void);
int     controls_count_rules(void); 
int     controls_count_premise(void); 
int     controls_count_thenAction(void); 
int     controls_count_elseAction(void); 
int     controls_get_premise_data(
                int* locationL,        int* locationR,
                int* linknodesimTypeL,          int* linknodesimTypeR,
                int* attributeL,       int* attributeR, 
                int* thisPremiseLevel, int rIdx);
int     controls_get_action_data(
                int* location, int* attribute, int* thisActionLevel, int rIdx, int isThen);

//=============================================================================
// END add_to_funcs.h
//=============================================================================