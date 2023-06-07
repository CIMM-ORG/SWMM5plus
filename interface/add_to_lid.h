// FILE add_to_lid.h
//==========================================================================
// SWMM5+ release, version 1.0.0
// 20230608
// Hydraulics engine that links with EPA SWMM-C
// June 8, 2023
//
// Description:
// FUNCTIONS FOR SWMM5+ ACCESS TO EPA-SWMM CONTROLS
// The following should be appended to the bottom of the lid.h file
// prior to compiling with SWMM5+
//
// Methods
// This function must be added for SWMM5+ to access the drain flow
// from LID that is computed in EPA-SWMM
//=============================================================================
void     lid_get_DrainInflow(int sIdx, double ff, int *nodeIdx, double *flowrate);
//=============================================================================
// END add_to_lid.h
//=============================================================================