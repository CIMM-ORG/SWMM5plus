//=============================================================================
// FUNCTIONS FOR SWMM5+ ACCESS TO EPA-SWMM LID
// The following should be appended to the bottom of the lid.h file
// 20221104
//=============================================================================
void     lid_get_DrainInflow(int sIdx, double ff, int *nodeIdx, double *flowrate);
//EOF