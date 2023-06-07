// FILE add_to_lid.c
//==========================================================================
// SWMM5+ release, version 1.0.0
// 20230608
// Hydraulics engine that links with EPA SWMM-C
// June 8, 2023
//
// Description:
// FUNCTIONS FOR SWMM5+ ACCESS TO EPA-SWMM CONTROLS
// The following should be appended to the bottom of the lid.c file
// prior to compiling with SWMM5+
//
// Methods
// This function must be added for SWMM5+ to access the drain flow
// from LID that is computed in EPA-SWMM
//=============================================================================
void  lid_get_DrainInflow(int sIdx, double ff, int *nodeIdx, double *flowrate)
//
//  Purpose: gets LID drain flow 
//      similar to lid_addDrainInflow but doesn't add flows. Used for SWMM5+
//  Author:  B R Hodges 20221104, derived from lid_addDrainInflow
//
//  Input:   sIdx = subcatchment index
//           ff = time interval weighting factor
//
//  Output:  nIdx = node index that LID is assigned to
//           flowrate = flowrate to conveyance nodes
//  
//  Note:    this function does NOT update anything in EPA SWMM
{        
    int nIdx;
    TLidUnit*  lidUnit;
    TLidList*  lidList;
    TLidGroup  lidGroup;

    // --- check if LID group exists
    lidGroup = LidGroups[sIdx];
    if ( lidGroup != NULL )
    {
        // --- zero the flowrate accumulator
        *flowrate = 0.0;
        *nodeIdx  = 0;
        
        // --- examine each LID in the group
        lidList = lidGroup->lidList;        
        while ( lidList )
        {
            // --- see if LID's drain discharges to conveyance system node
            lidUnit = lidList->lidUnit;
            nIdx    = lidUnit->drainNode;
            if ( nIdx >= 0 )
            {
                // ---increment node for fortran
                *nodeIdx = nIdx+1;
                // --- get drain flow and add it to the accumulator
                *flowrate += (1.0 - ff) * (lidUnit->oldDrainFlow) + ff * (lidUnit->newDrainFlow);
            }
            lidList = lidList->nextLidUnit;
        }
    }
}
//=============================================================================
// END add_to_lid.c
//=============================================================================