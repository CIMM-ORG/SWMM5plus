#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "swmm5.h"
#include "headers.h"
#include "api_error.h"
#include "lid.h"
#include "api.h"


//-----------------------------------------------------------------------------
//  Imported variables
//-----------------------------------------------------------------------------
#define REAL4 float

extern REAL4* NodeResults;             //  "
extern REAL4* LinkResults;             //  "
struct stat st = {0};

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------

Interface* api;
static char *Tok[MAXTOKS];             // String tokens from line of input
static int  Ntokens;                   // Number of tokens in line of input
char errmsg[1024];

//===============================================================================
// --- Simulation
//===============================================================================

// Initializes the EPA-SWMM simulation. It creates an Interface
// variable (details about Interface in interface.h). The
// function opens the SWMM input file and creates report and
// output files. Although the .inp is parsed within swmm_start,
// not every property is extracted, e.g., slopes of trapezoidal channels.
// In swmm.c
// to parse the .inp again using EPA-SWMM functionalities
// within api_load_vars
// The interface object is passed to many functions in interface.c
// but it is passed as a void pointer. This is because the object
// is also used by the Fortran engine. Treating Interface as void*
// facilitates interoperability
//===============================================================================
int DLLEXPORT api_controls_count(
    int* nRules, int* nPremise, int* nThenAction, int* nElseAction)
//===============================================================================
    ///
    /// Input: integers that will be overwritten
    /// Output: counts of various control rules
    /// Purpose: provide sizes for setting arrays in SWMM5+
{
    int ii;
    ii = controls_display();
    *nRules      = controls_count_rules();
    *nPremise    = controls_count_premise();
    *nThenAction = controls_count_thenAction();
    *nElseAction = controls_count_elseAction();

    return 0;
}

//===============================================================================
int DLLEXPORT api_controls_get_premise_data(
    int* locationL,        int* locationR,
    int* linknodesimTypeL, int* linknodesimTypeR,
    int* attributeL,       int* attributeR, 
    int* thisPremiseLevel, int rIdx)
//===============================================================================
    ///
    /// Input: values that will be ovewritten
    /// Output: values containging control premise information
    /// Purpose: gets data for the monitoring arrays for SWMM5+
{
    int success = 0;

    //printf(" rIDX here in api %d \n ",rIdx);

    success = controls_get_premise_data( 
                locationL,  locationR, 
                linknodesimTypeL,linknodesimTypeR, 
                attributeL, attributeR, 
                thisPremiseLevel,rIdx);         

    return success;
}

//===============================================================================
int DLLEXPORT api_controls_get_action_data(
    int* location,     
    int* attribute,
    int* thisActionLevel, int rIdx, int isThen)
//===============================================================================
    ///
    /// Input: values that will be ovewritten
    /// Output: values containging control action information
    /// Purpose: gets data for the action arrays for SWMM5+
{
    int success = 0;

    //printf(" rIDX here in api %d \n ",rIdx);

    success = controls_get_action_data( 
                location,  
                attribute,
                thisActionLevel,rIdx,isThen);

    return success;
}

//===============================================================================
int DLLEXPORT api_controls_transfer_monitor_data(     
    double Depth, double Volume, double Inflow, double Flow, 
    double StatusSetting, double TimeLastSet, int idx, int linknodesimType)
//===============================================================================
    ///
    /// Input: values from SWMM5+ finite-volume solution applied to 
    /// EPA-SWMM Link or Node with LinkNodeNum index.
    /// Purpose: set up EPA-SWMM link/nodes for control evaluation
{
    if (linknodesimType == 1)
    {
        Link[idx].newDepth = MTOFT(Depth);
        // Head is not valid for links
        // Volume is not valid for links
        // Inflow is not valid for links
        Link[idx].newFlow  = CMTOCFT(Flow);
        // note that for conduit and pumps, the following is the r_STATUS
        // while for orifice, weirs and outlets this is the r_SETTING
        // but they both get stored in the same place in the Link[].setting
        // HACK POSSIBLE PROBLEM IN EPA-SWMM, OUTLET not handled in controls.c/getVariableValue
        Link[idx].setting  = StatusSetting;
        Link[idx].timeLastSet = TimeLastSet;
    }
    else if (linknodesimType == 0)
    {
        Node[idx].newDepth = MTOFT(Depth);
        // Head is not stored for nodes in EPA-SWMM, so r_HEAD in controls
        // is handled by newDepth and stored invertElev
        Node[idx].newVolume  = CMTOCFT(Volume);
        Node[idx].newLatFlow = CMTOCFT(Inflow);
        // Flow is not valid for nodes
        // Status/Setting is not valid for nodes
        // TimeLastSet is not valid for nodes
    }
    else
    {
        // do nothing for simulation type
    }
    return 0;
}    

//===============================================================================
int DLLEXPORT api_controls_execute(
    double currentTimeEpoch, double ElapsedDays, double dtDays)
//===============================================================================
    /// Purpose: calls the controls_evaluate in EPA-SWMM
{
    int numactions = 0;

    numactions = controls_evaluate(currentTimeEpoch, ElapsedDays, dtDays); 

    return numactions; 
}

//===============================================================================
int DLLEXPORT api_teststuff()
//===============================================================================
{
    int ii, nPremise, nThenAction, nElseAction;

    //printf(" \n \n in api_teststuff \n \n ");

    ii = controls_display();

    nPremise    = controls_count_premise();
    nThenAction = controls_count_thenAction();
    nElseAction = controls_count_elseAction();

    //printf("\n %d \n ",nElseAction);

    //printf(" \n N premise, then, else locations: %d %d %d \n", nPremise, nThenAction, nElseAction);

    return 0;
}

//===============================================================================
int DLLEXPORT api_initialize(
    char* f1, char* f2, char* f3, int run_routing)
//===============================================================================
{
    int error;
    api = (Interface*) malloc(sizeof(Interface));
    error = swmm_open(f1, f2, f3);

    if (error) {
        api->IsInitialized = FALSE;
        return 0;
    }
    error = swmm_start(0);
    if (error) {
        api->IsInitialized = FALSE;
        return 0;
    }
    api->IsInitialized = TRUE;
    error = api_load_vars();
    if (!run_routing) IgnoreRouting = TRUE;

    return 0;
}

//===============================================================================
int DLLEXPORT api_finalize()
//===============================================================================
    //
    //  Input: f_api is an Interface object passed as a void*
    //  Output: None
    //  Purpose: Closes the link with the EPA-SWMM library
{
    int ii;

     //printf("\n in api finalize \n");

    swmm_end();

    //printf(" ------------------\n");
    //printf("Fout.mode %d, %d %d \n",Fout.mode,SCRATCH_FILE,SAVE_FILE);
    //printf(" ------------------\n");

    // NOTE -- because we always initialize with
    // an *.out file name, the Fout.mode is SAVE_FILE
    // so the swmm_report() is NOT called at the end
    //if ( Fout.mode == SCRATCH_FILE ) swmm_report();

    swmm_close();

    // frees double variables in API
    for (ii = 0; ii < NUM_API_DOUBLE_VARS; ii++)
    {
        if (api->double_vars[ii] != NULL)
            free(api->double_vars[ii]);
    }

    // // frees integer variables in API --- these do not exist brh20211217
    // for (i = 0; i < NUM_API_INT_VARS; i++)
    // {
    //     if (api->int_vars[i] != NULL)
    //         free(api->int_vars[i]);
    // }

    free(api);

    return 0;
}

//===============================================================================
double DLLEXPORT api_run_step()
//===============================================================================
{
    swmm_step(&(api->elapsedTime));

    return api->elapsedTime;
}

//===============================================================================
// --- Property-extraction
// * During Simulation
//===============================================================================

//===============================================================================
int DLLEXPORT api_get_node_results(
    char* node_name, float* inflow, float* overflow, float* depth, float* volume)
//===============================================================================
    //
    //  Input:    f_api = Interface object passed as a void*
    //            node_name = string identifier of node
    //            inflow, overflow, depth, volume =
    //  Output: None
    //  Purpose: Closes the link with the SWMM C library
    //
{
    int jj, error;

    error = check_api_is_initialized("api_get_node_results");
    if (error != 0) return error;

    jj = project_findObject(NODE, node_name);

    *inflow   = CFTOCM(Node[jj].inflow);
    *overflow = CFTOCM(Node[jj].overflow);
    *depth    = FTTOM (Node[jj].newDepth);
    *volume   = CFTOCM(Node[jj].newVolume);

    return 0;
}

//===============================================================================
int DLLEXPORT api_get_link_results(
    char* link_name, float* flow, float* depth, float* volume)
//===============================================================================
{
    int jj, error;

    error = check_api_is_initialized("api_get_link_results");
    if (error != 0) return error;

    jj = project_findObject(LINK, link_name);

    *flow   = CFTOCM(Link[jj].newFlow);
    *depth  = FTTOM (Link[jj].newDepth);
    *volume = CFTOCM(Link[jj].newVolume);

    return 0;
}

//===============================================================================
// --- Property-extraction
// * After Initialization
//===============================================================================

//===============================================================================
double DLLEXPORT api_get_start_datetime()
//===============================================================================
{
    return StartDateTime;
}

//===============================================================================
double DLLEXPORT api_get_end_datetime()
//===============================================================================
{
    return EndDateTime;
}

//===============================================================================
double DLLEXPORT api_get_TotalDuration()
//===============================================================================
// this is the SWMM-C Total Duration in milliseconds
{
    return TotalDuration;
}
//===============================================================================
int DLLEXPORT api_get_flowBC(
    int node_idx, double current_datetime, double* flowBC)
//===============================================================================
{
    int ptype, pcount, ii, jj, pp;
    int tyear, tmonth, tday;
    int thour, tmin, tsec, tweekday, attr;
    double val;
    double bline, sfactor;
    double total_factor = 1;
    double total_extinflow = 0;

   //printf("          starting api_get_flowBC-===========================- \n");
    //printf(" filename %s \n",Finflows.name);

    //jj = Node[node_idx].extInflow->tSeries; 
    //printf("          table-x1, table->x2 here %e %e \n", Tseries[jj].x1- 36526.0,Tseries[jj].x2- 36526.0);

    *flowBC = 0;
    datetime_decodeDate(current_datetime, &tyear, &tmonth, &tday);
    datetime_decodeTime(current_datetime, &thour, &tmin,   &tsec);
    tweekday = datetime_dayOfWeek(current_datetime);

    // handle dry weather inflows
    attr = nodef_has_dwfInflow;
    api_get_nodef_attribute(node_idx, attr, &val);
    
    //printf("          node idx %d \n ",node_idx);
    //printf("          dry weather flow val %g \n ",val);

    if (val == 1) { // node_has_dwfInflow 
        for(ii=0; ii<4; ii++)
        {
            // jj is the pattern #
            jj = Node[node_idx].dwfInflow->patterns[ii];
            if (jj > 0)
            {
                ptype = Pattern[jj].type;
                // if (ptype == MONTHLY_PATTERN)
                //     total_factor *= Pattern[jj].factor[tmonth-1];
                // else if (ptype == DAILY_PATTERN)
                //     total_factor *= Pattern[jj].factor[tweekday-1];
                // else if (ptype == HOURLY_PATTERN)
                //     total_factor *= Pattern[jj].factor[thour];
                // else if (ptype == WEEKEND_PATTERN)
                // {
                //     if ((tweekday == 1) || (tweekday == 7))
                //         total_factor *= Pattern[jj].factor[thour];
                // }
                switch(ptype) {
                    case MONTHLY_PATTERN :
                        total_factor *= Pattern[jj].factor[tmonth-1];
                        break;
                    case DAILY_PATTERN :
                        total_factor *= Pattern[jj].factor[tweekday-1];
                        break;
                    case HOURLY_PATTERN :
                        total_factor *= Pattern[jj].factor[thour];
                        break;
                    case WEEKEND_PATTERN :
                        if ((tweekday == 1) || (tweekday == 7))
                            total_factor *= Pattern[jj].factor[thour];
                        break;
                    default :
                        printf(" unexpected default reached in api_get_flowBC at A -- needs debugging");
                        return -1;
                }
            }
        }
        // convert EPA-SWMM cubic ft/s to cubic meters/s
        *flowBC += CFTOCM(total_factor * Node[node_idx].dwfInflow->avgValue);
        //printf("*flowBC %g \n",CFTOCM(total_factor * Node[node_idx].dwfInflow->avgValue));
    }

    // check to see if external inflow exists
    attr = nodef_has_extInflow;
    api_get_nodef_attribute(node_idx, attr, &val);
    //printf("          here node_idx %d \n",node_idx);
    //printf("          external inflow val = %f \n ",val);
    //printf("          MONTHLY %d \n ", MONTHLY_PATTERN);

    if (val == 1) { // node_has_extInflow
        // pp is the pattern #
        pp = Node[node_idx].extInflow->basePat; // pattern

        //printf("          base pat %d \n ",pp);
        //printf("          pattern numbers %d %d %d %d \n", MONTHLY_PATTERN, DAILY_PATTERN,HOURLY_PATTERN,WEEKEND_PATTERN);

        // convert EPA-SWMM cubic ft/s to cubic m/s
        bline = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->baseline); // baseline value

        //printf("          bline at top %f \n ",bline);

        if (pp >= 0)
        {
            ptype = Pattern[pp].type;
            // if (ptype == MONTHLY_PATTERN)
            // {
            //     //total_extinflow += Pattern[j].factor[mm-1] * bline;  //brh20211221
            //     total_extinflow += Pattern[pp].factor[tmonth-1] * bline;  //brh20211221
            // } 
            // else if (ptype == DAILY_PATTERN)
            //     //total_extinflow += Pattern[j].factor[dow-1] * bline; //brh20211221
            //     total_extinflow += Pattern[pp].factor[tweekday-1] * bline;//brh20211221
            // else if (ptype == HOURLY_PATTERN)
            //     //total_extinflow += Pattern[j].factor[h] * bline;  /brh20211221
            //     total_extinflow += Pattern[pp].factor[thour] * bline;   //brh20211221
            // else if (ptype == WEEKEND_PATTERN)
            // {
            //     if ((tweekday == 1) || (tweekday == 7))
            //         //total_extinflow += Pattern[j].factor[h] * bline;
            //         total_extinflow += Pattern[pp].factor[thour] * bline; //brh20211221
            // }
            switch(ptype) {
                case MONTHLY_PATTERN :
                    total_extinflow += Pattern[pp].factor[tmonth-1] * bline;
                    break;
                case DAILY_PATTERN :
                    total_extinflow += Pattern[pp].factor[tweekday-1] * bline;
                    break;
                case HOURLY_PATTERN :
                    total_extinflow += Pattern[pp].factor[thour] * bline;
                    break;
                case WEEKEND_PATTERN :
                    if ((tweekday == 1) || (tweekday == 7))
                        total_extinflow += Pattern[pp].factor[thour] * bline;
                    break;
                default :
                    printf(" unexpected default reached in api_get_flowBC at B -- needs debugging");
                    return -1;
            }
            //printf("total_extinflow at pattern %g \n ", total_extinflow);
        }
        else if (bline > 0)  // no pattern, but baseline inflow provided
        {
            total_extinflow += bline;  //already converted to cubic m/s
            //printf("total_extinflow (at bline) %g \n ", total_extinflow);
        }

        // external inflow from time series are added to baseline and pattern
        // jj is the time series
        jj = Node[node_idx].extInflow->tSeries; // tseries 
        
        //printf("          Node time series %d \n", jj);
        //printf("          totalextinflow = %g  \n ",total_extinflow);

        sfactor = Node[node_idx].extInflow->sFactor; // scale factor
        if (jj >= 0)
        {
            // Unit conversion for External Files
            // FlowUnitWords[]  = { w_CFS, w_GPM, w_MGD, w_CMS, w_LPS, w_MLD, NULL};
            // Corresponding FlowUnits are  0  ,    1 ,    2 ,    3 ,    4 ,    5 ,
            //printf("          totalextinflow = %g  \n ",total_extinflow);
            total_extinflow += (table_tseriesLookup(&Tseries[jj], current_datetime, FALSE) * sfactor) / QCF_SI[FlowUnits];
            //printf("          current datetime %g \n ",current_datetime - 36526.0);
            //printf("          CALLING table_tseriesLookup \n");
            //total_extinflow = table_tseriesLookup(&Tseries[jj], current_datetime, FALSE);
            //printf("          after table_tseriesLookup \n");
            //printf(" table lookup %g \n",table_tseriesLookup(&Tseries[jj], current_datetime, FALSE));
            //printf("          inflow = %g %g %g %d \n ",total_extinflow, sfactor,QCF_SI[FlowUnits], FlowUnits);
        }

        // add the external inflows to the dry weather flows stored in flowBC
         //printf("          flowBC = %g  \n ", total_extinflow);
        *flowBC += total_extinflow; // already converted to cubic m/s
        
    }
    //printf("          table-x1, table->x2 here %g %g \n", Tseries[jj].x1,Tseries[jj].x2);
    //printf("          leaving api_get_flowBC ==================================\n");
    //printf("\n");

    return 0;
}

//===============================================================================
int DLLEXPORT api_get_headBC(
    int node_idx, double current_datetime, double* headBC)
//===============================================================================
{
    double   xx, yy;                           // x,y values in a table
    int      kidx;                            // table index
    int      ii = Node[node_idx].subIndex;    // outfall index
    
    //     printf("  in api_get_headBC with outfall %d \n",Outfall[ii].type);
    //     printf("     node_idx = %d \n",node_idx);
    //     printf("     subIdx   = %d \n",ii);
    //     printf("     FIXED_OUTFALL = %d \n",FIXED_OUTFALL);
    //     printf("     NORMAL_OUTFALL = %d \n",NORMAL_OUTFALL);
    //     printf("     FREE_OUTFALL = %d \n",FREE_OUTFALL);
    //     printf("     TIDAL_OUTFALL = %d \n",TIDAL_OUTFALL);
    //    printf("     TIMESERIES_OUTFALL = %d \n",TIMESERIES_OUTFALL);
    
    // this mimics node.c/ outfall_setOutletDepth
    switch (Outfall[ii].type)
    { 
        case FREE_OUTFALL:
            // free outfall merely returns invert elevation, with yCrit and yNorm
            // computation made in BC code of SWMM5+
            *headBC = FTTOM(Node[node_idx].invertElev);
            return 0;

        case NORMAL_OUTFALL:
            // normal outfall merely returns invert elevation, with yNorm
            // computation made in BC code of SWMM5+
            *headBC = FTTOM(Node[node_idx].invertElev);
            return 0;    

        case FIXED_OUTFALL:
            *headBC = FTTOM(Outfall[ii].fixedStage);
            //printf(" here %d %g \n ", ii, Outfall[ii].fixedStage);
            return 0;

        case TIDAL_OUTFALL:
            // --- tidal outfall is always in terms of stage height
            kidx = Outfall[ii].tideCurve;
            table_getFirstEntry(&Curve[kidx], &xx, &yy);
            xx += (current_datetime - floor(current_datetime) ) * 24.0;
            *headBC = FTTOM(table_lookup(&Curve[kidx], xx )); 
            return 0;

        case TIMESERIES_OUTFALL:  
            // --- timeseries outfall is given as a stage (elevation) series 
            //     note that SWMM does NOT seem to automatically convert the data from its input 
            //     form into CFS; thus, we have to choose whether or not to convert based on 
            //     the FlowUnits, which stores the FLOW_UNITS in the *.inp file. 
            //     CFS, GPM and MGD have FlowUnits of 0,1,2;  CMS, LPS, MLD have FlowUnits of 3,4,5

            kidx = Outfall[ii].stageSeries;
            if (FlowUnits < MGD)
            {
                *headBC     = FTTOM(table_tseriesLookup(&Tseries[kidx], current_datetime, TRUE)); 
            }
            else
            {
                *headBC     = table_tseriesLookup(&Tseries[kidx], current_datetime, TRUE); 
            }
            return 0;
        
        default:
            *headBC = API_NULL_VALUE_I;
            sprintf(errmsg, "OUTFALL TYPE %d at NODE %s", Outfall[ii].type, Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_not_developed, errmsg);
            printf("Unexpected default case for Outfall[ii].type \n");
            return api_err_not_developed;
    }
}
//===============================================================================
int DLLEXPORT api_count_subobjects (int *N_groundwater)
//===============================================================================
{
    int error, ii;

    error = check_api_is_initialized("api_count_subobjects");
    if (error) return error;

    // number of subcatchments with groundwater
    *N_groundwater = 0;
    for (ii = 0; ii < Nobjects[SUBCATCH]; ii++)
    {
        if (Subcatch[ii].groundwater)
            *N_groundwater++ ;   
    }

    return 0;
}
//===============================================================================
int DLLEXPORT api_get_SWMM_setup(
    int*    flow_units,
    int*    infiltration_model_type,
    int*    flow_routing_model_type,
    int*    link_offset_type,
    int*    force_main_equation,
    int*    ignore_rainfall,
    int*    ignore_snowmelt,
    int*    ignore_groundwater,
    int*    ignore_rdii,
    int*    ignore_routing,
    int*    ignore_quality,
    int*    allow_ponding,
    int*    steadystate_skip,
    double* steadystate_system_flow_tolerance,
    double* steadystate_lateral_flow_tolerance,
    double* routing_step_lengthening_time,
    double* routing_step_Courant_factor,
    double* routing_step_minimum, 
    int*    inertial_damping_type,
    int*    normal_flow_limiter_type,
    double* minimum_surface_area,
    double* minimum_conduit_slope,
    int*    maximum_number_of_trials,
    double* head_convergence_tolerance,
    int*    number_parallel_threads,
    int*    tempdir_provided,
    int*    control_rule_step,
    int*    surcharge_method
    )

//===============================================================================
    // Gets all SWMM options from input file plus some derived options
{
    int error;

    error = check_api_is_initialized("api_get_SWMM_setup");
    if (error) return error;

    // OPTIONS from SWMM input file

    *flow_units = FlowUnits;

    *infiltration_model_type = InfilModel;

    *flow_routing_model_type = RouteModel;

    *link_offset_type = LinkOffsets;    

    *force_main_equation = ForceMainEqn;

    *ignore_rainfall    = IgnoreRainfall;
    *ignore_snowmelt    = IgnoreSnowmelt;
    *ignore_groundwater = IgnoreGwater;
    *ignore_rdii        = IgnoreRDII;
    *ignore_routing     = IgnoreRouting;
    *ignore_quality     = IgnoreQuality;

    *allow_ponding = AllowPonding;

    *steadystate_skip = SkipSteadyState;
    *steadystate_system_flow_tolerance  = CFTOCM(SysFlowTol);
    *steadystate_lateral_flow_tolerance = CFTOCM(LatFlowTol);

    *routing_step_lengthening_time = LengtheningStep;
    *routing_step_Courant_factor   = CourantFactor;
    *routing_step_minimum          = MinRouteStep;

    *inertial_damping_type    = InertDamping;
    *normal_flow_limiter_type = NormalFlowLtd;

    *minimum_surface_area  = FT2TOM2(MinSurfArea);
    *minimum_conduit_slope = MinSlope;

    *maximum_number_of_trials = MaxTrials;

    *head_convergence_tolerance = FTTOM(HeadTol);

    *number_parallel_threads = NumThreads;

    *tempdir_provided = 0;
    if (strlen(TempDir) >0) *tempdir_provided = 1;
    
    // OTHER SWMM setup values

    *control_rule_step = RuleStep;

    *surcharge_method = SurchargeMethod;

    //printf(" RouteModel = %d \n",RouteModel);

    //printf(" \n tempdir len is %ld \n ",strlen(TempDir));
    //printf(" courant factor %f \n", CourantFactor);
    //printf("\n testing variable %s \n",TempDir);

    //printf(" \n %s \n ",Title[0]);

    return 0;
}

//===============================================================================
int DLLEXPORT api_get_SWMM_times(
    double* starttime_epoch,
    double* endtime_epoch,
    double* report_start_datetime, 
    int*    report_step, 
    int*    hydrology_wet_step, 
    int*    hydrology_dry_step, 
    int*    sweep_start_dayofyear,
    int*    sweep_end_dayofyear,
    int*    dry_days,
    double* hydraulic_step,
    double* total_duration
    ) 
//===============================================================================    
{
    int error;

    error = check_api_is_initialized("api_get_SWMM_times");
    if (error) return error;

    *starttime_epoch       = StartDateTime;
    *endtime_epoch         = EndDateTime;
    *report_start_datetime = ReportStart;
    *report_step           = ReportStep;
    *hydrology_wet_step    = WetStep;
    *hydrology_dry_step    = DryStep;
    *sweep_start_dayofyear = SweepStart;
    *sweep_end_dayofyear   = SweepEnd;
    *dry_days              = StartDryDays;
    *hydraulic_step        = RouteStep;
    *total_duration        = TotalDuration / 1000.0;

    // printf(" report start datetime = %f \n",ReportStart);
    // printf(" report step = %d \n", ReportStep);
    // printf(" hydrology_step = %d \n",WetStep);
    // printf(" hydrology_dry_step = %d \n",DryStep);
    // printf(" hydraulic_step = %f \n",RouteStep);

    return 0;
}

//===============================================================================
double DLLEXPORT api_get_NewRunoffTime()
//===============================================================================
{
    return NewRunoffTime;
}

//===============================================================================
int DLLEXPORT api_get_nodef_attribute(
    int node_idx, int attr, double* value)
//===============================================================================
{
    int error, bpat, idx, tseries_idx;

    // printf("==== IN api_get_nodef_attribute %d \n",attr);
    // printf(" node index ",node_idx);

    error = check_api_is_initialized("api_get_nodef_attribute");
    if (error) return error;

    switch (attr) {

        case nodef_type :
            *value = Node[node_idx].type;
            break;

        case nodef_outfall_type  :
            switch (Node[node_idx].type) {
                case OUTFALL :
                    *value = Outfall[Node[node_idx].subIndex].type;
                    break;
                default :    
                    *value = API_NULL_VALUE_I;
                    sprintf(errmsg, "Extracting nodef_outfall_type for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                    api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                    return api_err_wrong_type;
            }
            break;   

        case nodef_outfall_idx :
            // DO NOT add 1 to this! SWMM5+ needs the EPA SWMM index
            switch (Node[node_idx].type) {
                case OUTFALL :             
                    *value = Node[node_idx].subIndex;
                    break;
                default :    
                    *value = API_NULL_VALUE_I;
                    sprintf(errmsg, "Extracting nodef_outfall_idx for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                    api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                    return api_err_wrong_type;
            }
            break;

        case nodef_hasFlapGate :
            switch (Node[node_idx].type) {
                case OUTFALL :             
                    *value = Outfall[Node[node_idx].subIndex].hasFlapGate;
                    break;
                default :    
                    *value = API_NULL_VALUE_I;
                    sprintf(errmsg, "Extracting nodef_hasFlapGate for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                    api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                    return api_err_wrong_type;
            }
            break;

        case nodef_RouteTo :  
            // NOTE +1 added to the *value in interface_get_nodef_attribute
            switch (Node[node_idx].type) {
                case OUTFALL :             
                    *value = Outfall[Node[node_idx].subIndex].routeTo;
                    //  printf("==== node index %d \n",node_idx);
                    //  printf("===== Node[node_idx].subIndex %d \n",Node[node_idx].subIndex);
                    //  printf("==== routeTo %e \n",*value);
                    break;
                default :    
                    *value = API_NULL_VALUE_I;
                    sprintf(errmsg, "Extracting nodef_RouteTo for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                    api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                    return api_err_wrong_type;
            }
            break;

        
        case nodef_head_tSeries :
            // NOTE +1 added to the *value in interface_get_nodef_attribute
            if (Node[node_idx].type == OUTFALL)   
            {
                // the outfall index
                idx = Node[node_idx].subIndex; 
                //printf(" here in API %d %d \n ",Outfall[idx].type, TIMESERIES_OUTFALL);
        
                switch (Outfall[idx].type) {
                    case TIMESERIES_OUTFALL:
                        //printf(" outfall type in timeseries %d \n ",Outfall[idx].type);
                         // the stage series
                        *value = Outfall[idx].stageSeries;
                        break;
                    case TIDAL_OUTFALL:
                        *value = Outfall[idx].tideCurve;
                        break;
                    default:
                        //printf(" outfall type in default %d \n ",Outfall[idx].type);
                        *value = API_NULL_VALUE_I;
                        sprintf(errmsg, "Attemptint to extract node_head_tSeries for NODE %s, which doesn't have an time series [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                        api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                        return api_err_wrong_type;
                }    
            }
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Attempting to extract node_head_tSeries for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_head_tSeries_x1  :
            if (Node[node_idx].type == OUTFALL)
            {
                // the outfall index
                idx = Node[node_idx].subIndex;
                tseries_idx = Outfall[idx].stageSeries;
                if (tseries_idx >= 0)        
                    *value = Tseries[tseries_idx].x1;
                else
                {
                    *value = API_NULL_VALUE_I;
                    sprintf(errmsg, "Extracting node_head_tSeries_x1 for NODE %s, which doesn't have a head timeseries [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                    api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                    return api_err_wrong_type;   
                }
                break;
            }    
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Attempting to extract node_head_tSeries_x1 for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_head_tSeries_x2  :
           if (Node[node_idx].type == OUTFALL)
            {
                // the outfall index
                idx = Node[node_idx].subIndex;
                tseries_idx = Outfall[idx].stageSeries;
                if (tseries_idx >= 0)        
                    *value = Tseries[tseries_idx].x2;
                else
                {
                    *value = API_NULL_VALUE_I;
                    sprintf(errmsg, "Extracting node_head_tSeries_x2 for NODE %s, which doesn't have a head timeseries [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                    api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                    return api_err_wrong_type;   
                }
                break;
            }    
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Attempting to extract node_head_tSeries_x2 for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;


        case nodef_invertElev  :
            *value = FTTOM(Node[node_idx].invertElev);
            break;

        case nodef_fullDepth  :
            *value = FTTOM(Node[node_idx].fullDepth);
            break;  

        case nodef_surDepth :
            *value = FTTOM(Node[node_idx].surDepth);
            break;        

        case nodef_initDepth  :
            switch (Node[node_idx].type) {
                case OUTFALL :
                    // get the head at the outfall
                    error = api_get_headBC(node_idx, StartDateTime, value);
                    // subtract the bottom elevation
                    *value -= FTTOM(Node[node_idx].invertElev);
                    // printf(" \n outfall head %e \n ",*value);
                    // printf(" \n %e \n ", FTTOM(Node[node_idx].invertElev));
                    // printf(" \n %e \n ",FTTOM(Node[node_idx].initDepth));
                    //printf(" \n return value %e \n ",*value);
                    if (error) return error;
                    break;
                default :
                    *value = FTTOM(Node[node_idx].initDepth);
            }
            break;

        case nodef_StorageConstant  :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].aConst;
                    break;
                default :
                    *value = -1;
            }
            break; 

        case nodef_StorageCoeff :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].aCoeff;
                    break;
                default :
                    *value = -1;
            }
            break;

        case nodef_StorageExponent  :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].aExpon;
                    break;
                default :
                    *value = -1;
            }
            break;   

        case nodef_StorageCurveID  :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].aCurve + 1;
                    break;
                default :
                    *value = -1;
            }
            break;

        case nodef_StorageFevap :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].fEvap;
                    break;
                default :
                    *value = -1;
            }
            break;

        case nodef_extInflow_tSeries :
            // NOTE +1 added to the *value in interface_get_nodef_attribute
            if (Node[node_idx].extInflow)
            {
                // printf(" in api_get_nodef_attribute node+idx, tseries %d %d \n", node_idx,Node[node_idx].extInflow->tSeries);
                *value = Node[node_idx].extInflow->tSeries;  
            }
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_extInflow_tSeries for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;   

        case nodef_extInflow_tSeries_x1  :
            tseries_idx = Node[node_idx].extInflow->tSeries;
            if (tseries_idx >= 0)
                *value = Tseries[tseries_idx].x1;
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_extInflow_tSeries_x1 for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_extInflow_tSeries_x2  :
            tseries_idx = Node[node_idx].extInflow->tSeries;
            if (tseries_idx >= 0)
                *value = Tseries[tseries_idx].x2;
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting tseries_idx for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;     

        case nodef_extInflow_basePat_idx  :
            // NOTE +1 added to the *value in interface_get_nodef_attribute
            if (Node[node_idx].extInflow)
            {
                *value = Node[node_idx].extInflow->basePat;
                //*value = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->basePat); // BRH20211221 HACK THIS IS WRONG
                //printf("%g \n",Node[node_idx].extInflow->cFactor);
                //printf("%d \n",Node[node_idx].extInflow->basePat);
            }    
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_extInflow_basePat for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_extInflow_basePat_type  :
            bpat = Node[node_idx].extInflow->basePat;
            //printf(" bpat %d",bpat);
            if (bpat >= 0) // baseline pattern exists
                *value = Pattern[bpat].type;
            else
            {
                *value = bpat;  // brh changed to bpat (-1) because API_NULL_VALUE_I does not have scope for where its needed
                //*value = API_NULL_VALUE_I;
                //printf("  bpat = %d \n", bpat);
                //printf("  location 3098705 problem with basePatType \n");
                // brh20211207s  commenting this so that it moves through with null result
                //sprintf(errmsg, "Extracting node_extInflow_basePat_type for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                //api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                //return api_err_wrong_type;
                // brh20211207e
            }
            break;   

        case nodef_extInflow_baseline :
            if (Node[node_idx].extInflow)
                *value = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->baseline);
            else
                *value = 0;
            break;

        case nodef_extInflow_sFactor  :
            if (Node[node_idx].extInflow)
                *value = Node[node_idx].extInflow->sFactor;
            else
                *value = 1;
            break;      

        case nodef_has_extInflow :
            if (Node[node_idx].extInflow)
                *value = 1;
            else
                *value = 0;
            break;   

        case nodef_dwfInflow_monthly_pattern  :
            if (Node[node_idx].dwfInflow)
                *value = Node[node_idx].dwfInflow->patterns[0];
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_dwfInflow_monthly_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_dwfInflow_daily_pattern :
            if (Node[node_idx].dwfInflow)
                *value = Node[node_idx].dwfInflow->patterns[1];
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_dwfInflow_daily_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;    

        case nodef_dwfInflow_hourly_pattern  :
            if (Node[node_idx].dwfInflow)
                *value = Node[node_idx].dwfInflow->patterns[2];
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_dwfInflow_hourly_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_dwfInflow_weekend_pattern  :
            if (Node[node_idx].dwfInflow)
                *value = Node[node_idx].dwfInflow->patterns[3];
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_dwfInflow_weekend_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;   

        case nodef_dwfInflow_avgvalue  :
            if (Node[node_idx].dwfInflow)
                *value = CFTOCM(Node[node_idx].dwfInflow->avgValue);
            else
                *value = 0;
            break;

        case nodef_has_dwfInflow  :
            if (Node[node_idx].dwfInflow)
                *value = 1;
            else
                *value = 0;
            break;      

        case nodef_newDepth  :
            *value = FTTOM(Node[node_idx].newDepth);
            break;

        case nodef_inflow  :
            *value = CFTOCM(Node[node_idx].inflow);
            break;   

        case nodef_volume  :
            *value = CFTOCM(Node[node_idx].newVolume);
            break;

        case nodef_overflow  :
            // THIS SHOULD NOT BE NEEDED
            *value = CFTOCM(Node[node_idx].overflow);
            break;     

        case nodef_pondedarea :
            *value = FT2TOM2(Node[node_idx].pondedArea);
            break;

        case nodef_rptFlag  :
            if (Node[node_idx].rptFlag)
                *value = 1;
            else
                *value = 0;
            break;

        default :
            printf(" ****** api_get_nodef_attribute called without supported attr at 3979874 %d ",attr);
            *value = API_NULL_VALUE_I;
    }
  
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_linkf_attribute(
    int link_idx, int attr, double* value)
//===============================================================================
{
    int error;

//   printf(" \n      ****** in api_get_linkf_attribute  %d \n  \n",attr);
  

    error = check_api_is_initialized("api_get_linkf_attribute");
    if (error) return error;

    // printf("\n in api_get_linkf_attribute %d \n \n",attr);

// the following are in the order of the enumeration in define_api_keys.f90 and api.h
    switch (attr) {

        case linkf_ID :
            *value = link_idx;
            break;

        case linkf_subIndex :
            *value = Link[link_idx].subIndex;
            break;
        
        case linkf_direction :
            *value = Link[link_idx].direction;
            break;

        case linkf_node1 :
            *value = Link[link_idx].node1;
            break;

        case linkf_node2 : 
            *value = Link[link_idx].node2;
            break;

        case linkf_offset1 :
            *value = FTTOM(Link[link_idx].offset1);
            break;

        case linkf_offset2 : 
            *value = FTTOM(Link[link_idx].offset2);
            break;

        case linkf_q0 :
            *value = CFTOCM(Link[link_idx].q0);
            break;    

        case linkf_qlimit :
            *value = CFTOCM(Link[link_idx].qLimit);
            break;            

        case linkf_flow :
            *value = CFTOCM(Link[link_idx].newFlow);
            break;

        case linkf_depth :
            *value = FTTOM(Link[link_idx].newDepth);
            break;

        case linkf_volume :
            *value = CFTOCM(Link[link_idx].newVolume);
            break;

        case linkf_froude :
            *value = Link[link_idx].froude;
            break;
        
        case linkf_setting :
            *value = Link[link_idx].setting;
            break;

        case linkf_targetsetting :
            *value = Link[link_idx].targetSetting;
            break;

        case linkf_timelastset :
            *value = Link[link_idx].timeLastSet;
            break;

        case linkf_left_slope :
            *value = api->double_vars[api_left_slope][link_idx];
            break;

        case linkf_right_slope :
            *value = api->double_vars[api_right_slope][link_idx];
            break;

        case linkf_weir_end_contractions :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].endCon;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_weir_side_slope :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].slope;
                    break;
                default :
                    *value = 0;
            }
            break;
        
        case linkf_weir_road_width :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = FTTOM(Weir[Link[link_idx].subIndex].roadWidth);
                    break;
                default :
                    *value = 0;
            }
            break;
        
        case linkf_weir_road_surface :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].roadSurface;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_curveid :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].cdCurve+1;
                    break;
                case PUMP :
                    *value = Pump[Link[link_idx].subIndex].pumpCurve+1;
                    break;
                case OUTLET :
                    *value = Outlet[Link[link_idx].subIndex].qCurve+1;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_discharge_coeff1 :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].cDisch1;
                    break;
                case ORIFICE :
                    *value = Orifice[Link[link_idx].subIndex].cDisch;
                    break;
                case OUTLET :
                    *value = Outlet[Link[link_idx].subIndex].qCoeff;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_discharge_coeff2 :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].cDisch2;
                    break;
                case OUTLET :
                    *value = Outlet[Link[link_idx].subIndex].qExpon;
                    break;
                case ORIFICE :
                    *value = Orifice[Link[link_idx].subIndex].orate;
                    break;
                default :
                    *value = 0;
            }
            break;
        
        case linkf_initSetting :
            switch (Link[link_idx].type) {
                case PUMP :
                    *value = Pump[Link[link_idx].subIndex].initSetting;
                    break;
                default :
                    *value = 0;
            }
            break;

       case linkf_yOn :
            switch (Link[link_idx].type) {
                case PUMP :
                    *value = Pump[Link[link_idx].subIndex].yOn;
                    break;
                default :
                    *value = 0;
            }
            break;
        
        case linkf_yOff :
            switch (Link[link_idx].type) {
                case PUMP :
                    *value = Pump[Link[link_idx].subIndex].yOff;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_conduit_roughness :
            switch (Link[link_idx].type) {
                case CONDUIT :
                    *value = Conduit[Link[link_idx].subIndex].roughness;
                    break;
                default :
                    *value = 0;
            }    
            break;

        case linkf_conduit_length : 
            switch (Link[link_idx].type) {
                case CONDUIT :
                    *value = FTTOM(Conduit[Link[link_idx].subIndex].length);
                    break;
                case ORIFICE :
                    *value = 0.01;
                    break;
                case WEIR :
                    *value = 0.01;
                case OUTLET :
                    *value = 0.01;
                    break;
                case PUMP :
                    *value = 0.01;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_conduit_barrels :
            *value = Conduit[Link[link_idx].subIndex].barrels;
            break;

        case linkf_culvertCode : 
            *value = Link[link_idx].xsect.culvertCode;
            break;   
    
        case linkf_rptFlag :
            if (Link[link_idx].rptFlag)
                *value = 1;
            else
                *value = 0; 
            break;  

        case linkf_hasFlapGate :
            if (Link[link_idx].hasFlapGate)    
                *value = 1;
            else
                *value = 0;
            break;

        case linkf_cLossInlet :
            *value = Link[link_idx].cLossInlet;
            break;

        case linkf_cLossOutlet :
            *value = Link[link_idx].cLossOutlet;
            break;

        case linkf_cLossAvg :
            *value = Link[link_idx].cLossAvg;
            break;

        case linkf_seepRate :
            *value = FTTOM(Link[link_idx].seepRate);
            //printf("\n ****** seep rate  %e \n ", FTTOM(Link[link_idx].seepRate));
            break;

        case linkf_commonBreak :
            // placeholder with no action
            *value = 0;
            break;

        case linkf_type : 
            *value = Link[link_idx].type;
            break;

        case linkf_sub_type :
            switch (Link[link_idx].type) {
                case CONDUIT :
                    *value = API_NULL_VALUE_I;
                    break;
                case ORIFICE :
                    *value = Orifice[Link[link_idx].subIndex].type;
                    break;
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].type;
                    break;
                case OUTLET :
                    *value = Outlet[Link[link_idx].subIndex].curveType;
                    break;
                case PUMP :
                    *value = Pump[Link[link_idx].subIndex].type;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_typeBreak :
            // placehoder with no action
            *value = 0;
            break;

        case linkf_xsect_type :
            *value = Link[link_idx].xsect.type;
            break;

        case linkf_geometry :
            printf(" ****** api_get_linkf_attribute called for unsupported attr = linkf_geometry at 2875 %d \n ",attr);
            break;

        case linkf_xsect_wMax :
            *value = FTTOM(Link[link_idx].xsect.wMax); 
            break;

        case linkf_xsect_yBot :
            *value = FTTOM(Link[link_idx].xsect.yBot);
            break;

        case linkf_xsect_yFull : 
            *value = FTTOM(Link[link_idx].xsect.yFull);
            break;
        
        case linkf_xsect_aFull : 
            *value = FT2TOM2(Link[link_idx].xsect.aFull);
            break;
        
        case linkf_xsect_rFull : 
            *value = FTTOM(Link[link_idx].xsect.rFull);
            break;
        
        case linkf_xsect_rBot : 
            *value = FTTOM(Link[link_idx].xsect.rBot);
            break; 

        case linkf_transectid :
            *value = Link[link_idx].xsect.transect;
            break;

        case linkf_forcemain_coef :
             switch ( ForceMainEqn )
            {
                case H_W:
                    // Link.xsect.rBot stores the input H-W coefficient
                    //printf(" ****** in api_get_linkf_attribute  %e \n ", Link[link_idx].xsect.rBot);
                    *value = Link[link_idx].xsect.rBot;
                    break;
                case D_W:
                    // Link.xsect.rBot stores the input D-W roughness divided by UCF(RAINDEPTH), either 12.0 (in) or 304.8 (mm)
                    // which is converting the input (inches or mm) into ft. Here we take the ft from EPA-SWM and covert to meters
                    //printf(" ****** in api_get_linkf_attribute  %e \n ", FTTOM(Link[link_idx].xsect.rBot));
                    *value = FTTOM(Link[link_idx].xsect.rBot);
                    break;
            }
            //printf(" ****** in api_get_linkf_attribute  %e \n ", Link[link_idx].xsect.rBot);
            //*value = Link[link_idx].xsect.rBot;  // works for H-W
            break;
 
        default :             
            printf(" ****** api_get_linkf_attribute called without supported attr at 837954 %d \n ",attr);
            *value = API_NULL_VALUE_I;              
    }

    // if (attr == linkf_subIndex)
    //     *value = Link[link_idx].subIndex;
    // else if (attr == linkf_type)
    //     *value = Link[link_idx].type;
    // else if (attr == linkf_node1)
    //     *value = Link[link_idx].node1;
    // else if (attr == linkf_node2)
    //     *value = Link[link_idx].node2;
    // else if (attr == linkf_offset1)
    //     *value = FTTOM(Link[link_idx].offset1);
    // else if (attr == linkf_offset2)
    //     *value = FTTOM(Link[link_idx].offset2);
    // else if (attr == linkf_xsect_type)
    //     *value = Link[link_idx].xsect.type;
    // else if (attr == linkf_xsect_wMax)
    //     *value = FTTOM(Link[link_idx].xsect.wMax);
    // else if (attr == linkf_xsect_yBot)
    //     *value = FTTOM(Link[link_idx].xsect.yBot);
    // else if (attr == linkf_xsect_yFull)
    //     *value = FTTOM(Link[link_idx].xsect.yFull);
    // else if (attr == linkf_q0)
    //     *value = CFTOCM(Link[link_idx].q0);
    // // brh20211207s duplicate ov above    
    // // rm else if (attr == linkf_type)
    // // rm    *value =  Link[link_idx].type;
    // // brh20211207e
    // else if (attr == linkf_pump_type)
    //     *value =  Pump[Link[link_idx].subIndex].type;
    // else if (attr == linkf_orifice_type)
    //     *value = Orifice[Link[link_idx].subIndex].type;
    // else if (attr == linkf_outlet_type)
    //     *value = Outlet[Link[link_idx].subIndex].curveType;
    // else if (attr == linkf_weir_type)
    //     *value = Weir[Link[link_idx].subIndex].type;
    // else if (attr == linkf_conduit_roughness)
    // {
    //     if (Link[link_idx].type == CONDUIT)
    //         *value = Conduit[Link[link_idx].subIndex].roughness;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_conduit_length)
    // {
    //     if (Link[link_idx].type == CONDUIT)
    //         *value = FTTOM(Conduit[Link[link_idx].subIndex].length);

    //     else if (Link[link_idx].type == ORIFICE)
    //         *value = 0.01;
    //     else if (Link[link_idx].type == WEIR)
    //         *value = 0.01;
    //     else if (Link[link_idx].type == OUTLET)
    //         *value = 0.01;
    //     else if (Link[link_idx].type == PUMP)
    //         *value = 0.01;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_weir_end_contractions)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].endCon;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_curveid)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].cdCurve+1;
    //     else if (Link[link_idx].type == PUMP)
    //         *value = Pump[Link[link_idx].subIndex].pumpCurve+1;
    //     else if (Link[link_idx].type == OUTLET)
    //         *value = Outlet[Link[link_idx].subIndex].qCurve+1;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_discharge_coeff1)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].cDisch1;
    //     else if (Link[link_idx].type == ORIFICE)
    //         *value = Orifice[Link[link_idx].subIndex].cDisch;
    //     else if (Link[link_idx].type == OUTLET)
    //         *value = Outlet[Link[link_idx].subIndex].qCoeff;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_discharge_coeff2)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].cDisch2;
    //     else if (Link[link_idx].type == OUTLET)
    //         *value = Outlet[Link[link_idx].subIndex].qExpon;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_weir_side_slope)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].slope;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_flow)
    //     *value = CFTOCM(Link[link_idx].newFlow);
    // else if (attr == linkf_depth)
    //     *value = FTTOM(Link[link_idx].newDepth);
    // else if (attr == linkf_volume)
    //     *value = CFTOCM(Link[link_idx].newVolume);
    // else if (attr == linkf_froude)
    //     *value = Link[link_idx].froude;
    // else if (attr == linkf_setting)
    //     *value = Link[link_idx].setting;
    // else if (attr == linkf_left_slope)
    //     *value = api->double_vars[api_left_slope][link_idx];
    // else if (attr == linkf_right_slope)
    //     *value = api->double_vars[api_right_slope][link_idx];
    // // brh 20211207s
    // else if (attr == linkf_geometry)  
    // {
    //     printf(" ****** api_get_linkf_attribute called for unsupported attr = linkf_geometry at 2875 %d ",attr); 
    // }    
    // else if (attr == linkf_rptFlag)
    // {
    //     if (Link[link_idx].rptFlag)
    //         *value = 1;
    //     else
    //         *value = 0;   
    // }        
    // // brh 20211207e    
    // else
    // {
    //     printf(" ****** api_get_linke_attribute called without supported attr at 837954 %d ",attr);
    //     *value = API_NULL_VALUE_I;
    // }    
    return 0;
}


//===============================================================================
int DLLEXPORT api_get_adjustments(
    int adj_len, double* adjTemperature, double* adjEvaporation, 
    double* adjRainfall, double* adjConductivity)
//===============================================================================
{
    int error;
    int ii;

    error = check_api_is_initialized("api_get_adjustments");
    if (error) return error;

    for(ii=0; ii< adj_len; ii++)
    {
        adjTemperature[ii]  = Adjust.temp[ii];
        adjEvaporation[ii]  = Adjust.evap[ii];
        adjRainfall[ii]     = Adjust.rain[ii];
        adjConductivity[ii] = Adjust.hydcon[ii];
    }

    return 0;

}

//===============================================================================
int DLLEXPORT api_get_transectf_attribute(
    int transect_idx, int attr, double* value)
//===============================================================================
{
    int error;

    error = check_api_is_initialized("api_get_transectf_attribute");
    if (error) return error;

    // the following are in the order of the enumeration in define_api_keys.f90 and api.h
    switch (attr) {
        case transectf_yFull :
            *value = FTTOM(Transect[transect_idx].yFull);
            break;
        case transectf_aFull :
            *value = FT2TOM2(Transect[transect_idx].aFull);
            break;
        case transectf_rFull :
            *value = FTTOM(Transect[transect_idx].rFull);    
            break;
        case transectf_wMax :
            *value = FTTOM(Transect[transect_idx].wMax);    
            break;
        case transectf_ywMax :
            // Note that EPA-SWMM does NOT store data in Transect[transect_idx].ywMax
            // See xsect_setIrreguXsectParams() where ywMax is computed for the cross
            // section but is not stored in the transect structure itself.
            // Unless this bug is fixed, we will overwrite this read init_transect_array()
            *value = FTTOM(Transect[transect_idx].ywMax);    
            break;
        case transectf_sMax :
            // units are ft^(4/3), so requires conversion to m^4/3
            *value = pow(FTTOM(pow(Transect[transect_idx].sMax,0.75)),4.0/3.0);    
            break;
        case transectf_aMax :
            *value = FT2TOM2(Transect[transect_idx].aMax);    
            break;
        case transectf_lengthFactor :
            *value = Transect[transect_idx].lengthFactor;    // non-dimensional
            break;
        case transectf_roughness :
            *value = Transect[transect_idx].roughness;    // non-dimensional
            break;
        default :
            printf(" ***** api_get_transect_attribute called without support attr at 239873 %d ",attr);   
    }
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_transect_table(
    int transect_idx, int table_len,
    double* tarea, double* twidth, double* thydradius)
//===============================================================================
    // obtains the area, depth, and hydraulic radius entries for a transect
{
    int error;
    int ii;

    error = check_api_is_initialized("api_get_transect_table");
    if (error) return error;

    // note that these do not need unit conversions because they
    // are normalized 0 to 1
    for(ii=0; ii<table_len; ii++)
         {
             tarea[ii]      = Transect[transect_idx].areaTbl [ii];
             twidth[ii]     = Transect[transect_idx].widthTbl[ii];
             thydradius[ii] = Transect[transect_idx].hradTbl [ii];
         }
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_N_TRANSECT_TBL()
//===============================================================================
    // this is the SWMM-C number of depth levels in transect table
{
    int error;

    error = check_api_is_initialized("api_get_N_TRANSECT_TBL");
    if (error) return error;

    return N_TRANSECT_TBL;
}

//===============================================================================
int DLLEXPORT api_get_num_objects(
    int object_type)
//===============================================================================
{
    int error;
    error = check_api_is_initialized("api_get_num_objects");
    if (error) return error;
    // if (object_type > API_START_INDEX) // Objects for API purposes
    //     return api->num_objects[object_type - API_START_INDEX];
    return Nobjects[object_type];
}

//===============================================================================
int DLLEXPORT api_get_object_name(
    int object_idx, char* object_name, int object_type)
//===============================================================================
{
    int error, ii;
    int obj_len = -1;

    // printf("in get_object_name object type %d \n ",object_type);

    error = check_api_is_initialized("api_get_object_name");
    if (error) return error;

    error = api_get_object_name_len(object_idx, object_type, &obj_len);
    if (error) return error;

    // switch (object_type) {
    //     case NODE :
    //         for(ii=0; ii<obj_len; ii++)
    //         {
    //             object_name[ii] = Node[object_idx].ID[ii];
    //         }
    //         break;
    //     case LINK :
    //         for(ii=0; ii<obj_len; ii++)
    //         {
    //             object_name[ii] = Link[object_idx].ID[ii];
    //         }
    //     default :
    //         strcpy(object_name, "");
    //         sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name]", object_type);
    //         api_report_writeErrorMsg(api_err_not_developed, errmsg);
    //         return api_err_not_developed;
    // }

   

    if (object_type == NODE)
    {
        for(ii=0; ii<obj_len; ii++)
        {
            object_name[ii] = Node[object_idx].ID[ii];
        }
    }
    else if (object_type == LINK)
    {
        for(ii=0; ii<obj_len; ii++)
        {
            object_name[ii] = Link[object_idx].ID[ii];
        }
    }
    else if (object_type == TRANSECT)
    {
        for(ii=0; ii<obj_len; ii++)
        {
            object_name[ii] = Transect[object_idx].ID[ii];
        }
    }    
    else
    {
        strcpy(object_name, "");
        sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name]", object_type);
        api_report_writeErrorMsg(api_err_not_developed, errmsg);
        return api_err_not_developed;
    }
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_object_name_len(
    int object_idx, int object_type, int* len)
//===============================================================================
{
    int error;

    error = check_api_is_initialized("api_get_object_name_len");
    if (error) {
        *len = API_NULL_VALUE_I;
        return error;
    }

    switch (object_type) {
        case NODE :
            *len = strlen(Node[object_idx].ID);
            return 0;
            break;
        case LINK :
            *len = strlen(Link[object_idx].ID);
            return 0;
            break;
        case TRANSECT :
            *len = strlen(Transect[object_idx].ID);
            return 0;
            break;
        default :
            *len = API_NULL_VALUE_I;
            sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name_len]", object_type);
            api_report_writeErrorMsg(api_err_not_developed, errmsg);
            return api_err_not_developed;
    }

    // if (object_type == NODE)
    // {
    //     *len = strlen(Node[object_idx].ID);
    //     return 0;
    // }
    // else if (object_type == LINK)
    // {
    //     *len = strlen(Link[object_idx].ID);
    //     return 0;
    // }
    // else
    // {
    //     *len = API_NULL_VALUE_I;
    //     sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name_len]", object_type);
    //     api_report_writeErrorMsg(api_err_not_developed, errmsg);
    //     return api_err_not_developed;
    // }
}

//===============================================================================
int DLLEXPORT api_get_num_table_entries(
    int table_idx, int table_type, int* num_entries)
//===============================================================================
{
    double xx, yy;
    int success;

    *num_entries = 0;
    // printf("1 Number of entries in Curve %d\n", *num_entries);

    switch (table_type) {
        case CURVE :
            // ERROR handling
            if (table_idx >= Nobjects[CURVE] || Nobjects[CURVE] == 0) return -1;
            success = table_getFirstEntry(&Curve[table_idx], &xx, &yy); // first values in the table
            (*num_entries)++;
            while (success)
            {
                success = table_getNextEntry(&(Curve[table_idx]), &xx, &yy);
                if (success) (*num_entries)++;
                // printf("0 Number of entries in Curve %d\n", *num_entries);
            }    
            break;
        default :
            return -1;
    }

    // if (table_type == CURVE)
    // {
    //     // ERROR handling
    //     if (table_idx >= Nobjects[CURVE] || Nobjects[CURVE] == 0) return -1;
    //     success = table_getFirstEntry(&Curve[table_idx], &x, &y); // first values in the table
    //     (*num_entries)++;
    //     while (success)
    //     {
    //         success = table_getNextEntry(&(Curve[table_idx]), &x, &y);
    //         if (success) (*num_entries)++;
    //         // printf("0 Number of entries in Curve %d\n", *num_entries);
    //     }
    // }
    // else
    // {
    //     return -1;
    // }
    // // printf("Number of entries in Curve %d\n", *num_entries);
    // // printf("SUCCESS %d\n", success);

    return 0;
}

//===============================================================================
int DLLEXPORT api_get_table_attribute(
    int table_idx, int attr, double* value)
//===============================================================================
{
    int error;

    error = check_api_is_initialized("api_get_table_attribute");
    if (error) return error;

    switch (attr) {
        case table_ID :
            *value = table_idx;
            break;
        case table_type :
            *value = Curve[table_idx].curveType;
            break;
        case table_refers_to :
            *value = Curve[table_idx].refersTo;
            break;
        default :
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "attr %d [api.c -> api_get_table_attribute]", attr);
            api_report_writeErrorMsg(api_err_not_developed, errmsg);
            return api_err_not_developed;
    }

    // if (attr == table_ID)
    //     *value = table_idx;
    // else if (attr == table_type)
    //     *value = Curve[table_idx].curveType;
    // else if (attr == table_refers_to)
    //     *value = Curve[table_idx].refersTo;
    // else
    // {
    //     *value = API_NULL_VALUE_I;
    //     sprintf(errmsg, "attr %d [api.c -> api_get_table_attribute]", attr);
    //     api_report_writeErrorMsg(api_err_not_developed, errmsg);
    //     return api_err_not_developed;
    // }
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_first_entry_table(
    int table_idx, int table_type, double* xx, double* yy)
//===============================================================================
{
    int success;

    switch (table_type) {
        case CURVE :
            success = table_getFirstEntry(&(Curve[table_idx]), xx, yy);
            // printf("...success, %d \n",success);
            // printf("...curveType, %d \n",Curve[table_idx].curveType);
            // unit conversion depending on the type of curve
            switch (Curve[table_idx].curveType) {
                case STORAGE_CURVE:
                    *xx /= SI_Unit_Conversion(LENGTH); 
                    *yy /= (SI_Unit_Conversion(LENGTH) * SI_Unit_Conversion(LENGTH));
                    break;
                case DIVERSION_CURVE:
                    *xx /= SI_Unit_Conversion(FLOW);
                    *yy /= SI_Unit_Conversion(FLOW);
                    break;
                case TIDAL_CURVE:
                    *yy /= SI_Unit_Conversion(LENGTH);
                    break;
                case RATING_CURVE:
                    *xx /= SI_Unit_Conversion(LENGTH);
                    *yy /= SI_Unit_Conversion(FLOW);
                    break;
                case SHAPE_CURVE:
                    *xx /= SI_Unit_Conversion(LENGTH);
                    *yy /= SI_Unit_Conversion(LENGTH);
                    break;
                case CONTROL_CURVE:
                    printf(" \n \n CONTROL CURVE CALLED FOR, NOT IMPLMEMENTED YET? \n \n");
                    break;
                case WEIR_CURVE:
                    *xx /= SI_Unit_Conversion(LENGTH);
                    break;
                case PUMP1_CURVE:
                    *xx /= SI_Unit_Conversion(VOLUME);
                    *yy /= SI_Unit_Conversion(FLOW);
                    break;
                case PUMP2_CURVE:
                    *xx /= SI_Unit_Conversion(LENGTH);
                    *yy /= SI_Unit_Conversion(FLOW);
                    break;
                case PUMP3_CURVE:
                    *xx /= SI_Unit_Conversion(LENGTH);
                    *yy /= SI_Unit_Conversion(FLOW);
                    break;
                case PUMP4_CURVE:
                    *xx /= SI_Unit_Conversion(LENGTH);
                    *yy /= SI_Unit_Conversion(FLOW);
                    break;
                default:
                    return 0;
            }
            break;
        case TSERIES :
            success = table_getFirstEntry(&(Tseries[table_idx]), xx, yy);
            break;          
        default :
            return 0;
    }
    // if (table_type == CURVE)
    //     success = table_getFirstEntry(&(Curve[table_idx]), x, y);
    // else if (table_type == TSERIES)
    //     success = table_getFirstEntry(&(Tseries[table_idx]), x, y);
    // else
    //     return 0;

    return success;
}

//===============================================================================
int DLLEXPORT api_get_next_entry_table(
    int table_idx, int table_type, double* xx, double* yy)
//===============================================================================
{
    int success;

    switch (table_type) {
        case TSERIES :
            success = table_getNextEntry(&(Tseries[table_idx]), &(Tseries[table_idx].x2), &(Tseries[table_idx].y2));
            if (success)
            {
                *xx = Tseries[table_idx].x2;
                *yy = Tseries[table_idx].y2;
            }
            break;
        case CURVE :
            success = table_getNextEntry(&(Curve[table_idx]), &(Curve[table_idx].x2), &(Curve[table_idx].y2));
            if (success)
            {
                    // unit conversion depending on the type of curve
            switch (Curve[table_idx].curveType) {
                case STORAGE_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(LENGTH);
                    *yy = Curve[table_idx].y2 * SI_Unit_Conversion(LENGTH) / SI_Unit_Conversion(VOLUME);
                    break;
                case DIVERSION_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(FLOW);
                    *yy = Curve[table_idx].y2 / SI_Unit_Conversion(FLOW);
                    break;
                case TIDAL_CURVE:
                    *xx = Curve[table_idx].x2;
                    *yy = Curve[table_idx].y2 / SI_Unit_Conversion(LENGTH);
                    break;
                case RATING_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(LENGTH);
                    *yy = Curve[table_idx].y2 / SI_Unit_Conversion(FLOW);
                    break;
                case CONTROL_CURVE:
                    *xx = Curve[table_idx].x2;
                    *yy = Curve[table_idx].y2;
                    break;
                case SHAPE_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(LENGTH);
                    *yy = Curve[table_idx].y2 / SI_Unit_Conversion(LENGTH);
                    break;
                case WEIR_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(LENGTH);
                    *yy = Curve[table_idx].y2;
                    break;
                case PUMP1_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(VOLUME);
                    *yy = Curve[table_idx].y2 / SI_Unit_Conversion(FLOW);
                    break;
                case PUMP2_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(LENGTH);
                    *yy = Curve[table_idx].y2 / SI_Unit_Conversion(FLOW);
                    break;break;
                case PUMP3_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(LENGTH);
                    *yy = Curve[table_idx].y2 / SI_Unit_Conversion(FLOW);
                    break;
                case PUMP4_CURVE:
                    *xx = Curve[table_idx].x2 / SI_Unit_Conversion(LENGTH);
                    *yy = Curve[table_idx].y2 / SI_Unit_Conversion(FLOW);
                    break;
                default:
                    *xx = Curve[table_idx].x2;
                    *yy = Curve[table_idx].y2;
                    break;
                }
            }
            break;
        default :
            return 0;
    }
    // if (table_type == TSERIES)
    // {
    //     success = table_getNextEntry(&(Tseries[table_idx]), &(Tseries[table_idx].x2), &(Tseries[table_idx].y2));
    //     if (success)
    //     {
    //         *x = Tseries[table_idx].x2;
    //         *y = Tseries[table_idx].y2;
    //     }
    // }
    // else if (table_type == CURVE)
    // {
    //     success = table_getNextEntry(&(Curve[table_idx]), &(Curve[table_idx].x2), &(Curve[table_idx].y2));
    //     if (success)
    //     {
    //         *x = Curve[table_idx].x2;
    //         *y = Curve[table_idx].y2;
    //     }
    // }

    return success;
}

//===============================================================================
int DLLEXPORT api_get_next_entry_tseries(
    int tseries_idx, double timemax)
//===============================================================================
{
    int success;
    double x2, y2;

    success = TRUE;
    // store the present upper values of time series x=time, y=value
    x2 = Tseries[tseries_idx].x2;
    y2 = Tseries[tseries_idx].y2;

    //printf("y2, x2, timemax %g %g %g \n ",y2, x2, timemax);

    // only get a new table entry if the present upper time (x2) is less than the maximum time.
    if (x2 < timemax)
    {
        success = table_getNextEntry(&(Tseries[tseries_idx]), &(Tseries[tseries_idx].x2), &(Tseries[tseries_idx].y2));
        // overwrite the x1,y1 with the old values for x2, y2
        // otherwise, no changes.
        if (success == TRUE)
        {
            Tseries[tseries_idx].x1 = x2;
            Tseries[tseries_idx].y1 = y2;
        }
    }

    // // only overwrite the time
    // if (success == TRUE)
    // {
    //     if (x2 < timemax) 
    //     {
    //         Tseries[tseries_idx].x1 = x2;
    //         Tseries[tseries_idx].y1 = y2;
    //     }      
    // }
    return success;
}

//===============================================================================
int DLLEXPORT api_reset_timeseries_to_start(
    int tseries_idx)
//===============================================================================
    //  
    //  Input: SWMM index of time series
    //  Output: TRUE if success in resetting to first entry
    //          FALSE if unsuccessful.   
{
    int success;
    double xx, yy;

    //printf("calling table_getFirstEntry \n");
    success = table_getFirstEntry(&Tseries[tseries_idx], &xx, &yy);
    Tseries[tseries_idx].x1 = xx;
    Tseries[tseries_idx].y1 = yy;
    //printf(" xx, yy %e %e  \n", xx - 36526.0,yy);
    //printf("Success %d \n",success);

    return success;
}

//===============================================================================
// --- Output Writing (Post Processing)
// * The follwing functions should only be executed after finishing
//   and writing SWMM5+ report files. The following functions are
//   meant to be called from Fortran in order to export .rpt and
//   .out files according to the SWMM 5.13 standard. Fortran-generated
//   report files are not manipulated here, the manipulation of
//   SWMM5+ report files is kept within the Fortran code to ensure
//   compatibility with future updates of the SWMM5+ standard
//===============================================================================

//===============================================================================
int DLLEXPORT api_write_output_line(
    double t)
//===============================================================================
    // t: elapsed time in seconds
{

    // --- check that simulation can proceed
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    // Update routing times to skip interpolation when saving results.
    OldRoutingTime = 0; NewRoutingTime = t*1000; // times in msec
    output_saveResults(t*1000);
    return 0;
}

//===============================================================================
int DLLEXPORT api_update_nodeResult(
    int node_idx, int resultType, double newNodeResult)
//===============================================================================
    // PRESENTLY NOT USED 20220521
{
    // WARNING -- this stores data from SWMM5+ into EPA-SWMM and
    // needs to have unit conversions added if it is to be used

    // --- check that simulation can proceed
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    if (resultType == output_node_depth)
        Node[node_idx].newDepth = newNodeResult;
    else if (resultType == output_node_volume)
        Node[node_idx].newVolume = newNodeResult;
    else if (resultType == output_node_latflow)
        Node[node_idx].newLatFlow = newNodeResult;
    else if (resultType == output_node_inflow)
        Node[node_idx].inflow = newNodeResult;
    else
    {
        sprintf(errmsg, "resultType %d [api.c -> api_update_nodeResult]", resultType);
        api_report_writeErrorMsg(api_err_not_developed, errmsg);
        return api_err_not_developed;
    }
    return 0;
}

//===============================================================================
int DLLEXPORT api_update_linkResult(
    int link_idx, int resultType, double newLinkResult)
//===============================================================================
    // PRESENTLY NOT USED 20220521
{

    // WARNING -- this stores data from SWMM5+ into EPA-SWMM and
    // needs to have unit conversions added if it is to be used

    // --- check that simulation can proceed
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    if (resultType == output_link_depth)
        Link[link_idx].newDepth = newLinkResult;
    else if (resultType == output_link_flow)
        Link[link_idx].newFlow = newLinkResult;
    else if (resultType == output_link_volume)
        Link[link_idx].newVolume = newLinkResult;
    else if (resultType == output_link_direction)
        Link[link_idx].direction = newLinkResult;
    else
    {
        sprintf(errmsg, "resultType %d [api.c -> api_update_linkResult]", resultType);
        api_report_writeErrorMsg(api_err_not_developed, errmsg);
        return api_err_not_developed;
    }
    return 0;
}

//===============================================================================
// --- Print-out
//===============================================================================

//===============================================================================
int DLLEXPORT api_export_linknode_properties(
    int units)
//===============================================================================
    // PRESENTLY NOT USED 20220521
{
    //  link
    int li_idx[Nobjects[LINK]];
    int li_link_type[Nobjects[LINK]];
    int li_geometry[Nobjects[LINK]];
    int li_Mnode_u[Nobjects[LINK]];
    int li_Mnode_d[Nobjects[LINK]];
    float lr_Length[Nobjects[LINK]];
    float lr_Slope[Nobjects[LINK]];
    float lr_Roughness[Nobjects[LINK]];
    float lr_FlowrateInitial[Nobjects[LINK]];
    float lr_InitialUpstreamDepth[Nobjects[LINK]];
    float lr_InitialDnstreamDepth[Nobjects[LINK]];
    int li_InitialDepthType[Nobjects[LINK]]; //
    float lr_BreadthScale[Nobjects[LINK]]; //
    //float lr_InitialDepth[Nobjects[LINK]]; //

    //  node
    int ni_idx[Nobjects[NODE]];
    int ni_node_type[Nobjects[NODE]];
    int ni_N_link_u[Nobjects[NODE]];
    int ni_N_link_d[Nobjects[NODE]];
    int ni_Mlink_u1[Nobjects[NODE]];
    int ni_Mlink_u2[Nobjects[NODE]];
    int ni_Mlink_u3[Nobjects[NODE]];
    int ni_Mlink_d1[Nobjects[NODE]];
    int ni_Mlink_d2[Nobjects[NODE]];
    int ni_Mlink_d3[Nobjects[NODE]];
    float nr_Zbottom[Nobjects[NODE]]; //

    int NNodes = Nobjects[NODE];
    int NLinks = Nobjects[LINK];

    float length_units;
    float manning_units;
    float flow_units;

    FILE *f_nodes;
    FILE *f_links;

    int i;
    int error;

    error = check_api_is_initialized("api_export_linknode_properties");
    if (error) return error;

    // Initialization
    for (i = 0; i<Nobjects[NODE]; i++) {
        ni_N_link_u[i] = 0;
        ni_N_link_d[i] = 0;
        ni_Mlink_u1[i] = API_NULL_VALUE_I;
        ni_Mlink_u2[i] = API_NULL_VALUE_I;
        ni_Mlink_u3[i] = API_NULL_VALUE_I;
        ni_Mlink_d1[i] = API_NULL_VALUE_I;
        ni_Mlink_d2[i] = API_NULL_VALUE_I;
        ni_Mlink_d3[i] = API_NULL_VALUE_I;
    }

    // Choosing unit system
    if (units == US)
    {
        flow_units = 1;
        manning_units = 1;
        length_units = 1;
    }
    else if (units == SI)
    {
        flow_units = M3perFT3;
        manning_units = pow(1/MperFT, 1/3);
        length_units = MperFT;
    }
    else
    {
        sprintf(errmsg, "Incorrect type of units %d [api.c -> api_export_linknode_properties]", units);
        api_report_writeErrorMsg(api_err_wrong_type, errmsg);
        return api_err_wrong_type;
    }

    // Links
    for (i=0; i<Nobjects[LINK]; i++) {
        int li_sub_idx;
        float h;

        li_idx[i] = i;
        li_link_type[i] = Link[i].type;
        li_geometry[i] = Link[i].xsect.type;

        li_Mnode_u[i] = Link[i].node1;
        error = add_link(i, li_Mnode_u[i], DOWNSTREAM, 
            ni_N_link_u, 
                ni_Mlink_u1, 
                ni_Mlink_u2, 
                ni_Mlink_u3, 
            ni_N_link_d, 
                ni_Mlink_d1, 
                ni_Mlink_d2, 
                ni_Mlink_d3);
        if (error) return error;

        li_Mnode_d[i] = Link[i].node2;
        error = add_link(i, li_Mnode_d[i], UPSTREAM, 
            ni_N_link_u, 
                ni_Mlink_u1, 
                ni_Mlink_u2, 
                ni_Mlink_u3, 
            ni_N_link_d, 
                ni_Mlink_d1, 
                ni_Mlink_d2, 
                ni_Mlink_d3);
        if (error) return error;

        li_sub_idx = Link[i].subIndex;
        // [li_InitialDepthType] This condition is associated to nodes in SWMM
        if (li_link_type[i] == CONDUIT) {
            lr_Length[i] = Conduit[li_sub_idx].length * length_units;
            lr_Roughness[i] = Conduit[li_sub_idx].roughness * manning_units;
            h = (Node[li_Mnode_u[i]].invertElev - Node[li_Mnode_d[i]].invertElev) * length_units;
            lr_Slope[i] = -SSIGN(h)*lr_Length[i]/(pow(lr_Length[i],2) - pow(fabsf(h),2));
        } else {
            lr_Length[i] = 0;
            lr_Roughness[i] = 0;
            lr_Slope[i] = 0;
        }

        lr_FlowrateInitial[i] = Link[i].q0 * flow_units;
        lr_InitialUpstreamDepth[i] = Node[li_Mnode_u[i]].initDepth * length_units;
        lr_InitialDnstreamDepth[i] = Node[li_Mnode_d[i]].initDepth * length_units; 

        // // check for a flap gate at an outfall node.
        // // if it exists, set the downstream depth to zero to be modified later
        // if (Node[li_Mnode_d[i]].type == OUTFALL) {
        //     lr_InitialDnstreamDepth[i] = 0.0;
        //         // if (Outfall[Node[li_Mnode_d[i]].subIndex].hasFlapGate) {
        //         //     lr_InitialDnstreamDepth[i] = 0.0;
        //         // } else {
        //         //     lr_InitialDnstreamDepth[i] = Node[li_Mnode_d[i]].initDepth * length_units;
        //         // }
        // } else {
        //     lr_InitialDnstreamDepth[i] = Node[li_Mnode_d[i]].initDepth * length_units;
        // }

        
    }

    // Nodes
    for (i=0; i<Nobjects[NODE]; i++) {
        ni_idx[i] = i;
        ni_node_type[i] = Node[i].type;
    }

    f_nodes = fopen("debug_input/node/nodes_info.csv", "w");
    f_links = fopen("debug_input/link/links_info.csv", "w");

    fprintf(f_nodes,
        "n_left,node_id,ni_idx,ni_node_type,ni_N_link_u,ni_N_link_d,ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,ni_Mlink_d2,ni_Mlink_d3\n");
    for (i=0; i<NNodes; i++) {
        fprintf(f_nodes, "%d,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
            NNodes-i,
            Node[i].ID,
            ni_idx[i],
            ni_node_type[i],
            ni_N_link_u[i],
            ni_N_link_d[i],
            ni_Mlink_u1[i],
            ni_Mlink_u2[i],
            ni_Mlink_u3[i],
            ni_Mlink_d1[i],
            ni_Mlink_d2[i],
            ni_Mlink_d3[i]);
    }
    fclose(f_nodes);

    fprintf(f_links,
        "l_left,link_id,li_idx,li_link_type,li_geometry,li_Mnode_u,li_Mnode_d,lr_Length,lr_Slope,lr_Roughness,lr_FlowrateInitial,lr_InitialUpstreamDepth,lr_InitialDnstreamDepth\n");
    for (i=0; i<NLinks; i++) {
        fprintf(f_links, "%d,%s,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
            NLinks-i,
            Link[i].ID,
            li_idx[i],
            li_link_type[i],
            li_geometry[i],
            li_Mnode_u[i],
            li_Mnode_d[i],
            lr_Length[i],
            lr_Slope[i],
            lr_Roughness[i],
            lr_FlowrateInitial[i],
            lr_InitialUpstreamDepth[i],
            lr_InitialDnstreamDepth[i]);
    }
    fclose(f_links);

    return 0;
}

//===============================================================================
int DLLEXPORT api_export_link_results(
    int link_idx)
//===============================================================================
    // PRESENTLY NOT USED 20220521
{
	FILE* tmp;
    DateTime days;
    int period;
    char theTime[20];
    char theDate[20];
	char path[50];
    int error;

    error = check_api_is_initialized("api_export_link_results");
    if (error) return error;

    /* File path writing */
    strcpy(path, "debug_output/swmm5/link/");
    strcat(path, Link[link_idx].ID); strcat(path, ".csv");
    tmp = fopen(path, "w");
    fprintf(tmp, "date,time,flow,velocity,depth,volume,capacity\n");

    for ( period = 1; period <= Nperiods; period++ )
    {
        output_readDateTime(period, &days);
        datetime_dateToStr(days, theDate);
        datetime_timeToStr(days, theTime);
        output_readLinkResults(period, link_idx);
        fprintf(tmp, "%10s,%8s,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            theDate,
            theTime,
            LinkResults[LINK_FLOW],
            LinkResults[LINK_VELOCITY],
            LinkResults[LINK_DEPTH],
            LinkResults[LINK_VOLUME],
            LinkResults[LINK_CAPACITY]);
    }
    fclose(tmp);

    return 0;
}

//===============================================================================
int DLLEXPORT api_export_node_results(
    int node_idx)
//===============================================================================
    // PRESENTLY NOT USED 20220521
{
	FILE* tmp;
    DateTime days;
    int period;
    char theTime[20];
    char theDate[20];
	char path[50];
    int error;

    error = check_api_is_initialized("api_export_node_results");
    if (error) return error;

    if (stat("NodeResults", &st) == -1) {
        mkdir("NodeResults", 0700);
    }

    /* File path writing */
    strcpy(path, "NodeResults/");
    strcat(path, Node[node_idx].ID);
    strcat(path, ".csv");
    tmp = fopen(path, "w");
    fprintf(tmp, "date,time,inflow,overflow,depth,volume\n");

    for ( period = 1; period <= Nperiods; period++ ) {
        output_readDateTime(period, &days);
        datetime_dateToStr(days, theDate);
        datetime_timeToStr(days, theTime);
        output_readNodeResults(period, node_idx);
        fprintf(tmp, "%10s,%8s,%.4f,%.4f,%.4f,%.4f\n",
            theDate,
            theTime,
            NodeResults[NODE_INFLOW],
            NodeResults[NODE_OVERFLOW],
            NodeResults[NODE_DEPTH],
            NodeResults[NODE_VOLUME]);
    }
    fclose(tmp);
    return 0;
}

//===============================================================================
// --- Utils
//===============================================================================

//===============================================================================
int DLLEXPORT api_find_object(
    int object_type, char *id)
//===============================================================================
{
    return project_findObject(object_type, id);
}


// -------------------------------------------------------------------------
// |
// |  Hydrology
// v
// -------------------------------------------------------------------------
//===============================================================================
int DLLEXPORT api_export_runon_volume(
    int outfall_idx, double volume)
//===============================================================================
    // exports the runon volume from SWMM5+ to EPA SWMM
{
    // printf(" \n in api_export_runon_volume \n");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    // convert the SWMM5+ m^3 flowrate to ft^3
    Outfall[outfall_idx].vRouted = CMTOCFT(volume);

    // printf("\n in api export_runon_volume \n ");
    // printf("   index %d  \n ",outfall_idx);
    // printf("   volume in cf  %e \n ", Outfall[outfall_idx].vRouted);
    // printf("   volume in cm  %e \n ", volume);
}
//===============================================================================
int DLLEXPORT api_call_runoff_execute()
//===============================================================================
    // calls the runoff_execute() procedure in SWMM-C
{
    // printf(" in api_call_runoff_execute");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    runoff_execute();
    
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_subcatch_runoff(
    int sc_idx, double *runoff)
//===============================================================================
{
    // printf(" in api_get_subcatch_runoff \n");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    // Get runoff and convert to cubic meters per second
    *runoff = CFTOCM(Subcatch[sc_idx].newRunoff);
    // printf("... sc_idx, newRunoff CMS %d , %f \n",sc_idx,Subcatch[sc_idx].newRunoff);
    
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_subcatch_runoff_nodeIdx(
    int sc_idx, int *node_idx)
//===============================================================================
{
    // printf(" in api_get_subcatch_runoff_nodeIdx \n");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    // Get node index
    *node_idx = Subcatch[sc_idx].outNode;

    //printf("... sc_idx, node_idx %d , %d \n",sc_idx,*node_idx);
    
    return 0;
}
//===============================================================================
int DLLEXPORT api_getNumRdiiFlows(
    double thisDateTime, int *nRDII)
//===============================================================================
    // calls the rdii_getNumRdiiFlows() procedure in SWMM-C
    // to get the count of the nodes with RDII inflows
{
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    *nRDII = rdii_getNumRdiiFlows(thisDateTime);

    return 0;
}
//===============================================================================
int DLLEXPORT api_getRdiiFlow(
    int rdiiIdx, int *nodeIdx, double *flowrate)
//===============================================================================
    // calls the rdii_getRdiiFlow() procedure in SWMM-C
{
    int nIdx;
    double fr;

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    // input RDII index is adjusted for C by -1
    //  note that function returns pointers
    rdii_getRdiiFlow(rdiiIdx-1, &nIdx, &fr);

    // adjust output node index for Fortran
    *nodeIdx = nIdx+1;

    // flowrate unit conversion from ft^3/s to m^3/s
    *flowrate = CFTOCM(fr);

    return 0;
}
//===============================================================================
int DLLEXPORT api_get_groundwaterFlow(
    double thisTime, double LastRunoffTime, double NextRunoffTime,
    int sIdx, int *nodeIdx, 
    double *flowrate)
//===============================================================================
    // gets the flowrate from groundwater at thisTime for the subcatchment
    // with sIdx that is into node at nodeIdx
{
    double ff;
    TGroundwater* gw;

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )

    // interpolate time between runoff
    ff = (thisTime - LastRunoffTime) / (NextRunoffTime - LastRunoffTime);
    if ( ff < 0.0 ) ff = 0.0;
    if ( ff > 1.0 ) ff = 1.0;

    gw = Subcatch[sIdx].groundwater;
    if ( gw )
    {
        
        *nodeIdx = gw->node;
        if (*nodeIdx >=0)
        {
            // increment node index for Fortran
            *nodeIdx++;
            *flowrate = CFTOCM( ((1.0-ff)*(gw->oldFlow) + ff*(gw->newFlow) )
                                 * Subcatch[sIdx].area );
        }
        
    }
    return 0;
}
//===============================================================================
int DLLEXPORT api_get_LID_DrainFlow(
    double thisTime, double LastRunoffTime, double NextRunoffTime,
    int sIdx, int *nodeIdx, double *flowrate)
//===============================================================================
    // gets the flowrate from LID at thisTime for the subcatchment
    // with sIdx that is into node at nodeIdx
    // Mimics routing.c/addLidDrainInflows and lid.c/lid_addDrainInflow
{
    double ff, qq;
    int nIdx;
    

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )

    // interpolate time between runoff
    ff = (thisTime - LastRunoffTime) / (NextRunoffTime - LastRunoffTime);
    if ( ff < 0.0 ) ff = 0.0;
    if ( ff > 1.0 ) ff = 1.0;

    // check that this subcatchment exists and has lid 
    if ( Subcatch[sIdx].area > 0.0 && Subcatch[sIdx].lidArea > 0.0 )
    {
        // get the inflow rate and node index
        // note that this function is an add-on to the lid.c functions
        // that is found in add_to_lid.c and requires add_to_lid.h
        lid_get_DrainInflow(sIdx, ff, &nIdx, &qq);

        if (nIdx >= 0)
        {
            *nodeIdx  = nIdx+1;
            *flowrate = CFTOCM(qq);
        }
        else
        {
            *nodeIdx  = 0;
            *flowrate = 0.0;
        }
    }
    return 0;
}
//===============================================================================
int DLLEXPORT api_call_climate_setState(double thisDate)
//===============================================================================
    // calls the climate_setState() procedure in SWMM-C
{
    //printf(" \n  in api_call_climate_setState \n");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    climate_setState(thisDate);
    
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_evaporation_rate(double *evapRate)
//===============================================================================
    // gets the Evap.rate (ft/s) from EPA SWMM
{
    //printf(" \n  in api_get_evaporation_rate \n");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    *evapRate = FTTOM(Evap.rate);

    //printf(" evap rate = %e \n ",Evap.rate);
    
    return 0;
}

// -------------------------------------------------------------------------
// |
// |  Private functionalities
// v
// -------------------------------------------------------------------------
//===============================================================================
int api_load_vars()
//===============================================================================
{
    char  line[MAXLINE+1];        // line from input data file
    char  wLine[MAXLINE+1];       // working copy of input line
    int sect, ii, jj, kk, error;
    int found = 0;
    double xx[4];

    error = check_api_is_initialized("api_load_vars");
    if (error) return error;

    for (ii = 0; ii < NUM_API_DOUBLE_VARS; ii++)
    {
        api->double_vars[ii] = (double*) calloc(Nobjects[LINK], sizeof(double));
    }

    rewind(Finp.file);
    while ( fgets(line, MAXLINE, Finp.file) != NULL )
    {
        // --- make copy of line and scan for tokens
        strcpy(wLine, line);
        Ntokens = getTokens(wLine);

        // --- skip blank lines and comments
        if ( Ntokens == 0 ) continue;
        if ( *Tok[0] == ';' ) continue;

        if (*Tok[0] == '[')
        {
            if (found) break;
            sect = findmatch(Tok[0], SectWords);
        }
        else
        {
            if (sect == s_XSECTION)
            {
                found = 1;
                jj = project_findObject(LINK, Tok[0]);
                kk = findmatch(Tok[1], XsectTypeWords);
                if ( kk == TRAPEZOIDAL )
                {
                    // --- parse and save geometric parameters
                    for (ii = 2; ii <= 5; ii++)
                        getDouble(Tok[ii], &xx[ii-2]);

                    // --- extract left and right slopes for trapezoidal channel
                    api->double_vars[api_left_slope] [jj] = xx[2];
                    api->double_vars[api_right_slope][jj] = xx[3];
                }
            }
        }
        continue;
    }
    return 0;
}
//===============================================================================
// int add_link_alt(
//     int li_idx,
//     int ni_idx,
//     int maxUp,
//     int maxDn,
//     int direction,
//     int* ni_N_link_u,
//     int* ni_N_link_d,
//     int* MlinkUp[3],
//     int*,MlinkDn[3])
// //===============================================================================
// {
//     if (direction == UPSTREAM) {
//         ni_N_link_u[ni_idx] ++;
//         if (ni_N_link_up[ni_idx] <= maxUp){
//             MlinkUp[ni_N_link_up[ni_idx]] = li_idx;
//         } else {
//             sprintf(errmsg, "incoming links for NODE %s > Max allowed [api.c -> add_link_alt]", Node[ni_idx].ID);
//             api_report_writeErrorMsg(api_err_model_junctions, errmsg);
//             return api_err_model_junctions;
//         }
//         return 0;
//     } else {
//         ni_N_link_d[ni_idx] ++;
//         if (ni_N_link_dn[ni_idx] <= maxDn){
//             MlinkDn[ni_N_link_dn[ni_idx]] = li_idx;
//         } else {
//             sprintf(errmsg, "outgoing links for NODE %s > 3 [api.c -> add_link_alt]", Node[ni_idx].ID);
//             api_report_writeErrorMsg(api_err_model_junctions, errmsg);
//             return api_err_model_junctions;
//         }
//         return 0;    
//     }
//     api_report_writeErrorMsg(api_err_internal, "[api.c -> add_link_alt]");
//     return api_err_internal;
// }

//===============================================================================
int add_link(
    int li_idx,
    int ni_idx,
    int direction,
    int* ni_N_link_u,
    int* ni_Mlink_u1,
    int* ni_Mlink_u2,
    int* ni_Mlink_u3,
    int* ni_N_link_d,
    int* ni_Mlink_d1,
    int* ni_Mlink_d2,
    int* ni_Mlink_d3)
//===============================================================================
{
    if (direction == UPSTREAM) {
        ni_N_link_u[ni_idx] ++;
        if (ni_N_link_u[ni_idx] <= 3) {
            if (ni_N_link_u[ni_idx] == 1) {
                ni_Mlink_u1[ni_idx] = li_idx;
            } else if (ni_N_link_u[ni_idx] == 2) {
                ni_Mlink_u2[ni_idx] = li_idx;
            } else if (ni_N_link_u[ni_idx] == 3) {
                ni_Mlink_u3[ni_idx] = li_idx;
            } else {
                sprintf(errmsg, "ni_Mlink_u3 == 0 at NODE %s [api.c -> add_link]", Node[ni_idx].ID);
                api_report_writeErrorMsg(api_err_internal, "errmsg");
                return api_err_internal;
            }
            return 0;
        } else {
            sprintf(errmsg, "incoming links for NODE %s > 3 [api.c -> add_link]", Node[ni_idx].ID);
            api_report_writeErrorMsg(api_err_model_junctions, errmsg);
            return api_err_model_junctions;
        }
    } else {
        ni_N_link_d[ni_idx] ++;
        if (ni_N_link_d[ni_idx] <= 3) {
            if (ni_N_link_d[ni_idx] == 1) {
                ni_Mlink_d1[ni_idx] = li_idx;
            } else if (ni_N_link_d[ni_idx] == 2) {
                ni_Mlink_d2[ni_idx] = li_idx;
            } else if (ni_N_link_d[ni_idx] == 3) {
                ni_Mlink_d3[ni_idx] = li_idx;
            } else {
                sprintf(errmsg, "ni_Mlink_d3 == 0 at NODE %s [api.c -> add_link]", Node[ni_idx].ID);
                api_report_writeErrorMsg(api_err_internal, "errmsg");
                return api_err_internal;
            }
            return 0;
        } else {
            sprintf(errmsg, "outgoing links for NODE %s > 3 [api.c -> add_link]", Node[ni_idx].ID);
            api_report_writeErrorMsg(api_err_model_junctions, errmsg);
            return api_err_model_junctions;
        }
    }
    api_report_writeErrorMsg(api_err_internal, "[api.c -> add_link]");
    return api_err_internal;
}

//===============================================================================
int check_api_is_initialized(
    char * function_name)
//===============================================================================
    //  provides error if api has not been initialized
{
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        sprintf(errmsg, "[api.c -> %s -> check_api_is_initialized]", function_name);
        api_report_writeErrorMsg(api_err_not_initialized, errmsg);
        return api_err_not_initialized;
    }
    return 0;
}

//===============================================================================
int getTokens(
    char *ss)
//===============================================================================
// Copy pasted getTokens from src/input.c to ensure independence
// from the original EPA-SWMM code. In the original code
// getTokens is not defined as an external API function
//
//  Input:   s = a character string
//  Output:  returns number of tokens found in s
//  Purpose: scans a string for tokens, saving pointers to them
//           in shared variable Tok[].
//
//  Notes:   Tokens can be separated by the characters listed in SEPSTR
//           (spaces, tabs, newline, carriage return) which is defined
//           in CONSTS.H. Text between quotes is treated as a single token.
//
{
    int  len, mm, nn;
    char *cc;

    // --- begin with no tokens
    for (nn = 0; nn < MAXTOKS; nn++) Tok[nn] = NULL;
    nn = 0;

    // --- truncate s at start of comment
    cc = strchr(ss,';');
    if (cc) *cc = '\0';
    len = strlen(ss);

    // --- scan s for tokens until nothing left
    while (len > 0 && nn < MAXTOKS)
    {
        mm = strcspn(ss,SEPSTR);              // find token length
        if (mm == 0) ss++;                    // no token found
        else
        {
            if (*ss == '"')                  // token begins with quote
            {
                ss++;                        // start token after quote
                len--;                      // reduce length of s
                mm = strcspn(ss,"\"\n");      // find end quote or new line
            }
            ss[mm] = '\0';                    // null-terminate the token
            Tok[nn] = ss;                     // save pointer to token
            nn++;                            // update token count
            ss += mm+1;                       // begin next token
        }
        len -= mm+1;                         // update length of s
    }
    return(nn);
}
//===============================================================================
// EOF
//===============================================================================