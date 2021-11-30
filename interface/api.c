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

// --- Simulation

// Initializes the EPA-SWMM simulation. It creates an Interface
// variable (details about Interface in interface.h). The
// function opens de SWMM input file and creates report and
// output files. Although the .inp is parsed within swmm_start,
// not every property is extracted, e.g., slopes of trapezoidal channels.
// In swmm.c
// to parse the .inp again using EPA-SWMM functionalities
// within api_load_vars
// The interface object is passed to many functions in interface.c
// but it is passed as a void pointer. This is because the object
// is also used by the Fortran engine. Treating Interface as void*
// facilitates interoperability

int DLLEXPORT api_initialize(char* f1, char* f2, char* f3, int run_routing)
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

void DLLEXPORT api_finalize()
//
//  Input: f_api is an Interface object passed as a void*
//  Output: None
//  Purpose: Closes the link with the EPA-SWMM library
//3
{
    int i;

    swmm_end();
    swmm_close();

    // frees double variables in API
    for (i = 0; i < NUM_API_DOUBLE_VARS; i++)
    {
        if (api->double_vars[i] != NULL)
            free(api->double_vars[i]);
    }

    // // frees integer variables in API
    // for (i = 0; i < NUM_API_INT_VARS; i++)
    // {
    //     if (api->int_vars[i] != NULL)
    //         free(api->int_vars[i]);
    // }

    free(api);
}

double DLLEXPORT api_run_step()
{
    swmm_step(&(api->elapsedTime));
    return api->elapsedTime;
}

// --- Property-extraction

// * During Simulation


int DLLEXPORT api_get_node_results(char* node_name, float* inflow, float* overflow, float* depth, float* volume)
//
//  Input:    f_api = Interface object passed as a void*
//            node_name = string identifier of node
//            inflow, overflow, depth, volume =
//  Output: None
//  Purpose: Closes the link with the SWMM C library
//
{
    int j, error;

    error = check_api_is_initialized("api_get_node_results");
    if (error != 0) return error;

    j = project_findObject(NODE, node_name);
    *inflow = Node[j].inflow;
    *overflow = Node[j].overflow;
    *depth = Node[j].newDepth;
    *volume = Node[j].newVolume;
    return 0;
}

int DLLEXPORT api_get_link_results(char* link_name, float* flow, float* depth, float* volume)
{
    int j, error;

    error = check_api_is_initialized("api_get_link_results");
    if (error != 0) return error;
    j = project_findObject(LINK, link_name);
    *flow = Link[j].newFlow;
    *depth = Link[j].newDepth;
    *volume = Link[j].newVolume;
    return 0;
}

// * After Initialization

double DLLEXPORT api_get_start_datetime()
{
    return StartDateTime;
}

double DLLEXPORT api_get_end_datetime()
{
    return EndDateTime;
}

int DLLEXPORT api_get_flowBC(int node_idx, double current_datetime, double* flowBC)
{

    int ptype, pcount, i, j, p;
    int yy, mm, dd;
    int h, m, s, dow, attr;
    double val;
    double x, y, next_datetime;
    double bline, sfactor;
    double total_factor = 1;
    double total_extinflow = 0;

    *flowBC = 0;
    datetime_decodeDate(current_datetime, &yy, &mm, &dd);
    datetime_decodeTime(current_datetime, &h, &m, &s);
    dow = datetime_dayOfWeek(current_datetime);
    attr = node_has_dwfInflow;
    api_get_node_attribute(node_idx, attr, &val);
    if (val == 1) { // node_has_dwfInflow
        for(i=0; i<4; i++)
        {
            j = Node[node_idx].dwfInflow->patterns[i];
            if (j > 0)
            {
                ptype = Pattern[j].type;
                if (ptype == MONTHLY_PATTERN)
                    total_factor *= Pattern[j].factor[mm-1];
                else if (ptype == DAILY_PATTERN)
                    total_factor *= Pattern[j].factor[dow-1];
                else if (ptype == HOURLY_PATTERN)
                    total_factor *= Pattern[j].factor[h];
                else if (ptype == WEEKEND_PATTERN)
                {
                    if ((dow == 1) || (dow == 7))
                        total_factor *= Pattern[j].factor[h];
                }
            }
        }
        *flowBC += total_factor * CFTOCM(Node[node_idx].dwfInflow->avgValue);
    }

    attr = node_has_extInflow;
    api_get_node_attribute(node_idx, attr, &val);
    if (val == 1) { // node_has_extInflow
        p = Node[node_idx].extInflow->basePat; // pattern
        bline = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->baseline); // baseline value
        if (p >= 0)
        {
            ptype = Pattern[p].type;
            if (ptype == MONTHLY_PATTERN)
                total_extinflow += Pattern[j].factor[mm-1] * bline;
            else if (ptype == DAILY_PATTERN)
                total_extinflow += Pattern[j].factor[dow-1] * bline;
            else if (ptype == HOURLY_PATTERN)
                total_extinflow += Pattern[j].factor[h] * bline;
            else if (ptype == WEEKEND_PATTERN)
            {
                if ((dow == 1) || (dow == 7))
                    total_extinflow += Pattern[j].factor[h] * bline;
            }
        }
        else if (bline > 0)
        {
            total_extinflow += bline;
        }

        j = Node[node_idx].extInflow->tSeries; // tseries
        sfactor = Node[node_idx].extInflow->sFactor; // scale factor
        if (j >= 0)
        {
            total_extinflow += table_tseriesLookup(&Tseries[j], current_datetime, FALSE) * sfactor;
        }
        *flowBC += total_extinflow;
    }
    return 0;
}

int DLLEXPORT api_get_headBC(int node_idx, double current_datetime, double* headBC)
{
    int i = Node[node_idx].subIndex;
    
    switch (Outfall[i].type)
    {
        case FIXED_OUTFALL:
            *headBC = FTTOM(Outfall[i].fixedStage);
            return 0;

        default:
            *headBC = API_NULL_VALUE_I;
            sprintf(errmsg, "OUTFALL TYPE %d at NODE %s", Outfall[i].type, Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_not_developed, errmsg);
            return api_err_not_developed;
    }
}

int DLLEXPORT api_get_report_times(double* report_start_datetime, int* report_step, int* hydrology_step) // WET_STEP in SWMM 5.13
{
    int error;

    error = check_api_is_initialized("api_get_report_times");
    if (error) return error;
    *report_start_datetime = ReportStart;
    *report_step = ReportStep;
    *hydrology_step = WetStep;
    return 0;
}

int DLLEXPORT api_get_node_attribute(int node_idx, int attr, double* value)
{
    int error, bpat, tseries_idx;

    error = check_api_is_initialized("api_get_node_attribute");
    if (error) return error;

    if (attr == node_type)
        *value = Node[node_idx].type;
    else if (attr == node_outfall_type)
    {
        if (Node[node_idx].type == OUTFALL)
        {
            *value = Outfall[Node[node_idx].subIndex].type;
        }
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_outfall_type for NODE %s, which is not an outfall [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_invertElev)
        *value = FTTOM(Node[node_idx].invertElev);

    else if (attr == node_fullDepth)
        *value = FTTOM(Node[node_idx].fullDepth);

    else if (attr == node_initDepth)
    {
        if (Node[node_idx].type == OUTFALL)
        {
            error = api_get_headBC(node_idx, StartDateTime, value);
            if (error) return error;
            *value -= FTTOM(Node[node_idx].invertElev);
        }
        else
            *value = FTTOM(Node[node_idx].initDepth);
    }
    else if (attr == node_StorageConstant)
    {
        if (Node[node_idx].type == STORAGE)
            *value = Storage[Node[node_idx].subIndex].aConst;
        else
            *value = -1;
    }
    else if (attr == node_StorageCoeff)
    {
        if (Node[node_idx].type == STORAGE)
            *value = Storage[Node[node_idx].subIndex].aCoeff;
        else
            *value = -1;
    }
    else if (attr == node_StorageExponent)
    {
        if (Node[node_idx].type == STORAGE)
            *value = Storage[Node[node_idx].subIndex].aExpon;
        else
            *value = -1;
    }
    else if (attr == node_StorageCurveID)
    {
        if (Node[node_idx].type == STORAGE)
            *value = Storage[Node[node_idx].subIndex].aCurve;
        else
            *value = -1;
    }
    else if (attr == node_extInflow_tSeries)
    {
        if (Node[node_idx].extInflow)
            *value = Node[node_idx].extInflow->tSeries;
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_extInflow_tSeries for NODE %s, which doesn't have an extInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_extInflow_tSeries_x1)
    {
        tseries_idx = Node[node_idx].extInflow->tSeries;
        if (tseries_idx >= 0)
            *value = Tseries[tseries_idx].x1;
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_extInflow_tSeries_x1 for NODE %s, which doesn't have an extInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_extInflow_tSeries_x2)
    {
        tseries_idx = Node[node_idx].extInflow->tSeries;
        if (tseries_idx >= 0)
            *value = Tseries[tseries_idx].x2;
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting tseries_idx for NODE %s, which doesn't have an extInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_extInflow_basePat)
    {
        if (Node[node_idx].extInflow)
            *value = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->basePat);
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_extInflow_basePat for NODE %s, which doesn't have an extInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_extInflow_basePat_type)
    {
        bpat = Node[node_idx].extInflow->basePat;
        if (bpat >= 0) // baseline pattern exists
            *value = Pattern[bpat].type;
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_extInflow_basePat_type for NODE %s, which doesn't have an extInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_extInflow_baseline)
    {
        if (Node[node_idx].extInflow)
        {
            *value = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->baseline);
        }
        else
            *value = 0;
    }
    else if (attr == node_extInflow_sFactor)
    {
        if (Node[node_idx].extInflow)
            *value = Node[node_idx].extInflow->sFactor;
        else
            *value = 1;
    }
    else if (attr == node_has_extInflow)
    {
        if (Node[node_idx].extInflow)
        {
            *value = 1;
        }
        else
            *value = 0;
    }
    else if (attr == node_dwfInflow_monthly_pattern)
    {
        if (Node[node_idx].dwfInflow)
            *value = Node[node_idx].dwfInflow->patterns[0];
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_dwfInflow_monthly_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_dwfInflow_daily_pattern)
    {
        if (Node[node_idx].dwfInflow)
        {
            *value = Node[node_idx].dwfInflow->patterns[1];
        }
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_dwfInflow_daily_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_dwfInflow_hourly_pattern)
    {
        if (Node[node_idx].dwfInflow)
            *value = Node[node_idx].dwfInflow->patterns[2];
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_dwfInflow_hourly_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_dwfInflow_weekend_pattern)
    {
        if (Node[node_idx].dwfInflow)
            *value = Node[node_idx].dwfInflow->patterns[3];
        else
        {
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "Extracting node_dwfInflow_weekend_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_node_attribute]", Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_wrong_type, errmsg);
            return api_err_wrong_type;
        }
    }
    else if (attr == node_dwfInflow_avgvalue)
    {
        if (Node[node_idx].dwfInflow)
            *value = CFTOCM(Node[node_idx].dwfInflow->avgValue);
        else
            *value = 0;
    }
    else if (attr == node_has_dwfInflow)
    {
        if (Node[node_idx].dwfInflow)
            *value = 1;
        else
            *value = 0;
    }
    else if (attr == node_depth)
        *value = FTTOM(Node[node_idx].newDepth);
    else if (attr == node_inflow)
        *value = CFTOCM(Node[node_idx].inflow);
    else if (attr == node_volume)
        *value = CFTOCM(Node[node_idx].newVolume);
    else if (attr == node_overflow)
        *value = CFTOCM(Node[node_idx].overflow);
    else
        *value = API_NULL_VALUE_I;
    return 0;
}

int DLLEXPORT api_get_link_attribute(int link_idx, int attr, double* value)
{
    int error;

    error = check_api_is_initialized("api_get_link_attribute");
    if (error) return error;

    if (attr == link_subIndex)
        *value = Link[link_idx].subIndex;
    else if (attr == link_type)
        *value = Link[link_idx].type;
    else if (attr == link_node1)
        *value = Link[link_idx].node1;
    else if (attr == link_node2)
        *value = Link[link_idx].node2;
    else if (attr == link_offset1)
        *value = FTTOM(Link[link_idx].offset1);
    else if (attr == link_offset2)
        *value = FTTOM(Link[link_idx].offset2);
    else if (attr == link_xsect_type)
        *value = Link[link_idx].xsect.type;
    else if (attr == link_xsect_wMax)
        *value = FTTOM(Link[link_idx].xsect.wMax);
    else if (attr == link_xsect_yBot)
        *value = FTTOM(Link[link_idx].xsect.yBot);
    else if (attr == link_xsect_yFull)
        *value = FTTOM(Link[link_idx].xsect.yFull);
    else if (attr == link_q0)
        *value = CFTOCM(Link[link_idx].q0);
    else if (attr == link_type)
        *value =  Link[link_idx].type;
    else if (attr == pump_type)
        *value =  Pump[Link[link_idx].subIndex].pumpCurve;
    else if (attr == orifice_type)
        *value = Orifice[Link[link_idx].subIndex].type;
    else if (attr == weir_type)
        *value = Weir[Link[link_idx].subIndex].type;
    else if (attr == conduit_roughness)
    {
        if (Link[link_idx].type == CONDUIT)
            *value = Conduit[Link[link_idx].subIndex].roughness;
        else
            *value = 0;
    }
    else if (attr == conduit_length)
    {
        if (Link[link_idx].type == CONDUIT)
            *value = FTTOM(Conduit[Link[link_idx].subIndex].length);

        else if (Link[link_idx].type == ORIFICE)
            *value = 0.01;
        else if (Link[link_idx].type == WEIR)
            *value = 0.01;
        else
            *value = 0;
    }
    else if (attr == weir_end_contractions)
    {
        if (Link[link_idx].type == WEIR)
            *value = Weir[Link[link_idx].subIndex].endCon;
        else
            *value = 0;
    }
    else if (attr == discharge_coeff1)
    {
        if (Link[link_idx].type == WEIR)
            *value = Weir[Link[link_idx].subIndex].cDisch1;
        else if (Link[link_idx].type == ORIFICE)
            *value = Orifice[Link[link_idx].subIndex].cDisch;
        else
            *value = 0;
    }
    else if (attr == discharge_coeff2)
    {
        if (Link[link_idx].type == WEIR)
            *value = Weir[Link[link_idx].subIndex].cDisch2;
        else
            *value = 0;
    }
    else if (attr == weir_side_slope)
    {
        if (Link[link_idx].type == WEIR)
            *value = Weir[Link[link_idx].subIndex].slope;
        else
            *value = 0;
    }
    else if (attr == link_flow)
        *value = CFTOCM(Link[link_idx].newFlow);
    else if (attr == link_depth)
        *value = FTTOM(Link[link_idx].newDepth);
    else if (attr == link_volume)
        *value = CFTOCM(Link[link_idx].newVolume);
    else if (attr == link_froude)
        *value = Link[link_idx].froude;
    else if (attr == link_setting)
        *value = Link[link_idx].setting;
    else if (attr == link_left_slope)
        *value = api->double_vars[api_left_slope][link_idx];
    else if (attr == link_right_slope)
        *value = api->double_vars[api_right_slope][link_idx];
    else
        *value = API_NULL_VALUE_I;
    return 0;
}

int DLLEXPORT api_get_num_objects(int object_type)
{
    int error;
    error = check_api_is_initialized("api_get_num_objects");
    if (error) return error;
    // if (object_type > API_START_INDEX) // Objects for API purposes
    //     return api->num_objects[object_type - API_START_INDEX];
    return Nobjects[object_type];
}

int DLLEXPORT api_get_object_name(int object_idx, char* object_name, int object_type)
{
    int error, i;
    int obj_len = -1;

    error = check_api_is_initialized("api_get_object_name");
    if (error) return error;
    error = api_get_object_name_len(object_idx, object_type, &obj_len);
    if (error) return error;

    if (object_type == NODE)
    {
        for(i=0; i<obj_len; i++)
        {
            object_name[i] = Node[object_idx].ID[i];
        }
    }
    else if (object_type == LINK)
    {
        for(i=0; i<obj_len; i++)
        {
            object_name[i] = Link[object_idx].ID[i];
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

int DLLEXPORT api_get_object_name_len(int object_idx, int object_type, int* len)
{
    int error;

    error = check_api_is_initialized("api_get_object_name_len");
    if (error) {
        *len = API_NULL_VALUE_I;
        return error;
    }

    if (object_type == NODE)
    {
        *len = strlen(Node[object_idx].ID);
        return 0;
    }
    else if (object_type == LINK)
    {
        *len = strlen(Link[object_idx].ID);
        return 0;
    }
    else
    {
        *len = API_NULL_VALUE_I;
        sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name_len]", object_type);
        api_report_writeErrorMsg(api_err_not_developed, errmsg);
        return api_err_not_developed;
    }
}

int DLLEXPORT api_get_num_table_entries(int table_idx, int table_type, int* num_entries)
{
    double x, y;
    int success;

    *num_entries = 0;
    // printf("1 Number of entries in Curve %d: %d\n", k, *num_entries);
    if (table_type == CURVE)
    {
        // ERROR handling
        if (table_idx >= Nobjects[CURVE] || Nobjects[CURVE] == 0) return -1;
        success = table_getFirstEntry(&Curve[table_idx], &x, &y); // first values in the table
        (*num_entries)++;
        while (success)
        {
            success = table_getNextEntry(&(Curve[table_idx]), &x, &y);
            if (success) (*num_entries)++;
            // printf("0 Number of entries in Curve %d: %d\n", k, *num_entries);
        }
    }
    else
    {
        return -1;
    }
    // printf("Number of entries in Curve %d: %d\n", k, *num_entries);
    // printf("SUCCESS %d\n", success);

    return 0;
}

int DLLEXPORT api_get_table_attribute(int table_idx, int attr, double* value)
{
    int error;

    error = check_api_is_initialized("api_get_table_attribute");
    if (error) return error;

    if (attr == table_ID)
        *value = table_idx;
    else if (attr == table_type)
        *value = Curve[table_idx].curveType;
    else if (attr == table_refers_to)
        *value = Curve[table_idx].refersTo;
    else
    {
        *value = API_NULL_VALUE_I;
        sprintf(errmsg, "attr %d [api.c -> api_get_table_attribute]", attr);
        api_report_writeErrorMsg(api_err_not_developed, errmsg);
        return api_err_not_developed;
    }
    return 0;
}

int DLLEXPORT api_get_first_entry_table(int table_idx, int table_type, double* x, double* y)
{
    int success;

    if (table_type == CURVE)
        success = table_getFirstEntry(&(Curve[table_idx]), x, y);
    else if (table_type == TSERIES)
        success = table_getFirstEntry(&(Tseries[table_idx]), x, y);
    else
        return 0;
    return success;
}

int DLLEXPORT api_get_next_entry_table(int table_idx, int table_type, double* x, double* y)
{
    int success;

    if (table_type == TSERIES)
    {
        success = table_getNextEntry(&(Tseries[table_idx]), &(Tseries[table_idx].x2), &(Tseries[table_idx].y2));
        if (success)
        {
            *x = Tseries[table_idx].x2;
            *y = Tseries[table_idx].y2;
        }
    }
    else if (table_type == CURVE)
    {
        success = table_getNextEntry(&(Curve[table_idx]), &(Curve[table_idx].x2), &(Curve[table_idx].y2));
        if (success)
        {
            *x = Curve[table_idx].x2;
            *y = Curve[table_idx].y2;
        }
    }

    return success;
}

int DLLEXPORT api_get_next_entry_tseries(int tseries_idx)
{
    int success;
    double x2, y2;

    x2 = Tseries[tseries_idx].x2;
    y2 = Tseries[tseries_idx].y2;
    success = table_getNextEntry(&(Tseries[tseries_idx]), &(Tseries[tseries_idx].x2), &(Tseries[tseries_idx].y2));
    if (success == TRUE)
    {
        Tseries[tseries_idx].x1 = x2;
        Tseries[tseries_idx].y1 = y2;
    }
    return success;
}

// --- Output Writing (Post Processing)
// * The follwing functions should only be executed after finishing
//   and writing SWMM5+ report files. The following functions are
//   meant to be called from Fortran in order to export .rpt and
//   .out files according to the SWMM 5.13 standard. Fortran-generated
//   report files are not manipulated here, the manipulation of
//   SWMM5+ report files is kept within the Fortran code to ensure
//   compatibility with future updates of the SWMM5+ standard

int DLLEXPORT api_write_output_line(double t)
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

int DLLEXPORT api_update_nodeResult(int node_idx, int resultType, double newNodeResult)
{
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

int DLLEXPORT api_update_linkResult(int link_idx, int resultType, double newLinkResult)
{
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

// --- Print-out

int DLLEXPORT api_export_linknode_properties(int units)
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
    float lr_InitialFlowrate[Nobjects[LINK]];
    float lr_InitialUpstreamDepth[Nobjects[LINK]];
    float lr_InitialDnstreamDepth[Nobjects[LINK]];
    int li_InitialDepthType[Nobjects[LINK]]; //
    float lr_BreadthScale[Nobjects[LINK]]; //
    float lr_InitialDepth[Nobjects[LINK]]; //

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
        error = add_link(i, li_Mnode_u[i], DOWNSTREAM, ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
        if (error) return error;

        li_Mnode_d[i] = Link[i].node2;
        error = add_link(i, li_Mnode_d[i], UPSTREAM, ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
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

        lr_InitialFlowrate[i] = Link[i].q0 * flow_units;
        lr_InitialUpstreamDepth[i] = Node[li_Mnode_u[i]].initDepth * length_units;
        lr_InitialDnstreamDepth[i] = Node[li_Mnode_d[i]].initDepth * length_units;
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
        "l_left,link_id,li_idx,li_link_type,li_geometry,li_Mnode_u,li_Mnode_d,lr_Length,lr_Slope,lr_Roughness,lr_InitialFlowrate,lr_InitialUpstreamDepth,lr_InitialDnstreamDepth\n");
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
            lr_InitialFlowrate[i],
            lr_InitialUpstreamDepth[i],
            lr_InitialDnstreamDepth[i]);
    }
    fclose(f_links);

    return 0;
}

int DLLEXPORT api_export_link_results(int link_idx)
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

int DLLEXPORT api_export_node_results(int node_idx)
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

// --- Utils

int DLLEXPORT api_find_object(int object_type, char *id)
{
    return project_findObject(object_type, id);
}

// -------------------------------------------------------------------------
// |
// |  Private functionalities
// v
// -------------------------------------------------------------------------

int api_load_vars()
{
    char  line[MAXLINE+1];        // line from input data file
    char  wLine[MAXLINE+1];       // working copy of input line
    int sect, i, j, k, error;
    int found = 0;
    double x[4];

    error = check_api_is_initialized("api_load_vars");
    if (error) return error;

    for (i = 0; i < NUM_API_DOUBLE_VARS; i++)
    {
        api->double_vars[i] = (double*) calloc(Nobjects[LINK], sizeof(double));
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
                j = project_findObject(LINK, Tok[0]);
                k = findmatch(Tok[1], XsectTypeWords);
                if ( k == TRAPEZOIDAL )
                {
                    // --- parse and save geometric parameters
                    for (i = 2; i <= 5; i++)
                        getDouble(Tok[i], &x[i-2]);

                    // --- extract left and right slopes for trapezoidal channel
                    api->double_vars[api_left_slope][j] = x[2];
                    api->double_vars[api_right_slope][j] = x[3];
                }
            }
        }
        continue;
    }
    return 0;
}

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

int check_api_is_initialized(char * function_name)
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

// Copy pasted getTokens from src/input.c to ensure independence
// from the original EPA-SWMM code. In the original code
// getTokens is not defined as an external API function

int getTokens(char *s)
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
    int  len, m, n;
    char *c;

    // --- begin with no tokens
    for (n = 0; n < MAXTOKS; n++) Tok[n] = NULL;
    n = 0;

    // --- truncate s at start of comment
    c = strchr(s,';');
    if (c) *c = '\0';
    len = strlen(s);

    // --- scan s for tokens until nothing left
    while (len > 0 && n < MAXTOKS)
    {
        m = strcspn(s,SEPSTR);              // find token length
        if (m == 0) s++;                    // no token found
        else
        {
            if (*s == '"')                  // token begins with quote
            {
                s++;                        // start token after quote
                len--;                      // reduce length of s
                m = strcspn(s,"\"\n");      // find end quote or new line
            }
            s[m] = '\0';                    // null-terminate the token
            Tok[n] = s;                     // save pointer to token
            n++;                            // update token count
            s += m+1;                       // begin next token
        }
        len -= m+1;                         // update length of s
    }
    return(n);
}