#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "swmm5.h"
#include "headers.h"
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

static char *Tok[MAXTOKS];             // String tokens from line of input
static int  Ntokens;                   // Number of tokens in line of input

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

void* DLLEXPORT api_initialize(char* f1, char* f2, char* f3)
{
    int error;
    Interface* api = (Interface*) malloc(sizeof(Interface));
    error = swmm_open(f1, f2, f3);
    if (error != 0) {
        api->IsInitialized = FALSE;
        return (void*) api;
    }
    error = swmm_start(0);
    if (error != 0) {
        api->IsInitialized = FALSE;
        return (void*) api;
    }
    api->IsInitialized = TRUE;
    error = api_load_vars((void*) api);
    return (void*) api;
}

void DLLEXPORT api_finalize(void* f_api)
//
//  Input: f_api is an Interface object passed as a void*
//  Output: None
//  Purpose: Closes the link with the EPA-SWMM library
//
{
    int i;
    Interface* api = (Interface*) f_api;

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

    free((Interface*) f_api);
}

// --- Property-extraction

// * After Initialization

double DLLEXPORT api_get_start_datetime()
{
    return StartDateTime;
}

double DLLEXPORT api_get_end_datetime()
{
    return EndDateTime;
}

int DLLEXPORT api_get_node_attribute(void* f_api, int k, int attr, double* value)
{
    int error;
    Interface * api = (Interface*) f_api;
    error = check_api_is_initialized(api);
    if (error != 0) return error;

    if (attr == node_type)
        *value = Node[k].type;
    else if (attr == node_invertElev)
        *value = FTTOM(Node[k].invertElev);
    else if (attr == node_initDepth)
        *value = FTTOM(Node[k].initDepth);
    else if (attr == node_extInflow_tSeries)
    {
        if (Node[k].extInflow)
            *value = Node[k].extInflow->tSeries;
        else
            *value = -1;
    }
    else if (attr == node_extInflow_basePat)
    {
        if (Node[k].extInflow)
            *value = Node[k].extInflow->basePat;
        else
            *value = -1;
    }
    else if (attr == node_extInflow_baseline)
    {
        if (Node[k].extInflow)
        {
            *value = CFTOCM(Node[k].extInflow->baseline);
        }
        else
            *value = 0;
    }
    else if (attr == node_extInflow_sFactor)
    {
        if (Node[k].extInflow)
            *value = Node[k].extInflow->sFactor;
        else
            *value = 1;
    }
    else if (attr == node_has_extInflow)
    {
        if (Node[k].extInflow)
        {
            *value = 1;
        }
        else
            *value = 0;
    }
    else if (attr == node_dwfInflow_monthly_pattern)
    {
        if (Node[k].dwfInflow)
            *value = Node[k].dwfInflow->patterns[0];
        else
            *value = -1;
    }
    else if (attr == node_dwfInflow_daily_pattern)
    {
        if (Node[k].dwfInflow)
        {
            *value = Node[k].dwfInflow->patterns[1];
        }
        else
            *value = -1;
    }
    else if (attr == node_dwfInflow_hourly_pattern)
    {
        if (Node[k].dwfInflow)
            *value = Node[k].dwfInflow->patterns[2];
        else
            *value = -1;
    }
    else if (attr == node_dwfInflow_weekend_pattern)
    {
        if (Node[k].dwfInflow)
            *value = Node[k].dwfInflow->patterns[3];
        else
            *value = -1;
    }
    else if (attr == node_dwfInflow_avgvalue)
    {
        if (Node[k].dwfInflow)
            *value = CFTOCM(Node[k].dwfInflow->avgValue);
        else
            *value = 0;
    }
    else if (attr == node_has_dwfInflow)
    {
        if (Node[k].dwfInflow)
            *value = 1;
        else
            *value = 0;
    }
    else if (attr == node_depth)
        *value = FTTOM(Node[k].newDepth);
    else if (attr == node_inflow)
        *value = CFTOCM(Node[k].inflow);
    else if (attr == node_volume)
        *value = CFTOCM(Node[k].newVolume);
    else if (attr == node_overflow)
        *value = CFTOCM(Node[k].overflow);
    else
        *value = nullvalue;
    return 0;
}

int DLLEXPORT api_get_link_attribute(void* f_api, int k, int attr, double* value)
{
    int error;
    Interface * api = (Interface*) f_api;

    error = check_api_is_initialized(api);
    if (error != 0) return error;

    if (attr == link_subIndex)
        *value = Link[k].subIndex;
    else if (attr == link_type)
        *value = Link[k].type;
    else if (attr == link_node1)
        *value = Link[k].node1;
    else if (attr == link_node2)
        *value = Link[k].node2;
    else if (attr == link_xsect_type)
        *value = Link[k].xsect.type;
    else if (attr == link_xsect_wMax)
        *value = FTTOM(Link[k].xsect.wMax);
    else if (attr == link_xsect_yBot)
        *value = FTTOM(Link[k].xsect.yBot);
    else if (attr == link_q0)
        *value = CFTOCM(Link[k].q0);
    else if (attr == conduit_roughness)
    {
        if (Link[k].type == CONDUIT)
            *value = Conduit[Link[k].subIndex].roughness;
        else
            *value = 0;
    }
    else if (attr == conduit_length)
    {
        if (Link[k].type == CONDUIT)
            *value = FTTOM(Conduit[Link[k].subIndex].length);
        else
            *value = 0;
    }
    else if (attr == link_flow)
        *value = CFTOCM(Link[k].newFlow);
    else if (attr == link_depth)
        *value = FTTOM(Link[k].newDepth);
    else if (attr == link_volume)
        *value = CFTOCM(Link[k].newVolume);
    else if (attr == link_froude)
        *value = Link[k].froude;
    else if (attr == link_setting)
        *value = Link[k].setting;
    else if (attr == link_left_slope)
    {
        *value = api->double_vars[api_left_slope][k];
    }
    else if (attr == link_right_slope)
    {
        *value = api->double_vars[api_right_slope][k];
    }
    else
        *value = nullvalue;
    return 0;
}

int DLLEXPORT api_get_num_objects(void* f_api, int object_type)
{
    int error;
    Interface * api = (Interface*) f_api;
    error = check_api_is_initialized(api);
    if (error != 0) return error;
    // if (object_type > API_START_INDEX) // Objects for API purposes
    //     return api->num_objects[object_type - API_START_INDEX];
    return Nobjects[object_type];
}

int DLLEXPORT api_get_object_name_len(void* f_api, int k, int object_type)
{
    int error;
    Interface * api = (Interface*) f_api;

    error = check_api_is_initialized(api);
    if (error != 0) return error;

    if (object_type == NODE)
    {
        return strlen(Node[k].ID);
    }
    else if (object_type == LINK)
    {
        return strlen(Link[k].ID);
    }
    else
        return -1;
}

int DLLEXPORT api_get_object_name(void* f_api, int k, char* object_name, int object_type)
{
    int error;
    Interface * api = (Interface*) f_api;

    error = check_api_is_initialized(api);
    if (error != 0) return error;
    if (object_type == NODE)
    {
        for (int i=0; i<api_get_object_name_len(f_api, k, object_type); i++)
        {
            object_name[i] = Node[k].ID[i];
        }
    }
    else if (object_type == LINK)
    {
        for (int i=0; i<api_get_object_name_len(f_api, k, object_type); i++)
        {
            object_name[i] = Link[k].ID[i];
        }
    }
    else
        return ERROR_FEATURE_NOT_COMPATIBLE;
    return 0;
}

int DLLEXPORT api_get_first_table_entry(int k, int table_type, double* x, double* y)
{
    int success;
    if (table_type == TSERIES)
    {
        success = table_getFirstEntry(&Tseries[k], x, y);
        if (Tseries[k].curveType == -1 && Tseries[k].refersTo == EXTERNAL_INFLOW)
            *y = CFTOCM(*y);
        return success;
    }
    else
    {
        return ERR_API_WRONG_TYPE;
    }
}

int DLLEXPORT api_get_next_table_entry(int k, int table_type, double* x, double* y)
{
    int success;
    if (table_type == TSERIES)
    {
        success = table_getNextEntry(&Tseries[k], x, y);
        if (Tseries[k].curveType == -1 && Tseries[k].refersTo == EXTERNAL_INFLOW)
            *y = CFTOCM(*y);
        return success;
    }
    else
    {
        return ERR_API_WRONG_TYPE;
    }
}

int DLLEXPORT api_get_pattern_count(int k)
{
    return Pattern[k].count;
}

double DLLEXPORT api_get_pattern_factor(int k, int j)
{
    return Pattern[k].factor[j];
}

int DLLEXPORT api_get_pattern_type(int k)
{
    return Pattern[k].type;
}

double DLLEXPORT api_get_next_inflow_bc(void* f_api, int node_idx, double current_datetime) {

    int ptype, pcount, i, j, p;
    int yy, mm, dd;
    int h, m, s, dow;
    double val;
    double x, y, next_datetime;
    double bline, sfactor;
    double total_factor = 1;
    double total_extinflow = 0;
    double total_inflow = 0;

    datetime_decodeDate(current_datetime, &yy, &mm, &dd);
    datetime_decodeTime(current_datetime, &h, &m, &s);
    dow = datetime_dayOfWeek(current_datetime);

    api_get_node_attribute(f_api, node_idx, node_has_dwfInflow, &val);
    if (val == 1) { // node_has_dwfInflow
        for (int i=0; i<4; i++)
        {
            j = Node[node_idx].dwfInflow->patterns[i];
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
        total_inflow += total_factor * Node[node_idx].dwfInflow->avgValue;
    }

    api_get_node_attribute(f_api, node_idx, node_has_extInflow, &val);
    if (val == 1) { // node_has_dwfInflow
        p = Node[node_idx].extInflow->basePat; // pattern
        if (p > 0)
        {
            ptype = Pattern[p].type;
            bline = Node[node_idx].extInflow->baseline; // baseline value
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
        j = Node[node_idx].extInflow->tSeries; // tseries
        sfactor = Node[node_idx].extInflow->sFactor; // scale factor
        if (j > 0)
        {
            total_extinflow += table_lookup(&Tseries[j], current_datetime) * sfactor;
        }
        total_inflow += total_extinflow;
    }

    return total_inflow;
}

double api_get_runoff(int node_idx, double current_datetime)
{
    // --- compute runoff until next routing time reached or exceeded
    while ( NewRunoffTime < current_datetime )
    {
        runoff_execute();
        if ( ErrorCode ) return;
    }
    return CFTOCM(Node[node_idx].newLatFlow);
}

// --- Print-out

void DLLEXPORT api_print_object_name(int k, int object_type)
{
    if (object_type == NODE) printf("%s\n", Node[k].ID);
    if (object_type == LINK) printf("%s\n", Link[k].ID);
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
                return -500;
            }
            return 0;
        } else {
            return -400;
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
                return -500;
            }
            return 0;
        } else {
            return -400;
        }
    }
    return -296;
}


int DLLEXPORT api_export_linknode_properties(void* f_api, int units)
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

    Interface * api = (Interface*) f_api;

    error = check_api_is_initialized(api);
    if (error != 0) return error;

    // Initialization
    for (i = 0; i<Nobjects[NODE]; i++) {
        ni_N_link_u[i] = 0;
        ni_N_link_d[i] = 0;
        ni_Mlink_u1[i] = nullvalueI;
        ni_Mlink_u2[i] = nullvalueI;
        ni_Mlink_u3[i] = nullvalueI;
        ni_Mlink_d1[i] = nullvalueI;
        ni_Mlink_d2[i] = nullvalueI;
        ni_Mlink_d3[i] = nullvalueI;
    }

    // Choosing unit system
    if (units == US) {
        flow_units = 1;
        manning_units = 1;
        length_units = 1;
    } else {
        flow_units = M3perFT3;
        manning_units = pow(1/MperFT, 1/3);
        length_units = MperFT;
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
        if (error != 0) return error;

        li_Mnode_d[i] = Link[i].node2;
        error = add_link(i, li_Mnode_d[i], UPSTREAM, ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
        if (error != 0) return error;

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

    f_nodes = fopen("debug/nodes_info.csv", "w");
    f_links = fopen("debug/links_info.csv", "w");

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

int DLLEXPORT api_export_link_results(void* f_api, char* link_name)
{
	FILE* tmp;
    DateTime days;
    int period, j;
    char theTime[20];
    char theDate[20];
	char path[50];

    int error;
    Interface * api = (Interface*) f_api;
    error = check_api_is_initialized(api);
    if (error != 0) return error;

    j = project_findObject(LINK, link_name);

    if (stat("LinkResults", &st) == -1) {
        mkdir("LinkResults", 0700);
    }

    /* File path writing */
    strcpy(path, "LinkResults/");
    strcat(path, Link[j].ID); strcat(path, ".csv");
    tmp = fopen(path, "w");
    fprintf(tmp, "date,time,flow,velocity,depth,volume,capacity\n");

    for ( period = 1; period <= Nperiods; period++ )
    {
        output_readDateTime(period, &days);
        datetime_dateToStr(days, theDate);
        datetime_timeToStr(days, theTime);
        output_readLinkResults(period, j);
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

int DLLEXPORT api_export_node_results(void* f_api, char* node_name)
{
	FILE* tmp;
    DateTime days;
    int period, j;
    char theTime[20];
    char theDate[20];
	char path[50];
    int error;

    Interface * api = (Interface*) f_api;
    error = check_api_is_initialized(api);
    if (error != 0) return error;

    j = project_findObject(NODE, node_name);

    if (stat("NodeResults", &st) == -1) {
        mkdir("NodeResults", 0700);
    }

    /* File path writing */
    strcpy(path, "NodeResults/");
    strcat(path, Node[j].ID);
    strcat(path, ".csv");
    tmp = fopen(path, "w");
    fprintf(tmp, "date,time,inflow,overflow,depth,volume\n");

    for ( period = 1; period <= Nperiods; period++ ) {
        output_readDateTime(period, &days);
        datetime_dateToStr(days, theDate);
        datetime_timeToStr(days, theTime);
        output_readNodeResults(period, j);
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
int check_api_is_initialized(Interface* api)
{
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    return 0;
}

// ---
int api_load_vars(void * f_api)
{
    Interface * api = (Interface*) f_api;
    char  line[MAXLINE+1];        // line from input data file
    char  wLine[MAXLINE+1];       // working copy of input line
    int sect, i, j, k, error;
    int found = 0;
    double x[4];

    error = check_api_is_initialized(api);
    if (error != 0) return error;

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


int api_findObject(int type, char *id)
{
    return project_findObject(type, id);
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