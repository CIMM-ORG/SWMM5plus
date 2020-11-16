#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "swmm5.h"
#include "headers.h"
#include "interface.h"
#include "arrays.h"

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
static char *Tok[MAXTOKS];             // String tokens from line of input
static int  Ntokens;                   // Number of tokens in line of input

void* DLLEXPORT api_initialize(char* f1, char* f2, char* f3)
{
    Interface* api = (Interface*) malloc(sizeof(Interface));
    swmm_open(f1, f2, f3);
    swmm_start(0);
    api_load_vars((void*) api);
    api->IsInitialized = TRUE;
    if ( Nobjects[SUBCATCH] > 0 ) api->DoRunoff = TRUE;
    else api->DoRunoff = FALSE;
    return (void*) api;
}

void DLLEXPORT api_finalize(void* f_api)
{
    int i;
    Interface* api = (Interface*) f_api;

    swmm_end();
    swmm_close();

    for (i = 0; i < NUM_API_VARS; i++)
        free(api->vars[i]);

    free((Interface*) f_api);
}

double DLLEXPORT api_get_node_attribute (void* f_api, int k, int attr)
{
    Interface * api = (Interface*) f_api;
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    switch (attr)
    {
        case node_type:
            return Node[k].type;
        case node_invertElev:
            return Node[k].invertElev;
        case node_initDepth:
            return Node[k].initDepth;
        case node_extInflow_tSeries:
            if (Node[k].extInflow)
                return Node[k].extInflow->tSeries;
            return -1;
        case node_extInflow_basePat:
            if (Node[k].extInflow)
                return Node[k].extInflow->basePat;
            return -1;
        case node_extInflow_baseline:
            if (Node[k].extInflow)
                return Node[k].extInflow->baseline;
            return -1;
        case node_depth:
            return Node[k].newDepth;
        case node_inflow:
            return Node[k].inflow;
        case node_volume:
            return Node[k].newVolume;
        case node_overflow:
            return Node[k].overflow;
        default:
            return nullvalue;
    }
}

double DLLEXPORT api_get_link_attribute (void* f_api, int k, int attr)
{
    Interface * api = (Interface*) f_api;
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    switch (attr)
    {
        case link_subIndex:
            return Link[k].subIndex;
        case link_type:
            return Link[k].type;
        case link_node1:
            return Link[k].node1;
        case link_node2:
            return Link[k].node2;
        case link_xsect_type:
            return Link[k].xsect.type;
        case link_xsect_wMax:
            return Link[k].xsect.wMax;
        case link_xsect_yBot:
            return Link[k].xsect.yBot;
        case link_q0:
            return Link[k].q0;
        case conduit_roughness:
            if (Link[k].type == CONDUIT)
                return Conduit[Link[k].subIndex].roughness;
            else
                return 0;
        case conduit_length:
            if (Link[k].type == CONDUIT)
                return Conduit[Link[k].subIndex].length;
            else
                return 0;
        case link_flow:
            return Link[k].newFlow;
        case link_depth:
            return Link[k].newDepth;
        case link_volume:
            return Link[k].newVolume;
        case link_froude:
            return Link[k].froude;
        case link_setting:
            return Link[k].setting;
        case link_left_slope:
            return api->vars[api_left_slope][k];
        case link_right_slope:
            return api->vars[api_right_slope][k];
        default:
            return nullvalue;
    }
}

int DLLEXPORT api_num_links ()
{
    return Nobjects[LINK];
}

int DLLEXPORT api_num_nodes ()
{
    return Nobjects[NODE];
}

int DLLEXPORT api_num_time_series ()
{
    return Nobjects[TSERIES];
}

int DLLEXPORT api_num_curves ()
{
    return Nobjects[CURVE];
}

int api_load_vars (void * f_api)
{
    Interface * api = (Interface*) f_api;
    char  line[MAXLINE+1];        // line from input data file
    char  wLine[MAXLINE+1];       // working copy of input line
    int sect, i, j, k;
    int found = 0;
    double x[4];

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    for (i = 0; i < NUM_API_VARS; i++)
    {
        api->vars[i] = (double*) calloc(Nobjects[LINK], sizeof(double));
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
                    api->vars[api_left_slope][j] = x[2];
                    api->vars[api_right_slope][j] = x[3];
                }
            }
        }
        continue;
    }
    return 0;
}

int api_runoff_step(void * f_api)
{
    Interface * api = (Interface*) f_api;

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    if (NewRunoffTime == TotalDuration) return -1;
    if ( api->DoRunoff && NewRunoffTime < TotalDuration )
    {
        runoff_execute();
        if ( ErrorCode ) return error_getCode(ErrorCode);
    }

    // --- update elapsed time (days)
    if ( NewRoutingTime < TotalDuration )
    {
        ElapsedTime = NewRoutingTime / MSECperDAY;
    }
    // --- otherwise end the simulation
    else ElapsedTime = 0.0;

    api->elapsedTime = ElapsedTime;
    return 0;
}

int api_get_node_inflow(void * f_api, int j, DateTime * current_time, double * inflow)
{
    // j : node id
    int i, k, success;
    int p_ext;
    int p_dry[4];
    int yy, mm, dd, h;
    double x, y;
    DateTime next_time = 1e20;
    Interface * api = (Interface*) f_api;

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    *current_time *= 86400;
    *inflow = 0;

    k = Node[j].extInflow->basePat;
    p_ext = -1;
    if (k >= 0) p_ext = Pattern[k].type;
    for (i = 0; i < 4; i++)
        p_dry[i] = Pattern[Node[j].dwfInflow->patterns[i]].type;

    datetime_decodeDate(*current_time, &yy, &mm, &dd);

    if (p_ext == HOURLY_PATTERN || p_dry[0] == HOURLY_PATTERN || p_dry[1] == HOURLY_PATTERN || p_dry[2] == HOURLY_PATTERN || p_dry[3] == HOURLY_PATTERN)
        next_time = fmin(next_time, (*current_time + 3600 - fmod(*current_time, 3600)) / 86400);
    if (p_ext == DAILY_PATTERN || p_dry[0] == DAILY_PATTERN || p_dry[1] == DAILY_PATTERN || p_dry[2] == DAILY_PATTERN || p_dry[3] == DAILY_PATTERN)
        next_time = fmin(next_time, (*current_time + 86400 - fmod(*current_time, 86400)) / 86400);
    if (p_ext == MONTHLY_PATTERN || p_dry[0] == MONTHLY_PATTERN || p_dry[1] == MONTHLY_PATTERN || p_dry[2] == MONTHLY_PATTERN || p_dry[3] == MONTHLY_PATTERN)
    {
        if (++mm > 12) { mm = 1; yy++; }
        next_time = fmin(next_time, datetime_encodeDate(yy, mm, dd));
    }
    if (p_ext == WEEKEND_PATTERN || p_dry[0] == WEEKEND_PATTERN || p_dry[1] == WEEKEND_PATTERN || p_dry[2] == WEEKEND_PATTERN || p_dry[3] == WEEKEND_PATTERN)
    {
        dd = datetime_dayOfYear(*current_time);
        if (dd == 1 || dd == 7) // Sunday or Saturday
        {
            if (dd == 1)
                next_time = fmin(next_time, floor(datetime_addDays(*current_time, 6)));
            else if (dd == 7)
                next_time = fmin(next_time, floor(datetime_addDays(*current_time, 7)));
        }
    }
    // i = Node[j].extInflow->tSeries;
    // if (i >= 0)
    // {
    //     if (*current_time/86400 == StartDateTime)
    //         success = table_getFirstEntry(&Tseries[i], &x, &y);
    //     else
    //         success = table_getNextEntry(&Tseries[i], &x, &y);
    //     next_time = fmin(next_time, x);
    // }

    *current_time = next_time; // Update current time
    *inflow += inflow_getExtInflow(Node[j].extInflow, next_time);
    datetime_decodeDate(next_time, &yy, &mm, &dd);
    h = datetime_hourOfDay(next_time);
    *inflow += inflow_getDwfInflow(Node[j].dwfInflow, mm, dd, h);
    return 0;
}

void DLLEXPORT api_print_pattern(void * f_api, int j)
{
    DateTime current_time = StartDateTime;
    double inflow;
    int yy, mm, dd, h, m, s;

    while (current_time <= EndDateTime)
    {
        api_get_node_inflow(f_api, j, &current_time, &inflow);
        datetime_decodeDate(current_time, &yy, &mm, &dd);
        datetime_decodeTime(current_time, &h, &m, &s);
        printf("%d/%d/%d %d:%d:%d %.2f\n", inflow, yy, mm, dd, h, m, s);
    }
}
//-----------------------------------------------------------------------------
//  External Local Functions
//-----------------------------------------------------------------------------

// --- Functions retrieved from input.c

int  getTokens(char *s)
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

int DLLEXPORT api_get_next_tseries_entry(void* f_api, int k, double* entries)
{
    // k : index of Tseries
    Interface * api = (Interface*) f_api;

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    entries[0] = Tseries[k].thisEntry->x;
    entries[1] = Tseries[k].thisEntry->y;

    return table_getNextEntry(&Tseries[k], &entries[2], &entries[3]);
}

int DLLEXPORT api_get_next_curve_entry(void* f_api, int k, double* entries)
{
    // k : index of Curves
    Interface * api = (Interface*) f_api;

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    entries[0] = Curve[k].thisEntry->x;
    entries[1] = Curve[k].thisEntry->y;

    return table_getNextEntry(&Curve[k], &entries[2], &entries[4]);
}