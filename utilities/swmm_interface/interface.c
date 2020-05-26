#include <stdlib.h>

#include "swmm5.h"
#include "headers.h"
#include "interface.h"


void* DLLEXPORT api_initialize(char* f1, char* f2, char* f3)
{
    Interface* api = (Interface*) malloc(sizeof(Interface));
    swmm_open(f1, f2, f3);
    swmm_start(0);
    api->IsInitialized = TRUE;
    return (void*) api;
}

void DLLEXPORT api_finalize(void* f_api)
{
    swmm_end();
    swmm_close();
    free((Interface*) f_api);
}

double DLLEXPORT api_get_node_attribute (void* f_api, int k, int attr)
{
    Interface * api = (void*) f_api;
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
    Interface * api = (void*) f_api;
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