#include "swmm5.h"
#include "headers.h"
#include "interface.h"
#include <math.h>
#include <string.h>

void* DLLEXPORT api_initialize (char* f1, char* f2, char* f3, int unit_system)
{
    int error;
    Interface* api = (Interface*) malloc(sizeof(Interface));

    // --- open the files & read input data
    swmm_open(f1, f2, f3);
    swmm_start(0);
    //  // Choosing unit system
    api->units.unit_system = unit_system;

	if (unit_system == US)
    {
        api->units.flow_units = 1;
		api->units.manning_units = 1;
		api->units.length_units = 1;
	}
	else
    {
		api->units.flow_units = M3perFT3;
		api->units.manning_units = pow(1 / MperFT, 1 / 3);
		api->units.length_units = MperFT;
	}

    api->num_links = Nobjects[LINK];
    api->num_nodes = Nobjects[NODE];
    api->node_attributes = (float **) malloc(sizeof(float*)*num_node_attributes);
    for (int i = 0; i < num_node_attributes; i++)
        api->node_attributes[i] = (float *) calloc(api->num_nodes, sizeof(float));

    // Initialization
    for (int i = 0; i < api->num_nodes; i++)
    {
        api->node_attributes[ni_N_link_u][i] = 0;
        api->node_attributes[ni_N_link_d][i] = 0;
        api->node_attributes[ni_Mlink_u1][i] = nullvalueI;
        api->node_attributes[ni_Mlink_u2][i] = nullvalueI;
        api->node_attributes[ni_Mlink_u3][i] = nullvalueI;
        api->node_attributes[ni_Mlink_d1][i] = nullvalueI;
        api->node_attributes[ni_Mlink_d2][i] = nullvalueI;
        api->node_attributes[ni_Mlink_d3][i] = nullvalueI;
    }

    // Count Links
    for (int i = 0; i < api->num_links; i++)
    {
        add_link(api, i, Link[i].node1, DOWNSTREAM);
        add_link(api, i, Link[i].node2, UPSTREAM);
    }

    return (void *) api;
}

int DLLEXPORT api_finalize (void* fapi)
{
    Interface* api = (Interface*) fapi;

    for (int i = 0; i < num_node_attributes; i++)
        free(api->node_attributes[i]);
    free(api->node_attributes);

    free((Interface*) fapi);

    swmm_end();
    swmm_close();

    return 0;
}

float DLLEXPORT api_get_node_attribute (void* fapi, int k, int attr)
{
    Interface* api = (Interface*) fapi;
    int i = k - FIDX;
    if (attr == ni_node_type)
    {
        float total_n_links = api->node_attributes[ni_N_link_u][i] + api->node_attributes[ni_N_link_d][i];
        return get_node_type(k, total_n_links);
    }
    else if (attr == ni_N_link_u)
        return api->node_attributes[ni_N_link_u][i];
    else if (attr == ni_N_link_d)
        return api->node_attributes[ni_N_link_d][i];
    else if (attr == ni_Mlink_u1)
        return api->node_attributes[ni_Mlink_u1][i] + FIDX;
    else if (attr == ni_Mlink_u2)
        return api->node_attributes[ni_Mlink_u2][i] + FIDX;
    else if (attr == ni_Mlink_u3)
        return api->node_attributes[ni_Mlink_u3][i] + FIDX;
    else if (attr == ni_Mlink_d1)
        return api->node_attributes[ni_Mlink_d1][i] + FIDX;
    else if (attr == ni_Mlink_d2)
        return api->node_attributes[ni_Mlink_d2][i] + FIDX;
    else if (attr == ni_Mlink_d3)
        return api->node_attributes[ni_Mlink_d3][i] + FIDX;
    else if (attr == nr_Zbottom)
        return Node[i].invertElev * api->units.length_units;
    return 0;
}

float DLLEXPORT api_get_link_attribute (void* fapi, int k, int attr)
{
    Interface* api = (Interface*) fapi;
    int i = k - FIDX;

    if (attr == li_link_type || attr == li_geometry || attr == lr_BreadthScale)
    {
        return get_link_xsect_attrs(k, attr, api->units.length_units);
    }
    else if (attr == li_roughness_type)
    {
        // Manning coefficient by default at the moment
        // We have to decide how to incorporate other types of coeffs
        return E_MANNINGS_N;
    }
    else if (attr == li_Mnode_u)
    {
		return Link[i].node1;
    }
    else if (attr == li_Mnode_d)
    {
		return Link[i].node2;
    }
    else if (attr == li_InitialDepthType)
    {
        // DEFAULT
        // Same issue as rou li_roughness_type
        return 1;
    }
    else if (attr == lr_Length)
    {
        if (Link[i].type == CONDUIT)
        {
            return Conduit[Link[i].subIndex].length * api->units.length_units;
        }
        return 0;
    }
    else if (attr == lr_Slope)
    {
        float h, sign_h, X;

        if (Link[i].type == CONDUIT)
        {
            h = (Node[Link[i].node1].invertElev - Node[Link[i].node2].invertElev) * api->units.length_units;
            sign_h = SSIGN(h);
            h = fabs(h);
            X = sqrt(pow(Conduit[Link[i].subIndex].length * api->units.length_units, 2) - pow(h, 2));
            return h / X;
        }

        return 0;
    }
    else if (attr == lr_Roughness)
    {
        if (Link[i].type == CONDUIT)
        {
            return Conduit[Link[i].subIndex].roughness * api->units.manning_units;
        }

        return 0;
    }
    else if (attr == lr_InitialFlowrate)
    {
        return Link[i].q0 * api->units.flow_units;
    }
    else if (attr == lr_InitialUpstreamDepth)
    {
        return Node[Link[i].node1].initDepth * api->units.length_units;
    }
    else if (attr == lr_InitialDnstreamDepth)
    {
        return Node[Link[i].node2].initDepth * api->units.length_units;
    }
    else if (attr == lr_InitialDepth)
    {
		return fabs(api_get_link_attribute(api, k, lr_InitialDnstreamDepth) - api_get_link_attribute(api, k, lr_InitialUpstreamDepth)) / 2;
    }
    else
    {
        return nullvalueI;
    }
}

void DLLEXPORT api_print_info (void* fapi)
{
    Interface* api = (Interface*) fapi;
	FILE *f_nodes;
	FILE *f_links;
    int error;

    #ifdef WINDOWS
	    fopen_s(&f_nodes, "nodes_info.csv", "w");
	    fopen_s(&f_links, "links_info.csv", "w");
    #else
	    f_nodes = fopen("nodes_info.csv", "w");
	    f_links = fopen("links_info.csv", "w");
    #endif

	fprintf(f_nodes,
		"n_left,node_id,ni_idx,ni_node_type,ni_N_link_u,ni_N_link_d,ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,ni_Mlink_d2,ni_Mlink_d3,nr_Zbottom\n");
	for (int i = 0; i < api->num_nodes; i++)
    {
		fprintf(f_nodes, "%d,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
			api->num_nodes - i,
			Node[i].ID,
			i + 1,
			api->node_attributes[ni_node_type][i],
			api->node_attributes[ni_N_link_u][i],
			api->node_attributes[ni_N_link_d][i],
			api->node_attributes[ni_Mlink_u1][i],
			api->node_attributes[ni_Mlink_u2][i],
			api->node_attributes[ni_Mlink_u3][i],
			api->node_attributes[ni_Mlink_d1][i],
			api->node_attributes[ni_Mlink_d2][i],
			api->node_attributes[ni_Mlink_d3][i],
			api->node_attributes[nr_Zbottom][i]);
	}
    fclose(f_nodes);

	fprintf(f_links,
		"l_left,link_id,li_idx,li_link_type,li_roughness_type,li_geometry,li_Mnode_u,li_Mnode_d,li_InitialDepthType,lr_Length,lr_BreadthScale,lr_Slope,lr_Roughness,lr_InitialFlowrate,lr_InitialDepth,lr_InitialUpstreamDepth,lr_InitialDnstreamDepth\n");
	for (int i = 0; i < api->num_links; i++) {
		fprintf(f_links, "%d,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.8f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
			api->num_links - i,
			Link[i].ID,
			i + 1,
			api_get_link_attribute(api, i+FIDX, li_link_type),
			api_get_link_attribute(api, i+FIDX, li_roughness_type),
			api_get_link_attribute(api, i+FIDX, li_geometry),
			api_get_link_attribute(api, i+FIDX, li_Mnode_u),
			api_get_link_attribute(api, i+FIDX, li_Mnode_d),
			api_get_link_attribute(api, i+FIDX, li_InitialDepthType),
			api_get_link_attribute(api, i+FIDX, lr_Length),
			api_get_link_attribute(api, i+FIDX, lr_BreadthScale),
			api_get_link_attribute(api, i+FIDX, lr_Slope),
			api_get_link_attribute(api, i+FIDX, lr_Roughness),
			api_get_link_attribute(api, i+FIDX, lr_InitialFlowrate),
			api_get_link_attribute(api, i+FIDX, lr_InitialDepth),
			api_get_link_attribute(api, i+FIDX, lr_InitialUpstreamDepth),
			api_get_link_attribute(api, i+FIDX, lr_InitialDnstreamDepth));
	} fclose(f_links);
}

int DLLEXPORT api_num_links (void* fapi)
{
    Interface* api = (Interface*) fapi;
    return api->num_links;
}

int DLLEXPORT api_num_nodes (void* fapi)
{
    Interface* api = (Interface*) fapi;
    return api->num_nodes;
}

int add_link (Interface* api, int li_idx, int ni_idx, int direction)
{
	if (direction == UPSTREAM)
    {
		api->node_attributes[ni_N_link_u][ni_idx] ++;
		if (api->node_attributes[ni_N_link_u][ni_idx] <= 3)
        {
			if (api->node_attributes[ni_N_link_u][ni_idx] == 1)
				api->node_attributes[ni_Mlink_u1][ni_idx] = li_idx;
			else if (api->node_attributes[ni_N_link_u][ni_idx] == 2)
				api->node_attributes[ni_Mlink_u2][ni_idx] = li_idx;
			else if (api->node_attributes[ni_N_link_u][ni_idx] == 3)
				api->node_attributes[ni_Mlink_u3][ni_idx] = li_idx;
			else
				return -500;
			return 0;
		}
		else
			return -400;
	}
	else
    {
		api->node_attributes[ni_N_link_d][ni_idx] ++;
		if (api->node_attributes[ni_N_link_d][ni_idx] <= 3)
        {
			if (api->node_attributes[ni_N_link_d][ni_idx] == 1)
				api->node_attributes[ni_Mlink_d1][ni_idx] = li_idx;
			else if (api->node_attributes[ni_N_link_d][ni_idx] == 2)
				api->node_attributes[ni_Mlink_d2][ni_idx] = li_idx;
			else if (api->node_attributes[ni_N_link_d][ni_idx] == 3)
				api->node_attributes[ni_Mlink_d3][ni_idx] = li_idx;
			else
				return -500;
			return 0;
		}
		else
			return -400;
	}
	return -296;
}

float get_node_type (int k, float total_n_links)
{
    int i = k - FIDX;

    if (Node[i].extInflow)
    {
        if (Node[i].extInflow ->tSeries != -1 || Node[i].extInflow->basePat != 1 || Node[i].extInflow->baseline > 0 )
        {
            return N_BC_UP;
        }
    }

    if (Node[i].type == OUTFALL)
    {
        return N_BC_DN;
    }
    else if (total_n_links == 2)
    {
        return NJ2;
    }
    else if (total_n_links > 2)
    {
        return NJM;
    }
    else
    {
        return nullvalueI;
    }
}

float get_link_xsect_attrs (int k, int attr, float length_units)
{
    int i = k - FIDX;

	if (Link[i].xsect.type == RECT_CLOSED)
	{
        if (attr == li_geometry)
        {
		    return L_RECTANGULAR;
        }
        else if (attr == li_link_type)
        {
		    return L_PIPE;
        }
        else if (attr == lr_BreadthScale)
        {
            if (Link[i].type == CONDUIT)
            {
    		    return Link[i].xsect.wMax * length_units;
            }
            return 0;
        }
        else
        {
            return nullvalueI;
        }
	}
	else if (Link[i].xsect.type == RECT_OPEN)
	{
        if (attr == li_geometry)
        {
		    return L_RECTANGULAR;
        }
        else if (attr == li_link_type)
        {
            return L_CHANNEL;
        }
        else if (attr == lr_BreadthScale)
        {
            return Link[i].xsect.wMax * length_units;
        }
        else
        {
            return nullvalueI;
        }
	}
	else if (Link[i].xsect.type == TRAPEZOIDAL)
	{
		if (attr == li_geometry)
        {
		    return L_TRAPEZOIDAL;
        }
		else if (attr == li_link_type)
        {
            return L_CHANNEL;
        }
		else if (attr == lr_BreadthScale)
        {
            return Link[i].xsect.yBot * length_units;
        }
        else
        {
            return nullvalueI;
        }
	}
	else if (Link[i].xsect.type == PARABOLIC)
	{
		if (attr == li_geometry)
        {
		    return L_PARABOLIC;
        }
		else if (attr == li_link_type)
        {
            return L_CHANNEL;
        }
		else if (attr == lr_BreadthScale)
        {
            return Link[i].xsect.wMax * length_units;
        }
        else
        {
            return nullvalueI;
        }
	}
    return nullvalueI;
}