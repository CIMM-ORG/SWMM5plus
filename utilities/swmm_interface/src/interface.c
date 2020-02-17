#include "headers.h"
#include "interface.h"
#include <math.h>
#include <string.h>

int i_add_link(
	int li_idx,
	int ni_idx,
	int direction,
	float * ni_N_link_u,
	float * ni_Mlink_u1,
	float * ni_Mlink_u2,
	float * ni_Mlink_u3,
	float * ni_N_link_d,
	float * ni_Mlink_d1,
	float * ni_Mlink_d2,
	float * ni_Mlink_d3)
{

	if (direction == UPSTREAM) {
		ni_N_link_u[ni_idx] ++;
		if (ni_N_link_u[ni_idx] <= 3) {
			if (ni_N_link_u[ni_idx] == 1) {
				ni_Mlink_u1[ni_idx] = li_idx;
			}
			else if (ni_N_link_u[ni_idx] == 2) {
				ni_Mlink_u2[ni_idx] = li_idx;
			}
			else if (ni_N_link_u[ni_idx] == 3) {
				ni_Mlink_u3[ni_idx] = li_idx;
			}
			else {
				return -500;
			}
			return 0;
		}
		else {
			return -400;
		}
	}
	else {
		ni_N_link_d[ni_idx] ++;
		if (ni_N_link_d[ni_idx] <= 3) {
			if (ni_N_link_d[ni_idx] == 1) {
				ni_Mlink_d1[ni_idx] = li_idx;
			}
			else if (ni_N_link_d[ni_idx] == 2) {
				ni_Mlink_d2[ni_idx] = li_idx;
			}
			else if (ni_N_link_d[ni_idx] == 3) {
				ni_Mlink_d3[ni_idx] = li_idx;
			}
			else {
				return -500;
			}
			return 0;
		}
		else {
			return -400;
		}
	}
	return -296;
}

float i_get_node_type (int k, float total_n_links)
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

float i_get_link_xsect_attrs (int k, int attr, float length_units)
{
    int i = k - FIDX;

	if (Link[i].xsect.type == RECT_CLOSED)
	{
        if (attr == e_li_geometry)
        {
		    return L_RECTANGULAR;
        }
        else if (attr == e_li_link_type)
        {
		    return L_PIPE;
        }
        else if (attr == e_lr_BreadthScale)
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
        if (attr == e_li_geometry)
        {
		    return L_RECTANGULAR;
        }
        else if (attr == e_li_link_type)
        {
            return L_CHANNEL;
        }
        else if (attr == e_lr_BreadthScale)
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
		if (attr == e_li_geometry)
        {
		    return L_TRAPEZOIDAL;
        }
		else if (attr == e_li_link_type)
        {
            return L_CHANNEL;
        }
		else if (attr == e_lr_BreadthScale)
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
		if (attr == e_li_geometry)
        {
		    return L_PARABOLIC;
        }
		else if (attr == e_li_link_type)
        {
            return L_CHANNEL;
        }
		else if (attr == e_lr_BreadthScale)
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

int i_init_node_tmp_table(
    float * ni_N_link_u,
    float * ni_Mlink_u1,
    float * ni_Mlink_u2,
    float * ni_Mlink_u3,
    float * ni_N_link_d,
    float * ni_Mlink_d1,
    float * ni_Mlink_d2,
    float * ni_Mlink_d3)
{
    int error;
    // Initialization
    for (int i = 0; i < Nobjects[NODE]; i++) {
        ni_N_link_u[i] = 0;
        ni_N_link_d[i] = 0;
        ni_Mlink_u1[i] = nullvalueI;
        ni_Mlink_u2[i] = nullvalueI;
        ni_Mlink_u3[i] = nullvalueI;
        ni_Mlink_d1[i] = nullvalueI;
        ni_Mlink_d2[i] = nullvalueI;
        ni_Mlink_d3[i] = nullvalueI;
    }

    // Count Links
    for (int i = 0; i < Nobjects[LINK]; i++)
    {
        error = i_add_link(i, Link[i].node1, DOWNSTREAM, ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
        if (error != 0) return error;

        error = i_add_link(i, Link[i].node2, UPSTREAM, ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
        if (error != 0) return error;
    }
    return 0;
}

// Fortran indexes start in FIDX, therefore
// k \in [FIDX, Nobjects[.]]
int i_get_node_attrs (
    int k,
    int length_units,
    float * attrs,
    float * ni_N_link_u,
    float * ni_Mlink_u1,
    float * ni_Mlink_u2,
    float * ni_Mlink_u3,
    float * ni_N_link_d,
    float * ni_Mlink_d1,
    float * ni_Mlink_d2,
    float * ni_Mlink_d3)
{
    int i = k - FIDX;
    float total_n_links = ni_N_link_u[i] + ni_N_link_d[i];
    attrs[e_ni_node_type] = i_get_node_type(k, total_n_links);
    attrs[e_ni_N_link_u] = ni_N_link_u[i];
    attrs[e_ni_N_link_d] = ni_N_link_d[i];
    attrs[e_ni_Mlink_u1] = ni_Mlink_u1[i] + FIDX;
    attrs[e_ni_Mlink_u2] = ni_Mlink_u2[i] + FIDX;
    attrs[e_ni_Mlink_u3] = ni_Mlink_u3[i] + FIDX;
    attrs[e_ni_Mlink_d1] = ni_Mlink_d1[i] + FIDX;
    attrs[e_ni_Mlink_d2] = ni_Mlink_d2[i] + FIDX;
    attrs[e_ni_Mlink_d3] = ni_Mlink_d3[i] + FIDX;
    attrs[e_nr_Zbottom] = Node[i].invertElev * length_units;
    return 0;
}

float i_get_link_attribute (int k, int attr, int units)
{
    float flow_units;
    float manning_units;
    float length_units;
    int i = k - FIDX;

    // Choosing unit system
	if (units == US) {
		flow_units = 1;
		manning_units = 1;
		length_units = 1;
	}
	else {
		flow_units = M3perFT3;
		manning_units = pow(1 / MperFT, 1 / 3);
		length_units = MperFT;
	}

    if (attr == e_li_link_type || attr == e_li_geometry || attr == e_lr_BreadthScale)
    {
        return i_get_link_xsect_attrs(k, attr, length_units);
    }
    else if (attr == e_li_roughness_type)
    {
        // Manning coefficient by default at the moment
        // We have to decide how to incorporate other types of coeffs
        return E_MANNINGS_N;
    }
    else if (attr == e_li_Mnode_u)
    {
		return Link[i].node1;
    }
    else if (attr == e_li_Mnode_d)
    {
		return Link[i].node2;
    }
    else if (attr == e_li_InitialDepthType)
    {
        // DEFAULT
        // Same issue as rou li_roughness_type
        return 1;
    }
    else if (attr == e_lr_Length)
    {
        if (Link[i].type == CONDUIT)
        {
            return Conduit[Link[i].subIndex].length * length_units;
        }
        return 0;
    }
    else if (attr == e_lr_Slope)
    {
        float h, sign_h, X;

        if (Link[i].type == CONDUIT)
        {
            h = (Node[Link[i].node1].invertElev - Node[Link[i].node2].invertElev) * length_units;
            sign_h = SSIGN(h);
            h = fabs(h);
            X = sqrt(pow(Conduit[Link[i].subIndex].length * length_units, 2) - pow(h, 2));
            return h / X;
        }

        return 0;
    }
    else if (attr == e_lr_Roughness)
    {
        if (Link[i].type == CONDUIT)
        {
            return Conduit[Link[i].subIndex].roughness * manning_units;
        }

        return 0;
    }
    else if (attr == e_lr_InitialFlowrate)
    {
        return Link[i].q0 * flow_units;
    }
    else if (attr == e_lr_InitialUpstreamDepth)
    {
        return Node[Link[i].node1].initDepth * length_units;
    }
    else if (attr == e_lr_InitialDnstreamDepth)
    {
        return Node[Link[i].node2].initDepth * length_units;
    }
    else if (attr == e_lr_InitialDepth)
    {
		return fabs(i_get_link_attribute(k, e_lr_InitialDnstreamDepth, units) - i_get_link_attribute(k, e_lr_InitialUpstreamDepth, units)) / 2;
    }
    else
    {
        return nullvalueI;
    }
}

int i_num_links ()
{
    return Nobjects[LINK];
}

int i_num_nodes ()
{
    return Nobjects[NODE];
}

void i_record_nodes_data(double * nodeMatrix, int units)
{
    float * ni_N_link_u = (float *) calloc(Nobjects[NODE], sizeof(float));
    float * ni_N_link_d = (float *) calloc(Nobjects[NODE], sizeof(float));
    float * ni_Mlink_u1 = (float *) calloc(Nobjects[NODE], sizeof(float));
    float * ni_Mlink_u2 = (float *) calloc(Nobjects[NODE], sizeof(float));
    float * ni_Mlink_u3 = (float *) calloc(Nobjects[NODE], sizeof(float));
    float * ni_Mlink_d1 = (float *) calloc(Nobjects[NODE], sizeof(float));
    float * ni_Mlink_d2 = (float *) calloc(Nobjects[NODE], sizeof(float));
    float * ni_Mlink_d3 = (float *) calloc(Nobjects[NODE], sizeof(float));
    float attrs[num_node_attributes];
    float flow_units, length_units, manning_units;
    int len = num_node_attributes * Nobjects[NODE];
    int x, y;

    // Choosing unit system
	if (units == US) {
		flow_units = 1;
		manning_units = 1;
		length_units = 1;
	}
	else {
		flow_units = M3perFT3;
		manning_units = pow(1 / MperFT, 1 / 3);
		length_units = MperFT;
	}

    i_init_node_tmp_table( ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
	for (int i = 0; i < len; i+=num_node_attributes)
    {
        x = i % num_node_attributes;
        y = floor(i / num_node_attributes);
        i_get_node_attrs(x+FIDX, length_units, attrs, ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
        for (int j = 0; j < num_node_attributes; j++)
        {
            nodeMatrix[i+j] = attrs[j];
        }
    }
    free(ni_N_link_u);
    free(ni_Mlink_u1);
    free(ni_Mlink_u2);
    free(ni_Mlink_u3);
    free(ni_N_link_d);
    free(ni_Mlink_d1);
    free(ni_Mlink_d2);
    free(ni_Mlink_d3);
}

void i_record_links_data(double * linkMatrix, int units)
{
    int len = num_link_attributes * Nobjects[LINK];
    int x, y;
	for (int i = 0; i < len; i++)
    {
        y = i % num_link_attributes;
        x = floor(i/num_link_attributes);
        linkMatrix[i] = i_get_link_attribute(x+FIDX, y, units);
    }
}

void i_print_info(int units)
{
	FILE *f_nodes;
	FILE *f_links;
    int NNodes = Nobjects[NODE];
    int NLinks = Nobjects[LINK];
    int error;
    float * ni_N_link_u = (float *) malloc(Nobjects[NODE]*sizeof(float));
    float * ni_N_link_d = (float *) malloc(Nobjects[NODE]*sizeof(float));
    float * ni_Mlink_u1 = (float *) malloc(Nobjects[NODE]*sizeof(float));
    float * ni_Mlink_u2 = (float *) malloc(Nobjects[NODE]*sizeof(float));
    float * ni_Mlink_u3 = (float *) malloc(Nobjects[NODE]*sizeof(float));
    float * ni_Mlink_d1 = (float *) malloc(Nobjects[NODE]*sizeof(float));
    float * ni_Mlink_d2 = (float *) malloc(Nobjects[NODE]*sizeof(float));
    float * ni_Mlink_d3 = (float *) malloc(Nobjects[NODE]*sizeof(float));
    float attrs[num_node_attributes];
    float flow_units, length_units, manning_units;

    // Choosing unit system
	if (units == US) {
		flow_units = 1;
		manning_units = 1;
		length_units = 1;
	}
	else {
		flow_units = M3perFT3;
		manning_units = pow(1 / MperFT, 1 / 3);
		length_units = MperFT;
	}

    #ifdef WINDOWS
	    fopen_s(&f_nodes, "nodes_info.csv", "w");
	    fopen_s(&f_links, "links_info.csv", "w");
    #else
	    f_nodes = fopen("nodes_info.csv", "w");
	    f_links = fopen("links_info.csv", "w");
    #endif

    error = i_init_node_tmp_table( ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
	fprintf(f_nodes,
		"n_left,node_id,ni_idx,ni_node_type,ni_N_link_u,ni_N_link_d,ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,ni_Mlink_d2,ni_Mlink_d3,nr_Zbottom\n");
	for (int i = 0; i < NNodes; i++)
    {
        i_get_node_attrs(i+FIDX, length_units, attrs, ni_N_link_u, ni_Mlink_u1, ni_Mlink_u2, ni_Mlink_u3, ni_N_link_d, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3);
		fprintf(f_nodes, "%d,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
			NNodes - i,
			Node[i].ID,
			i + 1,
			attrs[e_ni_node_type],
			attrs[e_ni_N_link_u],
			attrs[e_ni_N_link_d],
			attrs[e_ni_Mlink_u1],
			attrs[e_ni_Mlink_u2],
			attrs[e_ni_Mlink_u3],
			attrs[e_ni_Mlink_d1],
			attrs[e_ni_Mlink_d2],
			attrs[e_ni_Mlink_d3],
			attrs[e_nr_Zbottom]);
	}
    fclose(f_nodes);

    // free tmp arrays
    free(ni_N_link_u);
    free(ni_Mlink_u1);
    free(ni_Mlink_u2);
    free(ni_Mlink_u3);
    free(ni_N_link_d);
    free(ni_Mlink_d1);
    free(ni_Mlink_d2);
    free(ni_Mlink_d3);

	fprintf(f_links,
		"l_left,link_id,li_idx,li_link_type,li_roughness_type,li_geometry,li_Mnode_u,li_Mnode_d,li_InitialDepthType,lr_Length,lr_BreadthScale,lr_Slope,lr_Roughness,lr_InitialFlowrate,lr_InitialDepth,lr_InitialUpstreamDepth,lr_InitialDnstreamDepth\n");
	for (int i = 0; i < NLinks; i++) {
		fprintf(f_links, "%d,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.8f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
			NLinks - i,
			Link[i].ID,
			i + 1,
			i_get_link_attribute(i+FIDX, e_li_link_type, units),
			i_get_link_attribute(i+FIDX, e_li_roughness_type, units),
			i_get_link_attribute(i+FIDX, e_li_geometry, units),
			i_get_link_attribute(i+FIDX, e_li_Mnode_u, units),
			i_get_link_attribute(i+FIDX, e_li_Mnode_d, units),
			i_get_link_attribute(i+FIDX, e_li_InitialDepthType, units),
			i_get_link_attribute(i+FIDX, e_lr_Length, units),
			i_get_link_attribute(i+FIDX, e_lr_BreadthScale, units),
			i_get_link_attribute(i+FIDX, e_lr_Slope, units),
			i_get_link_attribute(i+FIDX, e_lr_Roughness, units),
			i_get_link_attribute(i+FIDX, e_lr_InitialFlowrate, units),
			i_get_link_attribute(i+FIDX, e_lr_InitialDepth, units),
			i_get_link_attribute(i+FIDX, e_lr_InitialUpstreamDepth, units),
			i_get_link_attribute(i+FIDX, e_lr_InitialDnstreamDepth, units));
	} fclose(f_links);
}