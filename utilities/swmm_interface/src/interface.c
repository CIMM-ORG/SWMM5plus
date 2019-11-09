#include "headers.h"
#include "interface.h"
#include <math.h>
#include <string.h>

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
	int* ni_Mlink_d3) {

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

int interface_print_info(int units) {
	int error;

	//  link
	int * li_idx = (int *) malloc(sizeof(int)*Nobjects[LINK]);
	int * li_link_type = (int *) malloc(sizeof(int)*Nobjects[LINK]);
	int * li_geometry = (int *) malloc(sizeof(int)*Nobjects[LINK]);
	int  * li_Mnode_u = (int *) malloc(sizeof(int)*Nobjects[LINK]);
	int * li_Mnode_d = (int *) malloc(sizeof(int)*Nobjects[LINK]);
	float * lr_Length = (float *) malloc(sizeof(float)*Nobjects[LINK]);
	float * lr_Slope = (float *) malloc(sizeof(float)*Nobjects[LINK]);
	float * lr_Roughness = (float *) malloc(sizeof(float)*Nobjects[LINK]);
	float * lr_InitialFlowrate = (float *) malloc(sizeof(float)*Nobjects[LINK]);
	float * lr_InitialUpstreamDepth = (float *) malloc(sizeof(float)*Nobjects[LINK]);
	float * lr_InitialDnstreamDepth = (float *) malloc(sizeof(float)*Nobjects[LINK]);
	int * li_InitialDepthType = (int *) malloc(sizeof(int)*Nobjects[LINK]); //
	float * lr_BreadthScale = (float *) malloc(sizeof(float)*Nobjects[LINK]); //
	float * lr_InitialDepth = (float *) malloc(sizeof(float)*Nobjects[LINK]); //

	//  node
	int * ni_idx = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_node_type = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_N_link_u = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_N_link_d = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_Mlink_u1 = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_Mlink_u2 = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_Mlink_u3 = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_Mlink_d1 = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_Mlink_d2 = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	int * ni_Mlink_d3 = (int *) malloc(sizeof(int)*Nobjects[NODE]);
	float * nr_Zbottom = (float *) malloc(sizeof(float)*Nobjects[NODE]); //

	int NNodes = Nobjects[NODE];
	int NLinks = Nobjects[LINK];

	float length_units;
	float manning_units;
	float flow_units;

	FILE *f_nodes;
	FILE *f_links;

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

	// Links
	for (int i = 0; i < Nobjects[LINK]; i++) {
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
			lr_Slope[i] = -SSIGN(h)*lr_Length[i] / (pow(lr_Length[i], 2) - pow(fabs(h), 2));
		}
		else {
			lr_Length[i] = 0;
			lr_Roughness[i] = 0;
			lr_Slope[i] = 0;
		}

		lr_InitialFlowrate[i] = Link[i].q0 * flow_units;
		lr_InitialUpstreamDepth[i] = Node[li_Mnode_u[i]].initDepth * length_units;
		lr_InitialDnstreamDepth[i] = Node[li_Mnode_d[i]].initDepth * length_units;
	}

	// Nodes
	for (int i = 0; i < Nobjects[NODE]; i++) {
		ni_idx[i] = i;
		ni_node_type[i] = Node[i].type;
	}


    #ifdef WINDOWS
	    fopen_s(&f_nodes, "nodes_info.csv", "w");
	    fopen_s(&f_links, "links_info.csv", "w");
    #else
	    f_nodes = fopen("nodes_info.csv", "w");
	    f_links = fopen("links_info.csv", "w");
    #endif

	fprintf(f_nodes,
		"n_left,node_id,ni_idx,ni_node_type,ni_N_link_u,ni_N_link_d,ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,ni_Mlink_d2,ni_Mlink_d3\n");
	for (int i = 0; i < NNodes; i++) {
		fprintf(f_nodes, "%d,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
			NNodes - i,
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
	for (int i = 0; i < NLinks; i++) {
		fprintf(f_links, "%d,%s,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
			NLinks - i,
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

	//  link
	free(li_idx);
	free(li_link_type);
	free(li_geometry);
	free(li_Mnode_u);
	free(li_Mnode_d);
	free(lr_Length);
	free(lr_Slope);
	free(lr_Roughness);
	free(lr_InitialFlowrate);
	free(lr_InitialUpstreamDepth);
	free(lr_InitialDnstreamDepth);
	free(li_InitialDepthType); //
	free(lr_BreadthScale); //
	free(lr_InitialDepth); //

	// //  node
	free(ni_idx);
	free(ni_node_type);
	free(ni_N_link_u);
	free(ni_N_link_d);
	free(ni_Mlink_u1);
	free(ni_Mlink_u2);
	free(ni_Mlink_u3);
	free(ni_Mlink_d1);
	free(ni_Mlink_d2);
	free(ni_Mlink_d3);
	free(nr_Zbottom); //
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
	int* ni_Mlink_d3);
