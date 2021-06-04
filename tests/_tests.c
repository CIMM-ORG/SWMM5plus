#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "swmm5.h"
#include "headers.h"
#include "_tests.h"

void DLLEXPORT interface_print_inflow(char* node_name)
{
    int j;
    j = project_findObject(NODE, node_name);
    printf("%f FLOW\n", Node[j].inflow);
}

void DLLEXPORT interface_print_tseries_props(int k)
{
    // k : Node index
    int j = -1;

    if (k > Nobjects[NODE])
        printf("There are %d nodes only\n", Nobjects[NODE]);

    if (Node[k].dwfInflow != NULL)
    {
        j = Node[k].dwfInflow->patterns[0];
        if (j == -1)
        {
            printf("\nNone %d\n", k);
        }
        else
        {
            printf("\nNode [%s]\n", Node[k].ID);
            printf("\nTseries[%d].curveType = %d\n", k, Tseries[j].curveType);
            printf("\nTseries[%d].refersTo = %d\n", k, Tseries[j].refersTo);
        }
    }

}

void DLLEXPORT print_tseries_name(int k)
{
    printf("TSERIES: %s\n", Tseries[k].file.name);
}