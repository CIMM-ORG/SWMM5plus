#include "utils.h"

int get_direct_inflow_times()
{
    int i;
    int pat;

    for(i=0; i<Nobjects[NODE]; i++)
    {
        pat = Node[i].extInflow->basePat;
    }
}

int get_dry_weather_inflow_times()
{

}

int get_wet_inflow_times()
{

}

int get_rdii_inflow_times()
{

}

int get_ifave_inflow_times()
{

}

int get_times_from_pattern(int pat)
{
    int pattern_type;
    int y0, mm0, d0;
    int y1, mm1, d1;
    int h0, m0, s0;
    int h1, m1, s1;
    int months;
    long total_time;

    total_time = datetime_timeDiff(StartDateTime, EndDateTime);
    pattern_type = Pattern[pat].type;

    datetime_decodeDate(StartDate, &y0, &mm0, &d0);
    datetime_decodeTime(StartTime, &h0, &m0, &s0);
    datetime_decodeDate(EndDate, &y1, &mm1, &d1);
    datetime_decodeTime(EndTime, &h1, &m1, &s1);

    if (pattern_type == MONTHLY_PATTERN)
    {
        months = m1 - m0;
    }
    else if (pattern_type == DAILY_PATTERN)
    {
        hours = 
    }
    else if (pattern_type == WEEKEND_PATTERN)
    {

    }
    else if (pattern_type == HOURLY_PATTERN)
    {

    }
}

//-----------------------------------------------------------------------------
//  Linked List Data Structure
//-----------------------------------------------------------------------------


/* Function to add a node at the beginning of Linked List.
   This function expects a pointer to the data to be added
   and size of the data type */
void push(struct UNode** head_ref, void *new_data, size_t data_size)
{
    // Allocate memory for node
    struct UNode* new_node = (struct UNode*)malloc(sizeof(struct UNode));

    new_node->data  = malloc(data_size);
    new_node->next = (*head_ref);

    // Copy contents of new_data to newly allocated memory.
    // Assumption: char takes 1 byte.
    int i;
    for (i=0; i<data_size; i++)
        *(char *)(new_node->data + i) = *(char *)(new_data + i);

    // Change head pointer as new node is added at the beginning
    (*head_ref)    = new_node;
}

/* Function to print nodes in a given linked list. fpitr is used
   to access the function to be used for printing current node data.
   Note that different data types need different specifier in printf() */
void printList(struct UNode *node, void (*fptr)(void *))
{
    while (node != NULL)
    {
        (*fptr)(node->data);
        node = node->next;
    }
}

// Function to print an integer
void printInt(void *n)
{
   printf(" %d", *(int *)n);
}

// Function to print a float
void printFloat(void *f)
{
   printf(" %f", *(float *)f);
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