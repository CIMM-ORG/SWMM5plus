//=============================================================================
// FUNCTIONS FOR SWMM5+ ACCESS TO EPA-SWMM CONTROLS
// The following should be appended to the bottom of the controls.c file
// 20220627
//=============================================================================
int controls_display(void) // brh20220617
    /// Purpose: displays controls for developing API
    /// Author B R Hodges  20220626
{
    int     rr;
    struct  TPremise *thisPremise;
    struct  TAction  *thisThenAction, *thisElseAction;

    for ( rr=0; rr<RuleCount; rr++ )
    {
        thisPremise = Rules[rr].firstPremise;
        while (thisPremise)
        {
            // printf(" Premises for Rule rr = %d, for Rules[rr].ID %s \n ",rr, Rules[rr].ID);

            // printf(" node lhs %d \n ",thisPremise->lhsVar.node);

            // printf(" link lhs %d \n ",thisPremise->lhsVar.link);

            // printf(" attr lhs %d \n ",thisPremise->lhsVar.attribute);

            // printf(" node rhs %d \n ",thisPremise->rhsVar.node);

            // printf(" link rhs %d \n ",thisPremise->rhsVar.link);

            // printf(" attr rhs %d \n ",thisPremise->rhsVar.attribute);

            thisPremise = thisPremise->next;
        }

        thisThenAction = Rules[rr].thenActions;

        while (thisThenAction)
        {
            // printf(" Then Actions for Rule rr = %d, for Rules[rr].ID %s \n ",rr,Rules[rr].ID);

            // printf(" rule %d \n ",thisThenAction->rule);

            // printf(" link %d \n ",thisThenAction->link);

            // printf(" attr %d \n ",thisThenAction->attribute);

            thisThenAction = thisThenAction->next;
        }

        thisElseAction = Rules[rr].elseActions;

        while (thisElseAction)
        {
            // printf(" Else Actions for Rule rr = %d, for Rules[rr].ID %s \n ",rr,Rules[rr].ID);

            // printf(" rule %d \n ",thisElseAction->rule);

            // printf(" link %d \n ",thisElseAction->link);

            // printf(" attr %d \n ",thisElseAction->attribute);

            thisElseAction = thisElseAction->next;
        }
    }

    //  printf(" node %d \n ",nextPremise->lhsVar.node);
    //  printf(" link %d \n ",nextPremise->lhsVar.link);
    //  printf(" attr %d \n ",nextPremise->lhsVar.attribute);

    return 0;
}

//=============================================================================
int controls_count_rules(void)
    /// Author B R Hodges  20220626
{
    return RuleCount;
}
//=============================================================================
int controls_count_premise(void)
    /// Output: returns the total number of premise locations (node or element)
    /// Purpose: Provides count used to initialize coarray in SWMM5+
    /// Author B R Hodges  20220626
    ///
    /// HACK -- should be combined with controls_count_thenAction and ... elseAction
    /// for a single function called with different parts of the rule structure.
{
    int     rr;
    int     ncount = 0;
    struct  TPremise *thisPremise;

    if (Rules == NULL) return 0;

    for ( rr=0; rr<RuleCount; rr++ )
    {

        thisPremise = Rules[rr].firstPremise;
        while (thisPremise)
        {
            ncount += 2; // each premise has lhs and rhs with possibly different locations
            thisPremise = thisPremise->next;
        }
    }
    return ncount;
}
//=============================================================================
int controls_count_thenAction(void)
    ///
    /// Output: returns the total number of "then" action locations (node or element)
    /// Purpose: Provides count used to initialize coarray in SWMM5+
    /// Author B R Hodges  20220626
    ///
{
    int    rr;
    int     ncount = 0;
    struct  TAction  *thisThenAction;

    if (Rules == NULL) return 0;

    for ( rr=0; rr<RuleCount; rr++ )
    {
        thisThenAction = Rules[rr].thenActions;
        while (thisThenAction)
        {
            ncount++;  // each action has only one location
            thisThenAction = thisThenAction->next;
        }
    }    
    return ncount;
}
//=============================================================================
int controls_count_elseAction(void)
    ///
    /// Output: returns the total number of "else" action locations (node or element)
    /// Purpose: Provides count used to initialize coarray in SWMM5+
    /// Author B R Hodges  20220626
    ///
{
    int     rr;
    int     ncount = 0;
    struct  TAction  *thisElseAction;

    //printf(" in here %d \n ",ncount);

    if (Rules == NULL) return 0;

    for ( rr=0; rr<RuleCount; rr++ )
    {

        thisElseAction = Rules[rr].elseActions;
        //printf(" in here %d \n ",ncount);
        while (thisElseAction)
        {
            ncount++; 
            thisElseAction = thisElseAction->next;
        }

        //printf(" in here %d \n ",ncount);
    }    
    return ncount;
}
//=============================================================================
int controls_get_premise_data(
    int* locationL,        int* locationR,
    int* islinkL,          int* islinkR,
    int* attributeL,       int* attributeR, 
    int* thisPremiseLevel, int rIdx)
    ///
    /// Input:   dummy values for location, islink, attribute for premise data
    ///          thisPremiseLevel is the Premise level we're trying to extract data from
    ///          rIdx is the rule index
    /// Output:  location (index) in link/node space, islink is whether a link or node
    ///          attribute is keyword (number) for monitor data type
    ///          thisPremiseLevel is incremented (if successful) to setup for next call  
    ///          Returns 0 if thisPremiseLevel does not exist in linked list for Rule(rIdx)
    ///          Returns 1 if thisPremiseLevel exists  
    /// Purpose: Gets data from premise linked list in control rule
    /// Author B R Hodges  20220626
{
    int pp=0;
    struct  TPremise *thisPremise;

    if (Rules == NULL) return 0;

    //printf(" AA thisPremiseLevel %d \n ",*thisPremiseLevel);

    // point at the first premise of the rule
    thisPremise = Rules[rIdx].firstPremise;

    // increment the linked list to get to the working premise level
    while ((thisPremise) && (pp < *thisPremiseLevel))
    {
        pp++;
        thisPremise = thisPremise->next;
    }

    //printf(" BB pp, thisPremiseLevel %d %d \n ",pp, *thisPremiseLevel);

    // Get the data for this premise level
    if (thisPremise)
    {
         // the above returns pp = thisPremiseLevel, so
         // increment so that next call handles next premise level
        *thisPremiseLevel = pp+1;

        // printf(" CC pp, thisPremiseLevel %d %d \n ",pp, *thisPremiseLevel);

        // get data for the LHS of this premise
        *attributeL = thisPremise->lhsVar.attribute;
        if (thisPremise->lhsVar.link > -1) // is link
        {
            *islinkL = 1;
            *locationL = thisPremise->lhsVar.link;
        }
        else
        {
            if (thisPremise->lhsVar.node > -1) //is node
            {
                *islinkL = 0;
                *locationL = thisPremise->lhsVar.node;
            }
            else // neither link nor node
            {
                *islinkL = -1;
                *locationL = -1;
            }
        }

        // get data for the RHS of this premise
        *attributeR = thisPremise->rhsVar.attribute;
        //printf(" RHS link and node %d %d \n ",thisPremise->rhsVar.link, thisPremise->rhsVar.node);
        if (thisPremise->rhsVar.link > -1) // is link
        {
            *islinkR = 1;
            *locationR = thisPremise->rhsVar.link;
        }
        else
        {
            if (thisPremise->rhsVar.node > -1) //is node
            {
                *islinkR = 0;
                *locationR = thisPremise->rhsVar.node;
            }
            else  // neither link nor node
            {
                *islinkR = -1;
                *locationR = -1;
            }
        }
    }    
    else
    {
        // printf(" DD pp, thisPremiseLevel %d %d \n ",pp, *thisPremiseLevel);
        return 0;  // no premise data found
    }
    return 1;  // premise data found
}
//=============================================================================
int controls_get_action_data(
    int* location,  int* attribute,  int* thisActionLevel, int rIdx, int isThen)
    /// Author B R Hodges  20220626
{
    int pp=0;
    struct  TAction *thisAction;

    if (Rules == NULL) return 0;

    // point at the first action of this rule
    if (isThen)
    {
        thisAction = Rules[rIdx].thenActions;
    }
    else
    {
        thisAction = Rules[rIdx].elseActions;
    }
    

    // increment linked list to get to the working action level
    while ((thisAction) && (pp < *thisActionLevel))
    {
        pp++;
        thisAction = thisAction->next;
    }
    
    // get the data for this action level
    if (thisAction)
    {
        // the above returns pp = thisAciton level, so
        // increment so that next call handles the next action level
        *thisActionLevel = pp+1;

        *attribute = thisAction->attribute;
        *location  = thisAction->link;
    }
    else
    {
        return 0; // no action data found
    }
    return 1; // action data found
}
//=============================================================================
// END ADDITIONAL FUNCTIONS FOR CONTROL.C
//=============================================================================