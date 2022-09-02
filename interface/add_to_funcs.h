//=============================================================================
// FUNCTIONS FOR SWMM5+ ACCESS TO EPA-SWMM CONTROLS
// The following should be appended to the bottom of the funcs.h file
// 20220627 
//=============================================================================
int     controls_display(void);
int     controls_count_rules(void); 
int     controls_count_premise(void); 
int     controls_count_thenAction(void); 
int     controls_count_elseAction(void); 
int     controls_get_premise_data(
                int* locationL,        int* locationR,
                int* linknodesimTypeL,          int* linknodesimTypeR,
                int* attributeL,       int* attributeR, 
                int* thisPremiseLevel, int rIdx);
int     controls_get_action_data(
                int* location, int* attribute, int* thisActionLevel, int rIdx, int isThen);

