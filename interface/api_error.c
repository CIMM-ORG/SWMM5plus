// FILE api_error.c
//==========================================================================
// SWMM5+ release, version 1.0.0
// 20230608
// Hydraulics engine that links with EPA SWMM-C
// June 8, 2023
//
// Description:
// error reporting for the SWMM5+/EPA-SWMM API
//==========================================================================
#define _CRT_SECURE_NO_DEPRECATE

#include <string.h>
#include <stdlib.h>

#include "headers.h"
#include "api_error.h"

////////////////////////////////////////////////////////////////////////////
//  NOTE: Need to update ApiErrorMsgs[] and ApiErrorCodes[]
//        (in api_error.h) whenever a new error message is added.
/////////////////////////////////////////////////////////////////////////////

#define WRITE(x) (report_writeLine((x)))
#define API_FIRST_ERROR_CODE 1000
#define api_err1001 "\n API Init Error: SWMM5 needs to be initialized before executing %s"
#define api_err1002 "\n API Not Developed Error: have not developed the logic to handle %s"
#define api_err1003 "\n API Type Error: Invalid property given object type. %s"
#define api_err1004 "\n API Parameter Error: Invalid parameter %s"
#define api_err1005 "\n API Model Error: The number of incoming or outgoing links is incompatible with SWMM5+ at junction %s"
#define api_err1006 "\n API Internal Error: Internal Error %s"

char* APIErrMsgs[] =
    { "",     api_err1001, api_err1002, api_err1003, api_err1004, api_err1005, api_err1006};

int APIErrCodes[] =
    { 0,      1001,        1002,        1003,        1004,        1005,        1006};

char  APIErrString[256];
int APIErrCode = 0;

//=============================================================================
char* api_error_getMsg(int i)
{
    if ( i >= 0 && i < api_max_err_msg ) return APIErrMsgs[i-API_FIRST_ERROR_CODE];
    else return APIErrMsgs[0];
};

//=============================================================================
int  api_error_getCode(int i)
{
    if ( i >= 0 && i < api_max_err_msg ) return APIErrCodes[i-API_FIRST_ERROR_CODE];
    else return 0;
}

//=============================================================================
int  api_error_setInpError(int api_errcode, char* s)
{
    strcpy(APIErrString, s);
    return api_errcode;
}

//=============================================================================
//      ERROR REPORTING
//=============================================================================

//=============================================================================
void DLLEXPORT api_display_ErrorMsg(int code, char* s)
{
    if (( code >= api_err_not_developed ) && ( code < api_max_err_msg ))
        printf(APIErrString, api_error_getMsg(code), s);
}

//=============================================================================
void DLLEXPORT api_display_ErrorCode(int code, char* s)
{
    if (( code >= api_err_not_developed ) && ( code < api_max_err_msg ))
        printf("%s", error_getMsg(code));
}

//=============================================================================
void api_report_writeErrorMsg(int code, char* s)
//
//  Input:   code = error code
//           s = error message text
//  Output:  none
//  Purpose: writes error message to report file.
//
{
    if ( Frpt.file )
    {
        WRITE("");
        fprintf(Frpt.file, api_error_getMsg(code), s);
    }

    APIErrCode = code;

    if ( APIErrCode >= api_err_not_developed )
    {
        sprintf(APIErrMsg, api_error_getMsg(APIErrCode), s);
    }
}

//=============================================================================
void api_report_writeErrorCode()
//
//  Input:   none
//  Output:  none
//  Purpose: writes error message to report file.
//
{
    if ( Frpt.file )
    {
        if  (( APIErrCode >= api_err_not_developed ) &&
             ( APIErrCode < api_max_err_msg ))
            fprintf(Frpt.file, "%s", error_getMsg(APIErrCode));
    }
}
//=============================================================================
// END api_error.c
//=============================================================================