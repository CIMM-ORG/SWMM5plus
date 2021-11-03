#ifndef API_ERROR_H
#define API_ERROR_H

#define MAXMSG 1024 // Max. # characters in message text

// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
  #define WINDOWS
#endif
#ifdef __WIN32__
  #define WINDOWS
#endif

// --- define DLLEXPORT

#ifdef WINDOWS
    #define DLLEXPORT __declspec(dllexport) __stdcall
#else
    #define DLLEXPORT
#endif

// Enums written in caps are extracted from native
// EPA-SWMM, whereas lower case vars are added to EPA-SWMM

// Interface error codes:
enum api_error_codes {
    api_err_not_initialized = 1001, //1001 1
    api_err_not_developed,          //1002 2
    api_err_wrong_type,             //1003 3
    api_err_wrong_parameter,        //1004 4
    api_err_model_junctions,        //1005 5
    api_err_internal,               //1006 6
    api_max_err_msg
};

char APIErrMsg[MAXMSG+1];

#ifdef __cplusplus
extern "C" {
#endif

char* api_error_getMsg(int i);
int  api_error_getCode(int i);
int  api_error_setInpError(int api_errcode, char* s);
void api_report_writeErrorMsg(int code, char* s);
void api_report_writeErrorCode();

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif //API_ERROR_H