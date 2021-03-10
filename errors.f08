module errors

    implicit none
    public

    ! Interface errors
    integer, parameter :: ERROR_FEATURE_NOT_COMPATIBLE = 100001
    character(40) :: MSG_FEATURE_NOT_COMPATIBLE = "ERROR [100001]: feature not compatible"
    integer, parameter :: ERROR_API_NOT_INITIALIZED = 100002
    character(40) :: MSG_API_NOT_INITIALIZED = "ERROR [100002]: API not initialized"
    integer, parameter :: ERROR_API_TSERIES_HANDLING_ERROR = 100003
    character(40) :: MSG_API_TSERIES_HANDLING_ERROR = "ERROR [100003]: API can't get Tseries"

    ! Parameter errors
    integer, parameter :: ERROR_INCORRECT_PARAMETER = 200001
    character(40) :: MSG_INCORRECT_PARAMETER = "ERROR [200001]: incorrect parameter"

end module errors