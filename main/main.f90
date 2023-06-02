program main
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Primary contact: Ben R. Hodges, email: benrhodges@gmail.com
    !% 
    !% Code authors (current):
    !% 2016-23  Prof. Ben Hodges -- lead developer and algorithm guru
    !% 2019-23  Sazzad Sharior -- SWMM hydraulic features and Preissmann Slot
    !% 2019-23  Eric Jenkins -- SWMM-C interface, HDF5 input and output
    !% 2016-22  Eddie Tiernan -- network partitioning
    !% 2021-23  Christopher Bashear -- Profiling and testing
    !% 2021-23  Abdulmuttalib Lokhandwala -- testing
    !% 
    !% Code authors (prior development)
    !% 2019-22  Gerardo Riano-Briceno -- SWMM-C interface, link/node translation
    !% 2020-22  Dr. Cheng-Wei (Justin) Yu -- large system test cases, compilation
    !% 2018-20  Dr. Ehsan Madidi-Kandjani -- preliminary development
    !%
    !% This code was developed under Cooperative Agreement No. 83595001 awarded 
    !% by the U.S. Environmental Protection Agency to The University of Texas at 
    !% Austin. It has not been formally reviewed by EPA. The views expressed in 
    !% this document are solely those of the authors and do not necessarily 
    !% reflect those of the Agency. EPA does not endorse any products or 
    !% commercial services mentioned in this publication.
    !% 
    !% This is free and unencumbered software released into the public domain.
    !%
    !% Anyone is free to copy, modify, publish, use, compile, sell, or
    !% distribute this software, either in source code form or as a compiled
    !% binary, for any purpose, commercial or non-commercial, and by any
    !% means.
    !%
    !% In jurisdictions that recognize copyright laws, the author or authors
    !% of this software dedicate any and all copyright interest in the
    !% software to the public domain. We make this dedication for the benefit
    !% of the public at large and to the detriment of our heirs and
    !% successors. We intend this dedication to be an overt act of
    !% relinquishment in perpetuity of all present and future rights to this
    !% software under copyright law.
    !% 
    !% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    !% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    !% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    !% IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
    !% OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
    !% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    !% OTHER DEALINGS IN THE SOFTWARE.
    !%
    !% For more information, please refer to <http://unlicense.org/>
    !%==========================================================================
    !
    use initialization
    use timeloop
    use finalization

    implicit none
!%
!%==========================================================================
!%==========================================================================
!%
!% MAIN PROGRAM
!%
    call initialize_toplevel()

    call timeloop_toplevel()

    call finalize_toplevel()

    print *, 'exiting main program; SWMM5+ completed'
!%
!%==========================================================================
!%==========================================================================
!%
end program main
