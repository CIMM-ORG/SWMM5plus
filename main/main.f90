program main
    !==========================================================================
    !% SWMM5+ Beta release, version 0.0.0-b1
    !% 20210831
    !% Hydraulics engine that links with EPA SWMM-C
    !% January XX 2022  
    !% see technical report at https://doi.org/10.18738/T8/WQZ5EX
    !% 
    !% Code authors:
    !% 2016-22  Prof. Ben Hodges -- lead developer and algorithm guru
    !% 2016-22  Eddie Tiernan -- network partitioning
    !% 2019-22  Sazzad Sharior -- SWMM hydraulic features
    !% 2019-22  Gerardo Riano-Briceno -- SWMM-C interface, intput and output
    !% 2019-22  Eric Jenkins -- SWMM-C interface, input and output
    !% 2020-22  Dr. Cheng-Wei (Justin) Yu -- large system test cases, compilation
    !% 2021-22  Christopher Bashear -- Profiling and testing
    !% 2021-22  Abdulmuttalib Lokhandwala -- testing
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
    !==========================================================================
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

    print *, 'exiting main'
!%
!%==========================================================================
!%==========================================================================
!%
end program main
