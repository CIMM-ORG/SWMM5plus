module interface_
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Application Programming Interface (API) for SWMM5+ to call functions
    !% in api.c
    !%
    !%==========================================================================
    use iso_c_binding
    use c_library
    use utility
    use utility_datetime
    use define_indexes
    use define_keys
    use define_api_keys
    use define_globals
    use define_settings, only: setting
    use utility_crash

    implicit none

    private

!%==========================================================================
!% PUBLIC
!%==========================================================================

    type(c_lib_type), public :: c_lib
    logical, public :: api_is_initialized = .false.

    public :: interface_teststuff
    public :: interface_controls_count
    public :: interface_controls_get_premise_data
    public :: interface_controls_get_action_data
    public :: interface_controls_transfer_monitor_data
    public :: interface_controls_execute
    public :: interface_controls_get_action_results

    !% --- Public subroutines/functions
    public :: interface_init
    public :: interface_finalize
    public :: interface_get_nodef_attribute
    public :: interface_get_linkf_attribute

    public :: interface_get_adjustments

    public :: interface_get_transectf_attribute
    public :: interface_get_transect_table
    public :: interface_get_N_TRANSECT_TBL

    public :: interface_get_table_attribute
    public :: interface_get_num_table_entries
    public :: interface_get_first_entry_table
    public :: interface_reset_timeseries_to_start
    public :: interface_get_next_entry_table
    public :: interface_get_obj_name_len
    public :: interface_update_linknode_names
    public :: interface_update_transectID_names
    public :: interface_get_BC_resolution
    public :: interface_get_next_inflow_time
    public :: interface_get_next_head_time
    public :: interface_get_flowBC
    public :: interface_get_headBC
   
    public :: interface_find_object
    public :: interface_update_nodeResult
    public :: interface_update_linkResult
    public :: interface_write_output_line
    public :: interface_export_link_results
    public :: interface_export_runon_volume
    public :: interface_call_runoff_execute
    public :: interface_get_subcatch_runoff
    public :: interface_get_subcatch_runoff_nodeIdx
    public :: interface_get_NewRunoffTime
    public :: interface_call_climate_setState
    public :: interface_get_evaporation_rate
    public :: interface_get_RDII_inflow
    public :: interface_get_groundwater_inflow
    public :: interface_get_LID_inflow
!%
!%==========================================================================
!% INTERFACE
!%==========================================================================
!%
    ! --- Interface with SWMM shared library
    interface
        !% -----------------------------------------------------------------
        !% dummy for testing
        !%------------------------------------------------------------------
        integer(c_int) function api_teststuff() &
            BIND(C, name='api_teststuff')
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_teststuff
        !%------------------------------------------------------------------
        !% Controls and monitoring
        !%------------------------------------------------------------------
        integer(c_int) function api_controls_count(nRules, nPremise, nThenAction, nElseAction) &
            BIND(C, name="api_controls_count")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), intent(inout) :: nRules
            integer(c_int), intent(inout) :: nPremise
            integer(c_int), intent(inout) :: nThenAction
            integer(c_int), intent(inout) :: nElseAction
        end function api_controls_count
        !%------------------------------------------------------------------
        integer(c_int) function api_controls_get_premise_data( &
                locationL,        locationR,                   &
                linknodesimTypeL, linknodesimTypeR,            &
                attributeL,       attributeR,                  & 
                thisPremiseLevel, rIdx)                        &
            BIND(C, name="api_controls_get_premise_data")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), intent(inout) :: locationL,  locationR
            integer(c_int), intent(inout) :: linknodesimTypeL, linknodesimTypeR
            integer(c_int), intent(inout) :: attributeL, attributeR
            integer(c_int), intent(inout) :: thisPremiseLevel
            integer(c_int), value, intent(in)    :: rIdx
        end function api_controls_get_premise_data
        !%------------------------------------------------------------------
        integer (c_int) function api_controls_get_action_data(              &
                location, attribute, thisActionLevel, rIdx, isThen)         &
            BIND(C, name="api_controls_get_action_data")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), intent(inout) :: location, attribute, thisActionLevel
            integer(c_int), value, intent(in) :: rIdx, isThen
        end function api_controls_get_action_data
       !%------------------------------------------------------------------
        integer(c_int) function api_controls_transfer_monitor_data ( &
                Depth, Head, Volume, Inflow, Flow, StatusSetting,    &
                TimeLastSet, LinkNodeIdx, linknodesimType)                    &
            BIND(C, name="api_controls_transfer_monitor_data")
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), value, intent(in) :: Depth, Head, Volume, Inflow, Flow
            real(c_double), value, intent(in) :: StatusSetting, TimeLastSet
            integer(c_int), value, intent(in) :: LinkNodeIdx, linknodesimType
        end function api_controls_transfer_monitor_data
        !%------------------------------------------------------------------
        integer (c_int) function api_controls_execute( &
                currentTimeEpoch, ElapsedDays, dtDays) &
            BIND(C, name="api_controls_execute")
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), value, intent(in) :: currentTimeEpoch, ElapsedDays, dtDays
        end function api_controls_execute
        !%------------------------------------------------------------------

        !%------------------------------------------------------------------
        !% Simulation
        !%------------------------------------------------------------------
        integer(c_int) function api_initialize(inp_file, report_file, out_file, run_routing) &
            BIND(C, name="api_initialize")
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char), intent(in) :: inp_file
            character(kind=c_char), intent(in) :: report_file
            character(kind=c_char), intent(in) :: out_file
            integer(c_int),  value, intent(in) :: run_routing
        end function api_initialize
        !%------------------------------------------------------------------
        integer(c_int) function api_finalize() &
            BIND(C, name="api_finalize")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_finalize
        !%------------------------------------------------------------------
        real(c_double) function api_run_step() &
            BIND(C, name="api_run_step")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_run_step

        !%------------------------------------------------------------------
        ! --- Property-extraction
        !%------------------------------------------------------------------
        !% * During Simulation
        integer(c_double) function api_get_node_results(node_name, inflow, overflow, depth, volume) &
            BIND(C, name="api_get_node_results")
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char), intent(in   ) :: node_name
            real(c_float),          intent(inout) :: inflow
            real(c_float),          intent(inout) :: overflow
            real(c_float),          intent(inout) :: depth
            real(c_float),          intent(inout) :: volume
        end function api_get_node_results
        !%------------------------------------------------------------------
        integer(c_double) function api_get_link_results(link_name, flow, depth, volume) &
            BIND(C, name="api_get_link_results")
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char), intent(in   ) :: link_name
            real(c_float),          intent(inout) :: flow
            real(c_float),          intent(inout) :: depth
            real(c_float),          intent(inout) :: volume
        end function api_get_link_results
        !%------------------------------------------------------------------
        !% * After Initialization
        real(c_double) function api_get_start_datetime() &
            BIND(C, name="api_get_start_datetime")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_get_start_datetime
        !%------------------------------------------------------------------
        real(c_double) function api_get_end_datetime() &
            BIND(C, name="api_get_end_datetime")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_get_end_datetime
        !%------------------------------------------------------------------
        integer(c_int) function api_get_flowBC(node_idx, current_datetime, flowBC) &
            BIND(C, name="api_get_flowBC")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: node_idx
            real(c_double), value, intent(in   ) :: current_datetime
            real(c_double),        intent(inout) :: flowBC
        end function api_get_flowBC
        !%------------------------------------------------------------------
        integer(c_int) function api_get_headBC(node_idx, current_datetime, headBC) &
            BIND(C, name="api_get_headBC")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: node_idx
            real(c_double), value, intent(in   ) :: current_datetime
            real(c_double),        intent(inout) :: headBC
        end function api_get_headBC
        !%------------------------------------------------------------------
        integer(c_int) function api_count_subobjects(N_groundwater) &
            BIND(C, name="api_count_subobjects")
            use, intrinsic :: iso_c_binding
            implicit none 
            integer(c_int),       intent(inout) :: N_groundwater
        end function api_count_subobjects
        !%------------------------------------------------------------------
        integer(c_int) function api_get_SWMM_setup( &
            flow_units, &
            infiltration_model_type, &
            flow_routing_model_type, &
            link_offset_type, &
            force_main_equation, &
            ignore_rainfall, &
            ignore_snowmelt, &
            ignore_groundwater, &
            ignore_rdii, &
            ignore_routing, &
            ignore_quality, &
            allow_ponding, &
            steadystate_skip, &
            steadystate_system_flow_tolerance, &
            steadystate_lateral_flow_tolerance, &
            routing_step_lengthening_time, &
            routing_step_Courant_factor, &
            routing_step_minimum,  &
            inertial_damping_type, &
            normal_flow_limiter_type, &
            minimum_surface_area, &
            minimum_conduit_slope, &
            maximum_number_of_trials, &
            head_convergence_tolerance, &
            number_parallel_threads, &
            tempdir_provided, &
            control_rule_step, &
            surcharge_method) &
            BIND(C, name="api_get_SWMM_setup")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), intent(inout) :: flow_units
            integer(c_int), intent(inout) :: infiltration_model_type
            integer(c_int), intent(inout) :: flow_routing_model_type
            integer(c_int), intent(inout) :: link_offset_type
            integer(c_int), intent(inout) :: force_main_equation
            integer(c_int), intent(inout) :: ignore_rainfall
            integer(c_int), intent(inout) :: ignore_snowmelt
            integer(c_int), intent(inout) :: ignore_groundwater
            integer(c_int), intent(inout) :: ignore_rdii
            integer(c_int), intent(inout) :: ignore_routing
            integer(c_int), intent(inout) :: ignore_quality
            integer(c_int), intent(inout) :: allow_ponding
            integer(c_int), intent(inout) :: steadystate_skip
            real(c_double), intent(inout) :: steadystate_system_flow_tolerance
            real(c_double), intent(inout) :: steadystate_lateral_flow_tolerance
            real(c_double), intent(inout) :: routing_step_lengthening_time
            real(c_double), intent(inout) :: routing_step_Courant_factor
            real(c_double), intent(inout) :: routing_step_minimum
            integer(c_int), intent(inout) :: inertial_damping_type
            integer(c_int), intent(inout) :: normal_flow_limiter_type
            real(c_double), intent(inout) :: minimum_surface_area
            real(c_double), intent(inout) :: minimum_conduit_slope
            integer(c_int), intent(inout) :: maximum_number_of_trials
            real(c_double), intent(inout) :: head_convergence_tolerance
            integer(c_int), intent(inout) :: number_parallel_threads
            integer(c_int), intent(inout) :: tempdir_provided
            integer(c_int), intent(inout) :: control_rule_step
            integer(c_int), intent(inout) :: surcharge_method
        end function api_get_SWMM_setup
        !%------------------------------------------------------------------
        integer(c_int) function api_get_SWMM_times &
            (starttime_epoch, endtime_epoch, report_start_datetime, report_step, &
             hydrology_wet_step, hydrology_dry_step, sweep_start_dayofyear, &
             sweep_end_dayofyear, dry_days, hydraulic_step, total_duration) &
            BIND(C, name="api_get_SWMM_times")
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), intent(inout) :: starttime_epoch
            real(c_double), intent(inout) :: endtime_epoch
            real(c_double), intent(inout) :: report_start_datetime
            integer(c_int), intent(inout) :: report_step
            integer(c_int), intent(inout) :: hydrology_wet_step
            integer(c_int), intent(inout) :: hydrology_dry_step
            integer(c_int), intent(inout) :: sweep_start_dayofyear
            integer(c_int), intent(inout) :: sweep_end_dayofyear
            integer(c_int), intent(inout) :: dry_days
            real(c_double), intent(inout) :: hydraulic_step
            real(c_double), intent(inout) :: total_duration
        end function api_get_SWMM_times
        !%------------------------------------------------------------------
        real(c_double) function api_get_NewRunoffTime() &
            BIND(C, name="api_get_NewRunoffTime")
            use, intrinsic :: iso_c_binding
            implicit none 
        end function api_get_NewRunoffTime
        !%------------------------------------------------------------------
        integer(c_int) function api_get_nodef_attribute(node_idx, attr, value) &
            BIND(C, name="api_get_nodef_attribute")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: node_idx
            integer(c_int), value, intent(in   ) :: attr
            real(c_double),        intent(inout) :: value
        end function api_get_nodef_attribute
        !%------------------------------------------------------------------
        integer(c_int) function api_get_linkf_attribute(link_idx, attr, value) &
            BIND(C, name="api_get_linkf_attribute")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: link_idx
            integer(c_int), value, intent(in   ) :: attr
            real(c_double),        intent(inout) :: value
        end function api_get_linkf_attribute
        !%------------------------------------------------------------------
        integer(c_int) function api_get_transectf_attribute(transect_idx, attr, value) &
            BIND(C, name="api_get_transectf_attribute")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: transect_idx
            integer(c_int), value, intent(in   ) :: attr
            real(c_double),        intent(inout) :: value
        end function api_get_transectf_attribute
        !%------------------------------------------------------------------
        integer (c_int) function api_get_adjustments &
        (adj_len, adjTemperature, adjEvaporation, adjRainfall, adjConductivity) &
            BIND(C, name="api_get_adjustments")
            use, intrinsic :: iso_c_binding
            implicit none 
            integer(c_int), value, intent(in   ) :: adj_len
            real(c_double),        intent(inout) :: adjTemperature(adj_len)
            real(c_double),        intent(inout) :: adjEvaporation(adj_len)
            real(c_double),        intent(inout) :: adjRainfall(adj_len)
            real(c_double),        intent(inout) :: adjConductivity(adj_len)
        end function api_get_adjustments
        !%------------------------------------------------------------------
        integer (c_int) function api_get_N_TRANSECT_TBL() &
            BIND(C, name="api_get_N_TRANSECT_TBL")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_get_N_TRANSECT_TBL
        !%------------------------------------------------------------------
        integer (c_int) function api_get_transect_table &
            (transect_idx, table_len, tarea, twidth, thydradius) &
            BIND(C, name="api_get_transect_table")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: transect_idx
            integer(c_int), value, intent(in   ) :: table_len
            real(c_double),        intent(inout) :: tarea(table_len)
            real(c_double),        intent(inout) :: twidth(table_len)
            real(c_double),        intent(inout) :: thydradius(table_len)
        end function api_get_transect_table
        !%------------------------------------------------------------------
        integer(c_int) function api_get_num_objects(obj_type) &
            BIND(C, name="api_get_num_objects")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: obj_type
        end function api_get_num_objects
        !%------------------------------------------------------------------
        integer(c_int) function api_get_object_name_len(object_idx, object_type, len_value) &
            BIND(C, name="api_get_object_name_len")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: object_idx
            integer(c_int), value, intent(in   ) :: object_type
            integer(c_int),        intent(inout) :: len_value
        end function api_get_object_name_len
        !%------------------------------------------------------------------
        integer(c_int) function api_get_object_name(object_idx, object_name, object_type) &
            BIND(C, name="api_get_object_name")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value,  intent(in   ) :: object_idx
            character(kind=c_char), intent(inout) :: object_name
            integer(c_int), value,  intent(in   ) :: object_type
        end function api_get_object_name
        !%------------------------------------------------------------------
        integer(c_int) function api_get_num_table_entries(table_idx, table_type, num_entries) &
            BIND(C, name="api_get_num_table_entries")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: table_idx
            integer(c_int), value, intent(in   ) :: table_type
            integer(c_int),        intent(inout) :: num_entries
        end function api_get_num_table_entries
        !%------------------------------------------------------------------
        integer(c_int) function api_get_table_attribute(table_idx, attr, value) &
            BIND(C, name="api_get_table_attribute")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: table_idx
            integer(c_int), value, intent(in   ) :: attr
            real(c_double),        intent(inout) :: value
        end function api_get_table_attribute
        !%------------------------------------------------------------------
        integer(c_int) function api_get_first_entry_table(table_idx, table_type, x, y) &
            BIND(C, name="api_get_first_entry_table")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: table_idx
            integer(c_int), value, intent(in   ) :: table_type
            real(c_double),        intent(inout) :: x
            real(c_double),        intent(inout) :: y
        end function api_get_first_entry_table
        !%------------------------------------------------------------------
        integer(c_int) function api_get_next_entry_table(table_idx, table_type, x, y) &
            BIND(C, name="api_get_next_entry_table")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: table_idx
            integer(c_int), value, intent(in   ) :: table_type
            real(c_double),        intent(inout) :: x
            real(c_double),        intent(inout) :: y
        end function api_get_next_entry_table
        !%------------------------------------------------------------------
        integer(c_int) function api_get_next_entry_tseries(tseries_idx, timemax) &
            BIND(C, name="api_get_next_entry_tseries")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: tseries_idx
            real(c_double), value, intent(in   ) :: timemax
        end function api_get_next_entry_tseries
        !%------------------------------------------------------------------
        integer(c_int) function api_reset_timeseries_to_start(tseries_idx) &
            BIND(C, name="reset_timeseries_to_start")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: tseries_idx
        end function api_reset_timeseries_to_start
        !%------------------------------------------------------------------
        ! --- Output Writing (Post Processing)
        !%------------------------------------------------------------------
        integer(c_int) function api_write_output_line(t) &
            BIND(C, name="api_write_output_line")
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), value, intent(in) :: t
        end function api_write_output_line
        !%------------------------------------------------------------------
        integer(c_int) function api_update_nodeResult(node_idx, resultType, newNodeResult) &
            BIND(C, name="api_update_nodeResult")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: node_idx
            integer(c_int), value, intent(in) :: resultType
            real(c_double), value, intent(in) :: newNodeResult
        end function api_update_nodeResult
        !%------------------------------------------------------------------
        integer(c_int) function api_update_linkResult(link_idx, resultType, newLinkResult) &
            BIND(C, name="api_update_linkResult")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: link_idx
            integer(c_int), value, intent(in) :: resultType
            real(c_double), value, intent(in) :: newLinkResult
        end function api_update_linkResult
        !%------------------------------------------------------------------
        ! --- Print-out
        !%------------------------------------------------------------------
        integer(c_int) function api_export_linknode_properties(units) &
            BIND(C, name="api_export_linknode_properties")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: units
        end function api_export_linknode_properties
        !%------------------------------------------------------------------
        integer(c_int) function api_export_link_results(link_idx) &
            BIND(C, name="api_export_link_results")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: link_idx
        end function api_export_link_results
        !%------------------------------------------------------------------
        integer(c_int) function api_export_node_results(node_idx) &
            BIND(C, name="api_export_node_results")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: node_idx
        end function api_export_node_results
        !%------------------------------------------------------------------
        ! --- Hydrology
        !%------------------------------------------------------------------   
        integer(c_int) function api_export_runon_volume(outfall_idx,flowrate) &
            BIND(C, name='api_export_runon_volume')
            use, intrinsic :: iso_c_binding
            implicit none 
            integer(c_int), value, intent(in) :: outfall_idx
            real(c_double), value, intent(in) :: flowrate
        end function api_export_runon_volume
        !%------------------------------------------------------------------   
        integer(c_int) function api_call_runoff_execute() &
            BIND(C, name='api_call_runoff_execute')
            use, intrinsic :: iso_c_binding
            implicit none 
        end function api_call_runoff_execute
        !%------------------------------------------------------------------
        integer(c_int) function api_get_subcatch_runoff(sc_idx,runoff) &
            BIND(C, name="api_get_subcatch_runoff")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in)    :: sc_idx
            real(c_double),        intent(inout) :: runoff
        end function api_get_subcatch_runoff
        !%------------------------------------------------------------------
        integer(c_int) function api_get_subcatch_runoff_nodeIdx(sc_idx, node_idx) &
            BIND(C, name="api_get_subcatch_runoff_nodeIdx")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in)    :: sc_idx
            integer(c_int),        intent(inout) :: node_idx
        end function api_get_subcatch_runoff_nodeIdx
        !%------------------------------------------------------------------
        integer(c_int) function api_getNumRdiiFlows(thisDateTime, nRDII) &
            BIND(C, name="api_getNumRdiiFlows")
            use, intrinsic :: iso_c_binding
            implicit none 
            real(c_double), value, intent(in)    :: thisDateTime
            integer(c_int),        intent(inout) :: nRDII
        end function api_getNumRdiiFlows
        !%------------------------------------------------------------------
        integer(c_int) function api_getRdiiFlow(rdiiIdx, nodeIdx, flowrate) &
            BIND(C, name="api_getRdiiFlow")
            use, intrinsic :: iso_c_binding
            implicit none 
            integer(c_int), value, intent(in)    :: rdiiIdx
            integer(c_int),        intent(inout) :: nodeIdx
            real(c_double),        intent(inout) :: flowrate
        end function api_getRdiiFlow
        !%------------------------------------------------------------------
        integer (c_int) function api_get_groundwaterFlow( &
            thisTime, LastRunoffTime, NextRunoffTime, sIdx, nodeIdx, flowrate) &
            BIND(C, name="api_get_groundwaterFlow")
            use, intrinsic :: iso_c_binding
            implicit none 
            real(c_double), value, intent(in)    :: thisTime
            real(c_double), value, intent(in)    :: LastRunoffTime
            real(c_double), value, intent(in)    :: NextRunoffTime
            integer(c_int), value, intent(in)    :: sIdx
            integer(c_int),        intent(inout) :: nodeIdx
            real(c_double),        intent(inout) :: flowrate
        end function api_get_groundwaterFlow
        !%------------------------------------------------------------------
        integer (c_int) function api_get_LID_DrainFlow( &
            thisTime, LastRunoffTime, NextRunoffTime, sIdx, nIdx, flowrate) &
            BIND(C, name="api_get_LID_DrainFlow")
            use, intrinsic :: iso_c_binding
            implicit none 
            real(c_double), value, intent(in)    :: thisTime
            real(c_double), value, intent(in)    :: LastRunoffTime
            real(c_double), value, intent(in)    :: NextRunoffTime
            integer(c_int), value, intent(in)    :: sIdx
            integer(c_int),        intent(inout) :: nIdx
            real(c_double),        intent(inout) :: flowrate
        end function api_get_LID_DrainFlow
        !%------------------------------------------------------------------
        integer(c_int) function api_call_climate_setState(thisDate) &
            BIND(C, name='api_call_climate_setState')
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), value, intent(in) :: thisDate
        end function api_call_climate_setState
        !%------------------------------------------------------------------
        integer(c_int) function api_get_evaporation_rate(evapRate) &
            BIND(C, name='api_get_evaporation_rate')
            use, intrinsic :: iso_c_binding
            implicit none 
            real(c_double),  intent(inout) :: evapRate
        end function api_get_evaporation_rate
        !%------------------------------------------------------------------
        ! --- Utils
        !%------------------------------------------------------------------
        integer(c_int) function api_find_object(object_type, object_name) &
            BIND(C, name="api_find_object")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value,  intent(in) :: object_type
            character(kind=c_char), intent(in) :: object_name
        end function api_find_object

    end interface
!%
!%==========================================================================
!% Procedures
!%==========================================================================
!%
    procedure(api_teststuff),                  pointer :: ptr_api_teststuff
    procedure(api_controls_count),             pointer :: ptr_api_controls_count
    procedure(api_controls_get_premise_data),  pointer :: ptr_api_controls_get_premise_data
    procedure(api_controls_get_action_data),   pointer :: ptr_api_controls_get_action_data
    procedure(api_controls_transfer_monitor_data), pointer :: ptr_api_controls_transfer_monitor_data
    procedure(api_controls_execute),           pointer :: ptr_api_controls_execute
    procedure(api_initialize),                 pointer :: ptr_api_initialize
    procedure(api_finalize),                   pointer :: ptr_api_finalize
    procedure(api_run_step),                   pointer :: ptr_api_run_step
    procedure(api_get_node_results),           pointer :: ptr_api_get_node_results
    procedure(api_get_link_results),           pointer :: ptr_api_get_link_results
    procedure(api_get_start_datetime),         pointer :: ptr_api_get_start_datetime
    procedure(api_get_end_datetime),           pointer :: ptr_api_get_end_datetime
    procedure(api_get_flowBC),                 pointer :: ptr_api_get_flowBC
    procedure(api_get_headBC),                 pointer :: ptr_api_get_headBC
    procedure(api_count_subobjects),           pointer :: ptr_api_count_subobjects
    procedure(api_get_SWMM_setup),             pointer :: ptr_api_get_SWMM_setup
    procedure(api_get_SWMM_times),             pointer :: ptr_api_get_SWMM_times
    procedure(api_get_NewRunoffTime),          pointer :: ptr_api_get_NewRunoffTime
    procedure(api_get_nodef_attribute),        pointer :: ptr_api_get_nodef_attribute
    procedure(api_get_linkf_attribute),        pointer :: ptr_api_get_linkf_attribute
    procedure(api_get_adjustments),            pointer :: ptr_api_get_adjustments
    procedure(api_get_transectf_attribute),    pointer :: ptr_api_get_transectf_attribute
    procedure(api_get_N_TRANSECT_TBL),         pointer :: ptr_api_get_N_TRANSECT_TBL
    procedure(api_get_transect_table),         pointer :: ptr_api_get_transect_table
    procedure(api_get_num_objects),            pointer :: ptr_api_get_num_objects
    procedure(api_get_object_name_len),        pointer :: ptr_api_get_object_name_len
    procedure(api_get_object_name),            pointer :: ptr_api_get_object_name
    procedure(api_get_num_table_entries),      pointer :: ptr_api_get_num_table_entries
    procedure(api_get_table_attribute),        pointer :: ptr_api_get_table_attribute
    procedure(api_get_first_entry_table),      pointer :: ptr_api_get_first_entry_table
    procedure(api_get_next_entry_table),       pointer :: ptr_api_get_next_entry_table
    procedure(api_get_next_entry_tseries),     pointer :: ptr_api_get_next_entry_tseries
    procedure(api_reset_timeseries_to_start),  pointer :: ptr_api_reset_timeseries_to_start
    procedure(api_write_output_line),          pointer :: ptr_api_write_output_line
    procedure(api_update_nodeResult),          pointer :: ptr_api_update_nodeResult
    procedure(api_update_linkResult),          pointer :: ptr_api_update_linkResult
    procedure(api_export_linknode_properties), pointer :: ptr_api_export_linknode_properties
    procedure(api_export_link_results),        pointer :: ptr_api_export_link_results
    procedure(api_export_node_results),        pointer :: ptr_api_export_node_results
    procedure(api_find_object),                pointer :: ptr_api_find_object
    procedure(api_export_runon_volume),        pointer :: ptr_api_export_runon_volume
    procedure(api_call_runoff_execute),        pointer :: ptr_api_call_runoff_execute
    procedure(api_get_subcatch_runoff),        pointer :: ptr_api_get_subcatch_runoff 
    procedure(api_get_subcatch_runoff_nodeIdx),pointer :: ptr_api_get_subcatch_runoff_nodeIdx 
    procedure(api_getNumRdiiFlows),            pointer :: ptr_api_getNumRdiiFlows
    procedure(api_getRdiiFlow),                pointer :: ptr_api_getRdiiFlow
    procedure(api_get_groundwaterFlow),        pointer :: ptr_api_get_groundwaterFlow
    procedure(api_get_LID_DrainFlow),          pointer :: ptr_api_get_LID_DrainFlow
    procedure(api_call_climate_setState),      pointer :: ptr_api_call_climate_setState
    procedure(api_get_evaporation_rate),       pointer :: ptr_api_get_evaporation_rate
    
    !% Error handling
    character(len = 1024) :: errmsg
    integer :: errstat

contains

!%=============================================================================
!% PUBLIC
!%=============================================================================
!%
    subroutine interface_teststuff()
        !%---------------------------------------------------------------------
        !% Description:
        !% stub routine used for testing ideas for api development
        !%---------------------------------------------------------------------
            character(64) :: subroutine_name = 'interface_teststuff'
        !%---------------------------------------------------------------------

        call load_api_procedure("api_teststuff")
        errstat = ptr_api_teststuff()
        call print_api_error(errstat, subroutine_name)

    end subroutine interface_teststuff
!%
!%=============================================================================
!% Controls and monitoring
!%=============================================================================
!%   
    subroutine interface_controls_count(nRules, nPremise, nThenAction, nElseAction)
        !%---------------------------------------------------------------------
        !% Description
        !% gets the data for controls and monitoring
        !%---------------------------------------------------------------------
            integer, intent (inout) :: nRules, nPremise, nThenAction, nElseAction
            character(64) :: subroutine_name = "interface_controls_count"
        !%---------------------------------------------------------------------
        
        !% --- count the number of various control/monitoring information
        call load_api_procedure("api_controls_count")
        errstat = ptr_api_controls_count(nRules, nPremise, nThenAction, nElseAction)
        call print_api_error(errstat, subroutine_name)

    end subroutine interface_controls_count
!%
!%=============================================================================
!%=============================================================================
!%    
    subroutine interface_controls_get_premise_data (     &
            locationL,        locationR,                 &
            linknodesimTypeL, linknodesimTypeR,          &
            attributeL,       attributeR,                &
            thisPremiseLevel, rIdx, success)
        !%---------------------------------------------------------------------
        !% Description
        !% Gets the monitoring locations for control premises 
        !%---------------------------------------------------------------------
            integer, intent(inout) :: locationL, locationR
            integer, intent(inout) :: linknodesimTypeL, linknodesimTypeR
            integer, intent(inout) :: attributeL, attributeR
            integer, intent(inout) :: thisPremiseLevel, success
            integer, intent(in)    :: rIdx
            character(65) :: subroutine_name = "interface_controls_get_premise_data"
        !%---------------------------------------------------------------------
        !%---------------------------------------------------------------------

        call load_api_procedure("api_controls_get_premise_data")

        success = ptr_api_controls_get_premise_data(    &
                    locationL,        locationR,        &
                    linknodesimTypeL, linknodesimTypeR, &
                    attributeL,       attributeR,       & 
                    thisPremiseLevel, rIdx)

        !% --- output data of -1 is not valid for location or attribute, 
        !%     so return nullvalue. Note that -1 for linknodesimType indicates
        !%     a simulation variable (e.g., time) rather than a monitor location
        if (locationL  == -1 ) locationL  = nullValueI  
        if (locationR  == -1 ) locationR  = nullValueI  
        if (attributeL == -1 ) attributeL = nullValueI 
        if (attributeR == -1 ) attributeR = nullValueI    

        !% --- increment the location by 1 since EPA-SWMM starts at 0 with indexes
        if (locationL .ne. nullvalueI) locationL = locationL + 1
        if (locationR .ne. nullvalueI) locationR = locationR + 1

    end subroutine interface_controls_get_premise_data
!%    
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_controls_get_action_data (     &
        location,                                            &
        attribute,                                           &
        thisActionLevel, rIdx, success, isThen)
        !%---------------------------------------------------------------------
        !% Description
        !% Gets the action locations for controls 
        !% isThen = 1 for a "then" action, 0 for an "else" action
        !%---------------------------------------------------------------------
            integer, intent(inout) :: location
            integer, intent(inout) :: attribute
            integer, intent(inout) :: thisActionLevel, success
            integer, intent(in)    :: rIdx, isThen
            character(65) :: subroutine_name = "interface_controls_get_action_data"
        !%---------------------------------------------------------------------
        !%---------------------------------------------------------------------

        call load_api_procedure("api_controls_get_action_data")

        success = ptr_api_controls_get_action_data( &
                    location,                       &
                    attribute,                      & 
                    thisActionLevel, rIdx, isThen)

        !% output data of -1 is not valid, so return nullvalue
        if (location  == -1 ) location  = nullValueI  
        if (attribute == -1 ) attribute = nullValueI              

        !% --- increment the output location by 1 since EPA-SWMM starts at 0 with indexes
        if (location .ne. nullvalueI) location = location + 1

    end subroutine interface_controls_get_action_data    
!%    
!%=============================================================================
!%=============================================================================
!%    
    subroutine interface_controls_transfer_monitor_data &
        (Depth, Head, Volume, Inflow, Flow, StatusSetting, TimeLastSet, &
         LinkNodeNum, linknodesimType)
        !%---------------------------------------------------------------------
        !% Description:
        !% transfers the monitoring data from SWMM5+ into EPA-SWMM so that it
        !% can be called for control actions using EPA-SWMM
        !%---------------------------------------------------------------------
        !% Declarations
            integer :: success
            integer, intent(in) :: LinkNodeNum, linknodesimType
            real(8), intent(in) :: Depth, Head, Volume, Inflow, Flow
            real(8), intent(in) :: StatusSetting, TimeLastSet
            real(8) :: TimeLastSetEpoch
            character(65) :: subroutine_name = 'interface_controls_transfer_monitor_data'
        !%---------------------------------------------------------------------
        !%---------------------------------------------------------------------
        call load_api_procedure("api_controls_transfer_monitor_data")

        !% ---convert timelast set to epoch days
        TimeLastSetEpoch = util_datetime_secs_to_epoch(TimeLastSet)

        !% --- send data into EPA-SWMM 
        !%     Note the link/node index is decremented by 1 to account for 
        !%     SWMM indexes starting a 0
        success = ptr_api_controls_transfer_monitor_data(             &
                    Depth, Head, Volume, Inflow, Flow, StatusSetting,  &
                    TimeLastSetEpoch, LinkNodeNum-1, linknodesimType)

    end subroutine interface_controls_transfer_monitor_data
!%    
!%=============================================================================
!%=============================================================================
!%    
    subroutine interface_controls_execute ()
        !%---------------------------------------------------------------------
        !% Description calls the api procedure that executes the controls_evaluate
        !% in EPA-SWMM
        !%---------------------------------------------------------------------
        !% Declarations
            integer :: number_of_actions
            real(8) :: currentTimeEpoch, ElapsedDays, dtDays
        !%---------------------------------------------------------------------
        !% -- convert elapsed seconds to SWMM time
        currentTimeEpoch = util_datetime_secs_to_epoch(setting%Time%Now)

        !% --- compute elapsed days since start of simulation
        ElapsedDays = setting%Time%Now / seconds_per_day

        !% --- compute time step in days
        dtDays = setting%Time%Hydraulics%Dt / seconds_per_day

        !% --- load the procedure
        call load_api_procedure("api_controls_execute")

        !% --- execute controls
        number_of_actions = ptr_api_controls_execute (currentTimeEpoch, ElapsedDays, dtDays )

    end subroutine interface_controls_execute
!%    
!%=============================================================================
!%=============================================================================
!%   
    subroutine interface_controls_get_action_results (targetsetting, timelastset, LinkIdx)
        !%---------------------------------------------------------------------
        !% Description
        !% gets the updated target setting and time last set associated with the 
        !% evaluation of EPA-SWMM controls
        !%---------------------------------------------------------------------
        !% Declarations
        integer :: error
        integer, intent(in) :: LinkIdx
        real(8), intent(inout) :: targetsetting, timelastset

        !% --- load the procedure
        call load_api_procedure("api_get_linkf_attribute")

        !% --- get the target setting
        error = ptr_api_get_linkf_attribute(LinkIdx-1, api_linkf_targetsetting, targetsetting)

        !% --- store in link
        link%R(LinkIdx,lr_TargetSetting) = targetsetting

        !% --- get the time last set
        error = ptr_api_get_linkf_attribute(LinkIdx-1, api_linkf_timelastset, timelastset)

        !% --- convert time last set to seconds and store in link
        timelastset = util_datetime_epoch_to_secs(timelastset)
        link%R(LinkIdx,lr_TimeLastSet) = timelastset

    end subroutine interface_controls_get_action_results   
!%    
!%=============================================================================
!%  Simulation subroutines/functions
!%=============================================================================
!%
    subroutine interface_init()
        !%---------------------------------------------------------------------
        !% Description:
        !%    initializes the EPA-SWMM shared library, creating input (.inp),
        !%    report (.rpt), and output (.out) files, necessary to run simulation with
        !%    EPA-SWMM. It also updates the number of objects in the SWMM model, i.e.,
        !%    number of links, nodes, and tables, and defines the start and end
        !%    simulation times.
        !%----------------------------------------------------------------------
            integer :: ppos, num_args, error
            character(64) :: subroutine_name = 'interface_init'
        
        !% Preliminaries:
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        setting%File%inp_file = trim(setting%File%inp_file) // c_null_char
        setting%File%rpt_file = trim(setting%File%rpt_file) // c_null_char
        setting%File%out_file = trim(setting%File%out_file) // c_null_char
        c_lib%filename = trim(setting%File%library_folder) // "/libswmm5.so"

        !% --- initialize the api between SWMM-C and SWMM5+
        !%     This returns and stores the SWMM-C input and output filenames
        call load_api_procedure("api_initialize")
        error = ptr_api_initialize( &
            setting%File%inp_file, &
            setting%File%rpt_file, &
            setting%File%out_file, &
            dummyI)
        call print_api_error(error, subroutine_name)
        api_is_initialized = .true.

        if (this_image() == 1) print *, ' ' !% needed because there is no \n after the EPA-SWMM printout of "Retrieving project data"

        !% --- get the time start, end, and interval data from SWMM-C input file
        call interface_get_SWMM_times()

        !% --- get the setup from the SWMM-C input file
        call interface_get_SWMM_setup()

        !% --- get counts of subobjects
        call interface_count_subobjects ()

        !% --- Set globals 
        N_link = setting%SWMMinput%N_link
        N_node = setting%SWMMinput%N_node
       
        if ((N_link == 200) .AND. (N_node == 200)) then
            print *, '********************************************************************'
            print *, '*                   USER CONFIGURATION ERROR                       *'
            print *, '* The EPA SWMM code has detected a parse error for the *.inp file. *'
            print *, '* Something appears to be misalligned or missing. This might be    *'
            print *, '* connections between nodes and links that are missing or not      *'
            print *, '* allowed, or something as simple as a missing [ around a keyword, *'
            print *, '* e.g., JUNCTIONS] instead of [JUNCTIONS]. This also occurs when   *'
            print *, '* a keyword is mispelled or the wrong words are used, e.g., if     *'
            print *, '* TRUE is used instead of YES for ALLOW_PONDING. Check the *.rpt   *'
            print *, '* file created by EPA-SWMM, which should be in the user output     *'
            print *, '* folder                                                           *'
            print *, '*                                                                  *'
            print *, '* Suggest you use the SWMM GUI to adjust and edit the *.inp file.  *'
            print *, '*                                                                  *'
            print *, '* Note that this error can be erroneously returned if you have a   *'
            print *, '* system with exactly 200 nodes and exactly 200 links.             *'
            print *, '********************************************************************'
            print *, ''
            call util_crashpoint( 309786)
            !% HACK -- developer's note:
            !% Unfortunately, the parse error returns a code of 200 in the get_num_objects()
            !% function, which is appears as setting%SWMMinput%N_link=200 and setting%SWMMinput%N_node=200.
            !% As it is relatively unlikely that a system will have exactly 200 of each, we are 
            !% simply calling the error condition when this happens.  We need to fix the
            !% API so that the error condition is correctly represented.
        end if


        !%----------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%interface) then
                print *, new_line("")
                print *, "setting%SWMMinput%N_link", setting%SWMMinput%N_link
                print *, "setting%SWMMinput%N_node", setting%SWMMinput%N_node
                print *, new_line("")
                print *, "SWMM start time", setting%Time%StartEpoch
                print *, "SWMM end time", setting%Time%EndEpoch
                print *, "setting%time%start", setting%Time%Start
                print *, "setting%time%end", setting%Time%End
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            end if

    end subroutine interface_init
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_finalize()
        !%---------------------------------------------------------------------
        !% Description:
        !% stops EPA SWMMC, deletes pointers and closese the shared library
        !%---------------------------------------------------------------------
        !% Declarations:
            integer :: errstat
            character(64) :: subroutine_name = 'interface_finalize'
        !%---------------------------------------------------------------------
        !% Preliminaries:
           if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        call load_api_procedure("api_finalize")
        errstat = ptr_api_finalize()
        call print_api_error(errstat, subroutine_name)

        !%---------------------------------------------------------------------
        !% Closing        
           if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_finalize
!%
!%=============================================================================
!%   Property-extraction functions (only run after initialization)
!%=============================================================================
!%
    subroutine interface_update_linknode_names()
        !%------------------------------------------------------------------
        !% Description:
        !%    Updates the link%Names and node%Names arrays which should've been
        !%    allocated by the time the subroutine is executed. Names are copied
        !%    from EPA-SWMM
        !%------------------------------------------------------------------
            integer :: ii
            character(64) :: subroutine_name = "interface_update_linknode_names"
        !%------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        do ii = 1, setting%SWMMinput%N_link
            call load_api_procedure("api_get_object_name")
            errstat = ptr_api_get_object_name(ii-1, link%Names(ii)%str, API_LINK)

            if (errstat /= 0) then
                write(*, "(A,i2,A)") "API ERROR : ", errstat, " [" // subroutine_name // "]"
                !stop 
                call util_crashpoint( 8673489)
                !return
            end if
        end do

        do ii = 1, setting%SWMMinput%N_node
            call load_api_procedure("api_get_object_name")
            errstat = ptr_api_get_object_name(ii-1, node%Names(ii)%str, API_NODE)

            if (errstat /= 0) then
                write(*, "(A,i2,A)") "API ERROR : ", errstat, " [" // subroutine_name // "]"
                !stop 
                call util_crashpoint( 48705)
                !return
            end if
        end do

        if (setting%Debug%File%interface) then
            print *, new_line("")
            print *, "List of Links"
            do ii = 1, setting%SWMMinput%N_link
                print *, "- ", link%Names(ii)%str
            end do
            print *, new_line("")
            print *, "List of Nodes"
            do ii = 1, setting%SWMMinput%N_node
                print *, "- ", node%Names(ii)%str
            end do
            print *, new_line("")
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end if

    end subroutine interface_update_linknode_names
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_update_transectID_names ()
        !%----------------------------------------------------------------------
        !% Description:
        !%    Updates the transectID arrays which should've been
        !%    allocated by the time the subroutine is executed. Names are copied
        !%    from EPA-SWMM
        !%------------------------------------------------------------------
        integer :: ii
        character(64) :: subroutine_name = "interface_update_transectID_names"
        !%------------------------------------------------------------------

        do ii = 1, setting%SWMMinput%N_transect

            !print *, ii, 'API_TRANSECT ',API_TRANSECT, trim(transectID(ii))
            call load_api_procedure("api_get_object_name")
            !errstat = ptr_api_get_object_name(ii-1, link%transectID(ii), API_TRANSECT)
            errstat = ptr_api_get_object_name(ii-1, link%transectID(ii)%str, API_TRANSECT)

            if (errstat /= 0) then
                write(*, "(A,i2,A)") "API ERROR : ", errstat, " [" // subroutine_name // "]"
                !stop 
                call util_crashpoint(498273)
                !return
            end if

        end do

    end subroutine interface_update_transectID_names    
!%
!%=============================================================================
!%=============================================================================
!%
    integer(c_int) function interface_get_obj_name_len(obj_idx, obj_type) result(len_value)
        !%------------------------------------------------------------------
        !% Description:
        !%    Returns the length of the name string associated to the EPA-SWMM object.
        !%    This function is necessary to allocate the entries of the link%Names and
        !%    node%Names arraysr. The function is currently compatible with NODE and
        !%    LINK types.
        !%------------------------------------------------------------------
            integer, intent(in) :: obj_idx  ! index of the EPA-SWMM object
            integer, intent(in) :: obj_type ! type of EPA-SWMM object (API_NODE, API_LINK)
            integer                :: error
            character(64)          :: subroutine_name = "interface_get_obj_name_len"
        !%------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_object_name_len")
        error = ptr_api_get_object_name_len(obj_idx-1, obj_type, len_value)
        call print_api_error(error, subroutine_name)
        
        if (setting%Debug%File%interface) then
            print *, obj_idx, obj_type, len_value
        end if

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function interface_get_obj_name_len
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_nodef_attribute(node_idx, attr) result(node_value)
        !%------------------------------------------------------------------
        !% Description:
        !%    Retrieves node attributes from EPA-SWMM. API node attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%------------------------------------------------------------------
            integer :: node_idx, attr, error
            real(c_double), target :: node_value
            character(64) :: subroutine_name = 'interface_get_nodef_attribute'
        !%------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(3(A,i5),A)") '*** enter ' // trim(subroutine_name) // &
            "(node_idx=", node_idx, ", attr=", attr, ")" // " [Processor ", this_image(), "]"

        if ((attr .ge. api_nodef_end) .or. (attr .le. api_nodef_start)) then
            print *, "error: unexpected node attribute value", attr
            print *, trim(reverseKey_api(attr))
            call util_crashpoint( 948705)
        end if

        if ((node_idx > N_node) .or. (node_idx < 1)) then
            print *, "error: unexpected node index value", node_idx
            print *, trim(reverseKey_api(attr))
            call util_crashpoint( 397904)
        end if
        
        !% --- Subtract 1 from every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_nodef_attribute")
        !% --- get the node value
        error = ptr_api_get_nodef_attribute(node_idx-1, attr, node_value)
        call print_api_error(error, subroutine_name)

        !% Adds 1 to every C index extracted from EPA-SWMM (it becomes a Fortran index)
        if (    (attr == api_nodef_extInflow_tSeries    )    &
           .or. (attr == api_nodef_extInflow_basePat_idx)    &
           .or. (attr == api_nodef_head_tSeries)             &
           .or. (attr == api_nodef_RouteTo)                  &
           ) then
            if (node_value /= -1) node_value = node_value + 1
        end if

        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface) &
                write(*,"(3(A,i5),A)") '*** leave ' // trim(subroutine_name) // &
                "(node_idx=", node_idx, ", attr=", attr, ")" // " [Processor ", this_image(), "]"
            
    end function interface_get_nodef_attribute
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_linkf_attribute(link_idx, attr, isInt) result(link_value)
        !%---------------------------------------------------------------------
        !% Description:
        !%    Retrieves link attributes from EPA-SWMM. API link attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%---------------------------------------------------------------------
        !% Decclarations
            integer :: link_idx, attr, error, ilink_value

            logical :: isInt !% true output integer is expected, false for real

            real(c_double), target :: link_value
            character(64) :: thisposition
            character(64) :: subroutine_name = 'interface_get_linkf_attribute'
        !%----------------------------------------------------------------------
        !% Preliminaries

            if (setting%Debug%File%interface) &
                write(*,"(3(A,i5),A)") '*** enter ' // trim(subroutine_name) // &
                "(link_idx=", link_idx, ", attr=", attr, ")" // " [Processor ", this_image(), "]"
        
            if ((link_idx > setting%SWMMinput%N_link) .or. (link_idx < 1)) then
                print *, "error: unexpected link index value", link_idx
                print *, trim(reverseKey_api(attr))
                call util_crashpoint(9987355)
            end if
        !%----------------------------------------------------------------------
        !% --- parse the link section
        if     (  attr .le. api_linkf_start) then    
            !% --- link attr number too small
            print *, "error: unexpected link attribute value", attr
            print *, trim(reverseKey_api(attr)) 
            call util_crashpoint( 498705)

        elseif     ( (attr  >   api_linkf_start) .and. (  attr <  api_linkf_commonbreak)) then
            !% --- for link attributes 1 to < api_linkf_commonBreak
            !%      we simply read in the link value for the attribute and exit this routine
            call load_api_procedure("api_get_linkf_attribute")

            ! ---  NOTE: use link index-1 because Fortran index starts in 1, whereas in C starts in 0
            error = ptr_api_get_linkf_attribute(link_idx-1, attr, link_value)
            thisposition = trim(subroutine_name)//'_A01'
            call print_api_error(error, thisposition)

        elseif (  attr == api_linkf_commonBreak) then
            !% --- this should never be called -- usually indicates mismatch between api.h
            !%     and define_api_keys
            print *, "error: unexpected link attribute value", attr
            print *, trim(reverseKey_api(attr)) 
            call util_crashpoint(27287422)

        elseif ( (attr >  api_linkf_commonbreak) .and. (attr < api_linkf_typeBreak) ) then
            !% --- for input link attributes in the range for for special elements   
            !%     the link value (integer) tells us other stuff to read in
            !%
            !% --- The "attr" input is one of the special element attributes, e.g., type, sub_type
            !%     First we use the api_linkf_type to get the overarching link type for this
            !%     link. This allows us to properly categorize the sub_type values to read         
            
            call load_api_procedure("api_get_linkf_attribute")
            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_type, link_value)
            thisposition = trim(subroutine_name)//'_B02'
            call print_api_error(error, thisposition)

            ilink_value = int(link_value) !% the linkf_type is always an integer

            !% --- handle the different linkf_type
            select case (ilink_value)

                case (API_CONDUIT)
                    select case (attr)
                        case (api_linkf_type)
                            link_value = lPipe
                        case (api_linkf_sub_type)
                            !% conduits do not have a sub type
                            link_value = undefinedKey
                        case default
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                            write(*,*)
                            write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at 45782987'
                            write(*,*) '   attr = ',attr
                            write(*,*) '   allowable are ',api_linkf_type, api_linkf_sub_type
                            write(*,*) '   skipping error condition!'
                            write(*,*) '******'
                    end select

                case (API_PUMP)
                    select case (attr)
                        case (api_linkf_type)
                            link_value = lPump
                        case (api_linkf_sub_type)
                            !% --- load subtype for pump
                            !print *, 'call CCC'
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_sub_type, link_value)
                            thisposition = trim(subroutine_name)//'_C04'

                            !% --- assign SWMM5+ key for pump subtype
                            select case (int(link_value))
                                case (API_TYPE1_PUMP)
                                    link_value = lType1Pump
                                case (API_TYPE2_PUMP)
                                    link_value = lType2Pump
                                case (API_TYPE3_PUMP)
                                    link_value = lType3Pump
                                case (API_TYPE4_PUMP)
                                    link_value = lType4Pump
                                case (API_IDEAL_PUMP)
                                    link_value = lTypeIdealPump

                                case default
                                    link_value = nullvalueI
                                    write(*,*)
                                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  442789'
                                    write(*,*) '   link_value = ',link_value
                                    write(*,*) '   allowable are ',API_TYPE1_PUMP, API_TYPE2_PUMP, API_TYPE3_PUMP, API_TYPE4_PUMP, API_IDEAL_PUMP
                                    write(*,*) '   skipping error condition!'
                                    write(*,*) '******'   
                            end select

                        case default
                            link_value = nullvalueI
                            write(*,*)
                            write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  8836785'
                            write(*,*) '   attr = ',attr
                            write(*,*) '   allowable are ',api_linkf_type, api_linkf_sub_type
                            write(*,*) '   skipping error condition!'
                            write(*,*) '******'      
                    end select

                case (API_ORIFICE)
                    select case (attr)
                        case (api_linkf_type)
                            link_value = lOrifice
                        case (api_linkf_sub_type)
                            !% --- load sub_type for orifice
                            !print *, 'call DDD'
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_sub_type, link_value)
                            thisposition = trim(subroutine_name)//'_D05'

                            !% --- assign SWMM5+ key for orifice subtype
                            select case (int(link_value))
                                case (API_SIDE_ORIFICE)
                                    link_value = lSideOrifice
                                case (API_BOTTOM_ORIFICE)
                                    link_value = lBottomOrifice 
                                case default
                                    link_value = nullvalueI
                                    write(*,*)
                                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  442789'
                                    write(*,*) '   link_value = ',link_value
                                    write(*,*) '   allowable are ',API_SIDE_ORIFICE, API_BOTTOM_ORIFICE
                                    write(*,*) '   skipping error condition!'
                                    write(*,*) '******'  
                            end select

                        case default
                            link_value = nullvalueI
                            write(*,*)
                            write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  9937685'
                            write(*,*) '   attr = ',attr
                            write(*,*) '   allowable are ',api_linkf_type, api_linkf_sub_type
                            write(*,*) '   skipping error condition!'
                            write(*,*) '******'   
                    end select

                case (API_WEIR)
                    select case (attr)
                        case (api_linkf_type)
                            link_value = lWeir
                        case (api_linkf_sub_type)
                            !% --- load subtype for weir
                            !print *, 'call EEE'
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_sub_type, link_value)
                            thisposition = trim(subroutine_name)//'_E06'
                            call print_api_error(error, thisposition)

                            !% --- assign SWMM5+ key for weir subtype
                            select case (int(link_value))
                                case (API_TRANSVERSE_WEIR)
                                    link_value = lTransverseWeir
                                case (API_SIDEFLOW_WEIR)
                                    link_value = lSideFlowWeir
                                case (API_VNOTCH_WEIR)
                                    link_value = lVnotchWeir
                                case (API_TRAPEZOIDAL_WEIR)
                                    link_value = lTrapezoidalWeir
                                case (API_ROADWAY_WEIR)
                                    link_value = lRoadWayWeir  

                                case default
                                    link_value = nullvalueI
                                    write(*,*)
                                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at 620873'
                                    write(*,*) '   link_value = ',link_value
                                    write(*,*) '   allowable are ',API_TRANSVERSE_WEIR, API_SIDEFLOW_WEIR, API_VNOTCH_WEIR, API_TRAPEZOIDAL_WEIR, API_ROADWAY_WEIR
                                    write(*,*) '   skipping error condition!'
                                    write(*,*) '******' 
                            end select

                        case default
                            link_value = nullvalueI
                            write(*,*)
                            write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  9937685'
                            write(*,*) '   attr = ',attr
                            write(*,*) '   allowable are ',api_linkf_type, api_linkf_sub_type
                            write(*,*) '   skipping error condition!'
                            write(*,*) '******'   
                    end select

                case (API_OUTLET)
                    select case (attr)
                        case (api_linkf_type)
                            link_value = lOutlet
                        case (api_linkf_sub_type)
                            !% --- load subtype for outlet
                            !print *, 'call FFF'
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_sub_type, link_value)
                            thisposition = trim(subroutine_name)//'_F06'
                            call print_api_error(error, thisposition)

                            !% --- assign SWMM5+ key for outlet type
                            select case (int(link_value))
                                case (API_NODE_DEPTH)
                                    link_value = lNodeDepth
                                case (API_NODE_HEAD)
                                    link_value = lNodeHead 

                                case default
                                    link_value = nullvalueI
                                    write(*,*)
                                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at 220455'
                                    write(*,*) '   link_value = ',link_value
                                    write(*,*) '   allowable are ',API_NODE_DEPTH, API_NODE_HEAD
                                    write(*,*) '   skipping error condition!'
                                    write(*,*) '******'   
                            end select

                        case default
                            link_value = nullvalueI
                            write(*,*)
                            write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  11947'
                            write(*,*) '   attr = ',attr
                            write(*,*) '   allowable are ',api_linkf_type, api_linkf_sub_type
                            write(*,*) '   skipping error condition!'
                            write(*,*) '******' 
                    end select

                case default
                    print *, 'in ',trim(subroutine_name)
                    print *, 'CODE ERROR this case default that should not be reached'
                    print *, 'Link index of ',link_idx
                    print *, 'input attr of ',attr, reverseKey_api(attr)
                    print *, 'unknown link value',ilink_value
                    print *, 'problem in call to ptr_api_get_linkf_attribute with api_linkf_type'
                    print *, 'A possible cause is the EPA-SWMM library needs to be recompiled'
                    print *, 'This usually requires a cmake .. followed by a make'
                    print *, 'It is also possible that an old version of the libswmm5.so file'
                    print *, 'is being used instead of the correct version. Check the library'
                    print *, 'folder being used in setting%File%library_folder.'
                    print *, 'Library folder is: ',trim(setting%File%library_folder)
                    print *, 'The other likely source is failure to update transect_ID'
                    print *, 'and linkf_ID in api.h, which causes enum of api keys to be offset'
                    print *, 'from define_api_keys.f90'
                    call util_crashpoint(6698733)
            end select
      
        elseif (  attr == api_linkf_typeBreak  ) then
            !% --- this should never be called
            print *, "CODE ERROR unexpected link attribute value", attr
            print *, trim(reverseKey_api(attr)) 
            call util_crashpoint(98273)        

        elseif ( (attr > api_linkf_typeBreak)    .and. (attr < api_linkf_end) ) then

            !% --- load the cross-section type no matter what the input attr is.
            call load_api_procedure("api_get_linkf_attribute")
            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_type, link_value)
            thisposition = trim(subroutine_name)//'_E05'
            call print_api_error(error, thisposition)

            ilink_value = int(link_value) !% these attributes should be integers
            select case (ilink_value)

                case (API_CIRCULAR)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lCircular
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_Q16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_R17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        ! case (api_linkf_xsect_culvertCode)
                        !     !print *, 'call III'
                        !     call load_api_procedure("api_get_linkf_attribute")
                        !     error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_culvertCode, link_value)
                        !     thisposition = trim(subroutine_name)//'_U17'
                        !     call print_api_error(error, thisposition)
                        case default
                            !% circular geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_FILLED_CIRCULAR)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lFilled_circular
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_Q16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_R16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yBot) 
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yBot, link_value)
                            thisposition = trim(subroutine_name)//'_T19'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- filled circular geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_RECT_CLOSED)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lRectangular_closed
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_F06'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_G07'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_H06'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_I06'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- rectangular geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_RECT_OPEN)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lRectangular
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_H08'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_I09'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_J09'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_K09'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- rectangular geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_TRAPEZOIDAL)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lTrapezoidal
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yBot, link_value)
                            thisposition = trim(subroutine_name)//'_J10'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_K11'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_L11'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_M11'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- trapezoidal geometry does not have certain geometric features (i.e. top-width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_TRIANGULAR)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lTriangular
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_M12'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_N13'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_O13'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_P13'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- triangular geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_PARABOLIC)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lParabolic
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_O14'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_P15'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_Q15'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_R15'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- parabolic geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_POWERFUNC)
                    print *, 'CODE ERROR API_POWERFUNC geometry not handled yet'
                    call util_crashpoint(398472)
                    select case (attr)
                        case (api_linkf_geometry)
                        case (api_linkf_xsect_wMax)
                        case (api_linkf_xsect_yFull)
                        case default
                    end select

                case (API_RECT_TRIANG)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lRect_triang
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_S18'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_T19'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_U19'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_V19'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yBot)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yBot, link_value)
                            thisposition = trim(subroutine_name)//'_W19'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- rectangular triangular geometry does not have certain geometric features (i.e. bottom width) 
                            !%     thus, set that link%R column to nullvalueR
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_RECT_ROUND)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lRect_round
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_S18'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_T19'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rBot)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rBot, link_value)
                            thisposition = trim(subroutine_name)//'_T19'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yBot) 
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yBot, link_value)
                            thisposition = trim(subroutine_name)//'_T19'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- Mod_basket geometry does not have certain geometric features (i.e. bottom width) 
                            !%     thus, set that link%R column to nullvalueR
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_MOD_BASKET)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lMod_basket
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_S18'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_T19'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rBot)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rBot, link_value)
                            thisposition = trim(subroutine_name)//'_T19'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- Mod_basket geometry does not have certain geometric features (i.e. bottom width) 
                            !%     thus, set that link%R column to nullvalueR
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_HORIZ_ELLIPSE)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lHoriz_ellipse
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- basket handle geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_VERT_ELLIPSE)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lVert_ellipse
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- basket handle geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_ARCH)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lArch
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- basket handle geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_EGGSHAPED)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lEggshaped
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- basket handle geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_HORSESHOE)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lHorseshoe
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- basket handle geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_GOTHIC)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lGothic
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- gothic geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_CATENARY)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lCatenary
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- catenary geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_SEMIELLIPTICAL)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lSemi_elliptical
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- catenary geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_BASKETHANDLE)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lBasket_handle
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- basket handle geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_SEMICIRCULAR)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lSemi_circular
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_U16'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_V17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                            thisposition = trim(subroutine_name)//'_S17'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                            thisposition = trim(subroutine_name)//'_T17'
                            call print_api_error(error, thisposition)
                        case default
                            !% --- basket handle geometry does not have certain geometric features (i.e. bottom width) 
                            if (isInt) then
                                link_value = nullvalueI
                            else
                                link_value = nullvalueR
                            end if
                    end select

                case (API_IRREGULAR)
                    select case (attr)
                    case (api_linkf_geometry)
                        link_value = lIrregular
                    case (api_linkf_transectidx)
                        !% --- load the transect ID
                        call load_api_procedure("api_get_linkf_attribute")
                        error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_transectidx, link_value)
                        thisposition = trim(subroutine_name)//'_Q16'
                        call print_api_error(error, thisposition)
                    case default
                        !% --- irregular geometry does not have certain geometric features (i.e. bottom width) 
                        if (isInt) then
                            link_value = nullvalueI
                        else
                            link_value = nullvalueR
                        end if
                    end select

                case (API_CUSTOM)
                    print *, 'CODE ERROR API_CUSTOM geometry not handled yet'
                    call util_crashpoint(298733)
                    select case (attr)
                        case (api_linkf_geometry)
                        case (api_linkf_xsect_wMax)
                        case (api_linkf_xsect_yFull)
                        case default
                    end select

                case (API_FORCE_MAIN)
                    select case (attr)
                        case (api_linkf_geometry)
                            link_value = lForce_main
                        case (api_linkf_xsect_wMax)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                            thisposition = trim(subroutine_name)//'_S26'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yFull)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                            thisposition = trim(subroutine_name)//'_T27'
                            call print_api_error(error, thisposition)
                        case (api_linkf_xsect_aFull)
                             call load_api_procedure("api_get_linkf_attribute")
                             error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_aFull, link_value)
                             thisposition = trim(subroutine_name)//'_T28'
                             call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rFull)
                             call load_api_procedure("api_get_linkf_attribute")
                             error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rFull, link_value)
                             thisposition = trim(subroutine_name)//'_T29'
                             call print_api_error(error, thisposition)
                        case (api_linkf_xsect_yBot)
                            call load_api_procedure("api_get_linkf_attribute")
                             error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yBot, link_value)
                             thisposition = trim(subroutine_name)//'_T30'
                             call print_api_error(error, thisposition)
                        case (api_linkf_xsect_rBot)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_rBot, link_value)
                            thisposition = trim(subroutine_name)//'_T31'
                            call print_api_error(error, thisposition)     
                        case (api_linkf_forcemain_coef)
                            call load_api_procedure("api_get_linkf_attribute")
                            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_forcemain_coef, link_value)
                            thisposition = trim(subroutine_name)//'_U28'  
                            call print_api_error(error, thisposition)

                        case default
                            print *, 'case default for API_FORCE_MAIN ',trim(reverseKey_api(attr)),' (could be error): ',attr,trim(reverseKey_api(attr))
                    end select

                case default
                        !% --- some links like pumps or outlets does not have any geometric features
                        !%     thus, link%R geometry columns (i.e. fulldepth, width) will be set to nullvalueR
                        link_value = nullvalueR
            end select
          
        else
            !% t--- his should never be reached
            print *, "CODE ERROR unexpected link attribute value", attr
            print *, trim(reverseKey_api(attr)) 
            call util_crashpoint(878293)    
        
        end if

        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  then
                write(*,"(3(A,i5),A)") '*** leave ' // trim(subroutine_name) // &
                "(link_idx=", link_idx, ", attr=", attr, ")" // " [Processor ", this_image(), "]"
            end if

    end function interface_get_linkf_attribute
!%
!%=============================================================================    
!%=============================================================================
!%
    subroutine interface_get_adjustments ()
        !%---------------------------------------------------------------------
        !% Description
        !% Loads data in month arrays from SWMM input [ADJUSTMENTS]
        !%---------------------------------------------------------------------
        !% Declarations
            integer :: error
        !%---------------------------------------------------------------------

        call load_api_procedure("api_get_adjustments")

        !% --- set default values
        setting%Adjust%Temperature  = zeroR
        setting%Adjust%Evaporation  = zeroR
        setting%Adjust%Rainfall     = oneR
        setting%Adjust%Conductivity = oneR
        
         error = ptr_api_get_adjustments(         &
                12,                               &
                setting%Adjust%Temperature(:),    &
                setting%Adjust%Evaporation(:),    &
                setting%Adjust%Rainfall(:),       &
                setting%Adjust%Conductivity(:))

    end subroutine interface_get_adjustments
!%
!%=============================================================================
!%=============================================================================
!% 
    function interface_get_transectf_attribute (transect_idx, attr) result(transect_value)
        !%---------------------------------------------------------------------
        !% Obtains the scalar attribute values for transects from SWMM
        !%---------------------------------------------------------------------
        !% Declarations
            integer :: transect_idx, attr, error

            real(c_double), target :: transect_value
            character(64) :: thisposition
            character(64) :: subroutine_name = 'interface_get_transectf_attribute'
        !%------------------------------------------------------------------

        if ((transect_idx > setting%SWMMinput%N_transect) .or. (transect_idx < 1)) then
            print *, "error: unexpected tranect index value", transect_idx
            print *, trim(reverseKey_api(attr))
            call util_crashpoint(992255)
        end if

        if     (  attr .le. api_transectf_start) then
            print *, "error: unexpected transectf attribute value", attr
            print *, trim(reverseKey_api(attr)) 
            call util_crashpoint(2333235)

        elseif (  attr == api_transectf_ID) then
            !% --- skip this -- handled elsewhere
            return

        elseif ( (attr   >  api_transectf_ID) .and. (attr < api_transectf_end) ) then

            call load_api_procedure("api_get_transectf_attribute")
            !% --- transect index-1 because Fortran index starts in 1, whereas in C starts in 0
            error = ptr_api_get_transectf_attribute(transect_idx-1, attr, transect_value)
            thisposition = trim(subroutine_name)//'_A01'
            call print_api_error(error, thisposition)

            !% --- return with transect_value for output
            return

        elseif (  attr .ge. api_transectf_end ) then
            print *, "error: unexpected transectf attribute value", attr
            print *, trim(reverseKey_api(attr)) 
            call util_crashpoint(883782)
            !return

        else
            print *, "error: unexpected transectf attribute value", attr
            print *, trim(reverseKey_api(attr)) 
            call util_crashpoint(273672)
            !return
        end if

    end function interface_get_transectf_attribute
!%
!%=============================================================================
!%=============================================================================
!%
    integer function interface_get_N_TRANSECT_TBL () result(n_transect_tbl)   
        !%---------------------------------------------------------------------
        !% Description -- retrieves length of the transect table (i.e. number)
        !% of depth slices) from EPA-SWMM
        !%---------------------------------------------------------------------
        call load_api_procedure("api_get_N_TRANSECT_TBL")
        n_transect_tbl = ptr_api_get_N_TRANSECT_TBL ()

    end function interface_get_N_TRANSECT_TBL 
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_get_transect_table ()
        !%---------------------------------------------------------------------
        !% Description -- retrieves transect width, area, hyd radius from
        !% EPA-SWMM
        !%---------------------------------------------------------------------
        !% Declarations
            integer :: error, ii, jj
        !%---------------------------------------------------------------------

        call load_api_procedure("api_get_transect_table")
        
        do ii=1,setting%SWMMinput%N_transect

            error = ptr_api_get_transect_table(&
                ii-1, setting%SWMMInput%N_transect_depth_items,    &
                link%transectTableDepthR(ii,:,tt_area),  &
                link%transectTableDepthR(ii,:,tt_width),  &
                link%transectTableDepthR(ii,:,tt_hydradius))

        end do

        !%---------------------------------------------------------------------
    end subroutine interface_get_transect_table    
!%
!%=============================================================================    
!%=============================================================================
!%
    function interface_get_table_attribute(table_idx, attr)
        !%------------------------------------------------------------------
        !% Description:
        !%    Retrieves table attributes from EPA-SWMM. API table attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%------------------------------------------------------------------
            integer :: table_idx, attr, error
            real(8) :: interface_get_table_attribute

            real(c_double), target :: table_value
            character(64) :: thisposition
            character(64) :: subroutine_name = 'interface_get_table_attribute'
        !%------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if ((attr .ge. api_table_end) .or. (attr < 1) ) then
            print *, "error: unexpected table attribute value", attr
            call util_crashpoint( 28704)
        end if

        if ((table_idx > setting%SWMMinput%N_curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            print *, trim(reverseKey_api(attr))
            call util_crashpoint( 498734)
        end if

        !% --- Subtracts 1 from every Fortran index (it becomes a C index)
        if (attr == api_table_type) then
            call load_api_procedure("api_get_table_attribute")
            error = ptr_api_get_table_attribute(table_idx-1, attr, table_value)
            if (table_value == API_STORAGE_CURVE) then
                interface_get_table_attribute = StorageCurve
            else if (table_value == API_DIVERSION_CURVE) then
                interface_get_table_attribute = DiversionCurve
            else if (table_value == API_TIDAL_CURVE) then
                interface_get_table_attribute = TidalCurve
            else if (table_value == API_RATING_CURVE) then
                interface_get_table_attribute = RatingCurve
            else if (table_value == API_CONTROL_CURVE) then
                interface_get_table_attribute = ControlCurve
            else if (table_value == API_SHAPE_CURVE) then
                interface_get_table_attribute = ShapeCurve
            else if (table_value == API_WEIR_CURVE) then
                interface_get_table_attribute = WeirCurve
            else if (table_value == API_PUMP1_CURVE) then
                interface_get_table_attribute = Pump1Curve
            else if (table_value == API_PUMP2_CURVE) then
                interface_get_table_attribute = Pump2Curve
            else if (table_value == API_PUMP3_CURVE) then
                interface_get_table_attribute = Pump3Curve
            else if (table_value == API_PUMP4_CURVE) then
                interface_get_table_attribute = Pump4Curve
            else
                write(*,*) '****** Not clear why this is null brh20211207 at 873984'
                interface_get_table_attribute = nullvalueI
            end if
        else
            call load_api_procedure("api_get_table_attribute")
            error = ptr_api_get_table_attribute(table_idx-1, attr, table_value)
            call print_api_error(error, subroutine_name)
            interface_get_table_attribute = table_value
        end if

        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            print *, "table", table_value, attr
        end if
    end function interface_get_table_attribute
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_num_table_entries(table_idx)
        !%------------------------------------------------------------------
        !% Description:
        !%    Retrieves table attributes from EPA-SWMM. API table attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%------------------------------------------------------------------
            integer :: table_idx, table_type, error
            integer :: interface_get_num_table_entries
            integer(c_int), target :: table_entries
            character(64) :: thisposition
            character(64) :: subroutine_name = 'interface_get_num_table_entries'
        !%------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if ((table_idx > setting%SWMMinput%N_curve) .or. (table_idx < 1)) then
            print *, "CODE ERROR unexpected table index value", table_idx
            call util_crashpoint( 835551)
        end if

        !% --- Subtracts 1 from every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_num_table_entries")
        error = ptr_api_get_num_table_entries(table_idx-1, API_CURVE, table_entries)
        call print_api_error(error, subroutine_name)
        interface_get_num_table_entries = table_entries


        !%---------------------------------------------------------------------
        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            print *, "table", table_entries
        end if
        
    end function interface_get_num_table_entries
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_first_entry_table(table_idx)
        !%------------------------------------------------------------------
        !% Description:
        !%    Retrieves first table entries from EPA-SWMM. API table attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%------------------------------------------------------------------
            integer :: table_idx, error, success
            real(8) :: interface_get_first_entry_table(2)

            real(c_double), target :: x_entry, y_entry
            character(64) :: subroutine_name = 'interface_get_first_entry_table'
        !%------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if ((table_idx > setting%SWMMinput%N_curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            call util_crashpoint( 837014)
        end if

        !% --- Subtracts 1 from every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_first_entry_table")
        success = ptr_api_get_first_entry_table(table_idx-1, API_CURVE, x_entry, y_entry)

        if (success == 0) then
            error = -1
        else
            error = 0
        end if

        call print_api_error(error, subroutine_name)
        interface_get_first_entry_table(1) = x_entry
        interface_get_first_entry_table(2) = y_entry

        !%---------------------------------------------------------------------
        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end if
    end function interface_get_first_entry_table
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_next_entry_table(table_idx, table_type)
        !%------------------------------------------------------------------
        !% Description:
        !%    Retrieves next table entries from EPA-SWMM. API table attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%------------------------------------------------------------------
            integer :: table_idx, table_type, error, success
            real(8) :: interface_get_next_entry_table(2)

            real(c_double), target :: x_entry, y_entry
            character(64) :: subroutine_name
        !%------------------------------------------------------------------
            subroutine_name = 'interface_get_next_entry_table'
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        if ((table_idx > setting%SWMMinput%N_curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            !stop 
            call util_crashpoint( 8367894)
            !return
        end if

        !% Subtracts 1 from every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_next_entry_table")
        success = ptr_api_get_next_entry_table(table_idx-1, table_type, x_entry, y_entry)

        if (success == 0) then
            error = -1
        else
            error = 0
        end if

        call print_api_error(error, subroutine_name)
        interface_get_next_entry_table(1) = x_entry
        interface_get_next_entry_table(2) = y_entry

        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  then
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            end if

    end function interface_get_next_entry_table
!%
!%=============================================================================
!%   Boundary Conditions (execute after initialization only)
!%=============================================================================
!%
    function interface_get_BC_resolution(node_idx) result(resolution)
        !%------------------------------------------------------------------
        !% Description:
        !%    Computes the finest pattern resolution associated with a node's BC.
        !%    The resulting pattern type is stored in node%I(:,ni_pattern_resolution)
        !%    and is reused to fetch inflow BCs such that the amount of inflow points
        !%    is minimized.
        !% Notes:
        !%    * Currently patterns are only associated with inflow BCs. It is necessary
        !%      to preserve the order of api_monthly, api_weekend, api_daily, and
        !%      api_hourly in define_api_index.f08 for the function to work.
        !%    * The function is called during the intialization of the node%I table
        !%------------------------------------------------------------------
            integer, intent(in) :: node_idx
            integer             :: p0, p1, p2, p3, p4
            integer             :: resolution
            real(8)             :: baseline
            character(64)       :: subroutine_name
        !%------------------------------------------------------------------
            subroutine_name = 'interface_get_BC_resolution'
            if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name)  //  &
                " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------
        resolution = nullvalueI

        if (node%YN(node_idx, nYN_has_inflow)) then ! Upstream/Lateral BC

            resolution = api_nopattern

            if (node%YN(node_idx, nYN_has_extInflow)) then

                p0 = interface_get_nodef_attribute(node_idx, api_nodef_extInflow_basePat_type)

                select case (p0)
                case (api_hourly_pattern)
                    print *, 'CODE NEEDS TESTING: extInflow hourly patterns not tested'
                    call util_crashpoint(2240871)
                    resolution = api_hourly
                case (api_weekend_pattern) 
                    print *, 'CODE NEEDS TESTING: extInflow weekend patterns not tested'
                    call util_crashpoint(2240872)
                    resolution = api_weekend
                case (api_daily_pattern)
                    print *, 'CODE NEEDS TESTING: extInflow daily patterns not tested'
                    call util_crashpoint(2240873)
                    resolution = api_daily
                case (api_monthly_pattern)
                    print *, 'CODE NEEDS TESTING: extInflow monthly patterns not tested'
                    call util_crashpoint(2240874)
                    baseline = interface_get_nodef_attribute(node_idx, api_nodef_extInflow_baseline)
                    if (baseline > 0) resolution = api_monthly
                case (api_nopattern)
                    !% --- no pattern invoked
                case default
                    write(*,*) 'CODE ERROR unexpected case in ',trim(subroutine_name)
                    write(*,*), '  pattern resolution is ',p0
                    write(*,*), '  allowed are ',api_hourly, api_weekend, api_daily, api_monthly, api_nopattern
                    write(*,*), ' see define_api_keys.f90 for api_monthly, etc.'
                    call util_crashpoint(552873)
                end select
            else
                !% ---no external inflows to this node
            end if

            if (node%YN(node_idx, nYN_has_dwfInflow)) then

                ! write(*,*) 'call api_nodef_dwfInflow_hourly_pattern'
                p1 = interface_get_nodef_attribute(node_idx, api_nodef_dwfInflow_hourly_pattern)
                ! write(*,*) '   p1 = ',p1

                ! write(*,*) 'api_nodef_dwfInflow_weekend_pattern'
                p2 = interface_get_nodef_attribute(node_idx, api_nodef_dwfInflow_weekend_pattern)
                ! write(*,*) '   p2 = ',p2

                ! write(*,*) 'api_nodef_dwfInflow_daily_pattern'
                p3 = interface_get_nodef_attribute(node_idx, api_nodef_dwfInflow_daily_pattern)
                ! write(*,*) '   p3 = ',p3

                ! write(*,*) 'api_nodef_dwfInflow_monthly_pattern'
                p4 = interface_get_nodef_attribute(node_idx, api_nodef_dwfInflow_monthly_pattern)
                ! write(*,*) '   p4 = ',p4

                if (p1 > 0) then
                    ! print *, 'CODE NEEDS TESTING: dwfInflow hourly patterns not tested'
                    ! call util_crashpoint(18820871)
                    resolution = max(api_hourly, resolution)
                else if (p2 > 0) then
                    ! print *, 'CODE NEEDS TESTING: dwfInflow weekend patterns not tested'
                    ! call util_crashpoint(18820872)
                    resolution = max(api_weekend, resolution)
                else if (p3 > 0) then
                    ! print *, 'CODE NEEDS TESTING: dwfInflow daily patterns not tested'
                    ! call util_crashpoint(18820873)
                    resolution = max(api_daily, resolution)
                else if (p4 > 0) then
                    ! print *, 'CODE NEEDS TESTING: dwfInflow monthly patterns not tested'
                    ! call util_crashpoint(18820874)
                    resolution = max(api_monthly, resolution)
                else
                    write(*,*) '***** unexpected else in ',trim(subroutine_name), 'at 3987044'
                    write(*,*), '   dwf resolution  read are ',p1, p2, p3, p4
                    write(*,*), '   allowed are ',api_hourly_pattern, api_weekend_pattern, api_daily_pattern, api_monthly_pattern
                    write(*,*), '   skipping error condition!'
                    write(*,*) '******'  
                end if
            else
                !% --- no dwfInflow to this node
            end if
        else
            !% ---  no dwf or ext inflows to this node
        end if

    end function interface_get_BC_resolution
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_reset_timeseries_to_start(bc_idx) result(tstart)    
        !%---------------------------------------------------------------------
        !% Description:
        !% resets the current entry point for the time series associated with
        !% the bc_idx
        !%---------------------------------------------------------------------
            integer, intent(in) :: bc_idx
            real(8)             :: tstart
            integer             :: tseries_idx, nidx, error
            real(8)             :: tdata(2)
        !%---------------------------------------------------------------------

        !% --- get the node index
        nidx = BC%flowI(bc_idx, bi_node_idx)
        !% --- get the time series index
        tseries_idx = interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries)
        print *, 'in interface_reset_timeseries_to_start, tseries_idx', tseries_idx
        !% --- reset the entry point
        call load_api_procedure("api_reset_timeseries_to_start")
        error = ptr_api_reset_timeseries_to_start(tseries_idx-1)
        print *, 'after api_reset_timeseries...'

    end function interface_reset_timeseries_to_start
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_next_inflow_time(bc_idx, tnow, timemaxEpoch) result(tnext)
        !%---------------------------------------------------------------------
        !% Description:
        !% Gets the next inflow time. If the next time is less than the maximum
        !% time (timemax) then the Tseries.x1 and .y1 stored values will be changed
        !% to the new value.
        !% NOTE: timemax is the "Epoch" time used in EPA-SWMM, but the
        !% output from this is local time with time=0 as the start of the simulation
        !%
        !%---------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: bc_idx
            real(8), intent(in) :: tnow, timemaxEpoch
            real(8)             :: tnext, t1, t2, tnextp
            integer             :: tseries_idx, success
            integer             :: year, month, day, hours, minutes, seconds
            integer, pointer    :: nidx, nres
            character(64) :: subroutine_name = 'interface_get_next_inflow_time'
        !%---------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------
        !% Aliases:
            !% --- node index for this BC
            nidx => BC%flowI(bc_idx, bi_node_idx)
            !% --- get node pattern resolution
            nres => node%I(nidx, ni_pattern_resolution)
        !%---------------------------------------------------------------------
        
            !% --- consistency checking
        if (.not. node%YN(nidx, nYN_has_inflow)) then
            print *, "Error, node " // node%Names(nidx)%str // " does not have an inflow"
        end if

        !% --- get the next time for pattern resolution
        !%     Note that nres=0 returns nullvalueR
        tnextp = util_datetime_get_next_time(tnow, nres) 

        if (node%YN(nidx, nYN_has_extInflow) .or. node%YN(nidx, nYN_has_dwfInflow)) then
            !% --- for external inflows (file), get the timeseries index
            tseries_idx = interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries)

            if (tseries_idx >= 0) then
                !% --- this gets the Tseries.x2 values
                !% --- gets time in days at what is now the x2 pointer 
                tnext = interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries_x2)

                tnext = util_datetime_epoch_to_secs(tnext)

                !% --- move the linked list to the next step (which sets up the next value for x2)
                success = get_next_entry_tseries(tseries_idx, timemaxEpoch)

                if (success .ne. oneI) then
                    !% --- failure to read time later than tnow from file
                    print *, ' '
                    write(*,"(A)") 'INPUT FILE FAILURE'
                    write(*,"(A,i8,A,A,A)") 'at node ',nidx, ' ',trim(node%Names(nidx)%str)
                    write(*,"(A,f12.0,A)") 'Input file reader cannot find time past ',tnow /3600.d0, ' hours'
                    tnext = util_datetime_secs_to_epoch(tnow)
                    call util_datetime_decodedate(tnext, year, month, day)
                    call util_datetime_decodetime(tnext, hours, minutes, seconds)
                    write(*,"(A,i4,a,i2,a,i2,a,i2,a,i2)") 'or date ',year,'-',month,'-',day,' at ',hours,':',minutes
                    tnext = util_datetime_epoch_to_secs(timemaxEpoch)
                    write(*,"(A,f12.0,A)")  'Note that simulation end time is ',tnext/3600.d0,' hours'
                    call util_datetime_decodedate(timemaxEpoch, year, month, day)
                    call util_datetime_decodetime(timemaxEpoch, hours, minutes, seconds)
                    write(*,"(A,i4,a,i2,a,i2,a,i2,a,i2)") 'or date ',year,'-',month,'-',day,' at ',hours,':',minutes
                    write(*,"(A)") 'The input file must have a data up through the end of the simulation period.'
                    print *, ' '
                    
                    call util_crashpoint(6698734)
                end if
            else
                !% --- if no external file, use the end time
                !%     HACK -- what are we doing here?
                tnext = setting%Time%End
            end if
        else
            !% --- continue, no action if there isn't an external inflow    
        end if
        !% --- the next time is the smaller of the value in the the timeseries or
        !%     the time associated with the pattern.
        !%     HACK -- NEED TO CHECK PATTERN OPERATION
        tnext = min(tnext, tnextp)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function interface_get_next_inflow_time
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_next_head_time(bc_idx, tnow, timemaxEpoch) result(tnext)
        !%---------------------------------------------------------------------
        !% Description
        !% Gets the next  headtime. If the next time is less than the maximum
        !% time (timemax) then the Tseries.x1 and .y1 stored values will be changed
        !% to the new value.
        !% NOTE: timemax is the "Epoch" time used in EPA-SWMM, but the
        !% output from this is local time with time=0 as the start of the simulation
        !%---------------------------------------------------------------------
            integer, intent(in) :: bc_idx
            real(8), intent(in) :: tnow, timemaxEpoch
            real(8)             :: tnext, t1, t2, tnextp
            integer             :: tseries_idx, success
            integer             :: year, month, day, hours, minutes, seconds
            integer, pointer    :: nidx
            character(64) :: subroutine_name = 'interface_get_next_head_time'
        !%---------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------    
        !% Aliases
            nidx => BC%headI(bc_idx, bi_node_idx)
        !%---------------------------------------------------------------------

        select case (BC%headI(bc_idx, bi_subcategory))  
            case (BCH_fixed, BCH_normal, BCH_free) 
                !% --- these cases should not be here!
                print *, 'CODE ERROR unexpected boundary condition case'
                call util_crashpoint(609874)
                return
            case (BCH_tseries)
                !% --- get the timeseries index
                tseries_idx = interface_get_nodef_attribute(nidx, api_nodef_head_tSeries)
                if (tseries_idx >= 0) then
                    !% --- this gets the Tseries.x2 values
                    !%     Note the Tseries.x1 values will be overwritten by the .x2 values
                    !%     only if the x2 value is less than timemax. This prepares for the
                    !%     the next step of storing for SWMM5+      
                    !% --- gets time in days at what is now the x2 pointer 
                    tnext = interface_get_nodef_attribute(nidx, api_nodef_head_tSeries_x2)
                    tnext = util_datetime_epoch_to_secs(tnext)
                    success = get_next_entry_tseries(tseries_idx, timemaxEpoch)
                    if (success .ne. oneI) then
                        !% --- failure to read time later than tnow from file
                        print *, ' '
                        write(*,"(A)") 'INPUT FILE FAILURE'
                        write(*,"(A,f12.0,A)") 'Input file reader cannot find time past ',tnow /3600.d0, ' hours'
                        tnext = util_datetime_secs_to_epoch(tnow)
                        call util_datetime_decodedate(tnext, year, month, day)
                        call util_datetime_decodetime(tnext, hours, minutes, seconds)
                        write(*,"(A,i4,a,i2,a,i2,a,i2,a,i2)") 'or date ',year,'-',month,'-',day,' at ',hours,':',minutes
                        tnext = util_datetime_epoch_to_secs(timemaxEpoch)
                        write(*,"(A,f12.0,A)")  'Note that simulation end time is ',tnext/3600.d0,' hours'
                        call util_datetime_decodedate(timemaxEpoch, year, month, day)
                        call util_datetime_decodetime(timemaxEpoch, hours, minutes, seconds)
                        write(*,"(A,i4,a,i2,a,i2,a,i2,a,i2)") 'or date ',year,'-',month,'-',day,' at ',hours,':',minutes
                        write(*,"(A)") 'The input file must have a data up through the end of the simulation period.'
                        print *, ' '
                        
                        call util_crashpoint(609834)
                    end if
                else
                    !% --- if no external file, use the end time
                    !% HACK -- what are we doing here? Should this be an error condition?
                    tnext = setting%Time%End
                end if
            case (BCH_tidal)
                print *, 'CODE ERROR tidal outfall not yet handled'
                call util_crashpoint(446929)
            case default
                print *, BC%headI(bc_idx, bi_subcategory), trim(reverseKey(BC%headI(bc_idx, bi_subcategory)))
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(44822)
        end select

        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_get_next_head_time
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_flowBC(bc_idx, tnow) result(bc_value)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
            integer, intent(in) :: bc_idx
            real(8), intent(in) :: tnow
            integer             :: error
            integer, pointer    :: nidx
            real(8)             :: epochNow, bc_value
            character(64) :: subroutine_name  = 'interface_get_flowBC'
        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        nidx => BC%flowI(bc_idx, bi_node_idx)

        epochNow = util_datetime_secs_to_epoch(tnow)
        call load_api_procedure("api_get_flowBC")
        error = ptr_api_get_flowBC(nidx-1, epochNow, bc_value)
        call print_api_error(error, subroutine_name)

        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_get_flowBC
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_headBC(bc_idx, tnow) result(bc_value)
        !%---------------------------------------------------------------------
        !% Description:
        !% Gets the head boundary condition for time now
        !%---------------------------------------------------------------------
            integer, intent(in) :: bc_idx
            real(8), intent(in) :: tnow
            integer             :: error
            integer, pointer    :: nidx
            real(8)             :: epochNow, bc_value
            character(64) :: subroutine_name = 'interface_get_headBC'
        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        nidx => BC%headI(bc_idx, bi_node_idx)

        epochNow = util_datetime_secs_to_epoch(tnow)
        call load_api_procedure("api_get_headBC")
        error = ptr_api_get_headBC(nidx-1, epochNow, bc_value)
        call print_api_error(error, subroutine_name)

        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_get_headBC

!%=============================================================================
!%    Write Outputs (execute after initialization only)
!%=============================================================================
!%
    subroutine interface_export_link_results(link_idx)
        !%---------------------------------------------------------------------
        !% Description:
        !% calls the api_export_link_results() function in api.c
        !%---------------------------------------------------------------------
            integer, intent(in) :: link_idx
            integer :: error
            character(64) :: subroutine_name = 'interface_export_link_results'
        !%---------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        call load_api_procedure("api_export_link_results")
        error = ptr_api_export_link_results(link_idx-1)
        call print_api_error(error, subroutine_name)

        !%---------------------------------------------------------------------
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine interface_export_link_results
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_update_nodeResult(node_idx, result_type, node_result)
        !%---------------------------------------------------------------------
        !% Description:
        !% calls the update_nodeResult() function in api.c
        !% OBSOLETE: Hold for future consideration
        !%---------------------------------------------------------------------
            integer, intent(in) :: node_idx, result_type
            real(8), intent(in) :: node_result
            integer             :: error
            character(64)       :: subroutine_name = "interface_update_nodeResult"
        !%----------------------------------------------------------------------

        print *, 'OBSOLETE/WRONG CODE (needs unit conversion)'
        stop 209873

        ! if (setting%Debug%File%interface)  &
        !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        ! call load_api_procedure("api_update_nodeResult")
        ! error = ptr_api_update_nodeResult(node_idx-1, result_type, node_result)
        ! call print_api_error(error, subroutine_name)

        ! if (setting%Debug%File%interface)  &
        !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_update_nodeResult
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_update_linkResult(link_idx, result_type, link_result)
        !%---------------------------------------------------------------------
        !% Description:
        !% calls the update_linkResult() function in api.c
        !% OBSOLETE: Hold for future consideration
        !%---------------------------------------------------------------------
            integer, intent(in) :: link_idx, result_type
            real(8), intent(in) :: link_result
            integer             :: error
            character(64)       :: subroutine_name = "interface_update_linkResult"
        !%----------------------------------------------------------------------

        print *, 'OBSOLETE/WRONG CODE (needs unit conversion)'
        stop 20987322

        ! if (setting%Debug%File%interface)  &
        !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        ! call load_api_procedure("api_update_linkResult")
        ! error = ptr_api_update_linkResult(link_idx-1, result_type, link_result)
        ! call print_api_error(error, subroutine_name)

        ! if (setting%Debug%File%interface)  &
        !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_update_linkResult
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_write_output_line(reportTime)
        !%------------------------------------------------------------------
        !% Description:
        !%    Writes .out file with SWMM5+ data
        !%------------------------------------------------------------------
            real(c_double), intent(in) :: reportTime ! time in seconds
            integer                    :: error
            character(64)              :: subroutine_name = "interface_write_output_line"
        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call load_api_procedure("api_write_output_line")
        error = ptr_api_write_output_line(reportTime)
        call print_api_error(error, subroutine_name)

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine interface_write_output_line
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_count_subobjects()
        !%---------------------------------------------------------------------
        !% Description:
        !% gets counts of subobjects of SWMM objects
        !%---------------------------------------------------------------------
        !% Declarations:
            integer :: error
            integer, pointer :: Ngroundwater
            character(64) :: subroutine_name = "interface_count_subobjects"
        !%---------------------------------------------------------------------

        Ngroundwater => setting%SWMMinput%N_groundwater
        call load_api_procedure("api_count_subobjects")
        error = ptr_api_count_subobjects(Ngroundwater)
        call print_api_error(error, subroutine_name)   

    end subroutine interface_count_subobjects    
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_get_SWMM_setup()
        !%---------------------------------------------------------------------
        !% Description
        !% gets control variables that have been input in the SWMM-C *.inp
        !% file and have been processed by SWMM-C. Stores data that will be
        !% used in SWMM5+ in the setting.SWMMinput... structure
        !%---------------------------------------------------------------------
            integer :: flow_units, infiltration_model_type, flow_routing_model_type
            integer :: link_offset_type, force_main_equation, ignore_rainfall
            integer :: ignore_snowmelt, ignore_groundwater, ignore_rdii
            integer :: ignore_routing, ignore_quality, allow_ponding
            integer :: steadystate_skip
            real(8) :: steadystate_system_flow_tolerance
            real(8) :: steadystate_lateral_flow_tolerance
            real(8) :: routing_step_lengthening_time, routing_step_Courant_factor
            real(8) :: routing_step_minimum
            integer :: inertial_damping_type, normal_flow_limiter_type
            real(8) :: minimum_surface_area, minimum_conduit_slope
            integer :: maximum_number_of_trials
            real(8) :: head_convergence_tolerance
            integer :: number_parallel_threads, tempdir_provided, control_rule_step
            integer :: surcharge_method
            integer :: error, ii
            integer, parameter  :: nset = 30
            logical       :: thisWarning(1:nset)
            logical       :: thisFailure(1:nset)
            character(64) :: thisProblem(1:nset)
            character(30) :: thisVariable(1:nset)
            character(64) :: subroutine_name = 'interface_get_SWMM_setup'
        !%----------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%----------------------------------------------------------------------

        call load_api_procedure("api_get_SWMM_setup")

        error = ptr_api_get_SWMM_setup( &
            flow_units,  &
            infiltration_model_type, &
            flow_routing_model_type, &
            link_offset_type, &
            force_main_equation, &
            ignore_rainfall, &
            ignore_snowmelt, &
            ignore_groundwater, &
            ignore_rdii, &
            ignore_routing, &
            ignore_quality, &
            allow_ponding, &
            steadystate_skip, &
            steadystate_system_flow_tolerance, &
            steadystate_lateral_flow_tolerance, &
            routing_step_lengthening_time, &
            routing_step_Courant_factor, &
            routing_step_minimum,  &
            inertial_damping_type, &
            normal_flow_limiter_type, &
            minimum_surface_area, &
            minimum_conduit_slope, &
            maximum_number_of_trials, &
            head_convergence_tolerance, &
            number_parallel_threads, &
            tempdir_provided, &
            control_rule_step, &
            surcharge_method &
            )

        call print_api_error(error, subroutine_name)

        !% --- count numbers of objects from define_api_keys (SWMM objects)
        setting%SWMMinput%N_gage        = get_num_objects(API_GAGE)
        setting%SWMMinput%N_subcatch    = get_num_objects(API_SUBCATCH)
        setting%SWMMinput%N_node        = get_num_objects(API_NODE)
        setting%SWMMinput%N_link        = get_num_objects(API_LINK)
        setting%SWMMinput%N_pollutant   = get_num_objects(API_POLLUT)
        setting%SWMMinput%N_landuse     = get_num_objects(API_LANDUSE)
        setting%SWMMinput%N_timepattern = get_num_objects(API_TIMEPATTERN)
        setting%SWMMinput%N_curve       = get_num_objects(API_CURVE)
        setting%SWMMinput%N_tseries     = get_num_objects(API_TSERIES)
        setting%SWMMinput%N_control     = get_num_objects(API_CONTROL)
        setting%SWMMinput%N_transect    = get_num_objects(API_TRANSECT)
        setting%SWMMinput%N_aquifer     = get_num_objects(API_AQUIFER)
        setting%SWMMinput%N_unithyd     = get_num_objects(API_UNITHYD)
        setting%SWMMinput%N_snowmelt    = get_num_objects(API_SNOWMELT)
        setting%SWMMinput%N_shape       = get_num_objects(API_SHAPE)
        setting%SWMMinput%N_LID         = get_num_objects(API_LID)

        !% --- count number of node types
        setting%SWMMinput%N_junction = get_num_objects(API_JUNCTION)
        setting%SWMMinput%N_outfall  = get_num_objects(API_OUTFALL)
        setting%SWMMinput%N_storage  = get_num_objects(API_STORAGE)
        setting%SWMMinput%N_divider  = get_num_objects(API_DIVIDER)

        !% --- seconds between control rule evaluations
        setting%SWMMinput%ControlRuleStep = control_rule_step
    
        thisFailure(:) = .false.
        thisWarning(:) = .false.
        thisVariable(:) = ''
        thisProblem(:) = 'none'

        !% --- Flow Units are always CMS for SWMM5+ (enum in enum.h)
        ii = 1
        select case(flow_units)
            case(0)
                setting%SWMMinput%FlowUnitsType = SWMM_FlowUnits_CFS
                thisWarning(ii) = .true.
                thisProblem(ii) = 'CFS input flow units are converted to CMS for all SWMM5+ computation and output.'
            case(1)
                setting%SWMMinput%FlowUnitsType = SWMM_FlowUnits_GPM
                thisWarning(ii) = .true.
                thisProblem(ii) = 'GPM input flow units are converted to CMS for all SWMM5+ computation and output.'
            case(2)
                setting%SWMMinput%FlowUnitsType = SWMM_FlowUnits_MGD
                thisWarning(ii) = .true.
                thisProblem(ii) = 'MGD input flow units are converted to CMS for all SWMM5+ computation and output.'
            case(3)
                setting%SWMMinput%FlowUnitsType = SWMM_FlowUnits_CMS
                !% continue with CMS units
            case(4)
                setting%SWMMinput%FlowUnitsType = SWMM_FlowUnits_LPS
                thisWarning(ii) = .true.
                thisProblem(ii) = 'LPS input flow units are converted to CMS for all SWMM5+ computation and output.'
            case(5)
                setting%SWMMinput%FlowUnitsType = SWMM_FlowUnits_MLD
                thisWarning(ii) = .true.
                thisProblem(ii) = 'MLD input flow units are converted to CMS for all SWMM5+ computation and output.'
            case default
                thisWarning(ii) = .true.
                thisFailure(ii) = .true.
                thisVariable(ii) = 'FLOW_UNITS'
                thisProblem(ii) = 'Unknown value for SWMM input: (CFS, GPM, NGD, CMS, LPS, MLD) supported'
        end select

        !% --- Infiltration method 
        !%     from enum InfilType in infil.h
        ii=ii+1
        select case (infiltration_model_type)
            case (0)
                setting%SWMMinput%InfiltrationType = SWMM_Infiltration_Horton
            case (1)
                setting%SWMMinput%InfiltrationType = SWMM_Infiltration_Mod_Horton
            case (2)
                setting%SWMMinput%InfiltrationType = SWMM_Infiltration_Green_Ampt
            case (3)
                setting%SWMMinput%InfiltrationType = SWMM_Infiltration_Mod_Green_Ampt
            case (4)
                setting%SWMMinput%InfiltrationType = SWMM_Infiltration_Curve_Number
            case default
                thisWarning(ii) = .true.
                thisFailure(ii) = .true.
                thisProblem(ii) = 'Unknown value for SWMM input (HORTON, MODIFIED_HORTON, GREEN_AMPT, MODIFIED_GREEN_AMPT, CURVE_NUMBER) supported'
        end select

           
        !% --- Routing model 
        !%     is always DYNWAVE for SWMM5+
        !%     from enum RouteModelType in enum.h
        ii=ii+1
        select case (flow_routing_model_type)
            case (0)
                setting%SWMMinput%FlowRoutingType = SWMM_FlowRouting_NoRouting
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'FLOW_ROUTING'
                thisProblem(ii)  = 'NO_ROUTING must be set to DYNWAVE in *.inp file'
            case (1)
                setting%SWMMinput%FlowRoutingType = SWMM_FlowRouting_SteadyFlow
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'FLOW_ROUTING'
                thisProblem(ii)  = 'STEADY must be set to DYNWAVE in *.inp file'
            case (2)
                setting%SWMMinput%FlowRoutingType = SWMM_FlowRouting_KinematicWave
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'FLOW_ROUTING'
                thisProblem(ii)  = 'KINWAVE must be set to DYNWAVE in *.inp file'
            case (3)
                setting%SWMMinput%FlowRoutingType = SWMM_FlowRouting_ExtendedKinematicWave
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'FLOW_ROUTING'
                thisProblem(ii)  = 'XKINWAVE must be set to DYNWAVE in *.inp file'
            case (4)
                setting%SWMMinput%FlowRoutingType = SWMM_FlowRouting_DynamicWave
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'FLOW_ROUTING'
                thisProblem(ii)  = 'Unknown value for SWMM input: SWMM5+ only supports DYNWAVE in *.inp file'
        end select

        !% --- link offset type
        !%     from enum OffsetType in enums.h
        ii=ii+1
        select case (link_offset_type)
            case(0)
                setting%SWMMinput%LinkOffsetsType = SWMM_LinkOffset_DepthOffset
            case(1)
                setting%SWMMinput%LinkOffsetsType = SWMM_LinkOffset_ElevOffset
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'LINK_OFFSETS'
                thisProblem(ii)  = 'Unknown value found: (DEPTH, ELEVATION) supported'
        end select

        !% --- Force mains in EPA SWMM use 0 = Hazen-Williams, 1 = Darcy-Weisbach 
        !%     from enum ForceMainType in enum.h
        ii=ii+1
        select case (force_main_equation)
            case (0)
                setting%SWMMinput%ForceMainEquationType = HazenWilliams
            case (1)
                setting%SWMMinput%ForceMainEquationType = DarcyWeisbach
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'FORCE_MAIN_EQUATION'
                thisProblem(ii)  = 'Unknown value for SWMM input: (H_W, D_W) supported.'
        end select    

        !% --- ignore rainfall (logical)
        ii=ii+1
        select case (ignore_rainfall)
            case (0)
                setting%SWMMinput%IgnoreRainfall = .false.
            case (1)
                setting%SWMMinput%IgnoreRainfall = .true.  
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'IGNORE_RAINFALL'
                thisProblem(ii)  = 'Unknown value for SWMM input: (YES, NO) supported'
        end select

        !% --- ignore snowmelt (logical)
        ii=ii+1
        select case (ignore_snowmelt)
            case (0)
                setting%SWMMinput%IgnoreSnowmelt = .false.
            case (1)
                setting%SWMMinput%IgnoreSnowmelt = .true.  
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'IGNORE_SNOWMELT'
                thisProblem(ii)  = 'Unknown value for SWMM input: (YES, NO) supported'
        end select

        !% --- ignore groundwater (logical)
        ii=ii+1
        select case (ignore_groundwater)
            case (0)
                setting%SWMMinput%IgnoreGroundwater = .false.
            case (1)
                setting%SWMMinput%IgnoreGroundwater = .true.  
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'IGNORE_GROUNDWATER'
                thisProblem(ii)  = 'Unknown value for SWMM input: (YES, NO) supported'
        end select

        !% --- ignore RDII (logical)
        ii=ii+1
        select case (ignore_RDII)
            case (0)
                setting%SWMMinput%IgnoreRDII = .false.
            case (1)
                setting%SWMMinput%IgnoreRDII = .true.  
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'IGNORE_RDII'
                thisProblem(ii)  = 'Unknown value for SWMM input: (YES, NO) supported'
        end select

        !% --- ignore routing (logical)
        ii=ii+1
        select case (ignore_routing)
            case (0)
                setting%SWMMinput%IgnoreRouting = .false.
            case (1)
                setting%SWMMinput%IgnoreRouting = .true.  
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'IGNORE_ROUTING'
                thisProblem(ii)  = 'SWMM5+ requires IGNORE_ROUTING = NO in SWMM input file'
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'IGNORE_ROUTING'
                thisProblem(ii)  = 'Unknown value for SWMM input: (YES, NO) supported'
        end select

        !% --- ignore quality (logical)
        ii=ii+1
        select case (ignore_quality)
            case (0)
                setting%SWMMinput%IgnoreQuality = .false.
                if (setting%SWMMinput%N_pollutant > 0) then
                    thisWarning(ii)  = .true.
                    thisFailure(ii)  = .true.
                    thisVariable(ii) = 'IGNORE_QUALITY'
                    thisProblem(ii)  = 'SWMM5+ requires IGNORE_QUALTIY = YES in SWMM input file'
                else 
                    !% -- no pollutants found, so this shouldn't matter
                end if
            case (1)
                setting%SWMMinput%IgnoreQuality = .true.  
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'IGNORE_QUALITY'
                thisProblem(ii)  = 'Unknown value for SWMM input: (YES, NO) supported'
        end select

        !% --- Ponding 
        ii=ii+1
        select case (allow_ponding)
            case (0)
                setting%SWMMinput%AllowPonding = .false.
            case (1)
                setting%SWMMinput%AllowPonding = .true.
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'ALLOW_PONDING'
                thisProblem(ii)  = 'Unknown value for SWMM input: (YES, NO) supported'
        end select

        !% --- Steady state skip routing
        ii=ii+1
        select case (steadystate_skip)
            case (0)
                setting%SWMMinput%SteadyState_Skip = .false.
            case (1)
                setting%SWMMinput%SteadyState_Skip = .true.
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'SKIP_STEADY_STATE'
                thisProblem(ii)  = 'SWMM5+ requires SKIP_STEADY_STATE = NO in SWMM input file'
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'SKIP_STEADY_STATE'
                thisProblem(ii)  = 'Unknown value for SWMM input: (YES, NO) supported'
        end select

        !% --- steady state system flow tolerance
        !%     only relevant if steadystate_skip = 1, which is a failure mode
        setting%SWMMinput%SteadyState_System_FlowrateTolerance  = steadystate_system_flow_tolerance
        setting%SWMMinput%SteadyState_Lateral_FlowrateTolerance = steadystate_lateral_flow_tolerance

        !% --- routing step lengthening (not used in SWMM5+)
        ii=ii+1
        setting%SWMMinput%RoutingStep_LengtheningTime = routing_step_lengthening_time
        if (routing_step_lengthening_time .ne. zeroR) then
            thisWarning(ii)  = .true.
            thisFailure(ii)  = .false.
            thisVariable(ii) = 'LENGTHENING_STEP'
            thisProblem(ii)  = 'input value is ignored in SWMM5+.'
        end if

        !% --- Courant factor for EPA SWMM variable time step
        ii=ii+1
        setting%SWMMinput%RoutingStep_CourantFactor  = routing_step_Courant_factor
        if (routing_step_Courant_factor .ne. zeroR) then
            thisWarning(ii)  = .true.
            thisFailure(ii)  = .false.
            thisVariable(ii) = 'VARIABLE_STEP'
            thisProblem(ii)  = 'input value is ignored in SWMM5+.'
        end if

        !% --- minimum hydraulics time step of EPA SWMM is ignored
        ii=ii+1
        setting%SWMMinput%RoutingStep_Minimum         = routing_step_minimum
        if (routing_step_minimum .ne. zeroR) then
            thisWarning(ii)  = .true.
            thisFailure(ii)  = .false.
            thisVariable(ii) = 'MINIMUM_STEP'
            thisProblem(ii)  = 'input value is ignored in SWMM5+.'
        end if

        !% --- Inertial damping is never used in SWMM5+
        !%     from enum InertialDampingType in enums.h
        ii=ii+1
        select case (inertial_damping_type)
            case (0)
                setting%SWMMinput%InertialDampingType = SWMM_InertialDamping_NoDamping
            case (1)
                setting%SWMMinput%InertialDampingType = SWMM_InertialDamping_PartialDamping
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'INERTIAL_DAMPING'
                thisProblem(ii)  = 'Only NONE is supported in SWMM5+'
            case (2)
                setting%SWMMinput%InertialDampingType = SWMM_InertialDamping_FullDamping
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'INERTIAL_DAMPING'
                thisProblem(ii)  = 'Only NONE is supported in SWMM5+'
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'INERTIAL_DAMPING'
                thisProblem(ii)  = 'Unknown value for SWMM input: Only NONE is supported in SWMM5+'
        end select

        !% --- Normal flow limiter
        !%     from enum NormalFlowType in enums.h
        ii=ii+1
        thisWarning(ii)  = .true.
        thisFailure(ii)  = .false.
        thisVariable(ii) = 'NORMAL_FLOW_LIMITED'
        thisProblem(ii)  = 'input value is ignored in SWMM5+.'
        select case (normal_flow_limiter_type)
        case(0)
            setting%SWMMinput%NormalFlowLimiterType = SWMM_NormalFlowLimiterType_Slope
        case(1)
            setting%SWMMinput%NormalFlowLimiterType = SWMM_NormalFlowLimiterType_Froude
        case(2)
            setting%SWMMinput%NormalFlowLimiterType = SWMM_NormalFlowLimiterType_Both
        case default
        end select

        !% --- Minimum surface area for nodes from EPA SWMM is not used
        setting%SWMMinput%SurfaceArea_Minimum = minimum_surface_area

        !% --- Minimum slope is never used (default of 0.0 is ignored)
        ii=ii+1
        setting%SWMMinput%ConduitSlope_Minimum = minimum_conduit_slope
        thisWarning(ii)  = .true.
        thisFailure(ii)  = .false.
        thisVariable(ii) = 'MIN_SLOPE'
        thisProblem(ii)  = 'input value is ignored in SWMM5+.'

        !% --- Max trials is irrelevant
        ii=ii+1
        setting%SWMMinput%NumberOfTrials_Maximum = maximum_number_of_trials
        thisWarning(ii)  = .true.
        thisFailure(ii)  = .false.
        thisVariable(ii) = 'MAX_TRIALS'
        thisProblem(ii)  = 'input value is ignored in SWMM5+.'

        !% --- Head Tolerance is irrelevant 
        ii=ii+1
        setting%SWMMinput%Head_ConvergenceTolerance = head_convergence_tolerance
        thisWarning(ii)  = .true.
        thisFailure(ii)  = .false.
        thisVariable(ii) = 'HEAD_TOL'
        thisProblem(ii)  = 'input value is ignored in SWMM5+.'

        !% --- The number of parallel threads cannot be set at runtime in SWMM5+ as of 20211223
        ii=ii+1
        setting%SWMMinput%NumberParallelThreads = number_parallel_threads
        select case (number_parallel_threads)
            case (1)
                !% --- continue
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .false.
                thisVariable(ii) = 'NUM_THREADS'
                thisProblem(ii)  = 'input value is ignored in SWMM5+ (set threads at compile time).'
        end select

        ii=ii+1        
        if (tempdir_provided == 1) then
            setting%SWMMinput%TempDirectory_Provided = .true.
            thisWarning(ii)  = .true.
            thisFailure(ii)  = .false.
            thisVariable(ii) = 'TEMPDIR'
            thisProblem(ii)  = 'input value is ignored for SWMM5+ output (available for some EPA-SWMM output).'
        else
            setting%SWMMinput%TempDirectory_Provided = .false.
        end if

        !%---------------------------------------------------------------------
        !% --- OTHER SWMM OPTIONS

        !% --- Control Rule Step is always OK
        ii=ii+1
        thisWarning(ii) = .false.
        thisFailure(ii) = .false.
        thisVariable(ii) = 'RULESTEP'

        !% --- Surcharge method
        !%     only Preissman SLOT is presently allowed for surcharge method -- handled by JSON file
        !%     from enum  SurchargeMethodType in enums.h
        ii=ii+1
        select case (surcharge_method)
            case (0)
                setting%SWMMinput%SurchargeMethod = SWMM_SurchargeMethod_Extran
                thisWarning(ii) = .true.
                thisVariable(ii) = 'SURCHARGE_METHOD'
                thisProblem(ii)  = 'input value ignored. SWMM5+ uses Preissmann slot.'
            case (1)
                setting%SWMMinput%SurchargeMethod = SWMM_SurchargeMethod_Slot
                !% Preissmann SLOT is specified
            case default
                thisWarning(ii)  = .true.
                thisFailure(ii)  = .true.
                thisVariable(ii) = 'SURCHARGE_METHOD'
                thisProblem(ii)  = 'Unknown value for SWMM input. SWMM5+ overwrites (EXTRAN,SLOT).'
        end select

        !% Pollutant transport in hydraulics not supported in SWMM5+ as of 20221103
        ii=ii+1
        if ((setting%SWMMinput%N_pollutant > 0) .and. (.not. setting%SWMMinput%IgnoreQuality))then
            thisWarning(ii)  = .true.
            thisFailure(ii)  = .true.
            thisVariable(ii) = '[POLLUTANTS]'
            thisProblem(ii)  = 'are not supported in SWMM5+, set IGNORE_QUALITY = YES.'
        end if

        !% Dividers are not used in SWMM5+ as we do not support Kinematic Wave
        ii=ii+1
        if (setting%SWMMinput%N_divider > 0) then
            thisWarning(ii)  = .true.
            thisFailure(ii)  = .false.
            thisVariable(ii) = '[DIVIDER]'
            thisProblem(ii)  = 'in SWMM input file are converted to junctions for dynamic wave.'
        end if

        !%----------------------------------------------------------------------
        !% --- printouts
        !%
        if ((any(thisWarning)) .and. (this_image() == 1) ) then
            write(*,'(A)') ' '
            write(*,'(A)') ' '
            write(*,'(A)') ' *******************************************************************'
            write(*,'(A)') ' **                          WARNING'
            write(*,'(A)') ' ** The following from the SWMM *.inp file values or code defaults  '
            write(*,'(A)') ' ** are ignored or changed in SWMM5+ due to present code limitations.'
            do ii=1,nset
                if (thisWarning(ii)) then
                    write(*,"(A,A,A,A)") ' **    ',trim(thisVariable(ii)),'--',trim(thisProblem(ii))
                end if
            end do
            write(*,'(A)') '*******************************************************************'
            write(*,*) ' '
        end if

        if ((any(thisFailure)) .and. (this_image() == 1) ) then
            write(*,'(A)') ' '
            write(*,'(A)') ' '
            write(*,'(A)') ' *******************************************************************'
            write(*,'(A)') ' **                   FATAL INPUT FILE FAILURE'
            write(*,'(A)') ' ** The following from the SWMM *.inp file values or code defaults  '
            write(*,'(A)') ' ** cannot be used in SWMM5+ due to present code limitations.'
            do ii=1,nset
                if (thisFailure(ii)) then
                    write(*,"(A,A,A,A)") ' **    ',trim(thisVariable(ii)),'--',trim(thisProblem(ii))
                end if
            end do
            write(*,'(A)') '*******************************************************************'
            write(*,*) ' '
            call util_crashpoint(559874)
        end if


        !%----------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine interface_get_SWMM_setup    
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_get_SWMM_times()
        !%---------------------------------------------------------------------
        !% Description:
        !% gets the time variables that have been input in the SWMM-C *.inp
        !% file and have been processed by SWMM-C
        !%---------------------------------------------------------------------
            integer                :: error
            real(c_double), target :: reportStart_epoch, hydraulics_route_step 
            real(c_double), target :: starttime_epoch, endtime_epoch
            real(c_double), target :: TotalDuration
            integer(c_int), target :: reportStep, hydrology_wet_step, hydrology_dry_step
            integer(c_int), target :: sweep_start_dayofyear, sweep_end_dayofyear, dry_days
            character(64)          :: subroutine_name = 'interface_get_SWMM_times'
        !%----------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%----------------------------------------------------------------------

        !% --- reportStart is given in epoch datetime
        !% --- reportStep and hydroStep are given in integer seconds
        !%
        !% --- Note that if reportStart is earlier than StartEpoch from get_start_datetime()
        !% ... then SWMM will make reportStart = start datetime.
        !%
        !% --- This generates EPA-SWMM Error Code: -1 if the reportStart is after the end time.
        call load_api_procedure("api_get_SWMM_times")
        error = ptr_api_get_SWMM_times( &
            starttime_epoch,  &
            endtime_epoch,    &
            reportStart_epoch, &
            reportStep,  &
            hydrology_wet_step,     &
            hydrology_dry_step,     &
            sweep_start_dayofyear,  &
            sweep_end_dayofyear, &
            dry_days,                &
            hydraulics_route_step,   &
            TotalDuration)
        call print_api_error(error, subroutine_name)

        setting%SWMMinput%StartEpoch            = starttime_epoch
        setting%SWMMinput%EndEpoch              = endtime_epoch                      
        setting%SWMMinput%ReportStartTimeEpoch  = reportStart_epoch
        setting%SWMMinput%ReportTimeInterval    = reportStep
        setting%SWMMinput%Hydrology_WetStep     = hydrology_wet_step
        setting%SWMMinput%Hydrology_DryStep     = hydrology_dry_step
        setting%SWMMinput%Sweep_Start_Day       = sweep_start_dayofyear
        setting%SWMMinput%Sweep_End_Day         = sweep_end_dayofyear
        setting%SWMMinput%DryDaysBeforeStart    = dry_days
        setting%SWMMinput%Hydraulics_RouteStep  = hydraulics_route_step
        setting%SWMMinput%TotalDuration         = TotalDuration

        !%----------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine interface_get_SWMM_times
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_find_object(object_type, object_name) result(object_idx)
        !%---------------------------------------------------------------------
        !% Description:
        !% Returns the index of the object, or 0 if the object couldn't be found
        !%---------------------------------------------------------------------
            character(*), intent(in) :: object_name
            integer, intent(in) :: object_type
            integer :: object_idx
            character(64) :: subroutine_name = 'interface_find_object'
        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call load_api_procedure("api_find_object")
        object_idx = ptr_api_find_object(object_type, trim(object_name)//c_null_char) + 1

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_find_object
!%
!%=============================================================================
!%=============================================================================
!%   
    subroutine interface_export_runon_volume()
        !%---------------------------------------------------------------------
        !% Description:
        !% takes the runon flowrate from hydraulics to subcatchments, stored
        !% in subcatchR(:,sr_RunOn_Volume), and transfers to the outfall[].vRouted
        !% in EPASWMM
        !%---------------------------------------------------------------------
        !% Declarations
            integer :: mm, kk, error
            integer, pointer :: outfallIdx
            real(8), pointer :: thisVolume
        !%---------------------------------------------------------------------

        !% --- cycle through the subcatchments
        do mm=1,setting%SWMMinput%N_subcatch
            call load_api_procedure("api_export_runon_volume")
            !% 
            if (subcatchYN(mm,sYN_hasRunOn)) then
                do kk=1,subcatchI(mm,si_RunOn_count)        
                    !% --- pointers for clarity
                    outfallIdx   => subcatchI(mm,si_RunOn_SWMMoutfallIdx(kk))
                    thisVolume   => subcatchR(mm,sr_RunOn_Volume(kk))
                    !% --- export the runon flowrate
                    error = ptr_api_export_runon_volume(outfallIdx, thisVolume)
                end do
            end if
        end do

    end subroutine interface_export_runon_volume
!%
!%=============================================================================
!%=============================================================================
!%   
    subroutine interface_call_runoff_execute()
        !%---------------------------------------------------------------------
        !% Description:
        !% calls the EPA SWMM runoff procedures
        !%---------------------------------------------------------------------
            integer       :: error
            character(64) :: subroutine_name = 'interface_call_runoff_execute'     
        !%---------------------------------------------------------------------   
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call load_api_procedure("api_call_runoff_execute")
        error = ptr_api_call_runoff_execute()   
        call print_api_error(error, subroutine_name)

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
        
    end subroutine interface_call_runoff_execute
!%
!%=============================================================================
!%=============================================================================
!%   
    real(8) function interface_get_subcatch_runoff (sc_idx) result(runoff)
        !%---------------------------------------------------------------------
        !% Description:
        !% Gets the SWMM-C "newRunoff" for j=sc_idx for the data
        !% Subcatch[j].newRunoff in subcatch.c/subcatch_getRunoff
        !%---------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: sc_idx  !% index of this subcatchment
            integer             :: error
            character(64) :: subroutine_name = 'interface_get_subcatch_runoff'     
        !%---------------------------------------------------------------------   
        !% Preliminaries
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        call load_api_procedure("api_get_subcatch_runoff")
        error = ptr_api_get_subcatch_runoff(sc_idx, runoff)   
        call print_api_error(error, subroutine_name)

        !%---------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%interface)  &
             write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
        
    end function interface_get_subcatch_runoff
!%
!%=============================================================================
!%=============================================================================
!%  
    integer function interface_get_subcatch_runoff_nodeIdx (sc_idx) result(node_idx)    
        !%---------------------------------------------------------------------
        !% Description:
        !% Gets the SWMM-C node index that a subcatchment runs off to
        !%---------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: sc_idx  !% index of this subcatchment
            integer             :: error
            character(64) :: subroutine_name = 'interface_get_subcatch_runoff_nodeIdx'     
        !%---------------------------------------------------------------------   
        !% Preliminaries
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        call load_api_procedure("api_get_subcatch_runoff_nodeIdx")
        error = ptr_api_get_subcatch_runoff_nodeIdx(sc_idx, node_idx)   
        call print_api_error(error, subroutine_name)
        
        !%---------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end function interface_get_subcatch_runoff_nodeIdx
!%
!%=============================================================================
!%=============================================================================
!%  
    subroutine interface_call_climate_setState ()
        !%---------------------------------------------------------------------
        !% Description:
        !% calls the EPA SWMM procedures to set the current climate state
        !%---------------------------------------------------------------------
            integer       :: error
            character(64) :: subroutine_name = 'interface_call_climate_setState'     
        !%---------------------------------------------------------------------  

        call load_api_procedure("api_call_climate_setState")
        error = ptr_api_call_climate_setState(setting%Time%Now)   
        call print_api_error(error, subroutine_name)

    end subroutine interface_call_climate_setState
!%
!%=============================================================================
!%=============================================================================
!% 
    real(8) function interface_get_evaporation_rate () result(evapRate)
        !%---------------------------------------------------------------------
        !% Description:
        !% Gets the SWMM-C evaporation rate, converted to m/s
        !%---------------------------------------------------------------------
        !% Declarations
            integer       :: error
            character(64) :: subroutine_name = 'interface_get_evaporation_rate'     
        !%---------------------------------------------------------------------   

        call load_api_procedure("api_get_evaporation_rate")
        error = ptr_api_get_evaporation_rate (evapRate)
        call print_api_error(error, subroutine_name)

    end function interface_get_evaporation_rate    
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_get_RDII_inflow ()
        !%---------------------------------------------------------------------
        !% Description:
        !% Mimics section in EPA SWMM routing.C where addRdiiInflows() is called
        !%---------------------------------------------------------------------
        !% Declarations
            integer             :: nRDII, ii, nIdx
            integer             :: error
            integer, pointer    :: eIdx
            real(8)             :: flowrate
            character(64)       :: subroutine_name = 'interface_get_RDII_inflow'
        !%---------------------------------------------------------------------
        !%---------------------------------------------------------------------

        !% --- get the number of RDII flows
        call load_api_procedure("api_getNumRdiiFlows")
        error = ptr_api_getNumRdiiFlows (setting%Time%Now,nRDII) 
        call print_api_error(error, subroutine_name)

        !% --- if there are no RDII
        if (nRDII < 1) return

        !% --- cycle through all the RDII nodes
        !%     HACK: in the future we should bring the RDII functions into SWMM5+ so that
        !%     the operations can be vectorized.
        do ii=1,nRDII
            !% --- get the RDII flows
            !%     input/output are Fortran indexes
            call load_api_procedure("api_getRdiiFlow")
            error = ptr_api_getRdiiFlow(ii, nIdx, flowrate)
            call print_api_error(error, subroutine_name)

            !% --- note that EPA SWMM uses a FLOW_TOL = 1e-5 cfs, not sure if we need to replicate
            !%     however it seems a good idea to prevent a negative value in case something goes
            !%     wonky in the RDII file.
            if (flowrate < zeroR) cycle

            !% --- check if node is on this processor -- if so then add RDII to lateral flowrate
            if (node%I(nIdx,ni_P_image) == this_image()) then
                !% --- get the element index
                eIdx => node%I(nIdx,ni_elem_idx)
                !% --- add RDII to lateral flowrate 
                elemR(eIdx,er_FlowrateLateral) = elemR(eIdx,er_FlowrateLateral) + flowrate
            else
                !% --- no action
            end if

        end do

    end subroutine interface_get_RDII_inflow
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_get_groundwater_inflow ()
        !%---------------------------------------------------------------------
        !% Description:
        !% mimics the call to routing.c/addGroundwaterInflows in EPA SWMM
        !%---------------------------------------------------------------------
        !% Declarations:
            integer             :: error, nIdx, ii
            integer, pointer    :: eIdx
            real(8)             :: flowrate
            real(8), pointer    :: NextRunoffTime, LastRunoffTime, thisTime
            character(64)       :: subroutine_name = 'interface_get_groundwater_inflow'
        !%---------------------------------------------------------------------
        !% Aliases:
            thisTime => setting%Time%Now
            NextRunoffTime => setting%Time%Hydrology%NextTime
            LastRunoffTime => setting%Time%Hydrology%LastTime
        !%---------------------------------------------------------------------
        do ii=1,setting%SWMMinput%N_subcatch
            !% --- get the node and flowrate for each groundwater inflow
            call load_api_procedure("api_get_groundwaterFlow")
            error = ptr_api_get_groundwaterFlow &
                (thisTime, LastRunoffTime, NextRunoffTime, ii, nIdx, flowrate)
            call print_api_error(error, subroutine_name)

            !% --- check for no node returned
            if (nIdx < 1) cycle
        
            !% --- check if node is on this processor -- if so then add groundwater to lateral flowrate
            if (node%I(nIdx,ni_P_image) == this_image()) then
                !% --- get the element index
                eIdx => node%I(nIdx,ni_elem_idx)
                !% --- add RDII to lateral flowrate 
                elemR(eIdx,er_FlowrateLateral) = elemR(eIdx,er_FlowrateLateral) + flowrate
            else
                !% --- no action
            end if
        end do

    end subroutine interface_get_groundwater_inflow
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_get_LID_inflow ()
        !%---------------------------------------------------------------------
        !% Description:
        !% gets the LID inflows into nodes
        !%---------------------------------------------------------------------
        !% Declarations:
            integer             :: error, nIdx, ii
            integer, pointer    :: eIdx
            real(8)             :: flowrate
            real(8), pointer    :: NextRunoffTime, LastRunoffTime, thisTime
            character(64)       :: subroutine_name = 'interface_get_LID_inflow'
        !%---------------------------------------------------------------------
        !% Aliases
            thisTime      => setting%Time%Now
            NextRunoffTime => setting%Time%Hydrology%NextTime
            LastRunoffTime => setting%Time%Hydrology%LastTime
        !%---------------------------------------------------------------------

        nIdx    = zeroI
        flowrate = zeroR

        do ii=1,setting%SWMMinput%N_subcatch
            !% --- get the node and flowrate for each groundwater inflow
            call load_api_procedure("api_get_LID_DrainFlow")
            error = ptr_api_get_LID_DrainFlow &
                (thisTime, LastRunoffTime, NextRunoffTime, ii, nIdx, flowrate)
            call print_api_error(error, subroutine_name)

            !% --- check for no node returned
            if (nIdx < 1) cycle
        
            !% --- check if node is on this processor -- if so then add groundwater to lateral flowrate
            if (node%I(nIdx,ni_P_image) == this_image()) then
                !% --- get the element index
                eIdx => node%I(nIdx,ni_elem_idx)
                !% --- add RDII to lateral flowrate 
                elemR(eIdx,er_FlowrateLateral) = elemR(eIdx,er_FlowrateLateral) + flowrate
            else
                !% --- no action
            end if
        end do

    end subroutine interface_get_LID_inflow    
!%
!%=============================================================================
!% PRIVATE
!%=============================================================================
!%
    subroutine load_api_procedure(api_procedure_name)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
            character(kind=c_char, len=*) :: api_procedure_name
            character(64) :: subroutine_name = 'load_api_procedure'
        !%---------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name)  // " (" // api_procedure_name // ")" // &
                    " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        !% --- Open the shared library
        call c_lib_open(c_lib, errstat, errmsg)
        if (crashI==1) return
        c_lib%procname = api_procedure_name
        call c_lib_load(c_lib, errstat, errmsg)
        if (crashI==1) return

        !% --- Loads shared library funcitonalities
        select case (api_procedure_name)
            case ("api_teststuff")
                call c_f_procpointer(c_lib%procaddr, ptr_api_teststuff)
            case ("api_controls_count")
                call c_f_procpointer(c_lib%procaddr, ptr_api_controls_count)
            case ("api_controls_get_premise_data")
                call c_f_procpointer(c_lib%procaddr, ptr_api_controls_get_premise_data)
            case ("api_controls_get_action_data")
                call c_f_procpointer(c_lib%procaddr, ptr_api_controls_get_action_data)
            case ("api_controls_transfer_monitor_data")
                call c_f_procpointer(c_lib%procaddr, ptr_api_controls_transfer_monitor_data)
            case ("api_controls_execute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_controls_execute)
            case ("api_initialize")
                call c_f_procpointer(c_lib%procaddr, ptr_api_initialize)
            case ("api_finalize")
                call c_f_procpointer(c_lib%procaddr, ptr_api_finalize)
            case ("api_get_nodef_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_nodef_attribute)
            case ("api_get_linkf_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_linkf_attribute)
            case ("api_get_adjustments")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_adjustments)
            case ("api_get_transectf_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_transectf_attribute)
            case ("api_get_N_TRANSECT_TBL")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_N_TRANSECT_TBL)
            case ("api_get_transect_table")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_transect_table)        
            case ("api_get_table_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_table_attribute)
            case ("api_get_num_objects")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_num_objects)
            case ("api_get_object_name_len")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_object_name_len)
            case ("api_get_object_name")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_object_name)
            case ("api_get_num_table_entries")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_num_table_entries)
            case ("api_get_first_entry_table")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_first_entry_table)
            case ("api_get_next_entry_table")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_next_entry_table)
            case ("api_get_start_datetime")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_start_datetime)
            case ("api_get_end_datetime")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_end_datetime)
            case ("api_get_flowBC")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_flowBC)
            case ("api_get_headBC")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_headBC)
            case ("api_count_subobjects")
                call c_f_procpointer(c_lib%procaddr, ptr_api_count_subobjects)
            case ("api_get_SWMM_setup")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_SWMM_setup)
            case ("api_get_SWMM_times")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_SWMM_times)
            case ("api_get_next_entry_tseries")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_next_entry_tseries)
            case ("api_reset_timeseries_to_start")
                call c_f_procpointer(c_lib%procaddr, ptr_api_reset_timeseries_to_start)    
            case ("api_find_object")
                call c_f_procpointer(c_lib%procaddr, ptr_api_find_object)
            case ("api_run_step")
                call c_f_procpointer(c_lib%procaddr, ptr_api_run_step)
            case ("api_export_link_results")
                call c_f_procpointer(c_lib%procaddr, ptr_api_export_link_results)
            case ("api_write_output_line")
                call c_f_procpointer(c_lib%procaddr, ptr_api_write_output_line)
            case ("api_update_nodeResult")
                call c_f_procpointer(c_lib%procaddr, ptr_api_update_nodeResult)
            case ("api_update_linkResult")
                call c_f_procpointer(c_lib%procaddr, ptr_api_update_linkResult)
            case ("api_export_runon_volume")    
                call c_f_procpointer(c_lib%procaddr, ptr_api_export_runon_volume)
            case ("api_call_runoff_execute")    
                call c_f_procpointer(c_lib%procaddr, ptr_api_call_runoff_execute)
            case ("api_get_subcatch_runoff")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_subcatch_runoff)    
            case ("api_get_subcatch_runoff_nodeIdx")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_subcatch_runoff_nodeIdx)        
            case ("api_get_NewRunoffTime")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_NewRunoffTime)
            case ("api_getNumRdiiFlows")
                call c_f_procpointer(c_lib%procaddr, ptr_api_getNumRdiiFlows)
            case ("api_getRdiiFlow")
                call c_f_procpointer(c_lib%procaddr, ptr_api_getRdiiFlow)
            case ("api_get_groundwaterFlow")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_groundwaterFlow)
            case ("api_get_LID_DrainFlow")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_LID_DrainFlow)
            case ("api_call_climate_setState")
                call c_f_procpointer(c_lib%procaddr, ptr_api_call_climate_setState)
            case ("api_get_evaporation_rate")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_evaporation_rate)
            case default
                write(*,"(A,A)") "Error, procedure " // api_procedure_name // &
                 " has not been handled in load_api_procedure in interface.f90"
                call util_crashpoint(420987)
        end select

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name)  // " (" // api_procedure_name // ")" // &
                    " [Processor ", this_image(), "]"

    end subroutine load_api_procedure
!%
!%=============================================================================
!%=============================================================================
!%
    function get_next_entry_tseries(tseries_idx,timemax) result(success)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
            integer, intent(in   ) :: tseries_idx
            real(8), intent(in   ) :: timemax
            integer                :: success
            character(64)          :: subroutine_name
        !%---------------------------------------------------------------------
            subroutine_name = 'get_next_entry_tseries'
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call load_api_procedure("api_get_next_entry_tseries")
        success = ptr_api_get_next_entry_tseries(tseries_idx-1,timemax) ! Fortran to C convention

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function get_next_entry_tseries
!%
!%=============================================================================
!%=============================================================================
!%
    function get_num_objects(obj_type)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
            integer :: obj_type
            integer :: get_num_objects
            character(64) :: subroutine_name
        !%---------------------------------------------------------------------
            subroutine_name = 'get_num_objects'
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call load_api_procedure("api_get_num_objects")
        get_num_objects = ptr_api_get_num_objects(obj_type)

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function get_num_objects
!%
!%=============================================================================
!%=============================================================================
!%
    function get_start_datetime()
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
            real(c_double) :: get_start_datetime
            character(64) :: subroutine_name
        !%---------------------------------------------------------------------
            subroutine_name = 'get_start_datetime'
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call load_api_procedure("api_get_start_datetime")
        get_start_datetime = ptr_api_get_start_datetime()

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function get_start_datetime
!%
!%=============================================================================
!%=============================================================================
!%
    function get_end_datetime()
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
            real(8) :: get_end_datetime
            character(64) ::  subroutine_name
        !%---------------------------------------------------------------------
            subroutine_name = 'get_end_datetime'
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call load_api_procedure("api_get_end_datetime")
        get_end_datetime = ptr_api_get_end_datetime()

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function get_end_datetime
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_NewRunoffTime()
        !%---------------------------------------------------------------------
        !% Description:
        !% gets the latest NewRunoffTime from SWMM-C (in milliseconds)
        !% and converts to seconds
        !%---------------------------------------------------------------------
            real(8)       :: interface_get_NewRunoffTime
            character(64) :: subroutine_name
        !%---------------------------------------------------------------------
            subroutine_name = 'interface_get_NewRunoffTime'
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call load_api_procedure("api_get_NewRunoffTime")
        interface_get_NewRunoffTime = ptr_api_get_NewRunoffTime() / onethousandR

        !%------------------------------------------------------------------
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_get_NewRunoffTime
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine print_api_error(error, subroutine_name)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
            integer, intent(in) :: error
            character(64), intent(in) :: subroutine_name
        !%---------------------------------------------------------------------

        if (error /= 0) then
            write(*,*)
            write(*,*)
            write(*,*) "************************************************************"
            write(*,*) "*   EPA-SWMM USER CONFIGURATION ERROR, SIMULATION STOPPED  *"
            write(*,*) "************************************************************"
            write(*,*) "EPA-SWMM Error Code: ", error
            write(*,*) "Failure occurred in subroutine: "// trim(subroutine_name)
            write(*,*) "See SWMM error list in report file, located at: "
            write(*,*) trim(setting%File%rpt_file)
            write(*,*) "SWMM User Manual (Appendix E) has details on error codes."
            write(*,*)
            call util_crashpoint( 63455)
        end if
    end subroutine print_api_error
!%
!%=============================================================================
!% END MODULE interface
!%=============================================================================
end module interface_