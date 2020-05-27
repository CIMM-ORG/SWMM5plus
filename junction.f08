! module junction
!
! Procedures associated with junctions that have more than 2 faces
!
!==========================================================================
!
module junction
   !
   use array_index
   use data_keys
   use globals
   use setting_definition
   use utility


   implicit none

   private

   public  :: junction_adjacent_element_average
   public  :: junction_adjacent_element_values_to_branches
   public  :: junction_branch_assigned_to_faces
   public  :: junction_branch_average_of_inflows_and_outflows
   public  :: junction_branch_velocity_and_flowrate_proportional_to_area
   public  :: junction_branch_sum_areas_by_direction
   public  :: junction_branch_velocities
   public  :: junction_geometry_setup
   !
   integer :: debuglevel = 0

contains
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_adjacent_element_average &
      (elem2R, elemMR, elemMI, faceI, e2r_data, eMr_out)
      !
      ! this computes the average of values for all the elements upstream and
      ! downstream. Note that this should ONLY be use
      ! in setup routines (i.e. when initializing junctions). This violates
      ! the "no-neighbor" rule and is time-consuming because it requires double
      ! mapping
      !
      character(64) :: subroutine_name = 'junction_adjacent_element_average'

      real,      target,     intent(in out)  :: elemMR(:,:)
      real,                  intent(in)      :: elem2R(:,:)
      integer,               intent(in)      :: elemMI(:,:), faceI(:,:)
      integer,               intent(in)      :: e2r_data, eMr_out

      real,      pointer :: Uvalue(:), Dvalue(:)
      integer :: eMr_tUp, eMr_tDn

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      eMr_tUp = eMr_Temp(next_eMr_temparray)
      Uvalue  => elemMR(:,eMr_tUp)
      next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

      eMr_tDn = eMr_Temp(next_eMr_temparray)
      Dvalue  => elemMR(:,eMr_tDn)
      next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

      Uvalue = zeroR
      call junction_summation_from_adjacent_elements_one_direction &
         (eMr_tUp, elem2R, elemMR, elemMI, faceI, &
         upstream_face_per_elemM, eMi_nfaces_u, eMi_MfaceUp, fi_Melem_u, e2r_data)

      Dvalue = zeroR
      call junction_summation_from_adjacent_elements_one_direction &
         (eMr_tDn, elem2R, elemMR, elemMI, faceI, &
         dnstream_face_per_elemM, eMi_nfaces_d, eMi_MfaceDn, fi_Melem_d, e2r_data)

      where (elemMI(:,eMi_elem_type) == eJunctionChannel)
         elemMR(:,eMr_out) = (Uvalue + Dvalue) / real( elemMI(:,eMi_nfaces_u) + elemMI(:,eMi_nfaces_d) )
      endwhere

      Dvalue = nullvalueR
      Uvalue = nullvalueR
      nullify(Dvalue,Uvalue)
      next_eMr_temparray = next_eMr_temparray-2

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_adjacent_element_average
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_adjacent_element_values_to_branches &
      (elem2R, elemMR, elemMI, faceI, e2r_data, eMr_outUp, eMr_outDn)
      !
      ! gets the values from the adjacent elements and stores in the branches
      !
      ! THIS SHOULD ONLY BE CALLED DURING INITIAL CONDITIONS OR SETUP ROUTINES
      !
      character(64) :: subroutine_name = 'junction_adjacent_element_values_to_branches'

      real,      intent(in out)  :: elemMR(:,:)
      real,      intent(in)      :: elem2R(:,:)
      integer,   intent(in)      :: elemMI(:,:), faceI(:,:), eMr_outUp(:), eMr_outDn(:), e2r_data


      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      call junction_adjacent_element_values_one_direction &
         (elem2R, elemMR, elemMI, faceI, &
         upstream_face_per_elemM, eMi_nfaces_u, eMi_MfaceUp, fi_Melem_u, &
         e2r_data, eMr_outUp)

      call junction_adjacent_element_values_one_direction &
         (elem2R, elemMR, elemMI, faceI, &
         dnstream_face_per_elemM, eMi_nfaces_d, eMi_MfaceDn, fi_Melem_d, &
         e2r_data, eMr_outDn)

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_adjacent_element_values_to_branches
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_assigned_to_faces &
      (faceI, elemMI)

      character(64) :: subroutine_name = 'junction_branch_assigned_to_faces'

      integer,   intent(in out)  :: faceI(:,:)
      integer,   intent(in)      :: elemMI(:,:)

      integer :: mm

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      do mm=1,upstream_face_per_elemM
         where ((elemMI(:,eMi_nfaces_u) >= mm) .and. &
            (elemMI(:,eMi_elem_type) == eJunctionChannel))
            faceI( elemMI(: , eMi_MfaceUp(mm)),fi_branch_d) = mm
         endwhere
      end do

      do mm=1,dnstream_face_per_elemM
         where ((elemMI(:,eMi_nfaces_d) >= mm) .and. &
            (elemMI(:,eMi_elem_type) == eJunctionChannel))
            faceI( elemMI(: , eMi_MfaceDn(mm)),fi_branch_u) = mm
         endwhere
      end do

      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
   end subroutine junction_branch_assigned_to_faces
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_average_of_inflows_and_outflows &
      (elemMR, elemMI)
      !
      ! Computes an average flowrate based on the branch inflows and branch outflows
      ! Includes effects of flow reverals
      !
      character(64) :: subroutine_name = 'junction_branch_average_of_inflows_and_outflows'

      real,      target, intent(in out)      :: elemMR(:,:)
      integer,           intent(in)          :: elemMI(:,:)

      integer        :: eMr_inflow, eMr_outflow
      real,  pointer :: inflow(:), outflow(:)

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      eMr_inflow = eMr_Temp(next_eMr_temparray)
      inflow  => elemMR(:,eMr_inflow)
      next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

      eMr_outflow = eMr_Temp(next_eMr_temparray)
      outflow  => elemMR(:,eMr_outflow)
      next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

      call junction_net_inflow_and_outflow &
         (eMr_inflow, eMr_outflow, elemMR, elemMI)

      !% simple average of inflow and outflow
      where (elemMI(:,eMi_elem_type) == eJunctionChannel)
         elemMR(:,eMr_Flowrate) = onehalfR * (inflow + outflow)
      endwhere

      inflow = nullvalueR
      outflow = nullvalueR
      nullify(inflow,outflow)
      next_eMr_temparray = next_eMr_temparray-2

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_branch_average_of_inflows_and_outflows
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_velocity_and_flowrate_proportional_to_area &
      (eMR_totalarea, eMR_totalflowrate, &
      this_face_per_elem, eMr_AreaThis, eMi_MfaceThis, eMi_nfaces_This, &
      eMr_FlowrateThis, eMr_VelocityThis, &
      rdir_face_per_elem, eMr_AreaRdir, eMi_MfaceRdir, eMi_nfaces_Rdir, &
      eMr_FlowrateRdir, eMr_VelocityRdir, &
      elemMR, elemMI, faceR)
      !
      ! Flowrates in each branch of a junction
      ! Total flowrate is distributed proportionally over the areas
      !
      ! HACK - this can be cleaned up with 2 calls to a function, but we have to be
      ! careful that we don't introduce pass-by-value in the call
      !
      character(64) :: subroutine_name = 'junction_branch_velocity_and_flowrate_proportional_to_area'

      integer, intent(in) :: eMi_nfaces_This,     eMi_nfaces_Rdir
      integer, intent(in) :: eMR_totalarea,       eMR_totalflowrate
      integer, intent(in) :: this_face_per_elem,  rdir_face_per_elem

      integer, intent(in) :: eMr_AreaThis(:),     eMi_MfaceThis(:)
      integer, intent(in) :: eMr_AreaRdir(:),     eMi_MfaceRdir(:)
      integer, intent(in) :: eMr_FlowrateThis(:), eMr_FlowrateRdir(:)
      integer, intent(in) :: eMr_VelocityThis(:), eMr_VelocityRdir(:)

      real,      target, intent(in out)  :: elemMR(:,:)
      real,              intent(in)      :: faceR(:,:)
      integer,   target, intent(in)      :: elemMI(:,:)

      real,      pointer :: area(:), totalarea(:), totalflowrate(:)
      real,      pointer :: flowrate(:), velocity(:)
      integer,   pointer :: fThis(:), fRdir(:)

      integer :: mm

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      totalarea     => elemMR(:,eMR_totalarea)
      totalflowrate => elemMR(:,eMR_totalflowrate)

      !%  distribute flow proportionally over the downstream outflow branches
      !%  or opposite call for the upstream inflow branches
      do mm=1,this_face_per_elem
         area     => elemMR(:, eMr_AreaThis(mm))
         flowrate => elemMR(:, eMr_FlowrateThis(mm))
         velocity => elemMR(:, eMr_VelocityThis(mm))
         fThis    => elemMI(:, eMi_MfaceThis(mm))
         where ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
            (elemMI(:,eMi_nfaces_This) >= mm) .and. &
            (faceR(fThis,fr_Flowrate) >= 0.0) )
            flowrate = totalflowrate * area / totalarea
            velocity = flowrate / area
         endwhere

      enddo

      !%  distribute flow proportionally over the upstream outflow branches
      !%  or opposite call for downstream inflow branches
      do mm=1,rdir_face_per_elem
         area     => elemMR(:, eMr_AreaRdir(mm))
         flowrate => elemMR(:, eMr_FlowrateRdir(mm))
         velocity => elemMR(:, eMr_VelocityRdir(mm))
         fRdir    => elemMI(:, eMi_MfaceRdir(mm))
         where ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
            (elemMI(:,eMi_nfaces_Rdir) >= mm) .and. &
            (faceR(fRdir,fr_Flowrate) < 0.0) )
            flowrate = -totalflowrate * area / totalarea
            velocity =  flowrate / area
         endwhere
      enddo

      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
   end subroutine junction_branch_velocity_and_flowrate_proportional_to_area
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_sum_areas_by_direction &
      (eMR_totalarea, &
      this_face_per_element, eMr_AreaThis, eMi_MfaceThis, eMi_nfaces_This, &
      rdir_face_per_element, eMr_AreaRdir, eMi_MfaceRdir, eMi_nfaces_Rdir, &
      elemMR, elemMI, faceR)
      !
      ! Sum of the areas in each branch of a junction for outflows or inflows
      ! Called for outflow areas with this = downstream, and rdir = upstream
      ! Called for inflow  areas with this = upstream    and rdir = downstream
      !
      ! HACK- This could be cleaned up with two calls to a separate function,
      ! but we have to be careful that we don'tend up introducing a pass-by-value
      ! into the function.
      !
      character(64) :: subroutine_name = 'junction_branch_sum_areas_by_direction'

      integer, intent(in) :: eMR_totalarea
      integer, intent(in) :: this_face_per_element, rdir_face_per_element

      integer, intent(in) :: eMr_AreaThis(:), eMi_MfaceThis(:)
      integer, intent(in) :: eMr_AreaRdir(:), eMi_MfaceRdir(:)
      integer, intent(in) :: eMi_nfaces_This, eMi_nfaces_Rdir

      real,      target, intent(in out)  :: elemMR(:,:)
      real,              intent(in)      :: faceR(:,:)
      integer,   target, intent(in)      :: elemMI(:,:)

      real,      pointer :: area(:), totalarea(:)
      integer,   pointer :: fThis(:), fRdir(:)

      integer :: mm

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


      totalarea => elemMR(:,eMR_totalarea)

      !print *, trim(subroutine_name)
      !print *, this_face_per_element
      !stop

      !%  get the outflow area for downstream branches (or inflow area for upstream)
      do mm=1,this_face_per_element
         area  => elemMR(:,eMr_AreaThis(mm))
         fThis => elemMI(:,eMi_MfaceThis(mm))

         where ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
            (elemMI(:,eMi_nfaces_This) >= mm))
            where (faceR(fThis,fr_Flowrate) >= 0.0)
               totalarea= totalarea + area
            endwhere
         endwhere
      enddo


      !%  add the area for any reversing upstream branches (or reversing downstream)
      do mm=1,rdir_face_per_element
         area  => elemMR(:,eMr_AreaRdir(mm))
         fRdir => elemMI(:,eMi_MfaceRdir(mm))
         where ( (elemMI(:,eMi_elem_type) == eJunctionChannel) .and. &
            (elemMI(:,eMi_nfaces_Rdir) >= mm) .and. &
            (faceR(fRdir,fr_Flowrate) < 0.0) )
            totalarea = totalarea + area
         endwhere
      enddo

      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
   end subroutine junction_branch_sum_areas_by_direction
   !
   !========================================================================== !==========================================================================
   !
   subroutine junction_branch_velocities &
      (elemMR, elemMI)
      !
      ! computes velocities from flowrates and areas in both upstream and downstream
      ! branches
      !
      character(64) :: subroutine_name = 'junction_branch_velocities'

      real,      intent(in out)  :: elemMR(:,:)
      integer,   intent(in)      :: elemMI(:,:)

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      call junction_branch_velocities_one_direction &
         (elemMR, elemMI, upstream_face_per_elemM, eMi_nfaces_u, &
         eMr_FlowrateUp, eMr_AreaUp, eMr_VelocityUp)

      call junction_branch_velocities_one_direction &
         (elemMR, elemMI, dnstream_face_per_elemM, eMi_nfaces_d, &
         eMr_FlowrateDn, eMr_AreaDn, eMr_VelocityDn)

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_branch_velocities
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_geometry_setup  &
      (elemMR, elemMI)
      !
      ! Get the setup for junction geometry
      ! This gets the breadthscale, topwidth, and length of the junction
      ! based on the branch values.
      !
      character(64) :: subroutine_name = 'junction_geometry_setup'

      real,      target,      intent(in out)  :: elemMR(:,:)
      integer,                intent(in)      :: elemMI(:,:)
      real,      pointer                      :: Uvalue(:), Dvalue(:)

      integer :: ii, eMr_tUp, eMr_tDn

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


      !% For the breadthscale and topwidth
      !% take the sum of the upstream breadths and the sum of the downstream breadths
      !% then compute their average

      call junction_branch_summation_and_updown_average &
         (elemMR, elemMI, eMr_BreadthScaleUp, eMr_BreadthScaleDn, eMr_BreadthScale)

      call junction_branch_summation_and_updown_average &
         (elemMR, elemMI, eMr_TopwidthUp, eMr_TopwidthDn, eMr_Topwidth)

      !% For length
      !% take the average lengths of the upstream and of the downstream branches
      !% and then sum these.

      call junction_branch_average_for_directions_then_sum &
         (elemMR, elemMI, eMr_LengthUp, eMr_LengthDn, eMr_Length)


      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
   end subroutine junction_geometry_setup
   !
   !==========================================================================
   !
   ! PRIVATE BELOW HERE
   !
   !==========================================================================
   !
   subroutine junction_adjacent_element_values_one_direction &
      (elem2R, elemMR, elemMI, faceI, &
      dir_face_per_elemM, eMi_nfaces_dir, eMi_MfaceDir, fi_Melem_dir, e2r_data, eMr_outDir)

      character(64) :: subroutine_name = 'junction_adjacent_element_values_one_direction'

      real,              intent(in out)  :: elemMR(:,:)
      real,              intent(in)      :: elem2R(:,:)
      integer,   target, intent(in)      :: elemMI(:,:), faceI(:,:)

      integer,   intent(in)  :: dir_face_per_elemM, eMi_nfaces_dir
      integer,   intent(in)  :: eMi_MfaceDir(:), eMr_outDir(:)
      integer,   intent(in)  :: fi_Melem_dir, e2r_data

      integer,   pointer :: tface, telem
      integer :: mm, ii

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      do ii=1,N_elemM
         do mm=1,dir_face_per_elemM
            if ((elemMI(ii,eMi_nfaces_dir) >= mm) .and. &
               (elemMI(ii,eMi_elem_type) == eJunctionChannel)) then
               !% the face on a branch in the direction specified by dir
               tface => elemMI(ii,eMi_MfaceDir(mm))
               !% the element upstream of the face
               telem => faceI(tface,fi_Melem_dir)
               !% the value at the element
               !print *, ii, mm
               !print *, tface
               !print *, telem
               !print *, elem2R(telem,e2r_data)
               elemMR(ii,eMr_outDir(mm)) = elem2R(telem,e2r_data)
            endif
         enddo
      enddo

      ! print *, eMr_outDir(1), eMr_Area_u1
      ! print *, elemMR(:,eMr_outDir(1)), elemMR(:,eMr_outDir(2))
      ! print *, trim(subroutine_name)
      ! stop
      !
      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_adjacent_element_values_one_direction
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_average &
      (eMr_avgvalue, elemMR, elemMI, eMr_column_up, eMr_column_dn)
      !
      ! compute the average of data type in column_up and column_dn over
      ! all the branches of the junciton
      !
      character(64) :: subroutine_name = 'junction_branch_average'

      integer,           intent(in)      :: eMr_avgvalue
      real,      target, intent(in out)  :: elemMR(:,:)
      integer,           intent(in)      :: elemMI(:,:)
      integer,           intent(in)      :: eMr_column_up(:), eMr_column_dn(:)
      real,              pointer         :: avgvalue(:)
      !--------------------------------------------------------------------------

      avgvalue => elemMR(:,eMr_avgvalue)
      avgvalue = zeroR

      call junction_branch_summation &
         (eMr_avgvalue, elemMR, elemMI, eMr_column_up, eMr_column_dn)

      where (elemMI(:,eMi_elem_type) == eJunctionChannel)
         avgvalue = avgvalue / real(elemMI(:,eMi_nfaces))
      endwhere

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_branch_average
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_average_for_directions_then_sum &
      (elemMR, elemMI, eMr_columnUp, eMr_columnDn, eMr_out)
      !
      ! separately compute the average of the branches in up and down directions
      ! and then sum the result. Applies to data in eMr_columnUp and eMr_columnDn
      ! with result in eMr_out. e.g. eMr_LengthUp, eMr_LengthDn -> eMr_Length
      ! to get a length scale for a junction.
      !
      character(64) :: subroutine_name = 'junction_branch_average_for_directions_then_sum'

      real,      target,      intent(in out)  :: elemMR(:,:)
      integer,                intent(in)      :: elemMI(:,:)
      integer,                intent(in)      :: eMr_columnUp(:), eMr_columnDn(:)
      integer,                intent(in)      :: eMr_out
      real,      pointer                      :: Uvalue(:), Dvalue(:)
      integer :: ii, eMr_tUp, eMr_tDn

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      eMr_tUp = eMr_Temp(next_eMr_temparray)
      Uvalue  => elemMR(:,eMr_tUp)
      next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

      eMr_tDn = eMr_Temp(next_eMr_temparray)
      Dvalue  => elemMR(:,eMr_tDn)
      next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

      Uvalue = zeroR
      call junction_branch_summation_one_direction &
         (eMr_tUp, elemMR, elemMI, upstream_face_per_elemM, eMi_nfaces_u, eMr_columnUp)

      Dvalue = zeroR
      call junction_branch_summation_one_direction &
         (eMr_tDn, elemMR, elemMI, dnstream_face_per_elemM, eMi_nfaces_d, eMr_columnDn)

      where (elemMI(:,eMi_elem_type) == eJunctionChannel)
         Uvalue = Uvalue / real(elemMI(:,eMi_nfaces_u))
         Dvalue = Dvalue / real(elemMI(:,eMi_nfaces_D))
         elemMR(:,eMr_out) = Uvalue + Dvalue
      endwhere

      Dvalue = nullvalueR
      Uvalue = nullvalueR
      nullify(Dvalue,Uvalue)
      next_eMr_temparray = next_eMr_temparray-2

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_branch_average_for_directions_then_sum
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_summation &
      (eMr_sumvalue, elemMR, elemMI, eMr_column_up, eMr_column_dn)
      !
      ! sum all the values in all the branches - both up and down
      !
      character(64) :: subroutine_name = 'junction_branch_summation'

      integer,           intent(in)      :: eMr_sumvalue
      real,      target, intent(in out)  :: elemMR(:,:)
      integer,           intent(in)      :: elemMI(:,:)

      integer,   intent(in)  :: eMr_column_up(:), eMr_column_dn(:)
      real,      pointer     :: sumvalue(:)

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      sumvalue => elemMR(:,eMr_sumvalue)
      sumvalue = zeroR

      call junction_branch_summation_one_direction &
         (eMr_sumvalue, elemMR, elemMI, upstream_face_per_elemM, eMi_nfaces_u, eMr_column_up)

      call junction_branch_summation_one_direction &
         (eMr_sumvalue, elemMR, elemMI, dnstream_face_per_elemM, eMi_nfaces_d, eMr_column_dn)

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_branch_summation
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_summation_and_updown_average &
      (elemMR, elemMI, eMr_columnUp, eMr_columnDn, eMr_out)
      !
      ! separately compute the sum of the upstream and the downstream branches
      ! for data in eMr_columnUp and eMr_columnDn (e.g. eMr_TopwidthUp)
      ! and then computes the average between the two sums. Output is stored in
      ! elemMR(:,eMr_out)
      ! NOTE THIS DOES NOT ACCOUNT FOR FLOW REVERSALS
      !
      character(64) :: subroutine_name = 'junction_branch_summation_and_updown_average'

      real,      target,      intent(in out)  :: elemMR(:,:)
      integer,                intent(in)      :: elemMI(:,:)
      integer,                intent(in)      :: eMr_columnUp(:), eMr_columnDn(:)
      integer,                intent(in)      :: eMr_out
      real,      pointer                      :: Uvalue(:), Dvalue(:)
      integer :: ii, eMr_tUp, eMr_tDn

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      eMr_tUp = eMr_Temp(next_eMr_temparray)
      Uvalue  => elemMR(:,eMr_tUp)
      next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

      eMr_tDn = eMr_Temp(next_eMr_temparray)
      Dvalue  => elemMR(:,eMr_tDn)
      next_eMr_temparray = utility_advance_temp_array (next_eMr_temparray,eMr_n_temp)

      Uvalue = zeroR
      call junction_branch_summation_one_direction &
         (eMr_tUp, elemMR, elemMI, upstream_face_per_elemM, eMi_nfaces_u, eMr_BreadthScaleUp)

      Dvalue = zeroR
      call junction_branch_summation_one_direction &
         (eMr_tDn, elemMR, elemMI, dnstream_face_per_elemM, eMi_nfaces_d, eMr_BreadthScaleDn)

      where (elemMI(:,eMi_elem_type) == eJunctionChannel)
         elemMR(:,eMr_out) = onehalfR * ( Uvalue + Dvalue)
      endwhere

      Dvalue = nullvalueR
      Uvalue = nullvalueR
      nullify(Dvalue,Uvalue)
      next_eMr_temparray = next_eMr_temparray-2

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_branch_summation_and_updown_average
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_summation_one_direction &
      (eMr_sumvalue, elemMR, elemMI, face_per_elemM, eMi_nfaces_dir, eMr_columnDir)
      !
      ! cycle through the branches in one direction (up or dn) to sum all their values
      ! NOTE THAT THIS DOES NOT RESET THE SUM TO ZERO BEFORE EXECUTION
      ! This behavior is required for accumulators.
      !
      character(64) :: subroutine_name = 'junction_branch_summation_one_direction'

      integer,   intent(in)      :: eMr_sumvalue
      real,      intent(in out)  :: elemMR(:,:)
      integer,   intent(in)      :: elemMI(:,:)
      integer,   intent(in)      :: face_per_elemM, eMi_nfaces_dir, eMr_columnDir(:)

      integer :: ii

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      do ii=1,face_per_elemM
         where ( elemMI(:,eMi_nfaces_dir) >= ii)
            elemMR(:,eMr_sumvalue) = elemMR(:,eMr_sumvalue) + elemMR(:,eMr_columnDir(ii))
         endwhere
      enddo

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_branch_summation_one_direction
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_branch_velocities_one_direction &
      (elemMR, elemMI, dir_face_per_elemM, eMi_nfaces_dir, &
      eMr_FlowrateDir, eMr_AreaDir, eMr_VelocityDir)
      !
      ! computes the velocities in each branch (with dir = upstream or downstream)
      ! from flowrate and area
      !
      character(64) :: subroutine_name = 'junction_branch_velocities_one_direction'

      real,      intent(in out)  :: elemMR(:,:)
      integer,   intent(in)      :: elemMI(:,:)
      integer,   intent(in)      :: eMr_FlowrateDir(:), eMr_AreaDir(:), eMr_VelocityDir(:)
      integer,   intent(in)      :: dir_face_per_elemM, eMi_nfaces_dir

      integer :: ii

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      do ii=1,dir_face_per_elemM
         where ((elemMI(:,eMi_nfaces_dir) >= ii) .and. &
            (elemMI(:,eMi_elem_type) == eJunctionChannel))
            elemMR(:,eMr_VelocityDir(ii)) = elemMR(:,eMr_FlowrateDir(ii)) / elemMR(:,eMr_AreaDir(ii))
         endwhere
      enddo

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_branch_velocities_one_direction
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_net_flow_in_or_out &
      (eMr_flow, elemMR, elemMI, dir_face_per_elemM, eMi_nfaces_dir, &
      eMr_FlowrateDir1, eMr_FlowrateDir2)

      character(64) :: subroutine_name = 'junction_net_flow_in_or_out'

      integer,           intent(in)      :: eMr_flow
      real,      target, intent(in out)  :: elemMR(:,:)
      integer,           intent(in)      :: elemMI(:,:)
      integer,           intent(in)      :: dir_face_per_elemM, eMi_nfaces_dir
      integer,           intent(in)      :: eMr_FlowrateDir1(:), eMr_flowrateDir2(:)

      real,      pointer         :: flow(:)
      integer :: mm
      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      flow => elemMR(:,eMr_flow)

      flow = zeroR
      do mm=1,dir_face_per_elemM
         where ((elemMI(:,eMi_nfaces_dir) >= mm) .and. &
            (elemMI(:,eMi_elem_type) == eJunctionChannel))
            where (elemMR(:,eMr_FlowrateDir1(mm)) > zeroR)
               flow = flow + elemMR(:,eMr_FlowrateDir1(mm))
            endwhere
            where (elemMR(:,eMr_FlowrateDir2(mm)) < zeroR)
               flow = flow  - elemMR(:,eMr_FlowrateDir2(mm))
            endwhere
         endwhere
      enddo

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_net_flow_in_or_out
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_net_inflow_and_outflow &
      (eMr_inflow, eMr_outflow, elemMR, elemMI)
      !
      ! Computes the net inflow and net outflow including effects of flow reversals
      ! Stores results in the columns eMr_inflow and eMr_outflow.
      !
      character(64) :: subroutine_name = 'junction_net_inflow_and_outflow'

      integer,   intent(in)      :: eMr_inflow, eMr_outflow
      real,      intent(in out)  :: elemMR(:,:)
      integer,   intent(in)      :: elemMI(:,:)

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

      call junction_net_flow_in_or_out &
         (eMr_inflow, elemMR, elemMI, upstream_face_per_elemM, eMi_nfaces_u, &
         eMr_FlowrateUp, eMr_flowrateDn)

      call junction_net_flow_in_or_out &
         (eMr_outflow, elemMR, elemMI, dnstream_face_per_elemM, eMi_nfaces_d, &
         eMr_FlowrateDn, eMr_FlowrateUp)

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_net_inflow_and_outflow
   !
   !==========================================================================
   !==========================================================================
   !
   subroutine junction_summation_from_adjacent_elements_one_direction &
      (eMr_sumvalue, elem2R, elemMR, elemMI, faceI, &
      dir_face_per_elemM, eMi_nfaces_dir, eMi_MfaceDir, fi_Melem_dir, e2r_data)

      character(64) :: subroutine_name = 'junction_summation_from_adjacent_elements_one_direction'
      !
      ! computes the sum of all the elements adjacent to a junction in either
      ! the upstream or downstream direction
      !
      ! THIS SHOULD ONLY BE USED IN SETUP AND INITIAL CONDITION ROUTINES
      !
      integer,           intent(in)      :: eMr_sumvalue
      real,      target, intent(in out)  :: elemMR(:,:)
      real,              intent(in)      :: elem2R(:,:)
      integer,   target, intent(in)      :: elemMI(:,:), faceI(:,:)

      integer,   intent(in)  :: dir_face_per_elemM, eMi_nfaces_dir, eMi_MfaceDir(:)
      integer,   intent(in)  :: fi_Melem_dir, e2r_data

      integer,   pointer :: tface, telem
      real   :: thisvalue(dir_face_per_elemM)
      integer :: mm, ii

      !--------------------------------------------------------------------------
      if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


      do ii=1,N_elemM
         thisvalue = zeroR
         do mm=1,dir_face_per_elemM
            if ((elemMI(ii,eMi_nfaces_dir) >= mm) .and. &
               (elemMI(ii,eMi_elem_type) == eJunctionChannel)) then
               !% the face on a branch in the direction specified by dir
               tface => elemMI(ii,eMi_MfaceDir(mm))
               !% the element upstream of the face
               telem => faceI(tface,fi_Melem_dir)
               !% the value at the element
               thisvalue(mm) = elem2R(telem,e2r_data)
            endif
         enddo
         elemMR(ii,eMr_sumvalue) = sum(thisvalue)
      enddo

      if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
   end subroutine junction_summation_from_adjacent_elements_one_direction
   !
   !==========================================================================

   ! END OF MODULE junction
   !==========================================================================
end module junction
