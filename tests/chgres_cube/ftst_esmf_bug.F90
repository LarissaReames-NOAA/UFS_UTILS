 program esmf_regrid_bug

! Unit test for surface routine interp that regrids surface 
! variables from input to target grid. 
!
! Author: Larissa Reames, OU CIMMS/NOAA NSSL

 use esmf
 
 implicit none

 integer, parameter           :: i_input=4
 integer, parameter           :: j_input=3
 integer, parameter           :: i_target=8
 integer, parameter           :: j_target=5

 real, parameter              :: EPSILON=0.0001
 real, parameter              :: veg_type_landice_input = 15.
 real, parameter              :: veg_type_landice_target = 15.
 real(esmf_kind_r8)           :: deltalon

 integer                      :: clb(4), cub(4)
 integer                      :: ierr, localpet, npets, rc
 integer                      :: i, j, k
 integer                      :: isrctermprocessing
 integer(esmf_kind_i4), pointer    :: unmapped_ptr(:)
 integer(esmf_kind_i4), pointer   :: mask_target_ptr(:,:),& 
                                     mask_input_ptr(:,:)
 real(esmf_kind_r8), allocatable  :: latitude(:,:), longitude(:,:)
 integer(esmf_kind_i4), allocatable  :: mask_input(:,:)
 integer(esmf_kind_i4), allocatable  :: mask_target(:,:)
 real(esmf_kind_r8), allocatable  :: sotyp_input(:,:), veg_type_input(:,:)
 real(esmf_kind_r8), allocatable  :: sotyp_target(:,:), veg_type_target(:,:), &
                                     sotyp_correct(:,:)
 real(esmf_kind_r8), pointer      :: lon_ptr(:,:), &
                                     lat_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lon_corner_ptr(:,:), &
                                     lat_corner_ptr(:,:)
 type(esmf_vm)                :: vm
 type(esmf_polekind_flag)     :: polekindflag(2)
 type(esmf_regridmethod_flag) :: method
 type(esmf_routehandle)       :: regrid_land
 type(esmf_grid)              :: input_grid, target_grid
 type(esmf_field)             :: soil_type_input_grid,&
                                 soil_type_target_grid

 print*,"Starting test of surface interp."

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 !--------------------------------------------------------------------!
 !----------------- Setup Input Grid & Coordinates -------------------!
 !--------------------------------------------------------------------!
 
 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 input_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_input,j_input/), &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))  
          
 ! This is a random regional grid. I tried a global grid here but it had an unstable
 ! solution.
 
 deltalon = 2.0_esmf_kind_r8
 do i = 1, i_input
   longitude(i,:) = 90+real((i-1),kind=esmf_kind_r8) * deltalon
 enddo

 do j = 1, j_input
   latitude(:,j) = 35.0-real((j-1),kind=esmf_kind_r8) * deltalon
 end do

 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     lon_ptr(i,j) = longitude(i,j)
     if (lon_ptr(i,j) > 360.0_esmf_kind_r8) lon_ptr(i,j) = lon_ptr(i,j) - 360.0_esmf_kind_r8
     lat_ptr(i,j) = latitude(i,j)
   enddo
 enddo
 nullify(lat_ptr,lon_ptr)
 
! Create staggered coordinates for conservative regridding
 print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

  print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_corner_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_corner_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (i == i_input+1) then
       lon_corner_ptr(i,j) = longitude(i_input,1) + (0.5_esmf_kind_r8*deltalon)
     else
       lon_corner_ptr(i,j) = longitude(i,1) - (0.5_esmf_kind_r8*deltalon)
     endif

     if (j == j_input+1) then
       lat_corner_ptr(i,j) = latitude(1,j_input) +(0.5_esmf_kind_r8*deltalon)
     else
       lat_corner_ptr(i,j) = latitude(1,j) -(0.5_esmf_kind_r8*deltalon)
     endif
   enddo
 enddo

 deallocate(latitude, longitude) 
 nullify(lat_corner_ptr, lon_corner_ptr)

 !Allocate and fill in the fields on the input grid that we need to create soil type
 call ESMF_GridAddItem(input_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddItem", rc)

 call ESMF_GridGetItem(input_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       farrayPtr=mask_input_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetItem", rc)

 soil_type_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)
 
 allocate(mask_input(i_input,j_input))
 allocate(sotyp_input(i_input,j_input))
 allocate(veg_type_input(i_input,j_input)) 
 mask_input = reshape((/0,1,1,1, 0,1,1,1, 0,1,1,0/),(/i_input,j_input/))
 sotyp_input = reshape((/0.,16.,4.,16.,  0.,3.,5.,16.,   0.,3.,5.,0./),(/i_input,j_input/))
 veg_type_input = reshape((/0.,15.,5.,15.,  0.,5.,6.,15.,  0.,5.,6.,0./),(/i_input,j_input/)) 
 call ESMF_FieldScatter(soil_type_input_grid,sotyp_input,rootpet=0,rc=rc)
 
 deallocate(sotyp_input)

 !--------------------------------------------------------------------!
 !---------------- Setup Target Grid & Coordinates -------------------!
 !--------------------------------------------------------------------!
 
 target_grid = ESMF_GridCreate1PeriDim(maxIndex=(/i_target,j_target/), &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

 allocate(latitude(i_target,j_target))
 allocate(longitude(i_target,j_target))  

 ! Regional grid that fits within the input regional grid but with smaller grid cells
 deltalon = 0.5
 do i = 1, i_target
   longitude(i,:) = 91.1_esmf_kind_r8 + real((i-1),kind=esmf_kind_r8) * deltalon
 enddo

 do i = 1, j_target
   latitude(:,i) = 34.1_esmf_kind_r8 - real((i-1),kind=esmf_kind_r8) * deltalon
 enddo            

  call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     lon_ptr(i,j) = longitude(i,j)
     if (lon_ptr(i,j) > 360.0_esmf_kind_r8) lon_ptr(i,j) = lon_ptr(i,j) -360.0_esmf_kind_r8
     lat_ptr(i,j) = latitude(i,j)
   enddo
 enddo
 nullify(lat_ptr,lon_ptr)

 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_corner_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_corner_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 ! Create staggered coordinates for regional target grid
 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (i == i_target+1) then
       lon_corner_ptr(i,j) = longitude(i_input,1) + (0.5_esmf_kind_r8*deltalon)
     else
       lon_corner_ptr(i,j) = longitude(i,1) - (0.5_esmf_kind_r8*deltalon)
     endif

     if (j == j_target+1) then
       lat_corner_ptr(i,j) = latitude(1,j_input) +(0.5_esmf_kind_r8*deltalon)
     else
       lat_corner_ptr(i,j) = latitude(1,j) -(0.5_esmf_kind_r8*deltalon)
     endif
   enddo
 enddo


 deallocate(latitude, longitude)
 nullify(lat_corner_ptr, lon_corner_ptr)

 ! Create the fields for the target grid land and seamask since these would normally
 ! be created in the appropriate model_grid  subroutine
 
 call ESMF_GridAddItem(target_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddItem", rc)

 call ESMF_GridGetItem(target_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       farrayPtr=mask_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetItem", rc)

 soil_type_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 ! Create masks on the target grid and the correcte (expected) soil type on the target grid
 ! to check against what returns from interp
 
 allocate(mask_target(i_target,j_target))
 allocate(sotyp_target(i_target,j_target))
 allocate(veg_type_target(i_target,j_target))
 allocate(sotyp_correct(i_target,j_target))

 mask_target(:,1) = (/0,0,1,1,1,1,1,1/) 
 mask_target(:,2) = (/0,0,1,1,1,1,1,1/)
 mask_target(:,3) = (/0,0,1,1,1,1,1,1/)
 mask_target(:,4) = (/0,0,1,1,1,1,0,0/)
 mask_target(:,5) = (/0,0,1,1,1,1,0,0/)
   
 veg_type_target = reshape((/0., 0., 15.,15.,5., 5., 5., 5., &  
                             0., 0., 5., 5., 6., 6., 6., 6., &
                             0., 0., 5., 5., 6., 6., 6., 6., &
                             0., 0., 5., 5., 6., 6., 0., 0., &
                             0., 0., 5., 5., 6., 6., 0., 0. /),(/i_target,j_target/))
 sotyp_correct = reshape((/0., 0., 0., 0., 4., 4., 4., 4., &
                           0., 0., 3., 3., 5., 5., 5., 5., &
                           0., 0., 3., 3., 5., 5., 5., 5., &
                           0., 0., 3., 3., 5., 5., 0., 0., &
                           0., 0., 3., 3., 5., 5., 0., 0. /),(/i_target,j_target/))


 !Call the routine to unit test.
  method=ESMF_REGRIDMETHOD_NEAREST_STOD

 isrctermprocessing = 1
 
 mask_input_ptr = 0
 where (mask_input == 1) mask_input_ptr = 1
 where (nint(veg_type_input) == veg_type_landice_input) mask_input_ptr = 0

 mask_target_ptr = 0
 where (mask_target == 1) mask_target_ptr = 1
 where (nint(veg_type_target) == veg_type_landice_target) mask_target_ptr =0 

 call ESMF_FieldRegridStore(soil_type_input_grid, &
                            soil_type_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_land, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldRegridStore", rc)

 call ESMF_FieldRegrid(soil_type_input_grid, &
                       soil_type_target_grid, &
                       routehandle=regrid_land, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ,  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldRegrid", rc)


 call ESMF_FieldGather(soil_type_target_grid, sotyp_target, rootPet=0, rc=rc)

 print*,"Check results."

 if (any((abs(sotyp_target - sotyp_correct)) > EPSILON)) then
   print*,'TEST FAILED '
   print*,'ARRAY SHOULD BE:', sotyp_correct
   print*,'ARRAY FROM TEST:', sotyp_target
   stop 2
 endif

 
 print*,"OK"

 deallocate(mask_target,mask_input,sotyp_target,sotyp_correct,veg_type_target,veg_type_input)

 call ESMF_FieldRegridRelease(routehandle=regrid_land, rc=rc)
 call ESMF_FieldDestroy(soil_type_input_grid,rc=rc)
 call ESMF_FieldDestroy(soil_type_target_grid,rc=rc)
 call ESMF_GridDestroy(input_grid,rc=rc)
 call ESMF_GridDestroy(target_grid,rc=rc)
 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)
 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program esmf_regrid_bug
