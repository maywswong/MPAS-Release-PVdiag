! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_todynamics
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_dmpar
 use mpas_timer              ! MC added for tendencies
 use mpas_atmphys_constants, only: R_d,R_v,degrad

 implicit none
 private
 public:: physics_get_tend


!Interface between the physics parameterizations and the non-hydrostatic dynamical core.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
!
! subroutines in mpas_atmphys_todynamics:
! ---------------------------------------
! physics_get_tend: add and mass-weigh tendencies before being added to dynamics tendencies.
! tend_toEdges   : interpolate wind-tendencies from centers to edges of grid-cells.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added calculation of the advective tendency of the potential temperature due to horizontal
!   and vertical advection, and horizontal mixing (diffusion).
!   Laura D. Fowler (birch.mmm.ucar.edu) / 2013-11-19.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * renamed "tiedtke" with "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * modified the sourcecode to accomodate the packages "cu_kain_fritsch_in" and "cu_ntiedtke_in".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-24.
! * added the option bl_mynn for the calculation of the tendency for the cloud ice number concentration.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-11.
! * in subroutine physics_get_tend_work, added the option cu_ntiedtke in the calculation of rucuten_Edge.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-28.
! * Added support for initial tendency diagnostics package variables and timers  
!   May Wong (mwong@ucar.edu) and Manda Chasteen (chasteen@ucar.edu) / 2024-05-29      

 !
 ! Abstract interface for routine used to communicate halos of fields
 ! in a named group
 !
 abstract interface
    subroutine halo_exchange_routine(domain, halo_group, ierr)

       use mpas_derived_types, only : domain_type

       type (domain_type), intent(inout) :: domain
       character(len=*), intent(in) :: halo_group
       integer, intent(out), optional :: ierr

    end subroutine halo_exchange_routine
 end interface


 contains

 
!=================================================================================================================
 subroutine physics_get_tend( block, mesh, state, diag, tend, tend_physics, diag_physics, configs, rk_step, dynamics_substep, &
                              tend_ru_physics, tend_rtheta_physics, tend_rho_physics, exchange_halo_group )
!=================================================================================================================
   
 use mpas_atm_dimensions

!input variables:
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: configs
 integer, intent(in):: rk_step
 integer, intent(in):: dynamics_substep
 procedure (halo_exchange_routine) :: exchange_halo_group

!inout variables:
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: tend
 type(mpas_pool_type),intent(inout):: tend_physics
 type(mpas_pool_type),intent(inout):: diag_physics                  ! MW: added for ITM

 real(kind=RKIND),dimension(:,:) :: tend_ru_physics, tend_rtheta_physics, tend_rho_physics

!local variables:
 character(len=StrKIND), pointer :: config_pbl_scheme, config_convection_scheme, &
                                    config_radt_lw_scheme, config_radt_sw_scheme, &
                                    config_gwdo_scheme              ! MW: for ITM 

 logical, pointer :: config_tend       ! MC -- adding for tendencies

 integer:: i,iCell,k,n
 integer,pointer:: index_qv, index_qc, index_qr, index_qi, index_qs, index_qg
 integer,pointer:: index_ni
 integer,pointer:: nCells,nCellsSolve,nEdges,nEdgesSolve

 real(kind=RKIND),dimension(:,:),pointer:: mass          ! time level 2 rho_zz
 real(kind=RKIND),dimension(:,:),pointer:: mass_edge     ! diag rho_edge
 real(kind=RKIND),dimension(:,:),pointer:: theta_m       ! time level 1
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rqvblten,rqcblten, &
                                           rqiblten,rqsblten,rublten,rvblten
 real(kind=RKIND),dimension(:,:),pointer:: rniblten
 real(kind=RKIND),dimension(:,:),pointer:: rubldiff, rvbldiff                ! MW on ITM: GWDO; MC added                                   
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten, &
                                           rqrcuten,rqicuten,rqscuten, &
                                           rucuten,rvcuten
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw                                    
 
 real(kind=RKIND),dimension(:,:),pointer:: tend_u_phys !nick           ! MC : no longer needed for PV; remove     

 real(kind=RKIND),dimension(:,:),pointer  :: tend_theta,tend_theta_euler,tend_u
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_scalars
 real(kind=RKIND):: coeff

 real(kind=RKIND):: tem
 real(kind=RKIND),dimension(:,:),pointer:: rublten_Edge,rucuten_Edge

 real(kind=RKIND),dimension(:,:),allocatable:: theta,tend_th

! MW: for ITM
 real(kind=RKIND),dimension(:,:),pointer:: rucuten_tend, rublten_tend, rugwdo_tend
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten_tend, rthblten_tend, rthratenlw_tend, rthratensw_tend
 real(kind=RKIND),dimension(:,:),pointer:: qvcuten_tend, qvblten_tend       

!=================================================================================================================
 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
 call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)

 call mpas_pool_get_config(configs, 'config_pbl_scheme', config_pbl_scheme)
 call mpas_pool_get_config(configs, 'config_convection_scheme', config_convection_scheme)
 call mpas_pool_get_config(configs, 'config_radt_lw_scheme', config_radt_lw_scheme)
 call mpas_pool_get_config(configs, 'config_radt_sw_scheme', config_radt_sw_scheme)
 call mpas_pool_get_config(configs, 'config_gwdo_scheme', config_gwdo_scheme)            ! MW on ITM: GWDO     

 call mpas_pool_get_config(configs, 'config_tend', config_tend)                          ! MC: for tendencies

 call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
 call mpas_pool_get_array(state, 'scalars', scalars, 1)
 call mpas_pool_get_array(state, 'rho_zz', mass, 2)
 call mpas_pool_get_array(diag , 'rho_edge', mass_edge)

 call mpas_pool_get_array(diag , 'tend_u_phys', tend_u_phys) !nick    ! MC -- need to remove

 call mpas_pool_get_dimension(state, 'index_qv', index_qv)
 call mpas_pool_get_dimension(state, 'index_qc', index_qc)
 call mpas_pool_get_dimension(state, 'index_qr', index_qr)
 call mpas_pool_get_dimension(state, 'index_qi', index_qi)
 call mpas_pool_get_dimension(state, 'index_qs', index_qs)
 call mpas_pool_get_dimension(state, 'index_qg', index_qg)
 call mpas_pool_get_dimension(state, 'index_ni', index_ni)

 call mpas_pool_get_array(tend_physics, 'rublten', rublten)
 call mpas_pool_get_array(tend_physics, 'rvblten', rvblten)
 call mpas_pool_get_array(tend_physics, 'rublten_Edge', rublten_Edge)
 call mpas_pool_get_array(tend_physics, 'rthblten', rthblten)
 call mpas_pool_get_array(tend_physics, 'rqvblten', rqvblten)
 call mpas_pool_get_array(tend_physics, 'rqcblten', rqcblten)
 call mpas_pool_get_array(tend_physics, 'rqiblten', rqiblten)
 call mpas_pool_get_array(tend_physics, 'rqsblten', rqsblten)
 call mpas_pool_get_array(tend_physics, 'rniblten', rniblten)

 call mpas_pool_get_array(diag_physics, 'rubldiff', rubldiff) ! MW on ITM: GWDO contrib. to rublten   
 call mpas_pool_get_array(diag_physics, 'rvbldiff', rvbldiff) ! MW on ITM

 call mpas_pool_get_array(tend_physics, 'rucuten', rucuten)
 call mpas_pool_get_array(tend_physics, 'rvcuten', rvcuten)
 call mpas_pool_get_array(tend_physics, 'rucuten_Edge', rucuten_Edge)
 call mpas_pool_get_array(tend_physics, 'rthcuten', rthcuten)
 call mpas_pool_get_array(tend_physics, 'rqvcuten', rqvcuten)
 call mpas_pool_get_array(tend_physics, 'rqccuten', rqccuten)
 call mpas_pool_get_array(tend_physics, 'rqrcuten', rqrcuten)
 call mpas_pool_get_array(tend_physics, 'rqicuten', rqicuten)
 call mpas_pool_get_array(tend_physics, 'rqscuten', rqscuten)

 call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
 call mpas_pool_get_array(tend_physics, 'rthratensw', rthratensw)

 call mpas_pool_get_array(tend,'u'           , tend_u         )
 call mpas_pool_get_array(tend,'theta_m'     , tend_theta     )
 call mpas_pool_get_array(tend,'theta_euler' ,tend_theta_euler)
 call mpas_pool_get_array(tend,'scalars_tend',tend_scalars    )

! MW on ITM 
 call mpas_pool_get_array(diag, 'rublten_tend', rublten_tend)
 call mpas_pool_get_array(diag, 'rugwdo_tend', rugwdo_tend)
 call mpas_pool_get_array(diag, 'rucuten_tend', rucuten_tend)
 call mpas_pool_get_array(diag, 'rthblten_tend', rthblten_tend)
 call mpas_pool_get_array(diag, 'rthcuten_tend', rthcuten_tend)
 call mpas_pool_get_array(diag, 'rthratenlw_tend', rthratenlw_tend)
 call mpas_pool_get_array(diag, 'rthratensw_tend', rthratensw_tend)
 call mpas_pool_get_array(diag, 'qvblten_tend', qvblten_tend)
 call mpas_pool_get_array(diag, 'qvcuten_tend', qvcuten_tend)


!initialize the tendency for the potential temperature and all scalars due to PBL, convection,
!and longwave and shortwave radiation:
! allocate(theta(nVertLevels,nCellsSolve)  )
 allocate(tend_th(nVertLevels,nCellsSolve))
 tend_th = 0._RKIND

 tend_scalars(:,:,:) = 0._RKIND

 tend_ru_physics(:,:) = 0._RKIND
 tend_rtheta_physics(:,:) = 0._RKIND
 tend_rho_physics(:,:) = 0._RKIND       ! NB: rho tendency is not currently supplied by physics, but this
                                        !     field may be later filled with IAU or other tendencies

 ! tend_u_phys(:,:) = 0._RKIND          ! MC: this fix was needed because tend_u_phys and tend_u_diff accumulate in PV code;
                                        ! will no longer be necessary once this var is removed   

 !
 ! In case some variables are not allocated due to their associated packages,
 ! we need to make their pointers associated here to avoid triggering run-time
 ! checks when calling physics_get_tend_work
 !
 if (.not. associated(rublten))  allocate(rublten(0,0) )
 if (.not. associated(rvblten))  allocate(rvblten(0,0) )
 if (.not. associated(rthblten)) allocate(rthblten(0,0))
 if (.not. associated(rqvblten)) allocate(rqvblten(0,0))
 if (.not. associated(rqcblten)) allocate(rqcblten(0,0))
 if (.not. associated(rqiblten)) allocate(rqiblten(0,0))
 if (.not. associated(rqsblten)) allocate(rqsblten(0,0))
 if (.not. associated(rniblten)) allocate(rniblten(0,0))
 if (.not. associated(rucuten))  allocate(rucuten(0,0) )
 if (.not. associated(rvcuten))  allocate(rvcuten(0,0) )
 if (.not. associated(rthcuten)) allocate(rthcuten(0,0))
 if (.not. associated(rqvcuten)) allocate(rqvcuten(0,0))
 if (.not. associated(rqccuten)) allocate(rqccuten(0,0))
 if (.not. associated(rqicuten)) allocate(rqicuten(0,0))
 if (.not. associated(rqrcuten)) allocate(rqrcuten(0,0))
 if (.not. associated(rqscuten)) allocate(rqscuten(0,0))

! MC -- adding below for tendency variables
 if (.not. associated(rublten_tend)) allocate(rublten_tend(nVertLevels,nEdges+1))
 if (.not. associated(rugwdo_tend)) allocate(rugwdo_tend(nVertLevels,nEdges+1))
 if (.not. associated(rucuten_tend)) allocate(rucuten_tend(nVertLevels,nEdges+1))
 if (.not. associated(rthblten_tend)) allocate(rthblten_tend(nVertLevels,nCells+1))
 if (.not. associated(rthcuten_tend)) allocate(rthcuten_tend(nVertLevels,nCells+1))
 if (.not. associated(rthratenlw_tend)) allocate(rthratenlw_tend(nVertLevels,nCells+1))
 if (.not. associated(rthratensw_tend)) allocate(rthratensw_tend(nVertLevels,nCells+1))
 if (.not. associated(qvblten_tend)) allocate(qvblten_tend(nVertLevels,nCells+1))
 if (.not. associated(qvcuten_tend)) allocate(qvcuten_tend(nVertLevels,nCells+1))


 call physics_get_tend_work(block, mesh, nCells, nEdges, nCellsSolve, nEdgesSolve, &
                           rk_step, dynamics_substep, &
                           config_pbl_scheme, config_convection_scheme, config_radt_lw_scheme, config_radt_sw_scheme, &
                           config_gwdo_scheme, config_tend, & ! MW on ITM: GWDO      
                           index_qv, index_qc, index_qr, index_qi, index_qs, index_ni, &
                           rublten, rvblten, mass_edge, rublten_Edge, &
                           tend_ru_physics, &
                           rucuten, rvcuten, rucuten_Edge, &
                           ! MW on ITM 
                           rubldiff, rvbldiff, & ! GWDO        
                           rublten_tend, rucuten_tend, rugwdo_tend, &                
                           rthblten_tend, rthcuten_tend, rthratenlw_tend, rthratensw_tend,&
                           qvblten_tend, qvcuten_tend, &
                           !
                           tend_th, tend_scalars, mass, rthblten, rqvblten, rqcblten, rqiblten, rqsblten, rniblten, &
                           rthcuten, rqvcuten, rqccuten, rqrcuten, rqicuten, rqscuten, &
                           rthratenlw, rthratensw, &
                           tend_u_phys, &           ! MC: tend_u_phys is still here but needs to eventually be removed
                           theta_m, scalars, &
                           tend_rtheta_physics, &
                           tend_theta_euler, &
                           exchange_halo_group &
                           )

 !
 ! Clean up any pointers that were allocated with zero size before the call to
 ! physics_get_tend_work
 !
 if (size(rublten)  == 0) deallocate(rublten )
 if (size(rvblten)  == 0) deallocate(rvblten )
 if (size(rthblten) == 0) deallocate(rthblten)
 if (size(rqvblten) == 0) deallocate(rqvblten)
 if (size(rqcblten) == 0) deallocate(rqcblten)
 if (size(rqiblten) == 0) deallocate(rqiblten)
 if (size(rqsblten) == 0) deallocate(rqsblten)
 if (size(rniblten) == 0) deallocate(rniblten)
 if (size(rucuten)  == 0) deallocate(rucuten )
 if (size(rvcuten)  == 0) deallocate(rvcuten )
 if (size(rthcuten) == 0) deallocate(rthcuten)
 if (size(rqvcuten) == 0) deallocate(rqvcuten)
 if (size(rqccuten) == 0) deallocate(rqccuten)
 if (size(rqicuten) == 0) deallocate(rqicuten)
 if (size(rqrcuten) == 0) deallocate(rqrcuten)
 if (size(rqscuten) == 0) deallocate(rqscuten)

! MC -- adding below for tendency variables. only deallocate if allocated above     
 if (.not. config_tend) then
    deallocate(rublten_tend)
    deallocate(rugwdo_tend)
    deallocate(rucuten_tend)   
    deallocate(rthblten_tend) 
    deallocate(rthcuten_tend)    
    deallocate(rthratenlw_tend)   
    deallocate(rthratensw_tend)     
    deallocate(qvblten_tend) 
    deallocate(qvcuten_tend) 
 end if

! deallocate(theta)
 deallocate(tend_th)

! if(rk_step .eq. 3) then
!    call mpas_log_write('')
!    call mpas_log_write('--- enter subroutine physics_get_tend:')
!    call mpas_log_write('max rthblten   = $r',realArgs=(/maxval(rthblten(:,1:nCellsSolve))/))
!    call mpas_log_write('min rthblten   = $r',realArgs=(/minval(rthblten(:,1:nCellsSolve))/))
!    call mpas_log_write('max rthcuten   = $r',realArgs=(/maxval(rthcuten(:,1:nCellsSolve))/))
!    call mpas_log_write('min rthcuten   = $r',realArgs=(/minval(rthcuten(:,1:nCellsSolve))/))
!    call mpas_log_write('max rthratenlw = $r',realArgs=(/maxval(rthratenlw(:,1:nCellsSolve))/))
!    call mpas_log_write('min rthratenlw = $r',realArgs=(/minval(rthratenlw(:,1:nCellsSolve))/))
!    call mpas_log_write('max rthratensw = $r',realArgs=(/maxval(rthratensw(:,1:nCellsSolve))/))
!    call mpas_log_write('min rthratensw = $r',realArgs=(/minval(rthratensw(:,1:nCellsSolve))/))
!    call mpas_log_write('--- end subroutine physics_get_tend')
!    call mpas_log_write('')
! endif

 end subroutine physics_get_tend

 !==================================================================================================
 subroutine physics_get_tend_work(block, mesh, nCells, nEdges, nCellsSolve, nEdgesSolve, &
                                 rk_step, dynamics_substep, &
                                 config_pbl_scheme, config_convection_scheme, config_radt_lw_scheme, config_radt_sw_scheme, &
                                 config_gwdo_scheme, config_tend, & ! MW on ITM 
                                 index_qv, index_qc, index_qr, index_qi, index_qs, index_ni, &
                                 rublten, rvblten, mass_edge, rublten_Edge, tend_u, &
                                 rucuten, rvcuten, rucuten_Edge, &
                                 ! MW on ITM   
                                 rubldiff, rvbldiff, & ! GWDO  
                                 rublten_tend, rucuten_tend, rugwdo_tend, &              
                                 rthblten_tend, rthcuten_tend, rthratenlw_tend, rthratensw_tend, &
                                 qvblten_tend, qvcuten_tend, &           
                                 !
                                 tend_th, tend_scalars, mass, rthblten, rqvblten, rqcblten, rqiblten, rqsblten, rniblten, &
                                 rthcuten, rqvcuten, rqccuten, rqrcuten, rqicuten, rqscuten, &
                                 rthratenlw, rthratensw, &
                                 tend_u_phys, &
                                 theta_m, scalars, tend_theta, tend_theta_euler, &
                                 exchange_halo_group &
                                )
!==================================================================================================

    use mpas_atm_dimensions
 
    implicit none

    type(block_type), intent(in) :: block
    type(mpas_pool_type), intent(in) :: mesh
    integer, intent(in) :: nCells, nEdges, nCellsSolve, nEdgesSolve
    integer, intent(in) :: rk_step, dynamics_substep
    character(len=StrKIND), intent(in) :: config_pbl_scheme
    character(len=StrKIND), intent(in) :: config_convection_scheme
    character(len=StrKIND), intent(in) :: config_radt_lw_scheme
    character(len=StrKIND), intent(in) :: config_radt_sw_scheme
    character(len=StrKIND), intent(in) :: config_gwdo_scheme ! MW on ITM: GWDO         
    logical, intent(in) :: config_tend  ! MC added for ITM                       

    integer, intent(in) :: index_qv, index_qc, index_qr, index_qi, index_qs, index_ni
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rublten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rvblten
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: mass_edge
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: rublten_Edge
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: tend_u
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rucuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rvcuten
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: rucuten_Edge
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: tend_th
    real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: tend_scalars
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: mass
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rthblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqvblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqcblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqiblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqsblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rniblten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rthcuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqvcuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqccuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqrcuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqicuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rqscuten
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rthratenlw
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rthratensw
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: tend_u_phys           ! MC - eventually remove
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: theta_m
    real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(in) :: scalars
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: tend_theta
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: tend_theta_euler
    procedure (halo_exchange_routine) :: exchange_halo_group

    ! MW on ITM  
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rubldiff
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rvbldiff
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: rucuten_tend
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: rublten_tend
    real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(inout) :: rugwdo_tend
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: rthcuten_tend
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: rthblten_tend
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: rthratenlw_tend
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: rthratensw_tend
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: qvblten_tend
    real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: qvcuten_tend

    real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: rubldiff_Edge !local 
    ! end MW

    integer :: i, k
    real (kind=RKIND) :: coeff

    !MW on ITM: diagnosing tendencies due to GWDO scheme: 
    call mpas_timer_start('Tendency and PV diagnostics')                   ! MC - added timer for dedicated diag calculation 
    if (config_tend .and. (config_gwdo_scheme .ne. 'off')) then            ! MC - added config_tend flag here. otherwise will break if config_tend is off due to halo group specs
       if (rk_step == 1 .and. dynamics_substep == 1) then
          ! note: uncoupled tendency doesn't change over dynamic full time step/
          !       physics is only called once per timestep
          ! MC adding halo group
          call exchange_halo_group(block % domain, 'physics:bldiff')
          call tend_toEdges(block,mesh,rubldiff,rvbldiff,rubldiff_Edge)
       end if

       do i = 1, nEdgesSolve
       do k = 1, nVertLevels
          rugwdo_tend(k,i) = rubldiff_Edge(k,i)*mass_edge(k,i)
       end do
       end do
    end if
    call mpas_timer_stop('Tendency and PV diagnostics')  
    ! end GWDO

    !add coupled tendencies due to PBL processes:
    if (config_pbl_scheme .ne. 'off') then
       if (rk_step == 1 .and. dynamics_substep == 1) then
          call exchange_halo_group(block % domain, 'physics:blten')
          call tend_toEdges(block,mesh,rublten,rvblten,rublten_Edge)

          !MGD for PV budget? should a similar line be in the cumulus section below?
          tend_u_phys(1:nVertLevels,1:nEdges) = rublten_Edge(1:nVertLevels,1:nEdges)         ! MC eventually needs to be removed 
       end if

       do i = 1, nEdgesSolve
       do k  = 1, nVertLevels
          tend_u(k,i)=tend_u(k,i)+rublten_Edge(k,i)*mass_edge(k,i)
          rublten_tend(k,i) = rublten_Edge(k,i)*mass_edge(k,i)  ! MW on ITM         
       enddo
       enddo

       do i = 1, nCellsSolve
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthblten(k,i)*mass(k,i)
          tend_scalars(index_qv,k,i) = tend_scalars(index_qv,k,i) + rqvblten(k,i)*mass(k,i)
          tend_scalars(index_qc,k,i) = tend_scalars(index_qc,k,i) + rqcblten(k,i)*mass(k,i)
          tend_scalars(index_qi,k,i) = tend_scalars(index_qi,k,i) + rqiblten(k,i)*mass(k,i)

          rthblten_tend(k,i) = rthblten(k,i)*mass(k,i)  ! MW on ITM
          qvblten_tend(k,i) = rqvblten(k,i)  ! MW on ITM; decoupled using mass at time t in advance_scalar_mono  
       enddo
       enddo

       pbl_select: select case (trim(config_pbl_scheme))

          case("bl_mynn")

             do i = 1, nCellsSolve
             do k = 1, nVertLevels
                tend_scalars(index_qs,k,i) = tend_scalars(index_qs,k,i) + rqsblten(k,i)*mass(k,i)
                tend_scalars(index_ni,k,i) = tend_scalars(index_ni,k,i) + rniblten(k,i)*mass(k,i)
             enddo
             enddo
   
          case default         

        end select pbl_select
    endif

    !add coupled tendencies due to convection:
    if (config_convection_scheme .ne. 'off') then

       do i = 1, nCellsSolve
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthcuten(k,i)*mass(k,i)
          tend_scalars(index_qv,k,i) = tend_scalars(index_qv,k,i) + rqvcuten(k,i)*mass(k,i)
          tend_scalars(index_qc,k,i) = tend_scalars(index_qc,k,i) + rqccuten(k,i)*mass(k,i)
          tend_scalars(index_qi,k,i) = tend_scalars(index_qi,k,i) + rqicuten(k,i)*mass(k,i)

          rthcuten_tend(k,i) = rthcuten(k,i)*mass(k,i)  ! MW on ITM   
          qvcuten_tend(k,i) = rqvcuten(k,i)             ! MW on ITM; decoupled using mass at time t in advance_scalar_mono
       enddo
       enddo

        convection_select: select case(config_convection_scheme)

           case('cu_kain_fritsch')
              do i = 1, nCellsSolve
              do k = 1, nVertLevels
                 tend_scalars(index_qr,k,i) = tend_scalars(index_qr,k,i) + rqrcuten(k,i)*mass(k,i)
                 tend_scalars(index_qs,k,i) = tend_scalars(index_qs,k,i) + rqscuten(k,i)*mass(k,i)
              enddo
              enddo
    
           case('cu_tiedtke','cu_ntiedtke')
              if (rk_step == 1 .and. dynamics_substep == 1) then
                 call exchange_halo_group(block % domain, 'physics:cuten')
                 call tend_toEdges(block,mesh,rucuten,rvcuten,rucuten_Edge)
                 
                 tend_u_phys(1:nVertLevels,1:nEdges) = tend_u_phys(1:nVertLevels,1:nEdges) &        ! MC -- needs to be removed 
                                                     + rucuten_Edge(1:nVertLevels,1:nEdges)
              end if
              do i = 1, nEdgesSolve
              do k  = 1, nVertLevels
                 tend_u(k,i)=tend_u(k,i)+rucuten_Edge(k,i)*mass_edge(k,i)
                 rucuten_tend(k,i) = rucuten_Edge(k,i)*mass_edge(k,i)  ! MW on ITM    
              enddo
              enddo

            case default
        end select convection_select
    endif

    !add coupled tendencies due to longwave radiation:
    if (config_radt_lw_scheme .ne. 'off') then
       do i = 1, nCellsSolve
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthratenlw(k,i)*mass(k,i)

          rthratenlw_tend(k,i) = rthratenlw(k,i)*mass(k,i)            ! MW on ITM   
       enddo
       enddo
    endif
    
    !add coupled tendencies due to shortwave radiation:
    if (config_radt_sw_scheme .ne. 'off') then
       do i = 1, nCellsSolve
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthratensw(k,i)*mass(k,i)

          rthratensw_tend(k,i) = rthratensw(k,i)*mass(k,i)            ! MW on ITM     
       enddo
       enddo
    endif

    !if non-hydrostatic core, convert the tendency for the potential temperature to a
    !tendency for the modified potential temperature:
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       coeff = (1. + R_v/R_d * scalars(index_qv,k,i))
       tend_th(k,i) = coeff * tend_th(k,i) + R_v/R_d * theta_m(k,i) * tend_scalars(index_qv,k,i) / coeff

       ! MW on ITM        
       rthblten_tend(k,i) = coeff* rthblten_tend(k,i) + R_v/R_d * theta_m(k,i) * mass(k,i)*qvblten_tend(k,i) / coeff
       rthcuten_tend(k,i) = coeff* rthcuten_tend(k,i) + R_v/R_d * theta_m(k,i) * mass(k,i)*qvcuten_tend(k,i) / coeff
       rthratenlw_tend(k,i) = coeff* rthratenlw_tend(k,i)
       rthratensw_tend(k,i) = coeff* rthratensw_tend(k,i)
       ! MW:end

       tend_theta(k,i) = tend_theta(k,i) + tend_th(k,i)
    enddo
    enddo

 end subroutine physics_get_tend_work

!=================================================================================================================
 subroutine tend_toEdges(block,mesh,Ux_tend,Uy_tend,U_tend)
!=================================================================================================================

 use mpas_atm_dimensions

!input arguments:
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: mesh
 real(kind=RKIND),intent(in),dimension(:,:),target:: Ux_tend,Uy_tend 

!output arguments:
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
 integer:: iCell,iEdge,k,j
 integer:: cell1, cell2
 integer,pointer:: nCells,nCellsSolve,nEdges
 integer,dimension(:,:),pointer:: cellsOnEdge

 real(kind=RKIND), dimension(:,:), pointer :: east, north, edgeNormalVectors

 
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)

 call mpas_pool_get_array(mesh, 'east', east)
 call mpas_pool_get_array(mesh, 'north', north)
 call mpas_pool_get_array(mesh, 'edgeNormalVectors', edgeNormalVectors)

 call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

 do iEdge = 1, nEdges
    cell1 = cellsOnEdge(1,iEdge)
    cell2 = cellsOnEdge(2,iEdge)
                           
    U_tend(:,iEdge) =  Ux_tend(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   &
                                              +  edgeNormalVectors(2,iEdge) * east(2,cell1)   &
                                              +  edgeNormalVectors(3,iEdge) * east(3,cell1))  &
                     + Uy_tend(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)   &
                                              +  edgeNormalVectors(2,iEdge) * north(2,cell1)   &
                                              +  edgeNormalVectors(3,iEdge) * north(3,cell1))  &
                     + Ux_tend(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   &
                                              +  edgeNormalVectors(2,iEdge) * east(2,cell2)   &
                                              +  edgeNormalVectors(3,iEdge) * east(3,cell2))  &
                     + Uy_tend(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)   &
                                              +  edgeNormalVectors(2,iEdge) * north(2,cell2)   &
                                              +  edgeNormalVectors(3,iEdge) * north(3,cell2))
 end do
 
 end subroutine tend_toEdges

!=================================================================================================================
 end module mpas_atmphys_todynamics
!=================================================================================================================
