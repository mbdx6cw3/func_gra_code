      MODULE global
!     -----------------------------------------------------------------------
!     This module provides the global variables for the functionalise_GO.f90
!     program.
!     -----------------------------------------------------------------------

      use iso_fortran_env, only: int64
      implicit none
      integer, parameter :: dp = SELECTED_REAL_KIND(15,99)

!     global variables
!     -----------------------------------------------------------------------
      integer :: index
      integer :: i_atm, j_atm, n_atm, i_gra, j_gra, n_gra, i_gro, j_gro, k_gro, l_gro, n_gro
      real(kind=dp) :: r_x, r_y, r_z, lx, ly, lz, r_ij
      real(kind=dp), allocatable, dimension (:) :: x, y, z
      real(kind=dp), allocatable, dimension (:) :: gra_x, gra_y, gra_z, gro_x, gro_y, gro_z, charge
      character(len=13) :: A, nam
      character(len=14) :: file_name
      character(len=2), allocatable, dimension (:) :: atm_type, gra_type, gro_type
      integer :: i_bond, i_loop, i_group
      integer, allocatable, dimension (:) :: n_gra_bond
      integer, parameter :: max_bond = 4
      integer :: bond(max_bond)
      real(kind=dp) :: bond_x, bond_y, bond_z, vec_x, vec_y, vec_z, vec_x0, vec_y0, vec_z0
      real(kind=dp) :: vec_x1, vec_y1, vec_z1
      real(kind=dp) :: theta, func_x_1, func_y_1, func_z_1, rads, func_xyz, func_x_2, func_y_2, func_z_2
      real(kind=dp), parameter :: pi = 3.141592654
      real(kind=dp) :: z_switch
      logical :: periodic

!     random number generation
!     -----------------------------------------------------------------------
      integer :: n_seed, i_seed, i_ran
      integer, allocatable, dimension (:) :: seed
      integer(int64) :: count
      real(kind=dp) :: ran

!     parameters fo intermolecular potential
!     -----------------------------------------------------------------------
      integer :: i_atm_type, i_bond_type, i_ang_type, i_dih_type, itp_ref, itp_count
      integer :: n_atm_type, n_bond_type, n_ang_type, n_dih_type
      integer :: n_bond, bond_count, angle_count, dihedral_count
      real(kind=dp), allocatable, dimension (:) :: mass, sigma, epsilon
      real(kind=dp), allocatable, dimension (:) :: bond_kr, bond_r
      real(kind=dp), allocatable, dimension (:) :: ang_theta, ang_k, ang_ub, ang_kub
      real(kind=dp), allocatable, dimension (:) :: dih_phi, dih_k
      integer, allocatable, dimension (:) :: dih_mult
      real(kind=dp), allocatable, dimension (:,:) :: gro_bond_type
      logical, allocatable, dimension (:,:) :: gro_bond_list
      real(kind=dp) :: qtot

      integer, parameter :: n_func_type = 4, max_func_loop = 10000
      integer, parameter :: max_func_atm = 4
      real(kind=dp) :: O_wp, gra_pos(3), edge_rat, plane_rat
      integer :: max_n_gro, n_func_atm, i_func, i_func_atm, n_func, func_ref
      integer :: n_hydroxy, n_epoxy
      integer :: max_carboxy, n_carboxy, n_carboxylate, i_carboxy, n_hydrogen
      logical, allocatable, dimension (:) :: carb_func	
      integer :: n_alcoholate, n_up_alcoholate, n_down_alcoholate
      real(kind=dp) :: mass_O, mass_tot, mass_edge
      character(len=2) :: func_type(n_func_type,max_func_atm)
      real(kind=dp) :: func_x(n_func_type,max_func_atm), func_y(n_func_type,max_func_atm), func_z(n_func_type,max_func_atm)
      character(len=14) :: func_name(n_func_type)
      real(kind=dp) :: edge_dp, plane_dp, dp_rat
      logical :: patchy, new_patch
      integer :: max_patch, patch_size

!     input deck
!     --------------------------------------------------------------------------
      namelist /input_deck/ n_atm,func_name,O_wp,max_carboxy,plane_rat,edge_dp,plane_dp,max_n_gro,periodic,patchy,max_patch
!     ------------------------------------------------------------------------
      END MODULE global
