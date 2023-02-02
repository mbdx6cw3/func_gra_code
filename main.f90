!     -------------------------------------------------------------------------------
!     Main program for Func_GO.
!     Version 2.0
!     Christopher D. Williams, Feb 2017      
!     -------------------------------------------------------------------------------
      program main
      use global
      use initialise
      use functionalise
      use finalise

      implicit none

!     --- read input file "input.dat"
      open(unit = 10, file = 'input.dat', status = 'old')
        read(unit = 10, nml = input_deck)
      close(unit = 10)

!     --- allocate arrays
      allocate(atm_type(n_atm),x(n_atm),y(n_atm),z(n_atm))
      allocate(gra_type(n_atm),gra_x(n_atm),gra_y(n_atm),gra_z(n_atm))
      allocate(gro_x(max_n_gro),gro_y(max_n_gro),gro_z(max_n_gro),gro_type(max_n_gro),charge(max_n_gro))
      allocate(gro_bond_list(max_n_gro,max_n_gro),gro_bond_type(max_n_gro,max_n_gro))
      allocate(carb_func(max_n_gro))
      carb_func(:) = .false.

!     --- read in force field parameters for the GO sheet
      call read_topol

!     --- read in coordinates for the graphene sheet
      call read_gra

!     --- read functionality reference coordinate files
      do i_func = 1, n_func_type
        call read_func(n_func_type,i_func,func_name(i_func),max_func_atm,func_type,func_x,func_y,func_z)
      end do

!     --- write graphene atom arrays and bond lists
      call bond_list

!     --- zero functionality counters
      n_hydrogen = 0
      n_hydroxy = 0
      n_epoxy = 0
      n_carboxy = 0
      n_carboxylate = 0
      n_alcoholate = 0
      n_up_alcoholate = 0
      n_down_alcoholate = 0

!     --- total number of functional groups to be added is defined by target wt% of O
      mass_O = 0.0d0
      mass_tot = real(n_gra) * 1.20110d1
      mass_edge = 0.0d0

!     --- initial random number generator
      call rand_init

!     --- if the sheet is not periodic then functionalise the edge
      if (periodic.eqv..false.) then
        call func_edge
      end if      

!     --- functionalise basal plane until correct oxygen content is reached
      call func_plane

!     --- write topology file for the GO sheet
      call write_itp

!     --- write coordinate file for the GO sheet
      call write_gro

!     --- write output file
      call write_output

!     --------------------------------------------------------------------------------------
      end program main

