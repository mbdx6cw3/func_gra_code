Module initialise

Contains

Subroutine read_topol

	use global

	open(unit=20,file='parameters.dat',status='unknown')
      	read(20,*) n_atm_type
      	allocate(mass(n_atm_type),sigma(n_atm_type),epsilon(n_atm_type))

      	do i_atm_type = 1, n_atm_type

        	read(20,*) mass(i_atm_type), sigma(i_atm_type), epsilon(i_atm_type)

      	end do

      	read(20,*) n_bond_type
      	allocate(bond_r(n_bond_type),bond_kr(n_bond_type))
      	do i_bond_type = 1, n_bond_type

       		read(20,*) bond_r(i_bond_type), bond_kr(i_bond_type)

     	end do

      	read(20,*) n_ang_type
      	allocate(ang_theta(n_ang_type),ang_k(n_ang_type),ang_ub(n_ang_type),ang_kub(n_ang_type))
      	do i_ang_type = 1, n_ang_type

       		read(20,*) ang_theta(i_ang_type), ang_k(i_ang_type), ang_ub(i_ang_type), ang_kub(i_ang_type)

      	end do

      	read(20,*) n_dih_type
      	allocate(dih_phi(n_dih_type),dih_k(n_dih_type),dih_mult(n_dih_type))
      	do i_dih_type = 1, n_dih_type

        	read(20,*) dih_phi(i_dih_type), dih_k(i_dih_type), dih_mult(i_dih_type)

      	end do  
      	close(unit=20)

End subroutine read_topol


Subroutine read_gra

	use global
        implicit none
      	100 format (a13, a2, i5, f8.3, f8.3, f8.3)
      	open(unit=20,file='graphene.gro')
      	read(20,*) A
      	read(20,*) n_atm
      	do i_atm = 1, n_atm

      		read(20,100) nam, atm_type(i_atm), index, x(i_atm), y(i_atm), z(i_atm)

     	end do
      	read(20,*) lx, ly, lz
      	close(unit=20)

End subroutine read_gra


Subroutine read_func(n_func_type,i_func,file_name,max_atm,type,x,y,z)

	implicit none
	integer, parameter :: dp = SELECTED_REAL_KIND(15,99)

	real(kind=dp) :: x(n_func_type,max_atm), y(n_func_type,max_atm), z(n_func_type,max_atm)
	character(len=2) :: type(n_func_type,max_atm)
	character(len=13) :: A, nam
	character(len=14) :: file_name
	integer :: index, i_atm, max_atm, n_atm, n_func_type, i_func

	100 format (a13, a2, i5, f8.3, f8.3, f8.3)
	open(unit=20,file=file_name)
	read(20,*) A
	read(20,*) n_atm
	do i_atm = 1, n_atm

        	read(20,100) nam, type(i_func,i_atm), index, x(i_func,i_atm), y(i_func,i_atm), z(i_func,i_atm)

	end do
      	close(unit=20)

End subroutine read_func


Subroutine bond_list
      	use global

!       --- remove undercoordinated carbons
	n_gra = 0
      	do i_atm = 1, n_atm

        	n_bond = 0

        	do j_atm = 1, n_atm

          		if (i_atm.ne.j_atm) then

            			r_x = x(i_atm) - x(j_atm)
            			r_y = y(i_atm) - y(j_atm)
            			r_z = z(i_atm) - z(j_atm)
            			r_x = r_x - lx * nint(r_x/lx)
           			r_y = r_y - ly * nint(r_y/ly) 
            			r_z = r_z - lz * nint(r_z/lz)
            			r_ij = sqrt(r_x*r_x + r_y*r_y + r_z*r_z)

!         			--- define neighbour atom distance
            			if (r_ij.le.1.43d-1) then

        	      			n_bond = n_bond + 1
   	
 			        end if

          		end if

        	end do
        
!               --- write new carbon atom array
!               --- remove atoms with only one bond
      		if (n_bond.eq.1) then

          		cycle

!               --- keep this atom because it has more than one bond
        	else

          		n_gra = n_gra + 1
          		gra_x(n_gra) = x(i_atm)
          		gra_y(n_gra) = y(i_atm)
          		gra_z(n_gra) = z(i_atm)
          		gra_type(n_gra) = atm_type(i_atm)

        	end if

      	end do

      	allocate(n_gra_bond(n_atm))

!       --- write bond lists for graphene carbon atoms
  	do i_gra = 1, n_gra

        	n_bond = 0

        	do j_gra = 1, n_gra

        	  	gro_bond_list(i_gra,j_gra) = .false.
	
          		if (i_gra.ne.j_gra) then

            			r_x = gra_x(i_gra) - gra_x(j_gra)
            			r_y = gra_y(i_gra) - gra_y(j_gra)
            			r_z = gra_z(i_gra) - gra_z(j_gra)
            			r_x = r_x - lx * nint(r_x/lx)
            			r_y = r_y - ly * nint(r_y/ly)
            			r_z = r_z - lz * nint(r_z/lz)
            			r_ij = sqrt(r_x*r_x + r_y*r_y + r_z*r_z)

!                               --- define neighbour atom distance
            			if (r_ij.le.1.43d-1) then

              				n_bond = n_bond + 1
              				gro_bond_list(i_gra,j_gra) = .true.
              				gro_bond_type(i_gra,j_gra) = 1

            			end if

          		end if

        	end do

        	n_gra_bond(i_gra) = n_bond

!		--- count number of edge sites
		if (n_gra_bond(i_gra).eq.2) then
			n_carboxy = n_carboxy + 1
		end if

      	end do

!	--- alert user if n_carboxy has had to be reduced
	if (max_carboxy.gt.n_carboxy) then
		max_carboxy = n_carboxy
		write(*,*) "Number of carboxylic acid groups requested is not possible"
		write(*,*) "Setting number of carboxylic acid groups to:", n_carboxy
	end if

!       --- copy graphene coordinates array to new GO coordinate array
      	n_gro = 0
      	do i_gra = 1, n_gra

        	n_gro = n_gro + 1
	        gro_x(n_gro) = gra_x(i_gra)
       		gro_y(n_gro) = gra_y(i_gra)
        	gro_z(n_gro) = gra_z(i_gra)
        	gro_type(n_gro) = gra_type(i_gra)
        	charge(n_gro) = 0.0d0

      	end do

End subroutine bond_list


Subroutine rand_init

	use global
      	call system_clock(count)
      	call random_seed(size=n_seed)
      	allocate (seed(n_seed))
      	do i_seed = 1, n_seed

	        seed(i_seed) = lcg(count)

	end do

	call random_seed(put=seed)

	Contains

      		function lcg(s)
        	use iso_fortran_env, only: int64
        	integer :: lcg
       	 	integer(int64) :: s
        	if (s.eq.0) then

        		s = 104729

     		else
	
        	  	s = mod(s,4294967296_int64)

        	end if
        	s = mod(s * 279470273_int64, 4294967291_int64)
        	lcg = int(mod(s,int(huge(0),int64)), kind(0))
        	end function lcg

End subroutine rand_init

End module initialise
