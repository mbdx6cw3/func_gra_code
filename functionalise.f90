Module functionalise

Contains

Subroutine func_edge

      	use global

!       --- loop over basal carbons and functionalise
        do i_carboxy = 1, max_carboxy

!		--- add carboxylic acid / carboxylate
		func_ref = 1
		n_bond = 0

!               --- keep randomly picking carbons at the edge of the flake
                do i_ran = 1, 10000

!                       --- choose one carbon at random
                        call random_number(ran)
			i_gra = ceiling(ran * n_gra)

!                       --- check that this carbon has two bonds (on the edge)
                        if (n_gra_bond(i_gra).ne.2) cycle

!			--- check that this carbon hasn't already been functionalised
			if (carb_func(i_gra).eqv..true.) cycle
			
			exit

		end do

!		--- this carbon has been functionalised
		carb_func(i_gra) = .true.

!               --- search for atoms bonded to this functional site
                do j_gra = 1, n_gra

                	if (gro_bond_list(i_gra,j_gra).eqv..true.) then
                        	n_bond = n_bond + 1
                                bond(n_bond) = j_gra
                        end if
                        
			if(n_bond.eq.n_gra_bond(i_gra)) exit

                end do


!               --- mid-point of the two atoms bonded to the functional site
                bond_x = (gra_x(bond(1)) + gra_x(bond(2))) / 2.0d0
!               --- temporary fix for periodic sheets
                if ((gra_x(i_gra))-bond_x.gt.1.0) then
                	bond_x = bond_x + (lx/2.0d0)
                else if ((gra_x(i_gra)-bond_x).lt.-1.0) then
                	bond_x = bond_x - (lx/2.0d0)
                end if

                bond_y = (gra_y(bond(1)) + gra_y(bond(2))) / 2.0d0
!               --- temporary fix for periodic sheets
                if ((gra_y(i_gra))-bond_y.gt.1.0) then
                	bond_y = bond_y + (ly/2.0d0)
                else if ((gra_y(i_gra)-bond_y).lt.-1.0) then
                	bond_y = bond_y - (ly/2.0d0)
                end if

                bond_z = (gra_z(bond(1)) + gra_z(bond(2))) / 2.0d0
!               --- temporary fix for periodic sheets
                if ((gra_z(i_gra))-bond_z.gt.1.0) then
                	bond_z = bond_z + (lz/2.0d0)
                else if ((gra_z(i_gra)-bond_z).lt.-1.0) then
                	bond_z = bond_z - (lz/2.0d0)
                end if

!               --- reference unit vector - from the mid-point to the functional site
                call vector(gra_x(i_gra),gra_y(i_gra),gra_z(i_gra),bond_x,bond_y,bond_z,vec_x0,vec_y0,vec_z0)

!               --- angle of the vector relative to the y axis (cosine rule)
                rads = acos(vec_y0/sqrt((vec_x0*vec_x0) + (vec_y0*vec_y0)))
                if (vec_x0.gt.0.0d0) then
                	theta = (2.0d0 * pi) - rads
                else
                        theta = rads
                end if

!               --- use the calculated vector for the first atom   
                func_x_2 = func_x(func_ref,1) * func_x(func_ref,1)
                func_y_2 = func_y(func_ref,1) * func_y(func_ref,1)
                func_z_2 = func_z(func_ref,1) * func_z(func_ref,1)
                func_xyz = sqrt(func_x_2 + func_y_2 + func_z_2)

!               --- add new atom
                n_gro = n_gro + 1
                gro_x(n_gro) = gra_x(i_gra) + func_xyz * vec_x0
                gro_y(n_gro) = gra_y(i_gra) + func_xyz * vec_y0
                gro_z(n_gro) = gra_z(i_gra) + func_xyz * vec_z0
                gro_type(n_gro) = func_type(func_ref,1)

                n_func_atm = 3
                mass_edge = mass_edge + 44.0098d0

!               --- add remaining atoms in the functional group
                do i_func_atm = 2, n_func_atm

                        call vector(func_x(func_ref,i_func_atm),func_y(func_ref,i_func_atm),&
                        func_z(func_ref,i_func_atm),func_x(func_ref,1),func_y(func_ref,1),&
                        func_z(func_ref,1),vec_x,vec_y,vec_z)

!                       --- rotate orginal vector by theta using rotation matrix
                        call rotate_z(vec_x,vec_y,vec_z,theta,vec_x1,vec_y1,vec_z1)

!                       --- calculate new coordinates
                        func_x_1 = func_x(func_ref,1) - func_x(func_ref,i_func_atm)
                        func_y_1 = func_y(func_ref,1) - func_y(func_ref,i_func_atm)
                        func_z_1 = func_z(func_ref,1) - func_z(func_ref,i_func_atm)
                        func_x_2 = func_x_1 * func_x_1
                        func_y_2 = func_y_1 * func_y_1
                        func_z_2 = func_z_1 * func_z_1
                        func_xyz = sqrt(func_x_2 + func_y_2 + func_z_2)

                	n_gro = n_gro + 1
                	gro_x(n_gro) = gro_x(n_gro-(i_func_atm-1)) + func_xyz * vec_x1
                	gro_y(n_gro) = gro_y(n_gro-(i_func_atm-1)) + func_xyz * vec_y1
                	gro_z(n_gro) = gro_z(n_gro-(i_func_atm-1)) + func_xyz * vec_z1
                	gro_type(n_gro) = func_type(func_ref,i_func_atm)

                end do

!		--- write CA - CC bond
                gro_bond_list(i_gra,n_gro-2) = .true.
                gro_bond_type(i_gra,n_gro-2) = 3
                gro_bond_list(n_gro-2,i_gra) = .true.

!		--- write CC - OC bond
                gro_bond_list(n_gro-2,n_gro-1) = .true.
                gro_bond_type(n_gro-2,n_gro-1) = 10
                gro_bond_list(n_gro-1,n_gro-2) = .true.
                gro_bond_type(n_gro-1,n_gro-2) = 10

!		--- write CC - OH bond
                gro_bond_list(n_gro-2,n_gro) = .true.
                gro_bond_type(n_gro-2,n_gro) = 11
                gro_bond_list(n_gro,n_gro-2) = .true.
                gro_bond_type(n_gro,n_gro-2) = 11

!               --- calculate fraction of carboxylic acids deprotonated
                dp_rat = real(n_carboxylate) / real(i_carboxy)

!               --- add carboxylate if fraction deprotonated is less than target ratio
                if (dp_rat.lt.edge_dp) then

!               	--- carboxylate counter
                        n_carboxylate = n_carboxylate + 1
			charge(i_gra) = -0.368d0
                        charge(n_gro-2) = 0.710d0
                        charge(n_gro-1) = -0.671d0
                        charge(n_gro) = -0.671d0

!               --- add proton
                else
			mass_edge = mass_edge + 1.0076d0
			i_func_atm = 4

                        call vector(func_x(func_ref,i_func_atm),func_y(func_ref,i_func_atm),&
                        func_z(func_ref,i_func_atm),func_x(func_ref,1),func_y(func_ref,1),&
                        func_z(func_ref,1),vec_x,vec_y,vec_z)

!                       --- rotate orginal vector by theta using rotation matrix
                        call rotate_z(vec_x,vec_y,vec_z,theta,vec_x1,vec_y1,vec_z1)

!                       --- calculate new coordinates
                        func_x_1 = func_x(func_ref,1) - func_x(func_ref,i_func_atm)
                        func_y_1 = func_y(func_ref,1) - func_y(func_ref,i_func_atm)
                        func_z_1 = func_z(func_ref,1) - func_z(func_ref,i_func_atm)
                        func_x_2 = func_x_1 * func_x_1
                        func_y_2 = func_y_1 * func_y_1
                        func_z_2 = func_z_1 * func_z_1
                        func_xyz = sqrt(func_x_2 + func_y_2 + func_z_2)

!			--- add proton
			n_gro = n_gro + 1
                        gro_x(n_gro) = gro_x(n_gro-(i_func_atm-1)) + func_xyz * vec_x1
                        gro_y(n_gro) = gro_y(n_gro-(i_func_atm-1)) + func_xyz * vec_y1
                        gro_z(n_gro) = gro_z(n_gro-(i_func_atm-1)) + func_xyz * vec_z1
                        gro_type(n_gro) = func_type(func_ref,i_func_atm)

!                       --- charges
                        charge(i_gra) = 0.066d0
                        charge(n_gro-3) = 0.723d0
                        charge(n_gro-2) = -0.585d0
                        charge(n_gro-1) = -0.640d0
                        charge(n_gro) = 0.436d0

!			--- write OH - HO bond
                        gro_bond_list(n_gro-1,n_gro) = .true.
                        gro_bond_type(n_gro-1,n_gro) = 12
                        gro_bond_list(n_gro,n_gro-1) = .true.
                        gro_bond_type(n_gro,n_gro-1) = 12

        	end if

	end do

!     	--- add hydrogens for those not functionalised with carboxylic acid    
	do i_gra = 1, n_gra

!		--- go to next carbon if this one has been functionalised with carboxy
		if (carb_func(i_gra).eqv..true.) cycle

!          	--- functionalise edge carbon
    	    	if (n_gra_bond(i_gra).eq.2) then

!    	        	--- update functionality counter
			n_bond = 0
			func_ref = 4
			n_hydrogen = n_hydrogen + 1

!         		--- search for atoms bonded to this functional site
          		do j_gra = 1, n_gra

	          		if (gro_bond_list(i_gra,j_gra).eqv..true.) then
              				n_bond = n_bond + 1
			        	bond(n_bond) = j_gra
		     	       	end if
	      		      	if(n_bond.eq.n_gra_bond(i_gra)) exit

         		end do
        
!         		--- mid-point of the two atoms bonded to the functional site
          		bond_x = (gra_x(bond(1)) + gra_x(bond(2))) / 2.0d0
!                       --- temporary fix for periodic sheets
                        if ((gra_x(i_gra))-bond_x.gt.1.0) then
                                bond_x = bond_x + (lx/2.0d0)
                        else if ((gra_x(i_gra)-bond_x).lt.-1.0) then
                                bond_x = bond_x - (lx/2.0d0)
                        end if

          		bond_y = (gra_y(bond(1)) + gra_y(bond(2))) / 2.0d0
!			--- temporary fix for periodic sheets
			if ((gra_y(i_gra))-bond_y.gt.1.0) then
                                bond_y = bond_y + (ly/2.0d0)
                        else if ((gra_y(i_gra)-bond_y).lt.-1.0) then
                                bond_y = bond_y - (ly/2.0d0)
                        end if

          		bond_z = (gra_z(bond(1)) + gra_z(bond(2))) / 2.0d0
!                       --- temporary fix for periodic sheets
                        if ((gra_z(i_gra))-bond_z.gt.1.0) then
                                bond_z = bond_z + (lz/2.0d0)
                        else if ((gra_z(i_gra)-bond_z).lt.-1.0) then
                                bond_z = bond_z - (lz/2.0d0)
                        end if

!         		--- reference unit vector - from the mid-point to the functional site
          		call vector(gra_x(i_gra),gra_y(i_gra),gra_z(i_gra),bond_x,bond_y,bond_z,vec_x0,vec_y0,vec_z0)

!         		--- angle of the vector relative to the y axis (cosine rule) 
          		rads = acos(vec_y0/sqrt((vec_x0*vec_x0) + (vec_y0*vec_y0)))
          		if (vec_x0.gt.0.0d0) then
				theta = (2.0d0 * pi) - rads
			else
				theta = rads
			end if

!         		--- use the calculated vector for the first atom        
         		func_x_2 = func_x(func_ref,1) * func_x(func_ref,1)
          		func_y_2 = func_y(func_ref,1) * func_y(func_ref,1)
          		func_z_2 = func_z(func_ref,1) * func_z(func_ref,1)
          		func_xyz = sqrt(func_x_2 + func_y_2 + func_z_2)

!        		--- add new atom
          		n_gro = n_gro + 1
          		gro_x(n_gro) = gra_x(i_gra) + func_xyz * vec_x0
          		gro_y(n_gro) = gra_y(i_gra) + func_xyz * vec_y0
          		gro_z(n_gro) = gra_z(i_gra) + func_xyz * vec_z0
          		gro_type(n_gro) = func_type(func_ref,1)

                        n_func_atm = 1
                        mass_edge = mass_edge + 1.0079d0

!			--- charges
            		charge(i_gra) = -0.158d0
           		charge(n_gro) = 0.158d0

! 			--- CA - HA bond
            		gro_bond_list(i_gra,n_gro) = .true.
            		gro_bond_type(i_gra,n_gro) = 4
            		gro_bond_list(n_gro,i_gra) = .true.
            		gro_bond_type(n_gro,i_gra) = 4

        	end if

      	end do 

	Contains

      	Subroutine vector (x1,y1,z1,x0,y0,z0,vec_x,vec_y,vec_z)

      		implicit none
	     	integer, parameter :: dp = SELECTED_REAL_KIND(15,99)
      		real(kind=dp) :: x0, y0, z0, x1, y1, z1
      		real(kind=dp) :: x, y, z, x_2, y_2, z_2, vec_xyz, vec_x, vec_y, vec_z

	      	x = x1 - x0
      		y = y1 - y0
      		z = z1 - z0
      		x_2 = x * x
      		y_2 = y * y 
      		z_2 = z * z
      		vec_xyz = sqrt(x_2 + y_2 + z_2)
      		vec_x = x / vec_xyz
      		vec_y = y / vec_xyz
      		vec_z = z / vec_xyz

	End subroutine vector

      	Subroutine rotate_z (x,y,z,theta,x1,y1,z1)

      		implicit none
    		integer, parameter :: dp = SELECTED_REAL_KIND(15,99)
      		real(kind=dp) :: x,y,z,x1,y1,z1,theta
		x1 = x * cos(theta) - y * sin(theta)
      		y1 = x * sin(theta) + y * cos(theta)
      		z1 = z

	End subroutine

End subroutine func_edge


Subroutine func_plane

Use global

!     	--- loop over basal carbons and functionalise
	do i_loop = 1, max_func_loop

!               --- exit loop if the target oxygen weight percent is zero
                if (O_wp.lt.1.0d-6) exit

! 		--- exit loop when the target weight percentage of O is reached
  	        if (((mass_O/mass_tot)*1.0d2).ge.O_wp) exit

!		--- first functional group is always added with random selection of carbon to functionalise
		if (i_loop.eq.1) then
			new_patch = .true.
		else
!			--- functional site will only be selected at random if patch size is equal to maximum
        		if (patchy.eqv..true.) then 
				if (patch_size.lt.max_patch) then
					new_patch = .false.
				else
					new_patch = .true.
					patch_size = 0
				end if
			end if
        	end if

!		--- keep randomly picking carbons until we find a pair suitable
        	do i_ran = 1, 10000

!			--- choose one carbon at random
         		call random_number(ran)
          		i_gra = ceiling(ran * n_gra)

! 			--- check that this is an aromatic carbon
          		if (gro_type(i_gra).ne.'CA') cycle

!	         	--- check that this carbon has three bonds (not on the edge)
        	  	if (n_gra_bond(i_gra).ne.3) cycle

! 			--- zero bond counter
          		n_bond = 0
          		j_gra = 0

! 			--- find all bonded carbons
		        do k_gra = 1, n_gra

            			! --- check that there is a bond
            			if (gro_bond_list(i_gra,k_gra).eqv..false.) cycle

            			! --- check that this is an aromatic carbon
            			if (gro_type(k_gra).ne.'CA') cycle

            			! --- check that this is a basal plane carbon
            			if (n_gra_bond(k_gra).ne.3) cycle
         
     		      		j_gra = k_gra
            			n_bond = n_bond + 1

          		end do

!                       --- if there are no possible bond sites select a new carbon
                        if (n_bond.eq.0) cycle

!			--- exit now if not looking for existing patch
			if (new_patch.eqv..true.) exit

!			--- if still here then we need to confirm new site is adjacent to previously functionalised one

!			--- zero counters
			n_bond = 0

!			--- check first atom
			do k_gra = 1, n_gra

!	 			--- check that there is a bond to first atom
                        	if (gro_bond_list(i_gra,k_gra).eqv..false.) cycle				

!				--- check if the bonded atom is CE or CH
                               	if (gro_type(k_gra).eq.'CA') cycle

				n_bond = n_bond + 1

			end do

!			--- if the first atom did not have adjacent functionalised carbon try second atom
			if (n_bond.eq.0) then

!                               --- check first atom
                               	do k_gra = 1, n_gra					

!                            		--- check that there is a bond to second atom
                                       	if (gro_bond_list(j_gra,k_gra).eqv..false.) cycle

!                                      	--- check if the bonded atom is CE or CH
                                       	if (gro_type(k_gra).eq.'CA') cycle

					n_bond = n_bond + 1

				end do

			end if

!                       --- if no sites are suitable then select a new one
                        if (n_bond.eq.0) then
                                cycle
                        else if (n_bond.gt.0) then
                                exit
                        end if				

        	end do

! 		--- end program if further functionalisation cannot be achieved
        	if (j_gra.eq.0) then
 
			write(*,*) "No more sites for functional groups"
		   	write(*,*) "Try decreasing O_wp"
		 	stop

	        end if

!               --- decide at random if the first group is to point up or down
                call random_number(ran)
		if (ran.lt.0.5d0) then
                      	z_switch = 1.0d0
                else
                       	z_switch = -1.0d0
                end if

!               --- decide whether to add one epoxide or two hydroxyl groups
                call random_number(ran)

! 		--- add hydroxyl groups
! 		--- factor of 2/3 accounts for 2 hydroxyl groups are added at a time
        	if (ran.lt.(plane_rat*0.66667d0)) then

          		func_ref = 2
			patch_size = patch_size + 2

! 			--- assign new carbon types
         		gro_type(i_gra) = 'CH'
          		gro_type(j_gra) = 'CH'

! 			--- add functional group
          		do i_group = 1, 2

				n_hydroxy = n_hydroxy + 1

!                       	--- calculate current fraction of hydroxyls deprotonated
                        	dp_rat = real(n_alcoholate) / real(n_hydroxy)

!				--- second functionality is added to opposite side to the first
            			if (i_group.eq.2) then

			              z_switch = -1.0d0 * z_switch
         			      i_gra = j_gra

			        end if

! 				--- add first atom
			        func_x_2 = func_x(func_ref,1) * func_x(func_ref,1)
          			func_y_2 = func_y(func_ref,1) * func_y(func_ref,1)
            			func_z_2 = func_z(func_ref,1) * func_z(func_ref,1)
            			func_xyz = sqrt(func_x_2 + func_y_2 + func_z_2)
            			n_gro = n_gro + 1
            			gro_x(n_gro) = gra_x(i_gra)
            			gro_y(n_gro) = gra_y(i_gra)
            			gro_z(n_gro) = gra_z(i_gra) + (func_xyz * z_switch)
            			gro_type(n_gro) = func_type(func_ref,1)

!                               --- add alcoholate
                                if (dp_rat.lt.plane_dp) then
                                        n_func_atm = 1
                                        mass_O = mass_O + 1.59994d1
                                        mass_tot = mass_tot + 1.59994d1

!                                       --- update charges
                                        charge(i_gra) = -0.123d0
                                        charge(n_gro) = -0.877d0

!                                       --- CH - OH bond
                                        gro_bond_list(i_gra,n_gro) = .true.
                                        gro_bond_type(i_gra,n_gro) = 7
                                        gro_bond_list(n_gro,i_gra) = .true.
                                        gro_bond_type(n_gro,i_gra) = 7

!					--- update alcoholate counter
					n_alcoholate = n_alcoholate + 1

!					--- count up and down alcoholates to ensure same charge on both sides
					if (z_switch.eq.1.0d0) then
						n_up_alcoholate = n_up_alcoholate + 1
					else if (z_switch.eq.-1.0d0) then
						n_down_alcoholate = n_down_alcoholate + 1
					end if

!                               --- add hydroxide
                                else
                                        n_func_atm = 2
                                        mass_O = mass_O + 1.59994d1
                                        mass_tot = mass_tot + 1.70068d1

					n_gro = n_gro + 1

!					--- determine coordinates
              				func_x_1 = func_x(func_ref,1) - func_x(func_ref,2)
              				func_y_1 = func_y(func_ref,1) - func_y(func_ref,2)
              				func_z_1 = func_z(func_ref,1) - func_z(func_ref,2)
             		 		func_x_2 = func_x_1 * func_x_1
             				func_y_2 = func_y_1 * func_y_1
              				func_z_2 = func_z_1 * func_z_1
              				func_xyz = sqrt(func_x_2 + func_y_2 + func_z_2)
              				gro_x(n_gro) = gro_x(n_gro-1)
              				gro_y(n_gro) = gro_y(n_gro-1)
             				gro_z(n_gro) = gro_z(n_gro-1) + (func_xyz * z_switch)
              				gro_type(n_gro) = func_type(func_ref,2)

!					--- update charges
              				charge(i_gra) = 0.263d0
              				charge(n_gro-1) = -0.661d0
              				charge(n_gro) = 0.398d0

! 					--- CH - OH bond
				        gro_bond_list(i_gra,n_gro-1) = .true.
              				gro_bond_type(i_gra,n_gro-1) = 7
              				gro_bond_list(n_gro-1,i_gra) = .true.
              				gro_bond_type(n_gro-1,i_gra) = 7

! 					--- OH - HO bond
              				gro_bond_list(n_gro-1,n_gro) = .true.
              				gro_bond_type(n_gro-1,n_gro) = 12
              				gro_bond_list(n_gro,n_gro-1) = .true.
              				gro_bond_type(n_gro,n_gro-1) = 12

            			end if

          		end do

!       	--- add epoxide group
        	else

			func_ref = 3
          		n_epoxy = n_epoxy + 1
			patch_size = patch_size + 1
          		mass_O = mass_O + 1.59994d1
         		mass_tot = mass_tot + 1.59994d1

! 			--- assign atom types
          		gro_type(i_gra) = 'CE'
         		gro_type(j_gra) = 'CE'
          
! 			--- add new atom
          		n_gro = n_gro + 1

!			--- calculate coordinates
			gro_x(n_gro) = (gra_x(i_gra) + gra_x(j_gra)) / 2.0d0
			if ((gra_x(i_gra) - gra_x(j_gra)).gt.(lx/2.0d0)) then
				gro_x(n_gro) = gro_x(n_gro) - (lx/2.0d0)
			else if ((gra_x(i_gra) - gra_x(j_gra)).lt.(-1.0d0*lx/2.0d0)) then
				gro_x(n_gro) = gro_x(n_gro) + (lx/2.0d0)
			end if
          		gro_y(n_gro) = (gra_y(i_gra) + gra_y(j_gra)) / 2.0d0
                        if ((gra_y(i_gra) - gra_y(j_gra)).gt.(ly/2.0d0)) then
                                gro_y(n_gro) = gro_y(n_gro) - (ly/2.0d0)
                        else if ((gra_y(i_gra) - gra_y(j_gra)).lt.(-1.0d0*ly/2.0d0)) then
                                gro_y(n_gro) = gro_y(n_gro) + (ly/2.0d0)
                        end if
          		gro_z(n_gro) = gra_z(i_gra) + (0.12406d0 * z_switch)
          		gro_type(n_gro) = func_type(func_ref,1)

!			--- update charges
          		charge(i_gra) = 0.184d0
          		charge(j_gra) = 0.184d0
          		charge(n_gro) = -0.368d0

! 			--- CE - OE bond
          		gro_bond_list(i_gra,n_gro) = .true.
          		gro_bond_type(i_gra,n_gro) = 9
          		gro_bond_list(n_gro,i_gra) = .true.
          		gro_bond_type(n_gro,i_gra) = 9

! 			--- CE - OE bond
          		gro_bond_list(j_gra,n_gro) = .true.
          		gro_bond_type(j_gra,n_gro) = 9
          		gro_bond_list(n_gro,j_gra) = .true.
          		gro_bond_type(n_gro,j_gra) = 9

        	end if

	end do

End subroutine func_plane

End module functionalise
