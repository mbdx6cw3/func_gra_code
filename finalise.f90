Module finalise

Contains

Subroutine write_itp

	use global

!       --- format statements
      	200 format (i6, a11, a7, a7, a7, i7, f12.3, f12.3)
      	300 format (i6, i6, a6, f12.5, f12.2)
      	400 format (i6, i6, i6, a6, f12.3, f12.3, f12.6, f12.3)
     	500 format (i6, i6, i6, i6, a6, f12.3, f12.6, i6)

!       --- modify bond types for existing graphene atoms
	do i_gra = 1, n_gra

        	do j_gra = 1, n_gra

          		if (gro_bond_list(i_gra,j_gra).eqv..true.) then

            			if ((gro_type(i_gra).eq.'CA').and.(gro_type(j_gra).eq.'CA')) then

              				gro_bond_type(i_gra,j_gra) = 1

            			else if ((gro_type(i_gra).eq.'CA').and.((gro_type(j_gra).eq.'CH').or.(gro_type(j_gra).eq.'CE'))) then

              				gro_bond_type(i_gra,j_gra) = 2

            			else if (((gro_type(i_gra).eq.'CH').or.(gro_type(i_gra).eq.'CE')).and.(gro_type(j_gra).eq.'CA')) then

             			 	gro_bond_type(i_gra,j_gra) = 2

		 	        else if ((gro_type(i_gra).eq.'CH').and.(gro_type(j_gra).eq.'CH')) then

			        	gro_bond_type(i_gra,j_gra) = 5

			        else if (((gro_type(i_gra).eq.'CH').and.(gro_type(j_gra).eq.'CE')).or.&

             				&((gro_type(i_gra).eq.'CH').and.(gro_type(j_gra).eq.'CE'))) then

              				gro_bond_type(i_gra,j_gra) = 6

            			else if ((gro_type(i_gra).eq.'CE').and.(gro_type(j_gra).eq.'CE')) then

              				gro_bond_type(i_gra,j_gra) = 8

            			end if

          		end if

        	end do

	end do
   
	open(unit=1,file='graphene_oxide.itp',status='unknown')

!     	--- write molecule section
      	write(1,*) '[ moleculetype ]'
      	write(1,*) 'GRO                3'
      	write(1,*)

!     	--- zero charge counter and topology reference
      	qtot = 0.0d0
      	itp_ref = 0

!     	--- write atom section
      	write(1,*) '[ atoms ]'
      	do i_gro = 1, n_gro

        	itp_count = itp_count + 1

	        if (gro_type(i_gro).eq.'CA') then

       			itp_ref = 1

        	else if (gro_type(i_gro).eq.'CH') then

          		itp_ref = 2

        	else if (gro_type(i_gro).eq.'CE') then

          		itp_ref = 3

        	else if (gro_type(i_gro).eq.'CC') then

         		itp_ref = 4

       		else if (gro_type(i_gro).eq.'OH') then

          		itp_ref = 5

        	else if (gro_type(i_gro).eq.'OE') then

          		itp_ref = 6

        	else if (gro_type(i_gro).eq.'OC') then

          		itp_ref = 7

        	else if (gro_type(i_gro).eq.'HA') then

          		itp_ref = 8

        	else if (gro_type(i_gro).eq.'HO') then

          		itp_ref = 9

        	end if

	        write(1,200) itp_count, gro_type(i_gro), '1', 'GRO', gro_type(i_gro), itp_count, charge(i_gro), mass(itp_ref)
       		qtot = qtot + charge(i_gro)

	end do

!     	--- check that total charge is zero
!      	if (abs(qtot).gt.1.0d-8) then

!		write(*,*) "WARNING: total charge = ", qtot
!		write(*,*) "Non-zero charge is only fine if pH > pKa(COOH)"

!	end if 

!     	--- copy new bond list for topology file
	write(1,*)
      	write(1,*) '[ bonds ]'
      	bond_count = 0
      	do i_gro = 1, n_gro - 1 

	        do j_gro = i_gro + 1, n_gro

          		if (gro_bond_list(i_gro,j_gro).eqv..true.) then

            			bond_count = bond_count + 1
            			itp_ref = gro_bond_type(i_gro,j_gro)
            			write(1,300) i_gro, j_gro, '1', bond_r(itp_ref), bond_kr(itp_ref)

          		end if

        	end do

	end do

	write(1,*)

!   	--- copy new angle list for topology file
      	write(1,*)
      	write(1,*) '[ angles ]'

!     	--- write angle list (j_gro - i_gro - k_gro)
      	angle_count = 0
      	do i_gro = 1, n_gro

	        do j_gro = 1, n_gro

			if (gro_bond_list(i_gro,j_gro).eqv..true.) then

			  	do k_gro = j_gro + 1, n_gro

			        	if (gro_bond_list(i_gro,k_gro).eqv..true.) then

				                angle_count = angle_count + 1

				                if (gro_type(i_gro).eq.'CA') then

					        	if ((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'CA')) then

					                	itp_ref = 1

					                else if (((gro_type(j_gro).eq.'CA').and.((gro_type(k_gro).eq.'CH').or.&
							& (gro_type(k_gro).eq.'CE'))).or.(((gro_type(j_gro).eq.'CH').or.&
							& (gro_type(j_gro).eq.'CE')).and.(gro_type(k_gro).eq.'CA'))) then
                    				
								itp_ref = 2

					                else if (((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'CC')).or.&
                 					& ((gro_type(j_gro).eq.'CC').and.(gro_type(k_gro).eq.'CA'))) then

					   			itp_ref = 3  

               						else if (((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'HA')).or.&
                 					& ((gro_type(j_gro).eq.'HA').and.(gro_type(k_gro).eq.'CA'))) then

                    						itp_ref = 4

                  					else if (((gro_type(j_gro).eq.'CH').and.(gro_type(k_gro).eq.'HA')).or.&
                 					& ((gro_type(j_gro).eq.'HA').and.(gro_type(k_gro).eq.'CH'))) then

								itp_ref = 5                

					                else if (((gro_type(j_gro).eq.'CE').and.(gro_type(k_gro).eq.'HA')).or.&
                 					& ((gro_type(j_gro).eq.'HA').and.(gro_type(k_gro).eq.'CE'))) then

                    						itp_ref = 5

                  					end if

                				else if (gro_type(i_gro).eq.'CH') then

                  					if (((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'CA')).or.&
							& ((gro_type(j_gro).eq.'CH').and.(gro_type(k_gro).eq.'CA')).or.&
                                  			& ((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'CH')).or.&
							& ((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'CE')).or.&
							& ((gro_type(j_gro).eq.'CE').and.(gro_type(j_gro).eq.'CA'))) then

                    						itp_ref = 6

                  					else if (((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'OH')).or.&
                 					& ((gro_type(j_gro).eq.'OH').and.(gro_type(k_gro).eq.'CA'))) then

                    						itp_ref = 7

                  					else if (((gro_type(j_gro).eq.'CH').and.(gro_type(k_gro).eq.'OH')).or.&
                 					& ((gro_type(j_gro).eq.'OH').and.(gro_type(k_gro).eq.'CH'))) then

                    						itp_ref = 8

                  					else if (((gro_type(j_gro).eq.'CE').and.(gro_type(k_gro).eq.'OH')).or.&
                 					& ((gro_type(j_gro).eq.'OH').and.(gro_type(k_gro).eq.'CE'))) then

                    						itp_ref = 8

                  					else if (((gro_type(j_gro).eq.'CH').and.((gro_type(k_gro).eq.'CH').or.&
							& (gro_type(k_gro).eq.'CE'))).or.((gro_type(j_gro).eq.'CE').and.&
							& ((gro_type(k_gro).eq.'CH').or.(gro_type(k_gro).eq.'CE')))) then 

                    						itp_ref = 9

                  					end if

                				else if (gro_type(i_gro).eq.'CE') then

                  					if (((gro_type(j_gro).eq.'CE').and.(gro_type(k_gro).eq.'OE')).or.&
                 					& ((gro_type(j_gro).eq.'OE').and.(gro_type(k_gro).eq.'CE'))) then

                    						itp_ref = 10

                  					else

                    						itp_ref = 11

                  					end if

                				else if (gro_type(i_gro).eq.'CC') then

                  					if (((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'OC')).or.&
                 					& ((gro_type(j_gro).eq.'OC').and.(gro_type(k_gro).eq.'CA'))) then 

                    						itp_ref = 12  

                  					else if (((gro_type(j_gro).eq.'CA').and.(gro_type(k_gro).eq.'OH')).or.&
                 					& ((gro_type(j_gro).eq.'OH').and.(gro_type(k_gro).eq.'CA'))) then

                    						itp_ref = 13

                  					else if (((gro_type(j_gro).eq.'OC').and.(gro_type(k_gro).eq.'OH')).or.&
                 					& ((gro_type(j_gro).eq.'OH').and.(gro_type(k_gro).eq.'OH'))) then

                    						itp_ref = 14

                  					end if

                				else if (gro_type(i_gro).eq.'OH') then

                  					if (((gro_type(j_gro).eq.'CH').and.(gro_type(k_gro).eq.'HO')).or.&
                 					& ((gro_type(j_gro).eq.'HO').and.(gro_type(k_gro).eq.'CH'))) then           

                    						itp_ref = 15

                  					else if (((gro_type(j_gro).eq.'CC').and.(gro_type(k_gro).eq.'HO')).or.&
                 					& ((gro_type(j_gro).eq.'HO').and.(gro_type(k_gro).eq.'CC'))) then

                    						itp_ref = 16

                  					end if

                				else if (gro_type(i_gro).eq.'OE') then

                  					itp_ref = 17

                				else

                  					itp_ref = 0
                  					angle_count = angle_count - 1

                				end if 

                				if (itp_ref.ne.0) then

                  					write(1,400) j_gro, i_gro, k_gro, '5', ang_theta(itp_ref), ang_k(itp_ref), ang_ub(itp_ref), ang_kub(itp_ref)

                				end if

              				end if

            			end do

          		end if

        	end do

	end do
      	write(1,*)

!     	--- copy new angle list for topology file
      	write(1,*)
      	write(1,*) '[ dihedrals ]'

!     	--- write dihedral list (k_gro - i_gro - j_gro - l_gro)
      	dihedral_count = 0
      	do i_gro = 1, n_gro

	        do j_gro = 1, n_gro

			itp_ref = 0
			if (gro_bond_list(i_gro,j_gro).eqv..true.) then
  
		        	do k_gro = 1, n_gro

		        		if (k_gro.eq.j_gro) cycle

					if (gro_bond_list(i_gro,k_gro).eqv..true.) then

						do l_gro = k_gro + 1, n_gro

							if (l_gro.eq.i_gro) cycle

							if (gro_bond_list(j_gro,l_gro).eqv..true.) then

								dihedral_count = dihedral_count + 1

                    						if ((gro_type(i_gro).eq.'CA').and.(gro_type(j_gro).eq.'CA')) then

									if ((gro_type(k_gro).eq.'HA').and.(gro_type(l_gro).eq.'HA')) then

										itp_ref = 3

						                        else if (((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'HA')).or.&
                     							& ((gro_type(k_gro).eq.'CH').and.(gro_type(l_gro).eq.'HA')).or.&
                     							& ((gro_type(k_gro).eq.'CE').and.(gro_type(l_gro).eq.'HA')).or.&
                     							& ((gro_type(k_gro).eq.'CC').and.(gro_type(l_gro).eq.'HA')).or.&   
                     							& ((gro_type(k_gro).eq.'HA').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'HA').and.(gro_type(l_gro).eq.'CH')).or.&
                     							& ((gro_type(k_gro).eq.'HA').and.(gro_type(l_gro).eq.'CE')).or.&
                     							& ((gro_type(k_gro).eq.'HA').and.(gro_type(l_gro).eq.'CC'))) then
                       		
										 itp_ref = 2

						                        else 

							                         itp_ref = 1

						                        end if

								else if (((gro_type(i_gro).eq.'CA').and.(gro_type(j_gro).eq.'CH')).or.&
                   						& ((gro_type(i_gro).eq.'CH').and.(gro_type(j_gro).eq.'CA'))) then

				  	           	                if (((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'CH')).or.&
                     		  				        & ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'CH').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'CE').and.(gro_type(l_gro).eq.'CA')).or.&         
                    							& ((gro_type(k_gro).eq.'CC').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'OH').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'CH')).or.&
                   							& ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'CE')).or.&
                     							& ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'CC')).or.&
                     							& ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'OH'))) then

							                        itp_ref = 4

						                        end if

						                else if (((gro_type(i_gro).eq.'CA').and.(gro_type(j_gro).eq.'CC')).or.&
                   						& ((gro_type(i_gro).eq.'CC').and.(gro_type(j_gro).eq.'CA'))) then
                      				
									if (((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'OC')).or.&
                     							& ((gro_type(k_gro).eq.'CH').and.(gro_type(l_gro).eq.'OC')).or.&
                     							& ((gro_type(k_gro).eq.'CE').and.(gro_type(l_gro).eq.'OC')).or.&
                     							& ((gro_type(k_gro).eq.'OC').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'OC').and.(gro_type(l_gro).eq.'CH')).or.&
                     							& ((gro_type(k_gro).eq.'OC').and.(gro_type(l_gro).eq.'CE'))) then
                        
										itp_ref = 5

                      							end if

					              		else if ((gro_type(i_gro).eq.'CH').and.(gro_type(j_gro).eq.'CH')) then

                       								itp_ref = 6

						            	else if (((gro_type(i_gro).eq.'CA').and.(gro_type(j_gro).eq.'CE')).or.&
						                & ((gro_type(i_gro).eq.'CH').and.(gro_type(j_gro).eq.'CE')).or.&
								& ((gro_type(i_gro).eq.'CE').and.(gro_type(j_gro).eq.'CA')).or.&
								& ((gro_type(i_gro).eq.'CE').and.(gro_type(j_gro).eq.'CH'))) then

						                	if (((gro_type(k_gro).eq.'OE').and.(gro_type(l_gro).eq.'HA')).or.&
                     							& ((gro_type(k_gro).eq.'HA').and.(gro_type(l_gro).eq.'OE'))) then

                        							itp_ref = 8

                      							else if (((gro_type(k_gro).eq.'OE').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'OE').and.(gro_type(l_gro).eq.'CH')).or.&
                     							& ((gro_type(k_gro).eq.'OE').and.(gro_type(l_gro).eq.'CE')).or.&
                     							& ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'OE')).or.&
                     							& ((gro_type(k_gro).eq.'CH').and.(gro_type(l_gro).eq.'OE')).or.&
                     							& ((gro_type(k_gro).eq.'CE').and.(gro_type(l_gro).eq.'OE'))) then

							                        itp_ref = 7

						                        end if

								else if ((gro_type(i_gro).eq.'CE').and.(gro_type(j_gro).eq.'CE')) then

						                        if (((gro_type(k_gro).eq.'OE').and.(gro_type(l_gro).eq.'CA')).or.&
                    							& ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'OE')).or.&
                     							& ((gro_type(k_gro).eq.'OE').and.(gro_type(l_gro).eq.'CH')).or.&
                     							& ((gro_type(k_gro).eq.'CH').and.(gro_type(l_gro).eq.'OE')).or.&
                     							& ((gro_type(k_gro).eq.'OE').and.(gro_type(l_gro).eq.'CE')).or.&
                     							& ((gro_type(k_gro).eq.'CE').and.(gro_type(l_gro).eq.'OE'))) then

 							                        itp_ref = 9

 						                  	end if

                    					        else if (((gro_type(i_gro).eq.'CH').and.(gro_type(j_gro).eq.'OH')).or.&
                   					        & ((gro_type(i_gro).eq.'OH').and.(gro_type(j_gro).eq.'CH'))) then

                      					       		if (((gro_type(k_gro).eq.'HO').and.(gro_type(l_gro).eq.'CA')).or.&
   							                & ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'HO')).or.&
                     							& ((gro_type(k_gro).eq.'HO').and.(gro_type(l_gro).eq.'CH')).or.&
                     							& ((gro_type(k_gro).eq.'CH').and.(gro_type(l_gro).eq.'HO')).or.&
                     							& ((gro_type(k_gro).eq.'HO').and.(gro_type(l_gro).eq.'CE')).or.&
                     							& ((gro_type(k_gro).eq.'CE').and.(gro_type(l_gro).eq.'HO'))) then

                        							itp_ref = 10

                      							end if

                    						else if (((gro_type(i_gro).eq.'CC').and.(gro_type(j_gro).eq.'OH')).or.&
                   						& ((gro_type(i_gro).eq.'OH').and.(gro_type(j_gro).eq.'CC'))) then

                      							if (((gro_type(k_gro).eq.'HO').and.(gro_type(l_gro).eq.'CA')).or.&
                     							& ((gro_type(k_gro).eq.'CA').and.(gro_type(l_gro).eq.'HO')).or.&
                     							& ((gro_type(k_gro).eq.'OC').and.(gro_type(l_gro).eq.'HO')).or.&
                     							& ((gro_type(k_gro).eq.'HO').and.(gro_type(l_gro).eq.'OC'))) then

                        							itp_ref = 11

                      							end if
                    
								else

! 									--- all other dihedral angles are set to zero
                      							itp_ref = 0
                     		 					dihedral_count = dihedral_count - 1

                    						end if

					                        if (itp_ref.ne.0) then

					   		       		write(1,500) k_gro, i_gro, j_gro, l_gro, '9', dih_phi(itp_ref), dih_k(itp_ref), dih_mult(itp_ref)

						                end if
							end if 
						end do
					end if
				end do
			end if
		end do
	end do

End subroutine write_itp

Subroutine write_gro
      	
	use global
      	100 format (a13, a2, i5, f8.3, f8.3, f8.3)
     	open(unit=20,file='graphene_oxide.gro')
      	write(20,*) 'Graphene Oxide'
      	write(20,*) n_gro
      	index = 0
      	do i_gro = 1, n_gro

	        index = index + 1
		write(20,100) '1GRO    ', gro_type(i_gro), index, gro_x(i_gro), gro_y(i_gro), gro_z(i_gro)
	end do
	write(20,*) lx, ly, lz
	close(unit=20)

End subroutine write_gro

Subroutine write_output

	use global

      	open(unit=30,file='output.dat')
     	write(30,*) "Total Number of Atoms in Topology = ", n_gro
      	write(30,*) "Total Number of Bonds in Topology = ", bond_count
      	write(30,*) "Total Number of Angles in Topology = ", angle_count
      	write(30,*) "Total Number of Dihedrals in Topology = ", dihedral_count
	write(30,*) "Total Charge = ", qtot
        write(30,*)
	write(30,*) "C : O : H atomic ratio =", real(n_gra + n_carboxy) / real(n_gra + n_carboxy),":",&
	& real(n_epoxy + n_hydroxy + n_carboxy * 2) / real(n_gra + n_carboxy),":",&
	& real(n_hydrogen + n_carboxy - n_carboxylate + n_hydroxy - n_alcoholate) / real(n_gra + n_carboxy)
        write(30,*) "Total Mass =", mass_tot + mass_edge
        write(30,*) "Total Weight Percent of Oxygen on Flake = ", real (n_epoxy + n_hydroxy + n_carboxy * 2) &
	& * 15.999d0 / (mass_tot + mass_edge)
	write(30,*)
	write(30,*) "Weight Percent of Oxygen on Basal Plane = ", (mass_O/mass_tot)*1.0d2, "%"
      	write(30,*) "Number of Hydroxyl Groups = ", n_hydroxy
	write(30,*) "Number of Deprotonated Hydroxyl Groups = ", n_alcoholate
	write(30,*) "Number of Deprotonated Hydroxyl Groups Per Side = ", n_up_alcoholate, n_down_alcoholate
      	write(30,*) "Number of Epoxide Groups = ", n_epoxy
      	write(30,*) "Fraction of Hydroxyl and Epoxide Groups = ", real(n_hydroxy) / real(n_hydroxy + n_epoxy), &
     	& real(n_epoxy) / real(n_hydroxy + n_epoxy)
	write(30,*)
      	write(30,*) "Number of Carboxylic Acid Groups = ", n_carboxy
	write(30,*) "Number of Deprotonated Carboxylic Acid Groups = ", n_carboxylate
      	write(30,*) "Number of Hydrogen Groups = ", n_hydrogen
      	write(30,*) "Fraction of Carboxylic Acid / Carboxylate and Hydrogen Groups = ",&
	&  real(n_carboxy+n_carboxylate) / real(n_carboxy+n_carboxylate+n_hydrogen), & 
     	& real(n_hydrogen) / real(n_carboxy+n_carboxylate+n_hydrogen)
	write(30,*)

End subroutine write_output

End module finalise
