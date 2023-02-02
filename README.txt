     Func_GO

     Christopher D. Williams (christopher.williams@manchester.ac.uk)
     Feb. 2017
     --------------------------------------------------------------------------

     Overview of this program:
     This program can be used to oxidise a pristine graphene flake to make 
     graphene oxide, generating the input files (.gro and .itp) required to 
     perform an MD simulation using GROMACS.

     --------------------------------------------------------------------------

     Compiling the program:
     	- simply use the "make" or "gmake" command to generate the program 
          executable, "Func_GO".
        - the default compiler is gfortran.
     	- the executable provided with this distribution was compiled using the
          gfortran compiler GNU Fortran (GCC) 4.9.2 on the Intel Core i5 
          processor, but I would strongly reccommend you compile your own version

     --------------------------------------------------------------------------

     Files provided with this distribution:
	- README.txt
	- Makefile
	- *.f90:	  source files
	- func_GO:	  executable
	- input.dat:	  contains input parameters for functionalisation 
	- parameters.dat: contains parameters for the intermolecular potential.
			  The example provided is for the CHARMM force field.
        - graphene.gro:   coordinate file.
		 	  The example provided is for a 3 nm x 3 nm flake.
			  Easy to generate using VMD nanotube builder.
        - carboyxl.gro, hydroxyl.gro, epoxide.gro, hydrogen.gro:
			  Example functional group coordinate files.
     Output files:
	- graphene_oxide.gro:	GO coordinate file (GROMACS format)
	- graphene_oxide.itp:	GO topology file (GROMACS format)
	- output.dat:		output data concerning the functionalisation, 
				e.g. oxygen content, functional group ratios, etc.

     --------------------------------------------------------------------------

     Notes:

	1) Carbon atoms are functionalised with oxygen containing groups taking a
       	user-specified weight percentage as oxygen to be the target. 

	2) The functionalisation routine considers two separate types of
	functional group: edge sites and surface sites

	3) The oxygen content is determined for surface sites ONLY and edge 
	functional groups are ignored in this calculation because we are 
	trying to mimic larger sheets. If surface sites were to be 
        included in this calculation then the oxygen content would be 
	dominated by edge sites and therefore unrealistic

      	4) Edge sites are functionalised with hydrogen and carboxylic acid groups
       	with a user-defined ratio and basal plane sites are functionalised with
       	hydroxy and epoxy groups with a user-defined ratio.

	5) Basal plane functionalisation is carried out in pairs, to preserve
	stoichiometry and prevent undercoordinated carbon atoms. This means
        that two adjacent hydroxyl groups are added simultaneously, but only 
	one epoxide group is added.

	6) Atom types:
         CA  - aromatic carbon
         CH  - carbon with OH group attached
         CE  - carbon in epoxide group
         CC  - carbon in carboxyl group
         OH  - oxygen in hydroxyl group
         OE  - oxygen in epoxide group
         OC  - oxygen in carboxyl group
         HA  - aromatic hydrogen
         HO  - hydrogen in hydroxyl group

        7) The input parameter "pH" can be used to determine whether -OH or -COOH
        are protonated or deprotonated. The pKa of each functional group should be 
	provided in the input file.

	8) Currently, surface functionalisation is entirely random (i.e.
	a group type is selected at random, added to a random carbon atom
 	and a direction is randomly chosen. It may be useful to adapt the
	program to perform "patchy" functionalisation. 

        9) To make graphene instead of graphene oxide set the oxygen weight
        percentage to zero and the carboxylic acid/hydrogen ratio to zero.

	10) Instead of a finite size flake, a periodic sheet can also be
	functionalised. In this case only surface sites are functionalised. 
	There is a switch in input.dat to allow for this.

     --------------------------------------------------------------------------

     Citation:
	Add this later.

