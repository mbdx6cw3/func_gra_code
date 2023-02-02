# Makefile: makes the func_GO program
# Use this Makefile with "gmake"

# Executable name
CMD = Func_GO

# Fortran compiler
FC = gfortran

# define general flags
GEN_FLAGS = -cpp -fcheck=bounds

# define optimisation flags
O_FLAGS = -O2

# targets and rules
%.o : %.f90
	$(FC) $(GEN_FLAGS) $(O_FLAGS) -c $<

OBJECTS = global.o initialise.o functionalise.o finalise.o main.o

$(CMD) : $(OBJECTS)
	$(FC) $(OFLAGS) $(OBJECTS) -o $(CMD)

clean:
	rm *.mod *.o

