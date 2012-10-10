FC=gfortran
LD=gfortran
FFLAGS=-ffree-line-length-0 -fbounds-check

.PHONY: all clean
all: neat

%.o: %.f90
	$(FC) $(FFLAGS) $< -c -o $@

neat: source/types.o source/extinction.o source/rec_lines.o source/helium.o source/equib_routines.o source/filereading.o source/abundances.o source/quicksort.o source/neat.o 
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

clean:
	rm -f neat source/*.o source/*.mod
