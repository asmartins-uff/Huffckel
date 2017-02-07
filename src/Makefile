FC = gfortran
FFLAGS = -O2 -ffast-math
LAPACK_LIBS = -L/usr/lib -llapack -lblas
OBJ = sto_ov.o supercell.o postproc.o diagonalize.o toten.o Huckel_TB.o

%.o: %.f90
	$(FC) -c -o $@ $(FFLAGS) $< 

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $< 

huffckel: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) $(LAPACK_LIBS)

clean: 
	rm -fr *.mod *.o huffckel
