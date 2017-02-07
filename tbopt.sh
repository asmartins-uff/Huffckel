#!/bin/bash
export LAPACK_LIBS="-L/usr/lib -llapack -lblas"
cd ./src
gfortran  -c -O2  -ffast-math  sto_ov.f90
gfortran  -c -O2  -ffast-math  supercell.f90
gfortran  -c -O2  -ffast-math  postproc.f90 
gfortran  -c -O2  -ffast-math  diagonalize.f90 ${LAPACK_LIBS}
gfortran  -c -O2  -ffast-math  toten.f90 ${LAPACK_LIBS}
gfortran  -c -O3  -ffast-math  Huckel_TB.f90 ${LAPACK_LIBS}
gfortran  -o  htb.x  Huckel_TB.o  sto_ov.o  supercell.o  postproc.o  diagonalize.o toten.o ${LAPACK_LIBS}
rm *.o *.mod
cd ..
mv ./src/htb.x .
