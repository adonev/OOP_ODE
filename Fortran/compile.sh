#!/bin/sh
# Set source directory:
export SRC=..
# Set compiler (usually gfortran)
export FC=gfortran
# If you have the Intel Fortran compiler
export FC=ifort

# Fortran compiler flags
#export FFLAGS=-I ${SRC} -I ${SRC}/FLib
# Recommended for debugging:
export FFLAGS="-g"
#export FFLAGS="-g -O1" # With -O1 the bug "disappears"

# Compile Fortran "library":
${FC} ${FFLAGS} -c ${SRC}/FLib/Precision.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/r8lib.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/Tridiagonal.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/MatrixExponential.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/MatrixInverse.f90

# Compile linear operators:
${FC} ${FFLAGS} -c ${SRC}/LinearOperator.f90
${FC} ${FFLAGS} -c ${SRC}/BandedMatrix.f90

# Make executable:
${FC} ${FFLAGS} -o TestBandedMatrix.x ${SRC}/TestBandedMatrix.f90 BandedMatrix.o LinearOperator.o Tridiagonal.o MatrixExponential.o MatrixInverse.o r8lib.o Precision.o -llapack
echo "Compiled executable TestBandedMatrix.x"

# Compile ODE solvers:
${FC} ${FFLAGS} -c ${SRC}/ODEIntegrand.f90

${FC} ${FFLAGS} -o TestODEIntegrator.x ${SRC}/TestODEIntegrator.f90 ODEIntegrand.o BandedMatrix.o LinearOperator.o Tridiagonal.o MatrixExponential.o MatrixInverse.o r8lib.o Precision.o -llapack
echo "Compiled executable TestODEIntegrator.x"
# ifort -c ../BandedMatrix.f90 2>&1 | head -6
