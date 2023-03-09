#!/bin/sh
# Set source directory:
export SRC=..
# Set compiler (usually gfortran)
export FC=gfortran
# If you have the Intel Fortran compiler
export FC=ifort

# Fortran compiler flags
#export FFLAGS=-I ${SRC} -I ${SRC}/FLib
export FFLAGS="-g -O1"

# Compile Fortran "library":
${FC} ${FFLAGS} -c ${SRC}/FLib/Precision.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/r8lib.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/Tridiagonal.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/MatrixExponential.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/MatrixInverse.f90

# Compile ODE solver:
${FC} ${FFLAGS} -c ${SRC}/LinearOperator.f90
${FC} ${FFLAGS} -c ${SRC}/BandedMatrix.f90

# Make executable:
${FC} ${FFLAGS} -o TestODEIntegrator.x ${SRC}/TestODEIntegrator.f90 BandedMatrix.o LinearOperator.o Tridiagonal.o MatrixExponential.o MatrixInverse.o r8lib.o Precision.o -llapack

# ifort -c ../BandedMatrix.f90 2>&1 | head -6
