cmake_minimum_required(VERSION 2.6)

project(fortran-lib Fortran)

set(lib-src ./constants_lib.f90 ./math_routines_lib.f90 ./stress_routines_lib.f90)

add_library(MechanicsFortran ${lib-src} )