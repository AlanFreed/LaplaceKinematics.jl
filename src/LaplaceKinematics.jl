#=
Created on Mon 22 Nov 2021
updated on Mon 01 Apr 2024
-------------------------------------------------------------------------------
References:
1) Freed, A.D., Erel, V. and Moreno, M.R., "Conjugate stress/strain base pairs
   for planar analysis of biological tissues", Journal of Mechanics of
   Materials and Structures, 12 (2017), 219-247.
   DOI: 10.2140/jomms.2017.12.219.
2) Freed A.D., Zamani, S., Szabo̕, L. and Clayton, J.D., "Laplace stretch:
   Eulerian and Lagrangian formulations", Zeitschrift fur angewandte Mathematik
   und Physik, 71 (2020), 157. DOI: 10.1007/s00033-020-01388-4.
3) Paul, S., Freed, A.D. and Clayton, J.D., "Coordinate indexing: On the use of
   Eulerian and Lagrangian Laplace stretches", Applications in Engineering
   Science, 5 (2021), 100029. DOI: 10.1016/j.apples.2020.100029.
-------------------------------------------------------------------------------
=#

"""
Module\n
    LaplaceKinematics\n
This module provides data structures for Lagrangian kinematic descriptions of motion, deformation, stretch and strain, and their rates in 1, 2 and 3 dimensions. These descriptions are based upon measures for stretch that are triangular in structure; specifically, an upper-triangular measure of stretch is used, which is referred to as the Lagrangian Laplace stretch.
"""
module LaplaceKinematics

using
    JSON3,
    PhysicalFields,
    StructTypes

export
    # Permutation Matrices
    P2D₁,
    P2D₂,
    P3D₁,
    P3D₂,
    P3D₃,
    P3D₄,
    P3D₅,
    P3D₆,

    # types
    FiberKinematics,
    MembraneKinematics,
    # Kinematics,

    # Methods
    copy,
    deepcopy,

    fromFile,
    toFile,

    advance!,
    update!

const DIMENSIONLESS = PhysicalFields.CGS_DIMENSIONLESS
const LENGTH        = PhysicalFields.CGS_LENGTH
const LENGTH_RATE   = PhysicalFields.PhysicalUnits("CGS", 1, 0, 0, -1, 0, 0, 0)
const AREA          = PhysicalFields.CGS_AREA
const AREA_RATE     = PhysicalFields.PhysicalUnits("CGS", 2, 0, 0, -1, 0, 0, 0)
const VOLUME        = PhysicalFields.CGS_VOLUME
const VOLUME_RATE   = PhysicalFields.PhysicalUnits("CGS", 3, 0, 0, -1, 0, 0, 0)
const TIME          = PhysicalFields.CGS_SECOND
const TIME_RATE     = PhysicalFields.PhysicalUnits("CGS", 0, 0, 0, -1, 0, 0, 0)

include("LaplaceKinematics1D.jl")

include("LaplaceKinematics2D.jl")

# include("LaplaceKinematics3D.jl")

end # LaplaceKinematics
