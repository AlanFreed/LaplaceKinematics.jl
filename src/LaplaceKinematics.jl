#=
Created on Mon 22 Nov 2021
updated on Sun 03 Nov 2024
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
4) Freed, A.D., "A Technical Note: Two-Step PECE Methods for Approximating
   Solutions To First- and Second-Order ODEs", arXiv/2056770, 1 Nov, 2017.
-------------------------------------------------------------------------------
=#

"""
# LaplaceKinematics

This module provides data structures for kinematic descriptions from a  Lagrangian perspective for motions, deformations, stretches and strains, and their rates in 1, 2 and 3 dimensions. In the 2D and 3D cases, these descriptions are based upon measures for stretch that are triangular in structure; specifically, an upper triangular measure of stretch is used, which is referred to as the Lagrangian Laplace stretch.

## Installation

To use this module, you will need to add the following Julia packages to your project:
```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/LaplaceKinematics.jl")
```
In the various documentations of this package, *PF* is used as an alias for *PhysicalFields* and *LK* is used as an alias for *LaplaceKinematics*.
"""
module LaplaceKinematics

using
    JSON3,
    PhysicalFields,
    StructTypes
    
import
    PhysicalFields as PF

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

const DIMENSIONLESS = PF.CGS_DIMENSIONLESS
const LENGTH        = PF.CGS_LENGTH
const LENGTH_RATE   = PF.PhysicalUnits("CGS", 1, 0, 0, -1, 0, 0, 0)
const AREA          = PF.CGS_AREA
const AREA_RATE     = PF.PhysicalUnits("CGS", 2, 0, 0, -1, 0, 0, 0)
const VOLUME        = PF.CGS_VOLUME
const VOLUME_RATE   = PF.PhysicalUnits("CGS", 3, 0, 0, -1, 0, 0, 0)
const TIME          = PF.CGS_SECOND
const TIME_RATE     = PF.PhysicalUnits("CGS", 0, 0, 0, -1, 0, 0, 0)

include("LaplaceKinematics1D.jl")

include("LaplaceKinematics2D.jl")

# include("LaplaceKinematics3D.jl")

end # LaplaceKinematics
