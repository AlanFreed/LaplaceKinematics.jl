To properly view this README file use, e.g., the Remarkable Markdown editor.

# LaplaceKinematics.jl

This module provides a **QR** or Gram-Schmidt factorization of a deformation gradient tensor denoted as **F** = Fᵢⱼ **E**ᵢ ⊗ **E**ⱼ, where *i*, *j* ∈ {1,2} for 2D analyses, or *i*, *j* ∈ {1,2,3} for 3D analyses, with Fᵢⱼ being its matrix of tensor components evaluated in basis (**E**₁, **E**₂) for 2D analyses, or in basis (**E**₁, **E**₂, **E**₃) for 3D analyses. Given **Q** = δᵢⱼ **e**ᵢ ⊗ **E**ⱼ = Qᵢⱼ **e**ᵢ ⊗ **e**ⱼ, matrix Qᵢⱼ from **QR** is a proper orthogonal matrix known as the Gram rotation tensor, while matrix Rᵢⱼ of **R** = Rᵢⱼ **e**ᵢ ⊗ **e**ⱼ from **QR** is an upper-triangular matrix with positive diagonal elements that in mechanics is called the Laplace stretch tensor, as Laplace introduced the mathematical technique employed by Gram in his Ph.D. thesis that Schmidt later made popular.

In the mechanics literature, when analyzed from a Lagrangian perspective, such decompositions are often denoted as **F** = **RU**, where **R** denotes a rotation and **U** denotes a stretch. There are an infinite number of possible **RU** matrix products that yield **F**. The deconstruction of a matrix into a product of two matrices of which one is orthogonal is not unique. We are interested in the one wherein **R** is a Gram rotation matrix in what is otherwise commonly known as a **QR** matrix decomposition in the linear algebra literature.

A deformation gradient is a mapping from one material configuration into another. It is a two-state property, independent of the path traversed between these two states. All deformation gradients considered herein satisfy the following conditions: 

1.  The deformation gradient **F** equates with the identity tensor **I** = δᵢⱼ **E**ᵢ ⊗ **E**ⱼ in its initial configuration κ₀, i.e., its motion is κ₀ ↦ κ₀. 

2.  All other deformation gradients associate with motions of κ₀ ↦ κₙ, where κₙ is the configuration at current time tₙ. This includes a strain-free reference configuration κᵣ, which can be attained through a motion of κ₀ ↦ κᵣ.

3.  The reference configuration κᵣ is, by definition, free from strain; therefore, the initial configuration κ₀ can be a strained configuration caused by motion κᵣ ↦ κ₀, even though the initial deformation gradient is **F**(t₀) = **I**. Although seldom adopted, this ought to be the norm for analyses of biologic tissue.

To use this module you will need to add the following Julia packages to your project:
```julia
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/LaplaceKinematics.jl")
```
and to run the test files, you'll also need to install the packages
```julia
Pkg.add(url = "https://github.com/AlanFreed/CubicSplines.jl")
Pkg.add(url = "https://github.com/AlanFreed/FijLung.jl")
```
In this documentation, *PF* is an alias for module *PhysicalFields* and *LK* is an alias for module *LaplaceKinematics*.

The following physical units are used in this module:
```julia
const DIMENSIONLESS = PF.CGS_DIMENSIONLESS
const RATE          = PF.PhysicalUnits("CGS", 0, 0, 0, -1, 0, 0, 0)
const RATEofRATE    = PF.PhysicalUnits("CGS", 0, 0, 0, -2, 0, 0, 0)

const TIME          = PF.CGS_SECOND
const LENGTH        = PF.CGS_LENGTH
const LENGTH_RATE   = PF.PhysicalUnits("CGS", 1, 0, 0, -1, 0, 0, 0)
const AREA          = PF.CGS_AREA
const AREA_RATE     = PF.PhysicalUnits("CGS", 2, 0, 0, -1, 0, 0, 0)
const VOLUME        = PF.CGS_VOLUME
const VOLUME_RATE   = PF.PhysicalUnits("CGS", 3, 0, 0, -1, 0, 0, 0)
```

1.  [Kinematics of 1D Fibers](./README_1D.md)

2.  [Kinematics of 2D Membranes](./README_2D.md)

3.  [Kinematics of 3D Materials](./README_3D.md)

## References

1) Freed, A.D., Erel, V. and Moreno, M.R., "Conjugate stress/strain base pairs for planar analysis of biological tissues", *Journal of Mechanics of Materials and Structures*, **12** (2017), 219-247. DOI: 10.2140/jomms.2017.12.219.
2) Freed A.D., Zamani, S., Szabo̕, L. and Clayton, J.D., "Laplace stretch: Eulerian and Lagrangian formulations", *Zeitschrift fur angewandte Mathematik und Physik*, **71** (2020), 157. DOI: 10.1007/s00033-020-01388-4.
3) Paul, S., Freed, A.D. and Clayton, J.D., "Coordinate indexing: On the use of Eulerian and Lagrangian Laplace stretches", *Applications in Engineering Science*, **5** (2021), 100029. DOI: 10.1016/j.apples.2020.100029.
4) Freed A.D. and Clayton, J.D., *Application of Laplace Stretch In Alveolar Mechanics: A Case Study of Blast and Blunt Trauma.* In development. (This software actualizes the theories being developed in this book.)
5) Freed, A.D., "A Technical Note: Two-Step PECE Methods for Approximating Solutions To First- and Second-Order ODEs", *arXiv/2056770*, 1 Nov, 2017.

[Next](./README_1D.md)

## Version History

Removed position, velocity and acceleration vectors for centroidal motions. All fields are point fields. Added second derivatives for stretch and strain.

### Version 0.2.0

Added the position, velocity and acceleration vectors for centroidal motions.  Changed the second argument of methods *advance!* and *update!* to be measurable values, not rates, e.g., length.

### Version 0.1.8

Changed the second argument of methods *advance!* and *update!* to be differentials instead of derivatives to make the user interface more intuitive--user friendly. Removed the mid-point option.

### Version 0.1.7

Bug fix in the BDF solvers.

### Version 0.1.5

The initial configuration is chosen to have a displacement gradient that is the identity tensor. Kinematic analyses are now considered to be at a mass or continuum point. Structural information is no longer contained within these data structures. For example, what were lengths are now dimensionless stretches.

### Version 0.1.4

Allow for distinction between mid-point and end-point quadrature rules. This is handled through their constructors.

### Version 0.1.3

Removed geometric fields that were in addition to the fundamental kinematic fields. This was done to simplify data structures built upon those exported here for Laplace kinematics.

### Version 0.1.2

Initial public release.
