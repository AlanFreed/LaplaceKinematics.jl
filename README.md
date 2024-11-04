# LaplaceKinematics.jl

This module provides a **QR** or Gram-Schmidt factorization of a deformation gradient tensor denoted as **F** = *Fᵢⱼ* **e**ᵢ ⊗ **e**ⱼ, where *i*,*j* ∈ {1,2} for 2D analyses, or *i*,*j* ∈ {1,2,3} for 3D analyses, with *Fᵢⱼ* being its matrix of tensor components evaluated in basis (**e**₁, **e**₂) for 2D analyses, or in basis (**e**₁, **e**₂, **e**₃) for 3D analyses. Here **Q** = Qᵢⱼ **e**ᵢ ⊗ **e**ⱼ is an orthogonal matrix known as the Gram rotation tensor, and **R** = *Rᵢⱼ* **e**ᵢ ⊗ **e**ⱼ is an upper-triangular matrix (the right matrix in product **QR**) called the Laplace stretch tensor, as Laplace introduced the technique employed by Gram in his Ph.D. thesis that Schmidt made popular.

In the mechanics literature, such decompositions are often denoted as **F** = **RU**, where **R** denotes a rotation and **U** denotes its stretch from a Lagrangian perspective. There are an infinite number of such matrix products. The deconstruction of a matrix into a product of two matrices of which one is orthogonal is not unique. We are interested in the one wherein **R** is a Gram rotation matrix in what is otherwise commonly known as a **QR** matrix decomposition in the linear algebra literature.

A deformation gradient is a mapping from one material configuration into another. It is a two-state property, independent of the path traversed between these two states. All deformation gradients considered herein satisfy the following conditions: 

1) The deformation gradient **F** equates with the identity tensor **I** = *δᵢⱼ* **e**ᵢ ⊗ **e**ⱼ in its initial configuration κ₀, i.e., its motion is κ₀ ↦ κ₀. 
2) All other deformation gradients associate with motions of κ₀ ↦ κₙ, where κₙ is the configuration at current time tₙ. This includes a strain-free reference configuration κᵣ, which can be attained through a motion of κ₀ ↦ κᵣ.
3) The reference configuration κᵣ is, by definition, free from strain; therefore, the initial configuration κ₀ can be a strained configuration caused by motion κᵣ ↦ κ₀, even though the initial deformation gradient is **F**(t₀) = **I**. Although seldom adopted, this ought to be the norm for analyses of biologic tissue.

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

## Kinematics of 1D Fibers

For one dimensional continua, their kinematic histories are stored in the data structure.
```julia
struct FiberKinematics
    # Properties of the arrays.
    dt::PhysicalScalar          # time increment separating neighboring nodes
    N::Int                      # number of intervals along a solution path
    n::MInteger                 # a counter that ratchets from 1 to N+1

    # Reference (strain free) length of a 1D fiber element.
    λᵣ::PhysicalScalar          # reference stretch, λᵣ = Lᵣ / L₀, L is length

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t::ArrayOfPhysicalScalars   # time at the solution nodes

    # Arrays for the fiber stretch and its rate.
    λ::ArrayOfPhysicalScalars   # stretch at the solution nodes
    λ′::ArrayOfPhysicalScalars  # stretch rate at the solution nodes

    # Arrays for the thermodynamic (true) strain and its rate.
    ε::ArrayOfPhysicalScalars   # strain at the solution nodes
    ε′::ArrayOfPhysicalScalars  # strain rate at the solution nodes
end
```
where types MInteger, PhysicalScalar and ArrayOfPhysicalScalars are exported by module PhysicalFields. The fields comprising this type are self explanatory.

### Internal Constructors

The constructor most likely to be used by a programmer is
```julia
k = FiberKinematics(dt::PhysicalScalar, N::Int,  λᵣ::PhysicalScalar))
```
which returns a new data structure *k* of type *FiberKinematics* that holds kinematic fields pertinent for the modeling of a 1D fiber. Arguments are: 

1) A differential step in time *dt*, with units of time, that separates neighboring nodes, which themselves are taken to be uniformly spaced over time. 
2) The total number of grid points or nodes *N* where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays, e.g., t[1] = 0. 
3) The reference (or strain free) stretch *λᵣ* of a fiber against which strains are to be measured, viz., ε(λᵣ) = 0. The fiber's initial stretch λ₀, associated with some initial configuration κ₀, is assigned a value of 1, with an outcome being that ε(λ₀) need not equal 0.

The constructor used by JSON3 and other external constructors is
```julia
k = FiberKinematics(dt::PhysicalScalar, N::Int, n::MInteger, λᵣ::PhysicalScalar, t::ArrayOfPhysicalScalars, λ::ArrayOfPhysicalScalars, λ′::ArrayOfPhysicalScalars, ε::ArrayOfPhysicalScalars, ε′::ArrayOfPhysicalScalars)
```
which is a serialization of the fields comprising type *FiberKinematics*.

### Methods

#### Copy

For making a copy *cc* of an object *k* of type *FiberKinematics*, use
```julia
cc = copy(k::FiberKinematics)
```

#### Persistence

To write an object *k* of type *FiberKinematics* to a JSON file, one can call
```julia
toFile(k::FiberKinematics, json_stream::IOStream)
```
while reading in such an object from a JSON file can be accomplished by calling
```julia
k = fromFile(FiberKinematics, json_stream::IOStream)
```
wherein a *json_stream* for writing to a JSON file can be created from a call to
```julia
json_stream =  PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)
```
while a *json_stream* for reading from a JSON file can be created by calling
```julia
json_stream = PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)
```
with
```julia
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```
closing a *json_stream*.

#### Solver

To advance a solution along its path, step by step, call the method
```julia
advance!(k::FiberKinematics, dλ::PhysicalScalar)
```
Method *advance!* moves a solution from previous step n-1 to current step n along a solution path of N solution nodes by advancing `λ(tₙ) = λ(tₙ₋₁) + dλ` over the time interval [tₙ₋₁, tₙ] subject to an initial condition of `λ(t₀) = 1.` Derivatives of this stretch are then approximated via third-order finite difference formula. From this differential change in stretch, the stretch λ, stretch rate λ′ = dλ/dt, strain ε and strain rate ε′ = dε/dt are all evaluated at time tₙ, and then stored in object *k*.

**Note**: It is not until step n=6 that all tabulated rates are third-order accurate. Given that node n=1 holds the initial condition, then at the first step, i.e., n=2, rates λ′ affiliated with nodes n=1,2 will be first-order accurate. At the second step, i.e., n=3, rates λ′ affiliated with nodes n=1,2,3 will become second-order accurate. It is not until the fifth step, i.e., n=6, and thereafter that the rates λ′ affiliated with nodes n=1,2,…,6,… will all become third-order accurate.

A solution at current node *k.n* can be refined by calling the method
```julia
update!(k::FiberKinematics, dλ::PhysicalScalar)
```
Method *update!* is to be called after *advance!*. It refines a solution at step n for an updated advancement in stretch of *dλ = λ(tₙ) - λ(tₙ₋₁)* over the previous time interval [tₙ₋₁, tₙ]. There is no need to call *update!* unless *dλ* is being iteratively refined at step n, e.g., by some external optimization process like a finite element engine.

## Laplace Kinematics for 2D Membranes

Membranes are planar structures that do not support an out-of-plane bending moment. The user's base vectors are denoted as (𝕚, 𝕛) that when pivoted to ensure a physical interpretation for shears become (**e**₁, **e**₂). A data structure that holds kinematic fields for a membrane described by a Gram-Schmidt deconstruction of the deformation gradient **F** = *Fᵢⱼ* **e**ᵢ ⊗ **e**ⱼ is given by
```julia
struct MembraneKinematics
    # Properties of the arrays.
    dt::PhysicalScalar           # time step separating neighboring entries
    N::Int                       # total number of steps or grid points
    n::MInteger                  # a counter that ratchets from 1 to N+1

    # 2D Laplace stretch attributes for a reference deformation of κ₀ ↦ κᵣ.
    aᵣ::PhysicalScalar           # reference elongation (stretch) in 𝕚 direction
    bᵣ::PhysicalScalar           # reference elongation (stretch) in 𝕛 direction
    γᵣ::PhysicalScalar           # reference in-plane shear in (𝕚,𝕛) plane

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t::ArrayOfPhysicalScalars    # time at the solution nodes, i.e., at the tₙ

    # Unpivoted 2D deformation gradients for a deformation of κ₀ ↦ κₙ in (𝕚, 𝕛).
    F::ArrayOfPhysicalTensors    # deformation gradients at tₙ: Fₙ κ₀ ↦ κₙ
    F′::ArrayOfPhysicalTensors   # deformation gradient rates at tₙ: dFₙ/dtₙ
    motion::Vector{Int64}        # the motion case that applies at time tₙ:
                                 # 1) with pure shear, no co-ordinate pivoting
                                 # 2) with pure shear and co-ordinate pivoting
                                 # 3) with rigid-body rotation, no pivoting
                                 # 4) with rigid-body rotation and pivoting

    # Gram angles of rotation and their rates at tₙ, mapped to (𝕚, 𝕛)
    ωₙ::ArrayOfPhysicalScalars   # angular rotations at tₙ: ωₙ
                                 # (𝕖₁, 𝕖₂) out of (𝕚, 𝕛) whenever P = P₁
                                 # (𝕖₂, 𝕖₁) out of (𝕚, 𝕛) whenever P = P₂
    ω′ₙ::ArrayOfPhysicalScalars  # angular rates of rotation (spin) at tₙ: dωₙ/dtₙ

    # 2D Laplace stretch attributes for deformation κᵣ ↦ κₙ, mapped to (𝕚, 𝕛)
    aₙ::ArrayOfPhysicalScalars   # elongations in 𝕚 direction at tₙ
    bₙ::ArrayOfPhysicalScalars   # elongations in 𝕛 direction at tₙ
    γₙ::ArrayOfPhysicalScalars   # in-plane shears in (𝕚, 𝕛) plane at tₙ

    # 2D Laplace stretch-rate attributes at configuration κₙ, mapped to (𝕚, 𝕛)
    a′ₙ::ArrayOfPhysicalScalars  # elongation rates in 𝕚 direction at tₙ: daₙ/dt
    b′ₙ::ArrayOfPhysicalScalars  # elongation rates in 𝕛 direction at tₙ: dbₙ/dt
    γ′ₙ::ArrayOfPhysicalScalars  # in-plane shear rates at tₙ: dγₙ/dt

    # 2D Laplace strain attributes for deformation κᵣ ↦ κₙ, mapped to (𝕚, 𝕛)
    δ::ArrayOfPhysicalScalars    # strains of dilation at tₙ: δ
    ε::ArrayOfPhysicalScalars    # strains of squeeze at tₙ: ε
    γ::ArrayOfPhysicalScalars    # strains of shear at tₙ: γ

    # 2D Laplace strain-rate attributes at configuration κₙ, mapped to (𝕚, 𝕛)
    δ′::ArrayOfPhysicalScalars   # strain rates of dilation at tₙ: dδ/dt
    ε′::ArrayOfPhysicalScalars   # strain rates of squeeze at tₙ: dε/dt
    γ′::ArrayOfPhysicalScalars   # strain rates of shear at tₙ: dγ/dt
end
```
where types *MInteger*, *PhysicalScalar*, PhysicalTensor, *ArrayOfPhysicalScalars* and *ArrayOfPhysicalTensors* are all exported by module *PhysicalFields*.

There are four stretch attributes that describe a planar Laplace stretch at nodal time *tₙ*: elongations *aₙ* and *bₙ*, an in-plane shear *γₙ*, and a Gram rotation *ωₙ*. Gram rotations are physical. They either associate with a rigid-body rotation or a pure shear. This physicality is ensured through a co-ordinate permutation. Permutation *P₁* associates with a right-handed co-ordinate frame, while permutation *P₂* associates with a left-handed co-ordinate frame. From the above Laplace stretch attributes come three thermodynamic strains: dilation *δ*, squeeze *ε*, and simple shear *γ*. It is in terms of these three strains and their rates that constitutive equations are to be constructed.

### Internal Constructors

The constructor most likely to be used by a programmer is
```julia
k = MembraneKinematics(dt::PhysicalScalar, N::Int, aᵣ::PhysicalScalar, bᵣ::PhysicalScalar, γᵣ::PhysicalScalar, Pᵣ::Int)
```
which returns a new data structure *k* of type *MembraneKinematics* that holds a variety of kinematic fields. Arguments include: 

1) A differential step in time *dt* that separates neighboring nodes, which themselves are taken to be uniformly spaced over time. 
2) The number of grid points or nodes *N* where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays, e.g., `t[1] = 0`. 
3) The reference Laplace stretch attributes, viz., *aᵣ*, *bᵣ* and *γᵣ*, against which isochoric strains are to be established so that `ε(aᵣ, bᵣ, γᵣ) = 0`. The initial deformation gradient **F**₀ is associated with some initial configuration κ₀. It is assigned the identity matrix **I** with an outcome being that strain ε(a₀, b₀, γ₀) need not equal 0. 
4) If γᵣ is to be a shearing in the 𝕚 direction then *Pᵣ* is to equal 1, else if γᵣ is to be a shearing in the 𝕛 direction then *Pᵣ* is to equal 2, where *Pᵣ* denotes which permutation matrix it to be applied in the reference configuration.

The constructor used by JSON3 and other external constructors is
```julia
function MembraneKinematics(dt::PhysicalScalar, N::Int, n::MInteger, aᵣ::PhysicalScalar, bᵣ::PhysicalScalar, γᵣ::PhysicalScalar, t::ArrayOfPhysicalScalars, F::ArrayOfPhysicalTensors, F′::ArrayOfPhysicalTensors, motion::Vector{Int}, ωₙ::ArrayOfPhysicalScalars, ω′ₙ::ArrayOfPhysicalScalars, aₙ::ArrayOfPhysicalScalars, bₙ::ArrayOfPhysicalScalars, γₙ::ArrayOfPhysicalScalars, a′ₙ::ArrayOfPhysicalScalars, b′ₙ::ArrayOfPhysicalScalars, γ′ₙ::ArrayOfPhysicalScalars, δ::ArrayOfPhysicalScalars, ε::ArrayOfPhysicalScalars, γ::ArrayOfPhysicalScalars, δ′::ArrayOfPhysicalScalars, ε′::ArrayOfPhysicalScalars, γ′::ArrayOfPhysicalScalars)
```
which is a serialization of the fields comprising an instance of type `MembraneKinematics.`

### Methods

#### Copy

For making a copy *cc* of an object *k* of type *MembraneKinematics*:
```julia
cc = copy(k::MembraneKinematics)
```


#### Persistence

To write an object of type `MembraneKinematics` to a JSON file, one can call
```julia
toFile(k::MembraneKinematics, json_stream::IOStream)
```
while reading in such an object from a JSON file can be accomplished by calling
```julia
k = fromFile(MembraneKinematics, json_stream::IOStream)
```
wherein a `json_stream` for writing to a JSON file can be created from a call to
```julia
json_stream = PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)
```
while a `json_stream` for reading from a JSON file can be created by calling
```julia
json_stream = PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)
```
with
```julia
PhysicalFields.closeJSONStream(json_stream::IOStream)
```
closing a `json_stream.`

#### Solver

To advance a solution along its path, step by step, call the procedure
```julia
advance!(k::MembraneKinematics, dF::PhysicalTensor)
```
Method *advance!* moves a solution from previous step n-1 to current step n along a solution path with N solution nodes by advancing `F(tₙ) = F(tₙ₋₁) + dF` over the time interval [tₙ₋₁, tₙ] subject to an initial condition of `F(t₀) = I.` Derivatives of the deformation gradient are then approximated via third-order finite difference formula. From these quantities, this method determines the stretch attributes a, b, γ and Gram rotation ω of a membrane, along with their rates a′ = da/dt, b′ = db/dt,  γ′ = dγ/dt and ω′ = dω/dt, plus the strain attributes δ, ε, γ and their rates δ′ = dδ/dt, ε′ = dε/dt, γ′ = dγ/dt. These are all evaluated at time tₙ, and then stored in object *k*. All tabulated fields have been mapped to the user's co-ordinate system whose base vectors are denoted as (𝕚, 𝕛). 

**Note**: The stretch attribute for shear is distinguished from its strain attribute in that the former have subscripts, e.g., *k.γₙ*, while the latter does not, viz., *k.γ*.

**Note**: It is not until step n=6 that all tabulated rates are third-order accurate. Given that node n=1 holds the initial condition, then at the first step, i.e., n=2, rates affiliated with nodes n=1,2 will be first-order accurate. At the second step, i.e., n=3, rates affiliated with nodes n=1,2,3 will become second-order accurate. It is not until the fifth step, i.e., n=6, and thereafter that the rates affiliated with nodes n=1,2,…,6,… will all become third-order accurate.

A solution at current node `k.n` can be refined by calling the method
```julia
update!(k::MembraneKinematics, dF::PhysicalTensor)
```
Such a refinement is accomplished by re-solving all the kinematic fields in *k* associated with step n according to an updated expression for the deformation gradient difference *dF*, thereby allowing for iterative improvements to be made on the deformation change *dF* from an external algorithm, e.g., a finite element engine. There is no need to call *update!* unless *dF* is being iteratively refined at step n by some external optimization process.

## Laplace Kinematics for 3D Materials

Not implemented yet

## References

1) Freed, A.D., Erel, V. and Moreno, M.R., "Conjugate stress/strain base pairs for planar analysis of biological tissues", *Journal of Mechanics of Materials and Structures*, **12** (2017), 219-247. DOI: 10.2140/jomms.2017.12.219.
2) Freed A.D., Zamani, S., Szabo̕, L. and Clayton, J.D., "Laplace stretch: Eulerian and Lagrangian formulations", *Zeitschrift fur angewandte Mathematik und Physik*, **71** (2020), 157. DOI: 10.1007/s00033-020-01388-4.
3) Paul, S., Freed, A.D. and Clayton, J.D., "Coordinate indexing: On the use of Eulerian and Lagrangian Laplace stretches", *Applications in Engineering Science*, **5** (2021), 100029. DOI: 10.1016/j.apples.2020.100029.
4) Freed A.D. and Clayton, J.D., *Application of Laplace Stretch In Alveolar Mechanics: A Case Study of Blast and Blunt Trauma.* In development. (This software actualizes the theories being developed in this book.)
5) Freed, A.D., "A Technical Note: Two-Step PECE Methods for Approximating Solutions To First- and Second-Order ODEs", *arXiv/2056770*, 1 Nov, 2017.

## Version History

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