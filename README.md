# LaplaceKinematics.jl

This module provides a **QR** or Gram-Schmidt factorization of a deformation gradient tensor **F** = *Fᵢⱼ* **e**ᵢ ⊗ **e**ⱼ, where *i*,*j* ∈ {1,2} for 2D analyses or *i*,*j* ∈ {1,2,3} for 3D analyses, with *Fᵢⱼ* being its matrix of tensor components evaluated in basis (**e**₁, **e**₂) for 2D analyses or in basis (**e**₁, **e**₂, **e**₃) for 3D analyses. Here **Q** = Qᵢⱼ **e**ᵢ ⊗ **e**ⱼ is an orthogonal matrix known as the Gram rotation tensor, and **R** = *Rᵢⱼ* **e**ᵢ ⊗ **e**ⱼ is an upper-triangular matrix known as Laplace stretch tensor. In the mechanics literature, this decomposition is often denoted as **F** = **RU**, where **R** denotes a rotation and **U** denotes a stretch.

To use this module you will need to add the following Julia packages to yours:
```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/LaplaceKinematics.jl")
```
and to run the test files, you'll also need to install the package
```
Pkg.add(url = "https://github.com/AlanFreed/FijLung.jl")
```

## Kinematics of 1D Fibers

```
struct FiberKinematics
    # Properties of the arrays.
    dt::PhysicalScalar          # time increment separating neighboring nodes
    N::Integer                  # number of nodes to traverse a solution path
    n::MInteger                 # a counter that ratchets from 1 to N+1

    # Reference (strain free) length of a 1D fiber element.
    Lᵣ::PhysicalScalar          # reference length

    # History arrays of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t::ArrayOfPhysicalScalars   # time at the solution nodes

    # Arrays for the fiber length and its rate.
    L::ArrayOfPhysicalScalars   # length at the solution nodes
    L′::ArrayOfPhysicalScalars  # length rate at the solution nodes

    # Thermodynamic (true) strains and their rates.
    ϵ::ArrayOfPhysicalScalars   # strain at the solution nodes
    ϵ′::ArrayOfPhysicalScalars  # strain rate at the solution nodes
end
```
where types `MInteger,` `PhysicalScalar` and `ArrayOfPhysicalScalars` are exported by module `PhysicalFields.` The fields of this type are self explanatory.

### Internal Constructors

The constructor most likely to be used by a programmer is
```
function FiberKinematics(dTime::PhysicalScalar, N::Integer, midPtQuad::Bool, Lᵣ::PhysicalScalar, L₀::PhysicalScalar)
```
which returns a new data structure of type `FiberKinematics` that holds kinematic fields pertinent to the modeling of a 1D fiber. Arguments are: (*i*) A differential step in time `dTime` that separates neighboring nodes, which themselves are taken to be uniformly spaced in time. (*ii*) The total number of grid points or nodes `N` where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays. (*iii*) A boolean flag `midPtQuad` that, if true, implies the nodal spacing is for a mid-point quadrature rule and, if false, implies the nodal spacing is for an end-point quadrature rule. This determines how the array of independent times is populated. (*iv*) The reference (or strain free) length `Lᵣ` of a fiber against which strains are to be measured. And (*v*) a fiber's initial length `L₀` in some initial configuration selected for analysis κ₀ where, typically, L₀ ≥ Lᵣ.

The constructor used by JSON3 and other external constructors is
```
function FiberKinematics(dt::PhysicalScalar, N::Integer, n::MInteger, Lᵣ::PhysicalScalar, t::ArrayOfPhysicalScalars, L::ArrayOfPhysicalScalars, L′::ArrayOfPhysicalScalars,ϵ::ArrayOfPhysicalScalars, ϵ′::ArrayOfPhysicalScalars)
```
which is a serialization of the fields comprising type `FiberKinematics.`

## Methods

### Copies

For making a shallow copy of an object of type `FiberKinematics,` use
```
function Base.:(copy)(k::FiberKinematics)::FiberKinematics
```
and for making a deep copy of an object of type `FiberKinematics,`, use
```
function Base.:(deepcopy)(k::FiberKinematics)::FiberKinematics
```

### Persistence

To write an object of type `FiberKinematics` to a JSON file, one can call
```
function toFile(k::FiberKinematics, json_stream::IOStream)
```
while reading in such an object from a JSON file can be accomplished by calling
```
function fromFile(::Type{FiberKinematics}, json_stream::IOStream)::FiberKinematics
```
wherein a `json_stream` for writing to a JSON file can be created from a call to
```
function PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)::IOStream
```
while a `json_stream` for reading from a JSON file can be created by calling
```
function PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)::IOStream
```
with
```
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```
closing a `json_stream.`

### Solver

To advance a solution along its path, step by step, call the method
```
function advance!(k::FiberKinematics, L′::PhysicalScalar)
```
Method `advance!` moves a solution from previous step `n-1` to current step `n` along a solution path of N solution nodes by integrating its governing differential equation for length using a backward difference formula (BDF) when given the fiber's time rate-of-change in length `L′`. For a time-step interval of [tₙ₋₁, tₙ], L′ = dL/dt associates with either time tₙ when using end-point quadrature, or with time (tₙ₋₁ + tₙ)/2 when using mid-point quadrature.

This method updates counter `k.n` and entries to its history arrays at the nᵗʰ array location in the `k` data structure; specifically: length `k.L[n]` and its rate `k.L′[n]`, plus strain `k.ϵ[n]` and its rate `k.ϵ′[n].` These fields are evaluated at either the end-point, i.e. at time tₙ, or at the mid-point, i.e. at time (tₙ₋₁ + tₙ)/2, according to the argument `midPtQuad` supplied to its constructor.


A solution at current node `k.n` can be refined by calling the method
```
function update!(k::FiberKinematics, L′::PhysicalScalar)
```
Method `update!` refines a solution at step `n` by re-integrating its governing differential equation for fiber length, thereby allowing for iterative improvements to be made on length rate `L′` from an external algorithm, e.g., a finite element engine. There is no need to call `update!` unless `L′` is being iteratively refined at step `n`, e.g., by some external optimization process. Here L′ = dL(k.t[n])/dt.

## Laplace Kinematics for 2D Membranes

Membranes are planar structures that do not support an out-of-plane bending moment. A data structure that holds kinematic fields for a membrane described by a Gram-Schmidt deconstruction of the deformation gradient is given by
```
struct MembraneKinematics
    # Properties of the arrays.
    dt::PhysicalScalar           # time step separating neighboring entries
    N::Integer                   # total number of steps or grid points
    n::MInteger                  # a counter that ratchets from 1 to N+1

    # 2D Laplace stretch attributes for a reference deformation of κ₀ ↦ κᵣ.
    aᵣ::PhysicalScalar           # reference elongation (stretch) in 𝕚 direction
    bᵣ::PhysicalScalar           # reference elongation (stretch) in 𝕛 direction
    γᵣ::PhysicalScalar           # reference in-plane shear in (𝕚,𝕛) plane

    # History arrays of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t::ArrayOfPhysicalScalars    # times at the solution nodes, i.e., the tₙ

    # Unpivoted 2D deformation gradients for a deformation of κ₀ ↦ κₙ in (𝕚, 𝕛).
    F::ArrayOfPhysicalTensors    # deformation gradients at tₙ: Fₙ κ₀ ↦ κₙ
    F′::ArrayOfPhysicalTensors   # deformation gradient rates at tₙ: dFₙ/dtₙ
    P::Vector                    # permutation case at tₙ: i.e., i in Pᵢ,
                                 # where {𝕖₁ 𝕖₂} = {𝕚 𝕛}[Pᵢ], i ∈ {1, 2}

    # Gram angles of rotation and their rates at tₙ mapped to (𝕚, 𝕛)\n
    ωₙ::ArrayOfPhysicalScalars   # angular rotations at tₙ: ωₙ
                                 # (𝕖₁, 𝕖₂) out of (𝕚, 𝕛) whenever P = P₁
                                 # (𝕖₂, 𝕖₁) out of (𝕚, 𝕛) whenever P = P₂
    ω′ₙ::ArrayOfPhysicalScalars  # angular rates of rotation at tₙ: dωₙ/dtₙ

    # 2D Laplace stretch attributes for deformation κᵣ ↦ κₙ mapped to (𝕚, 𝕛)
    aₙ::ArrayOfPhysicalScalars   # elongations in 𝕚 direction at tₙ
    bₙ::ArrayOfPhysicalScalars   # elongations in 𝕛 direction at tₙ
    γₙ::ArrayOfPhysicalScalars   # in-plane shears in (𝕚, 𝕛) plane at tₙ

    # 2D Laplace stretch-rate attributes at configuration κₙ mapped to (𝕚, 𝕛)
    a′ₙ::ArrayOfPhysicalScalars  # elongation rates in 𝕚 direction at tₙ: daₙ/dt
    b′ₙ::ArrayOfPhysicalScalars  # elongation rates in 𝕛 direction at tₙ: dbₙ/dt
    γ′ₙ::ArrayOfPhysicalScalars  # in-plane shear rates at tₙ: dγₙ/dt

    # 2D Laplace strain attributes for deformation κᵣ ↦ κₙ mapped to (𝕚, 𝕛)
    δ::ArrayOfPhysicalScalars    # strains of dilation at tₙ: δ
    ϵ::ArrayOfPhysicalScalars    # strains of squeeze at tₙ: ϵ
    γ::ArrayOfPhysicalScalars    # strains of shear at tₙ: γ

    # 2D Laplace strain-rate attributes at configuration κₙ mapped to (𝕚, 𝕛)
    δ′::ArrayOfPhysicalScalars   # strain rates of dilation at tₙ: dδ/dt
    ϵ′::ArrayOfPhysicalScalars   # strain rates of squeeze at tₙ: dϵ/dt
    γ′::ArrayOfPhysicalScalars   # strain rates of shear at tₙ: dγ/dt
end
```
where types `MInteger,` `PhysicalScalar,` `PhysicalTensor,` `ArrayOfPhysicalScalars` and `ArrayOfPhysicalTensors` are all exported by module `PhysicalFields.`

There are four stretch attributes that describe a planar Laplace stretch at nodal time `tₙ:` elongations `aₙ` and `bₙ,` an in-plane shear `γₙ,` and a Gram rotation `ωₙ.` Gram rotations are physical. They either associate with a rigid-body rotation or a pure shear. This physicality is ensured through a co-ordinate permutation. Permutation `P₁` associates with a right-handed co-ordinate frame, while permutation `P₂` associates with a left-handed co-ordinate frame. From the above Laplace stretch attributes come three thermodynamic strains: dilation `δ,` squeeze `ϵ,` and simple shear `γ.` It is in terms of these three strains and their rates that constitutive equations are to be constructed.

### Internal Constructors

The constructor most likely to be used by a programmer is
```
function MembraneKinematics(dTime::PhysicalScalar, N::Integer, midPtQuad::Bool, aᵣ::PhysicalScalar, bᵣ::PhysicalScalar, γᵣ::PhysicalScalar, F₀::PhysicalTensor)
```
which returns a new data structure of type `MembraneKinematics` that holds kinematic fields pertinent to the modeling of a 2D membrane. Arguments are: (*i*) A differential step in time `dTime` that separates neighboring nodes, which themselves are taken to be uniformly spaced in time. (*ii*) The total number of grid points or nodes `N` where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays. (*iii*) A boolean flag `midPtQuad` that, if true, implies the nodal spacing is for a mid-point quadrature rule and, if false, implies the nodal spacing is for an end-point quadrature rule. This determines how the array of independent times is populated. (*iv*) The reference (or strain free) elongations are: `aᵣ` and `bᵣ,` along with a simple shear `γᵣ.` They associate with a membrane deformation of `κ₀ ↦ κᵣ,` with configuration `κᵣ` being that configuration against which strains are to be measured. And (*v*) a membrane's 2x2 deformation gradient `F₀` for an initial configuration selected for analysis `κ₀` whose elements are evaluated in an user specified co-ordinate system with base vectors (𝕚, 𝕛). Laplace tensors are evaluated in a frame-indifferent co-ordinate system (𝕖₁, 𝕖₂), which are then mapped to the user's co-ordinate system (𝕚, 𝕛). It is in the (𝕚, 𝕛) co-ordinate system that the kinematic fields of this data structure are quantified in.

The constructor used by JSON3 and other external constructors is
```
function MembraneKinematics(dt::PhysicalScalar, N::Integer, n::MInteger, aᵣ::PhysicalScalar, bᵣ::PhysicalScalar, γᵣ::PhysicalScalar, t::ArrayOfPhysicalScalars, F::ArrayOfPhysicalTensors, F′::ArrayOfPhysicalTensors, P::Vector, ωₙ::ArrayOfPhysicalScalars, ω′ₙ::ArrayOfPhysicalScalars, aₙ::ArrayOfPhysicalScalars, bₙ::ArrayOfPhysicalScalars, γₙ::ArrayOfPhysicalScalars, a′ₙ::ArrayOfPhysicalScalars, b′ₙ::ArrayOfPhysicalScalars, γ′ₙ::ArrayOfPhysicalScalars, δ::ArrayOfPhysicalScalars, ϵ::ArrayOfPhysicalScalars, γ::ArrayOfPhysicalScalars, δ′::ArrayOfPhysicalScalars, ϵ′::ArrayOfPhysicalScalars, γ′::ArrayOfPhysicalScalars)
```
which is a serialization of the fields comprising type `MembraneKinematics.`

## Methods

### Copies

For making a shallow copy of an object of type `MembraneKinematics,` use
```
function Base.:(copy)(k::MembraneKinematics)::MembraneKinematics
```
and for making a deep copy of an object of type `MembraneKinematics,`, use
```
function Base.:(deepcopy)(k::MembraneKinematics)::MembraneKinematics
```

### Persistence

To write an object of type `MembraneKinematics` to a JSON file, one can call
```
function toFile(k::MembraneKinematics, json_stream::IOStream)
```
while reading in such an object from a JSON file can be accomplished by calling
```
function fromFile(::Type{MembraneKinematics}, json_stream::IOStream)::MembraneKinematics
```
wherein a `json_stream` for writing to a JSON file can be created from a call to
```
function PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)::IOStream
```
while a `json_stream` for reading from a JSON file can be created by calling
```
function PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)::IOStream
```
with
```
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```
closing a `json_stream.`

### Solver

To advance a solution along its path, step by step, call the method
```
function advance!(k::MembraneKinematics, F′::PhysicalTensor)
```
Method `advance!` moves a solution from a previous step `n-1` to its current step `n` along some solution path containing N+1 nodes. This is accomplished by integrating the supplied rate of deformation gradient `F′` using a backward difference formula (BDF) from which all other fields are then derived and assigned to their appropriate arrays. For a time-step interval of [tₙ₋₁, tₙ], F′ = dF/dt associates with either time tₙ when using end-point quadrature, or with time (tₙ₋₁ + tₙ)/2 when using mid-point quadrature.

This method advances counter `k.n` plus it assigns values to its history arrays located at the nᵗʰ array entries in the `k` data structure. Specifically, assignments are made to: the deformation gradient `k.F[n]` and its rate `k.F′[n],` the pivot case `k.P[n],` the Laplace attributes for stretch `k.aₙ[n],` `k.bₙ[n],` `k.γₙ[n]` and `k.ωₙ[n]` and their rates `k.a′ₙ[n],` `k.b′ₙ[n],` `k.γ′ₙ[n]` and `k.ω′ₙ[n],` plus the Laplace attributes for strain `k.δ[n],` `k.ϵ[n]` and `k.γ[n]` and their rates `k.δ′[n],` `k.ϵ′[n]` and `k.γ′[n].` All of these fields have values that associate with the user's co-ordinate system, whose base vectors are denoted as (𝕚, 𝕛). These fields are evaluated at either the end-point, i.e. at time tₙ, or at the mid-point, i.e. at time (tₙ₋₁ + tₙ)/2, according to the argument `midPtQuad` supplied to its constructor.

A solution at current node `k.n` can be refined by calling the method
```
function update!(k::MembraneKinematics, F′::PhysicalTensor)
```
Such a refinement is accomplished by re-integrating its governing differential equation for the 2D deformation gradient, thereby allowing for iterative improvements to be made throughout the data structure. There is no need to call `update!` unless `F′` is being iteratively refined at step `n,` e.g., through some external optimization process like a finite element engine. Here F′ = dF(k.t[n])/dt.

## Laplace Kinematics for 3D Materials

Not implemented yet

## References

1) Freed, A.D., Erel, V. and Moreno, M.R., "Conjugate stress/strain base pairs for planar analysis of biological tissues", Journal of Mechanics of
   Materials and Structures, 12 (2017), 219-247. DOI: 10.2140/jomms.2017.12.219.
2) Freed A.D., Zamani, S., Szabo̕, L. and Clayton, J.D., "Laplace stretch:
   Eulerian and Lagrangian formulations", Zeitschrift fur angewandte Mathematik und Physik, 71 (2020), 157. DOI: 10.1007/s00033-020-01388-4.
3) Paul, S., Freed, A.D. and Clayton, J.D., "Coordinate indexing: On the use of Eulerian and Lagrangian Laplace stretches", Applications in Engineering
   Science, 5 (2021), 100029. DOI: 10.1016/j.apples.2020.100029.

## Version History

0.1.4: Allow for distinction between mid-point and end-point quadrature rules. This is handled through their constructors.

0.1.3: Removed geometric fields that were in addition to the fundamental kinematic fields. This was done to simplify data structures built upon those exported here for Laplace kinematics.

0.1.2: Initial public release.