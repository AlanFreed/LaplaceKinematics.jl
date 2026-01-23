[Previous](./README.md)  [Next](./README_2D.md)

# Kinematics of 1D Fibers

For one-dimensional continua, their kinematic histories are stored in the data structure.
```julia
struct FiberKinematics
    # Properties of the arrays.
    dt::PF.PhysicalScalar          # time between neighboring solution nodes
    N::Int                         # number of intervals along a solution path
    n::PF.MInteger                 # a counter that ratchets from 1 to N+1

    # Reference stretch (stretch at zero strain) of a 1D fiber element.
    λᵣ::PF.PhysicalScalar          # reference stretch, λᵣ = Lᵣ/L₀, L is length

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array holding the independent variable, i.e., array of nodal times.
    t::PF.ArrayOfPhysicalScalars   # time at solution nodes

    # Arrays for position, velocity and acceleration at the fiber centroid
    # relative to an end of the fiber.
    x::PF.ArrayOfPhysicalScalars   # position at solution nodes
    v::PF.ArrayOfPhysicalScalars   # velocity at solution nodes
    a::PF.ArrayOfPhysicalScalars   # acceleration at solution nodes
    
    # Arrays for the fiber stretch and its rate.
    λ::PF.ArrayOfPhysicalScalars   # stretch at solution nodes
    λ′::PF.ArrayOfPhysicalScalars  # stretch rate at solution nodes

    # Arrays for the thermodynamic (true) strain and its rate.
    ε::PF.ArrayOfPhysicalScalars   # strain at solution nodes
    ε′::PF.ArrayOfPhysicalScalars  # strain rate at solution nodes
end
```
where types *MInteger*, *PhysicalScalar* and *ArrayOfPhysicalScalars* are exported by module *PhysicalFields*. The fields comprising this type are self explanatory.

## Internal Constructors

The constructor most likely to be used by a programmer is
```julia
k = function FiberKinematics(N::Int, 
                             dt::PF.PhysicalScalar,
                             Lᵣ::PF.PhysicalScalar,
                             L₀::PF.PhysicalScalar)
```
which returns a new data structure *k* of type *FiberKinematics* that holds kinematic fields pertinent for the modeling of a 1D fiber. Arguments are: 

1. `N` is the total number of grid points or nodes where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays, e.g., k.t[1] = 0 s. 
2. `dt` is a differential step in time. It is the time that separates neighboring nodes, which is taken to be uniformly spaced over time. dt must have physical units of *TIME*.
3. `Lᵣ` is the length of a fiber in its reference (strain-free) configuration κᵣ. Length Lᵣ must have physical units of *LENGTH*.
4. `L₀` is the initial length of a fiber in its initial configuration κ₀, i.e., the configuration at the start of an analysis. Length L₀ must have physical units of *LENGTH*.  Typically, Lᵣ < L₀ implicating that the fiber is prestressed in its initial configuration.

The constructor used by JSON3 and other external constructors is
```julia
k = FiberKinematics(dt::PF.PhysicalScalar, 
                    N::Int, 
                    n::PF.MInteger, 
                    Lᵣ::PF.PhysicalScalar, 
                    L₀::PF.PhysicalScalar,
                    t::PF.ArrayOfPhysicalScalars,
                    x::PF.ArrayOfPhysicalScalars,
                    v::PF.ArrayOfPhysicalScalars,
                    a::PF.ArrayOfPhysicalScalars,
                    λ::PF.ArrayOfPhysicalScalars,
                    λ′::PF.ArrayOfPhysicalScalars,
                    ε::PF.ArrayOfPhysicalScalars,
                    ε′::PF.ArrayOfPhysicalScalars)
```
which is a serialization of the fields comprising type *FiberKinematics*.

## Methods

### Copy

For making a copy *cc* of an object *k* of type *FiberKinematics*, use
```julia
cc = copy(k::FiberKinematics)
```

### Persistence

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
json_stream = PhysicalFields.openJSONWriter(my_dir_path::String, 
											my_file_name::String)
```
while a *json_stream* for reading from a JSON file can be created by calling
```julia
json_stream = PhysicalFields.openJSONReader(my_dir_path::String, 
											my_file_name::String)
```
where `<my_dir_path>` is the path to your working directory wherein the file to be written, i.e., `<my_file_name>`, will either exist or be created. This file name must have a `.json` extension.  Paired with these methods, a call to
```julia
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```
will close a *json_stream*.

### Solver

To advance a solution along its path, step by step, call the method
```julia
advance!(k::FiberKinematics, Lₙ::PhysicalScalar)
```
Method `advance!` moves a solution from previous step *n-1* to current step *n* along a solution path with *N* solution nodes. It is assumed that fiber length is controlled (and is therefore known as a function of time).  Argument `Lₙ` denotes fiber length in the current configuration κₙ, which must have physical units of *LENGTH*.  From this length all fields are updated, with rates and accelerations being approximated using third-order accurate, finite-difference formulæ.

This method updates counter *k.n*, plus it assigns entries to its history arrays at their nᵗʰ array location in this *k* data structure; specifically: position *k.x[n]*, velocity *k.v[n]*, acceleration *k.a[n]*, along with stretch *k.λ[n]* and its rate *k.λ′[n]*, plus strain *k.ε[n]* and its rate *k.ε′[n]*.

**Note**: It is not until step *n=6* that all tabulated derivatives in time become third-order accurate, or until step *n=8* that all tabulated second derivatives in time become third-order accurate. Derivatives at the early steps are updated with more accurate approximations as more data becomes available, until they become third-order accurate.

A solution at current node *k.n* can be refined by calling the method
```julia
update!(k::FiberKinematics, Lₙ::PhysicalScalar)
```
Method `update!` is to be called after `advance!`. It refines a solution at step *n* for an updated value of ` Lₙ` that is to be inserted into data structure `k`, e.g., updated by some external optimization process like a finite element engine.

[Previous](./README.md)  [Next](./README_2D.md)
