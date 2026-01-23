[Previous](./README_1D.md)  [Next](./README_3D.md)

# Laplace Kinematics for 2D Membranes

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

[Previous](./README_1D.md)  [Next](./README_3D.md)
