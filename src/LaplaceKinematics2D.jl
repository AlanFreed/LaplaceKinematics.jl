# Matrices P2D₁ and P2D₂ are the two possible permutation matrices in 2-space.

one  = PF.PhysicalScalar(1.0, DIMENSIONLESS)
P2D₁ = PF.PhysicalTensor(2, 2, DIMENSIONLESS)  # (𝕚, 𝕛) ↦ (𝕖₁, 𝕖₂)
P2D₁[1,1] = one
P2D₁[2,2] = one
P2D₂ = PF.PhysicalTensor(2, 2, DIMENSIONLESS)  # (𝕛, 𝕚) ↦ (𝕖₁, 𝕖₂)
P2D₂[1,2] = one
P2D₂[2,1] = one

# -----------------------------------------------------------------------------

"""
# MembraneKinematics

This documentation uses aliases PF for PhysicalFields and LK for LaplaceKinematics.

The fields for an object of type *MembraneKinematics* are:

    # Properties of the arrays.
    dt      time increment separating neighboring nodes
    N       total node count for traversing a solution path
    n       a counter that ratchets from 1 to N+1

    # 2D Laplace stretch attributes for a reference deformation of κ₀ ↦ κᵣ.
    aᵣ      reference elongation (stretch) in the 𝕚 direction
    bᵣ      reference elongation (stretch) in the 𝕛 direction
    γᵣ      reference in-plane shear in (𝕚, 𝕛) plane in the 𝕚 direction

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of nodal times.
    t       times at the solution nodes, i.e., the tₙ

    # Unpivoted 2D deformation gradients for deformation κ₀ ↦ κₙ in (𝕚, 𝕛).
    F       deformation gradients at tₙ: Fₙ, κ₀ ↦ κₙ in (𝕚, 𝕛)
    F′      deformation gradient rates at tₙ: dFₙ/dtₙ, κₙ in (𝕚, 𝕛)
    motion  the motion case that applies at time tₙ:
                1) with pure shear, no co-ordinate pivoting
                2) with pure shear and co-ordinate pivoting
                3) with rigid-body rotation, no pivoting
                4) with rigid-body rotation and pivoting

    # Gram angles of rotation and their rates at tₙ, mapped to (𝕚, 𝕛)
    ωₙ      angular rotations ωₙ at tₙ:
                (𝕖₁, 𝕖₂) out of (𝕚, 𝕛) whenever P = P₁
                (𝕖₂, 𝕖₁) out of (𝕚, 𝕛) whenever P = P₂
    ω′ₙ     angular rates of rotation at tₙ, i.e., dωₙ/dtₙ

    # 2D Laplace stretch attributes for κ₀ ↦ κₙ, mapped to (𝕚, 𝕛)
    aₙ      elongations in the 𝕚 direction at tₙ
    bₙ      elongations in the 𝕛 direction at tₙ
    γₙ      in-plane shears in (𝕚, 𝕛) plane in the 𝕚 direction at tₙ

    # 2D Laplace stretch-rate attributes at κₙ, mapped to (𝕚, 𝕛)
    a′ₙ     elongation rates in the 𝕚 direction at tₙ: daₙ/dt
    b′ₙ     elongation rates in the 𝕛 direction at tₙ: dbₙ/dt
    γ′ₙ     in-plane shear rates at tₙ in (𝕚, 𝕛) plane in 𝕚 direction: dγₙ/dt

    # 2D Laplace strain attributes for deformation κᵣ ↦ κₙ
    δ       strains of dilation at tₙ
    ε       strains of squeeze at tₙ
    γ       strains of shear at tₙ

    # 2D Laplace strain-rate attributes at configuration κₙ
    δ′      strain rates of dilation at tₙ, viz., dδ/dt
    ε′      strain rates of squeeze at tₙ, viz., dε/dt
    γ′      strain rates of shear at tₙ, viz., dγ/dt
    
*MembraneKinematics* is a data structure that contains the various Lagrangian fields associated with a Laplace (upper triangular) measure for stretch in an isochoric two-space. The arrays that comprise this data structure allow for a history of these kinematic variables to be compiled for later retrieval and use, e.g., for constitutive analysis, for graphing, etc. Fields of this type are evaluated in the user's co-ordinate frame (𝕚, 𝕛).

**Beware**: The prime character is used to denote a rate is the Unicode codepoint U+2032, i.e., ′, not the apostropy ', which is Unicode codepoint U+0027. They look the same in some font sets.

## Constructors

The constructor most likely to be used.
```julia
k = MembraneKinematics(dt, N, aᵣ, bᵣ, γᵣ, Pᵣ)
```
Returns a new data structure *k* of type *MembraneKinematics* that holds a variety of kinematic fields. Arguments include: 

1) A differential step in time *dt* that separates neighboring nodes, which themselves are taken to be uniformly spaced over time.
2) The number of grid points or nodes *N* where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays, e.g., `t[1] = 0`.
3) The reference Laplace stretch attributes, viz., *aᵣ*, *bᵣ* and *γᵣ*, against which isochoric strains are to be established so that `ε(aᵣ, bᵣ, γᵣ) = 0`, with the membrane's initial deformation gradient **F**₀, associated with some initial configuration κ₀, being assigned the identity matrix **I** with an outcome being that ε(a₀, b₀, γ₀) need not equal 0. 
4) If γᵣ is to be a shearing in the 𝕚 direction then `Pᵣ` is to equal 1; otherwise, if γᵣ is to be a shearing in the 𝕛 direction then `Pᵣ` is to equal 2. Argument Pᵣ denotes which permutation matrix is to be applied in the reference configuration.

There is also a constructor that is used by JSON3, it being
```julia
k = MembraneKinematics(dt, N, n, aᵣ, bᵣ, γᵣ, t, F, F′, motion, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ε, γ, δ′, ε′, γ′)
```

## Methods

```julia
cc = copy(k::LK.MembraneKinematics)
```
returns a copy `cc` of object `k` of type `MembraneKinematics`.

```julia
toFile(k::LK.MembraneKinematics, json_stream::IOStream)
```
writes data structure `k` to an IOStream `json_stream`.

```julia 
k = fromFile(LK.MembraneKinematics, json_stream::IOStream)
```
reads data structure `k` from an IOStream `json_stream`.

To manage a `json_stream` for writing, consider the code fragment:

    json_stream = PF.openJSONWriter(<my_dir_path>, <my_file_name>)
    LK.toFile(k, json_stream)
    PF.closeJSONStream(json_stream)
    
while to manage a `json_stream` for reading, consider the code fragment:

    json_stream = PF.openJSONReader(<my_dir_path>, <my_file_name>)
    k = LK.fromFile(LK.MembraneKinematics, json_stream)
    PF.closeJSONStream(json_stream)

where `<my_dir_path>` is the path to your working directory wherein the file to be written, i.e., `<my_file_name>`, either exists or will be created. This file must have a .json extension.

```julia
advance!(k::MembraneKinematics, dF::PhysicalTensor)
```
Method `advance!` moves a solution from previous step `n-1` to current step `n` along a solution path with `N` solution nodes. It is assumed that the deformation gradient is controlled (and is therefore known as a function of time). Argument `dF` denotes a differential, not a derivative. From these differentials, deformation gradient rates are computed via third-order, finite-difference formulæ from which strain rates are then established.

This method updates counter `k.n`, plus those entries to its history arrays that are at the nᵗʰ array location in this `k` data structure; specifically: deformation gradient `k.F[n]` and its rate `k.F′[n]`, Gram rotation `k.ωₙ[n]` and its rate `k.ω'ₙ[n]`, Laplace stretch attributes `k.aₙ[n]`, `k.bₙ[n]` and `k.γₙ[n]` plus their rates `k.a'ₙ[n]`, `k.b'ₙ[n]` and `k.γ'ₙ[n]`, and Laplace strains `k.δ[n]`, `k.ε[n]` and `k.γ[n]` plus their rates `k.δ'[n]`, `k.ε'[n]` and `k.γ'[n]`. All are mapped into the user's co-ordinate system whose base vectors are denoted as (𝕚, 𝕛).

**Note**: Deformation gradients are assigned at the end points of each solution interval; consequently, one is to send `dF = F(tₙ) - F(tₙ₋₁)` for `n=1,2,…,N`.
    
```julia
update!(k::MembraneKinematics, dF::PhysicalTensor)
```
Method `update!` refines a solution at step `n` whenever an improvement can be made for the deformation gradient differential `dF` through an external algorithm, e.g., a finite element engine. There is no need to call `update!` unless `dF` is being iteratively refined at a global step `n` via, say, some external optimization process.
"""
struct MembraneKinematics
    # Properties of the arrays.
    dt::PF.PhysicalScalar           # time step separating neighboring entries
    N::Int                          # total number of steps or grid points
    n::PF.MInteger                  # a counter that ratchets from 1 to N+1

    # 2D Laplace stretch attributes for a reference deformation of κ₀ ↦ κᵣ.
    aᵣ::PF.PhysicalScalar           # reference elongation in the 𝕚 direction
    bᵣ::PF.PhysicalScalar           # reference elongation in the 𝕛 direction
    γᵣ::PF.PhysicalScalar           # reference in-plane shear in (𝕚,𝕛) plane

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t::PF.ArrayOfPhysicalScalars    # time at the solution nodes, i.e., at tₙ

    # Unpivoted 2D deformation gradients for a deformation of κ₀ ↦ κₙ in (𝕚, 𝕛).
    F::PF.ArrayOfPhysicalTensors    # deformation gradients at tₙ: Fₙ κ₀ ↦ κₙ
    F′::PF.ArrayOfPhysicalTensors   # deformation gradient rates at tₙ: dFₙ/dtₙ
    motion::Vector{Int}             # the motion case that applies at time tₙ:
                                    # 1) pure shear, no co-ordinate pivoting
                                    # 2) pure shear and co-ordinate pivoting
                                    # 3) rigid-body rotation, no pivoting
                                    # 4) rigid-body rotation and pivoting

    # Gram angles of rotation and their rates at tₙ, mapped to (𝕚, 𝕛)\n
    ωₙ::PF.ArrayOfPhysicalScalars   # angular rotations at tₙ: ωₙ
                                    # (𝕖₁, 𝕖₂) out of (𝕚, 𝕛) whenever P = P₁
                                    # (𝕖₂, 𝕖₁) out of (𝕚, 𝕛) whenever P = P₂
    ω′ₙ::PF.ArrayOfPhysicalScalars  # angular rates of rotation at tₙ: dωₙ/dtₙ

    # 2D Laplace stretch attributes for deformation κᵣ ↦ κₙ, mapped to (𝕚, 𝕛)
    aₙ::PF.ArrayOfPhysicalScalars   # elongations in the 𝕚 direction at tₙ
    bₙ::PF.ArrayOfPhysicalScalars   # elongations in the 𝕛 direction at tₙ
    γₙ::PF.ArrayOfPhysicalScalars   # in-plane shears in (𝕚, 𝕛) plane at tₙ

    # 2D Laplace stretch-rate attributes at configuration κₙ, mapped to (𝕚, 𝕛)
    a′ₙ::PF.ArrayOfPhysicalScalars  # elongation rates in 𝕚 dir at tₙ: daₙ/dt
    b′ₙ::PF.ArrayOfPhysicalScalars  # elongation rates in 𝕛 dir at tₙ: dbₙ/dt
    γ′ₙ::PF.ArrayOfPhysicalScalars  # in-plane shear rates at tₙ: dγₙ/dt

    # 2D Laplace strain attributes for deformation κᵣ ↦ κₙ, mapped to (𝕚, 𝕛)
    δ::PF.ArrayOfPhysicalScalars    # strains of dilation at tₙ
    ε::PF.ArrayOfPhysicalScalars    # strains of squeeze at tₙ
    γ::PF.ArrayOfPhysicalScalars    # strains of shear at tₙ

    # 2D Laplace strain-rate attributes at configuration κₙ, mapped to (𝕚, 𝕛)
    δ′::PF.ArrayOfPhysicalScalars   # strain rates of dilation at tₙ: dδ/dt
    ε′::PF.ArrayOfPhysicalScalars   # strain rates of squeeze at tₙ: dε/dt
    γ′::PF.ArrayOfPhysicalScalars   # strain rates of shear at tₙ: dγ/dt

    # Internal constructors.

    function MembraneKinematics(dt::PF.PhysicalScalar, N::Int, aᵣ::PF.PhysicalScalar, bᵣ::PF.PhysicalScalar, γᵣ::PF.PhysicalScalar, Pᵣ::Int)

        # Verify inputs.
        if (Pᵣ < 1) || (Pᵣ > 2)
            error("Permutation case for reference configuration must be 1 or 2.")
        end

        # Convert all passed variables to CGS units.
        dt = PF.toCGS(dt)
        if Pᵣ == 1
            aᵣ = PF.toCGS(aᵣ)
            bᵣ = PF.toCGS(bᵣ)
        else
            ar = aᵣ
            br = bᵣ
            aᵣ = PF.toCGS(br)
            bᵣ = PF.toCGS(ar)
        end
        gᵣ = PF.toCGS(γᵣ)

        # Continue verification.
        if dt.units ≠ TIME
            error("The supplied time increment dt does not have units of time.")
        end
        dtₘᵢₙ = PF.PhysicalScalar(Float64(eps(Float64)), TIME)
        if dt < dtₘᵢₙ
            error("The supplied time increment dt must be positive valued.")
        end
        if N < 1
            error("Solution arrays must have a positive length.")
        end
        if !isDimensionless(aᵣ)
            error("The supplied reference stretch aᵣ is not dimensionless.")
        end
        λₘᵢₙ = PF.PhysicalScalar(Float64(eps(Float16)), DIMENSIONLESS)
        if aᵣ < λₘᵢₙ
            error("The supplied reference stretch aᵣ must be positive valued.")
        end
        if !isDimensionless(bᵣ)
            error("The supplied reference stretch bᵣ is not dimensionless.")
        end
        if bᵣ < λₘᵢₙ
            error("The supplied reference stretch bᵣ must be positive valued.")
        end
        if !isDimensionless(gᵣ)
            error("Supplied reference in-plane shear γᵣ is not dimensionless.")
        end

        # Establish the counter.
        n  = PF.MInteger(1)

        # Create and populate an array for nodal times.
        t  = PF.ArrayOfPhysicalScalars(N+1, TIME)
        for n in 1:N
            t[n+1] = n * dt
        end

        # Create data arrays for the independent kinematic fields.
        F  = PF.ArrayOfPhysicalTensors(N+1, 2, 2, DIMENSIONLESS)
        F′ = PF.ArrayOfPhysicalTensors(N+1, 2, 2, TIME_RATE)

        # Assign values to deformation gradient in its initial configuration κ₀.
        F₀ = PF.PhysicalTensor(2, 2, DIMENSIONLESS)
        F₀[1,1] = PF.PhysicalScalar(1.0, DIMENSIONLESS)
        F₀[2,2] = PF.PhysicalScalar(1.0, DIMENSIONLESS)
        F[1]  = F₀
        F′[1] = PF.PhysicalTensor(2, 2, TIME_RATE)

        # Assign initial conditions to the Laplace stretch attributes.
        a₀ = PF.PhysicalScalar(1.0, DIMENSIONLESS)
        b₀ = PF.PhysicalScalar(1.0, DIMENSIONLESS)
        γ₀ = PF.PhysicalScalar(DIMENSIONLESS)
        ω₀ = PF.PhysicalScalar(DIMENSIONLESS)

        # Assign initial conditions to the Laplace stretch rate attributes.
        a′₀ = PF.PhysicalScalar(TIME_RATE)
        b′₀ = PF.PhysicalScalar(TIME_RATE)
        γ′₀ = PF.PhysicalScalar(TIME_RATE)
        ω′₀ = PF.PhysicalScalar(TIME_RATE)

        # Data array that holds the various cases of motion.
        motion = zeros(Int, N+1)
        motion[1] = Pᵣ

        # Data arrays that hold the Gram rotations and their rates: κ₀ ↦ κₙ.
        ωₙ  = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ω′ₙ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)
        ωₙ[1]  = ω₀
        ω′ₙ[1] = ω′₀

        # Data arrays for Laplace stretch attributes and their rates: κᵣ ↦ κₙ.
        aₙ = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        bₙ = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        γₙ = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        aₙ[1] = a₀/aᵣ
        bₙ[1] = b₀/bᵣ
        γₙ[1] = (aᵣ/bᵣ)*(γ₀ - gᵣ)
        a′ₙ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)
        b′ₙ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)
        γ′ₙ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)
        a′ₙ[1] = a′₀/aᵣ
        b′ₙ[1] = b′₀/bᵣ
        γ′ₙ[1] = (aᵣ/bᵣ)*γ′₀

        # Data arrays for the thermodynamic strains and their rates: κᵣ ↦ κₙ.
        δ = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ε = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        γ = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        δ[1] = PF.PhysicalScalar(0.5log((a₀/aᵣ)*(b₀/bᵣ)), DIMENSIONLESS)
        ε[1] = PF.PhysicalScalar(0.5log((a₀/aᵣ)*(bᵣ/b₀)), DIMENSIONLESS)
        γ[1] = γ₀ - gᵣ
        δ′ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)
        ε′ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)
        γ′ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)
        δ′[1] = 0.5(a′₀/a₀ + b′₀/b₀)
        ε′[1] = 0.5(a′₀/a₀ - b′₀/b₀)
        γ′[1] = γ′₀

        # Create and return a new data structure for Laplace kinematics in 2D.
        new(dt, N, n, aᵣ, bᵣ, gᵣ, t, F, F′, motion, ωₙ, ω′ₙ, 
            aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ε, γ, δ′, ε′, γ′)::MembraneKinematics
    end

    # Internal constructor used by JSON3.

    function MembraneKinematics(dt::PF.PhysicalScalar, N::Int, n::PF.MInteger, aᵣ::PF.PhysicalScalar, bᵣ::PF.PhysicalScalar, γᵣ::PF.PhysicalScalar, t::PF.ArrayOfPhysicalScalars, F::PF.ArrayOfPhysicalTensors, F′::PF.ArrayOfPhysicalTensors, motion::Vector{Int}, ωₙ::PF.ArrayOfPhysicalScalars, ω′ₙ::PF.ArrayOfPhysicalScalars, aₙ::PF.ArrayOfPhysicalScalars, bₙ::PF.ArrayOfPhysicalScalars, γₙ::PF.ArrayOfPhysicalScalars, a′ₙ::PF.ArrayOfPhysicalScalars, b′ₙ::PF.ArrayOfPhysicalScalars, γ′ₙ::PF.ArrayOfPhysicalScalars, δ::PF.ArrayOfPhysicalScalars, ε::PF.ArrayOfPhysicalScalars, γ::PF.ArrayOfPhysicalScalars, δ′::PF.ArrayOfPhysicalScalars, ε′::PF.ArrayOfPhysicalScalars, γ′::PF.ArrayOfPhysicalScalars)

        new(dt, N, n, aᵣ, bᵣ, γᵣ, t, F, F′, motion, ωₙ, ω′ₙ, 
            aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ε, γ, δ′, ε′, γ′)::MembraneKinematics
    end
end # MembraneKinematics

# Methods

function Base.:(copy)(k::MembraneKinematics)::MembraneKinematics
    dt  = copy(k.dt)
    N   = copy(k.N)
    n   = copy(k.n)
    aᵣ  = copy(k.aᵣ)
    bᵣ  = copy(k.bᵣ)
    γᵣ  = copy(k.γᵣ)
    t   = copy(k.t)
    F   = copy(k.F)
    F′  = copy(k.F′)
    motion = copy(k.motion)
    ωₙ  = copy(k.ωₙ)
    ω′ₙ = copy(k.ω′ₙ)
    aₙ  = copy(k.aₙ)
    a′ₙ = copy(k.a′ₙ)
    bₙ  = copy(k.bₙ)
    b′ₙ = copy(k.b′ₙ)
    γₙ  = copy(k.γₙ)
    γ′ₙ = copy(k.γ′ₙ)
    δ   = copy(k.δ)
    δ′  = copy(k.δ′)
    ε   = copy(k.ε)
    ε′  = copy(k.ε′)
    γ   = copy(k.γ)
    γ′  = copy(k.γ′)
    return MembraneKinematics(dt, N, n, aᵣ, bᵣ, γᵣ, t, F, F′, motion, ωₙ, ω′ₙ,
                              aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ε, γ, δ′, ε′, γ′)
end

# The histories of MembraneKinematics are to be graphed, not printed,
# so a toString method is not provided for objects of this type.

# Methods for storing and retrieving a MembraneKinematics data structure 
# to and from a file.

StructTypes.StructType(::Type{MembraneKinematics}) = StructTypes.Struct()

function toFile(k::MembraneKinematics, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, k)
        write(json_stream, '\n')
    else
        error("The supplied JSON stream is not open.")
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{MembraneKinematics}, json_stream::IOStream)::MembraneKinematics
    if isopen(json_stream)
        k = JSON3.read(readline(json_stream), MembraneKinematics)
    else
        error("The supplied JSON stream is not open.")
    end
    return k
end

# Methods that serve as a solver for objects of type MembraneKinematics.

function _advance!(k::MembraneKinematics, m::Int)
    # Solve Laplace attributes and their rates from k.F and k.F′ at step m.
    Fₘ = k.F[m]
    F′ = k.F′[m]
    
    x  = Fₘ[1,1]
    y  = Fₘ[2,2]
    x′ = F′[1,1]
    y′ = F′[2,2]

    # Establish the Gram and Laplace attributes, and their rates.
    if Fₘ[2,1] ≈ 0.0 || sign(Fₘ[1,2]) == sign(Fₘ[2,1])
        # The deformation has an attribute that is a pure shear.
        # g is this pure shear.
        if (abs(Fₘ[1,2])/Fₘ[2,2] > abs(Fₘ[2,1])/Fₘ[1,1] ||
            abs(Fₘ[1,2])/Fₘ[2,2] ≈ abs(Fₘ[2,1])/Fₘ[1,1])
            # Co-ordinates are indexed as supplied--the default condition.
            case = 1
            # Pure shear contributions.
            g  = Fₘ[2,1] / Fₘ[1,1]
            g′ = (Fₘ[1,1]*F′[2,1] - Fₘ[2,1]*F′[1,1]) / (Fₘ[1,1]*Fₘ[1,1])
            # Simple shear contributions.
            G  = (Fₘ[1,1]*Fₘ[1,2] - Fₘ[2,2]*Fₘ[2,1]) / (Fₘ[1,1]*Fₘ[2,2])
            G′ = (Fₘ[2,2]*F′[1,2] - Fₘ[1,2]*F′[2,2]) / (Fₘ[2,2]*Fₘ[2,2]) - g′
        else
            # The Gram co-ordinate system is left handed.
            case = 2
            # Pure shear contributions.
            g  = Fₘ[1,2] / Fₘ[2,2]
            g′ = (Fₘ[2,2]*F′[1,2] - Fₘ[1,2]*F′[2,2]) / (Fₘ[2,2]*Fₘ[2,2])
            # Simple shear contributions.
            G  = -(Fₘ[1,1]*Fₘ[1,2] - Fₘ[2,2]*Fₘ[2,1]) / (Fₘ[1,1]*Fₘ[2,2])
            G′ =  (Fₘ[1,1]*F′[2,1] - Fₘ[2,1]*F′[1,1]) / (Fₘ[1,1]*Fₘ[1,1]) - g′
        end
        # Laplace stretch attributes in the (𝕚, 𝕛) co-ordinate frame.
        aₘ = x * sqrt(1+g*g)
        bₘ = y * (1 - g*(g+G)) / sqrt(1+g*g)
        γₘ = (y/x) * (2g+G) / (1+g*g)
        ωₘ = PF.PhysicalScalar(atan(g), DIMENSIONLESS)

        # Rates of Laplace attributes in the (𝕚, 𝕛) co-ordinate frame.
        a′ₘ = aₘ*(x′/x + g*g′/(1+g*g))
        b′ₘ = bₘ*(y′/y - g*g′/(1+g*g)) - y*((2g+G)*g′ + g*G′)/sqrt(1+g*g)
        γ′ₘ = γₘ*(y′/y - x′/x - 2g*g′/(1+g*g)) + (y/x)*(2g′+G′)/(1+g*g)
        ω′ₘ = g′/(1+g*g)
    else
        # The deformation has an attribute that is a rigid-body rotation.
        # g is this rigid-body rotation.
        if (abs(Fₘ[1,2])/Fₘ[2,2] > abs(Fₘ[2,1])/Fₘ[1,1] ||
            abs(Fₘ[1,2])/Fₘ[2,2] ≈ abs(Fₘ[2,1])/Fₘ[1,1])
            # The Gram co-ordinate system is right handed.
            case = 3
            # Rigid-body rotation contributions.
            g  = -Fₘ[2,1] / Fₘ[1,1]
            g′ = -(Fₘ[1,1]*F′[2,1] - Fₘ[2,1]*F′[1,1]) / (Fₘ[1,1]*Fₘ[1,1])
            # Angle of rigid-body rotation in the (𝕚, 𝕛) co-ordinate frame.
            ωₘ  = PF.PhysicalScalar(-atan(g), DIMENSIONLESS)
            ω′ₘ = -g′/(1+g*g)
        else
            # The Gram co-ordinate system is left handed.
            case = 4
            # Rigid-body rotation contributions.
            g  = -Fₘ[1,2] / Fₘ[2,2]
            g′ = -(Fₘ[2,2]*F′[1,2] - Fₘ[1,2]*F′[2,2]) / (Fₘ[2,2]*Fₘ[2,2])
            # Angle of rigid-body rotation in the (𝕚, 𝕛) co-ordinate frame.
            ωₘ  = PF.PhysicalScalar(atan(g), DIMENSIONLESS)
            ω′ₘ = g′/(1+g*g)
        end
        # G is a simple shear.
        G  = (Fₘ[1,1]*Fₘ[1,2] + Fₘ[2,2]*Fₘ[2,1]) / (Fₘ[1,1]*Fₘ[2,2])
        G′ = ((Fₘ[2,2]*F′[1,2] - Fₘ[1,2]*F′[2,2]) / (Fₘ[2,2]*Fₘ[2,2]) +
              (Fₘ[1,1]*F′[2,1] - Fₘ[2,1]*F′[1,1]) / (Fₘ[1,1]*Fₘ[1,1]))

        # Laplace attributes in the (𝕚, 𝕛) co-ordinate frame.
        aₘ = x * sqrt(1+g*g)
        bₘ = y * (1 + g*(g+G)) / sqrt(1+g*g)
        γₘ = (y/x) * G / (1+g*g)

        # Rates of Laplace attributes in the (𝕚, 𝕛) co-ordinate frame.
        a′ₘ = aₘ*(x′/x + g*g′/(1+g*g))
        b′ₘ = bₘ*(y′/y - g*g′/(1+g*g)) + y*((2g+G)*g′ + g*G′)/sqrt(1+g*g)
        γ′ₘ = γₘ*(y′/y - x′/x - 2g*g′/(1+g*g)) + (y/x)*G′/(1+g*g)
    end

    # Advance the data array that holds pivoting cases.
    k.motion[m] = case

    # Advance the Laplace stretch attributes and their rates for κᵣ ↦ κₙ.
    k.aₙ[m]  = aₘ/k.aᵣ
    k.bₙ[m]  = bₘ/k.bᵣ
    k.γₙ[m]  = (k.aᵣ/k.bᵣ)*(γₘ - k.γᵣ)
    k.ωₙ[m]  = ωₘ
    k.a′ₙ[m] = a′ₘ/k.aᵣ
    k.b′ₙ[m] = b′ₘ/k.bᵣ
    k.γ′ₙ[m] = (k.aᵣ/k.bᵣ)*γ′ₘ
    k.ω′ₙ[m] = ω′ₘ

    # Advance the thermodynamic Laplace strains and their rates for κᵣ ↦ κₙ.
    k.δ[m]  = PF.PhysicalScalar(0.5log(k.aₙ[m]*k.bₙ[m]), DIMENSIONLESS)
    k.ε[m]  = PF.PhysicalScalar(0.5log(k.aₙ[m]/k.bₙ[m]), DIMENSIONLESS)
    k.γ[m]  = γₘ - k.γᵣ
    k.δ′[m] = 0.5(a′ₘ/aₘ + b′ₘ/bₘ)
    k.ε′[m] = 0.5(a′ₘ/aₘ - b′ₘ/bₘ)
    k.γ′[m] = γ′ₘ

    return nothing
end # _advance!

function advance!(k::MembraneKinematics, dF::PF.PhysicalTensor)
    # Advance the counter.
    if k.n < k.N+1
        PF.set!(k.n, PF.get(k.n)+1)
    else
        println("The data structure is full and cannot accept further data.")
        return nothing
    end
    n = PF.get(k.n)

    # Convert the passed variable to CGS units.
    dFₙ = PF.toCGS(dF)

    # Verify inputs.
    if dF.matrix.rows ≠ 2 || dF.matrix.cols ≠ 2
        error("Deformation gradient differential dF must be a 2x2 matrix.")
    end
    if dF.units ≠ DIMENSIONLESS
        error("Deformation gradient differential dF must be dimensionless.")
    end
    
    # Update the deformation gradient, noting that dFₙ = Fₙ - Fₙ₋₁,
    # with F[1] = I, and given that n = 2,3,…,N+1, then
    k.F[n] = k.F[n-1] + dFₙ
    
    # Approximate the derivative of the deformation gradient using the
    # following finite difference formulæ.
    if n == 2
        # Euler's first-order forward and backward difference formulæ.
        k.F′[1] = (-k.F[1] + k.F[2]) / k.dt
        k.F′[2] = ( k.F[2] - k.F[1]) / k.dt
    elseif n == 3
        # Second-order forward, central and backward difference formulæ.
        k.F′[1] = (-3k.F[1] + 4k.F[2] - k.F[3]) / (2k.dt)
        k.F′[2] = (k.F[3] - k.F[1]) / (2k.dt)
        k.F′[3] = ( 3k.F[3] - 4k.F[2] + k.F[1]) / (2k.dt)
    elseif n == 4
        # Third-order forward and backward difference formulæ.
        k.F′[1] = (-11k.F[1] + 18k.F[2] - 9k.F[3] + 2k.F[4]) / (6k.dt)
        k.F′[4] = ( 11k.F[4] - 18k.F[3] + 9k.F[2] - 2k.F[1]) / (6k.dt)
    elseif n == 5
        # Third-order forward and backward difference formulæ.
        k.F′[2] = (-11k.F[2] + 18k.F[3] - 9k.F[4] + 2k.F[5]) / (6k.dt)
        k.F′[5] = ( 11k.F[5] - 18k.F[4] + 9k.F[3] - 2k.F[2]) / (6k.dt)
    elseif n == 6
        # Third-order forward and backward difference formulæ.
        k.F′[3] = (-11k.F[3] + 18k.F[4] - 9k.F[5] + 2k.F[6]) / (6k.dt)
        k.F′[6] = ( 11k.F[6] - 18k.F[5] + 9k.F[4] - 2k.F[3]) / (6k.dt)
    else
        # Third-order backward difference formula.
        k.F′[n] = (11k.F[n] - 18k.F[n-1] + 9k.F[n-2] - 2k.F[n-3]) / (6k.dt)
    end

    # Advance values for the Laplace attributes and their rates.
    # First update the history, if necessary, to improve overall accuracy.
    if n == 2
        _advance!(k, 1)
    elseif n == 3
        _advance!(k, 1)
        _advance!(k, 2)
    elseif n == 4
        _advance!(k, 1)
    elseif n == 5
        _advance!(k, 2)
    elseif n == 6
        _advance!(k, 3)
    else
        nothing
    end
    _advance!(k, n)

    return nothing
end # advance!

function update!(k::MembraneKinematics, dF::PF.PhysicalTensor)
    if k.n > 1
        PF.set!(k.n, PF.get(k.n)-1)
        advance!(k, dF)
    end
    return nothing
end # update!

