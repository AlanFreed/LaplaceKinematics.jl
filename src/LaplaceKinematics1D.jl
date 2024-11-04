"""
The fields for an object of type `FiberKinematics` are:

    # Properties of the arrays.
    dt   # time increment separating neighboring nodes
    N    # number of intervals along a solution path
    n    # a nodal counter that ratchets from 1 to N+1

    # Reference stretch (stretch at zero strain) of a 1D fiber element.
    λᵣ   # reference stretch, λᵣ = Lᵣ/L₀, where L denotes its length

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t    # time at the solution nodes

    # Arrays for the fiber stretch and its rate.
    λ    # stretch at the solution nodes
    λ′   # stretch rate at the solution nodes
        
    # Arrays for the thermodynamic (true) strains and their rates.
    ε    # strain at the solution nodes
    ε′   # strain rate at the solution nodes
    
*FiberKinematics* is a data structure that contains the physical fields needed to describe the kinematics of a 1D fiber. The arrays in this data structure allow for a history of these kinematic fields to be used, e.g., in a constitutive analysis, for graphing, etc.

## Constructors

The constructor most likely to be used.
```julia
k = FiberKinematics(dt::PF.PhysicalScalar, N::Int, λᵣ::PF.PhysicalScalar)
```
where
    
1. `dt` is a differential step in time that separates neighboring nodes, which themselves are taken to be uniformly spaced over time. 
2. `N` is the total number of grid points or nodes where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays, e.g., t[1] = 0. 
3. `λᵣ` is the reference (or strain free) stretch of a fiber against which strains are to be measured, viz., ε(λᵣ) = 0, with the fiber's initial stretch λ₀ associating with some initial configuration κ₀. Stretch λ₀ is assigned a value of 1, with an outcome being that strain ε(λ₀) need not equal 0.

The general constructor used by JSON3, etc.
```julia
k = FiberKinematics(dt::PF.PhysicalScalar, N::Int, n::PF.MInteger, λᵣ::PF.PhysicalScalar, t::PF.ArrayOfPhysicalScalars, λ::PF.ArrayOfPhysicalScalars, λ′::PF.ArrayOfPhysicalScalars, ε::PF.ArrayOfPhysicalScalars, ε′::PF.ArrayOfPhysicalScalars)
```

## Methods

```julia
cc = copy(k::LK.FiberKinematics)
```
returns a copy `cc` of object `k` of type `FiberKinematics`.

```julia
toFile(k::LK.FiberKinematics, json_stream::IOStream)
```
writes data structure `k` to an IOStream `json_stream`.

```julia 
k = fromFile(LK.FiberKinematics, json_stream::IOStream)
```
reads data structure `k` from an IOStream `json_stream`.

To manage a `json_stream` for writing, consider the code fragment:

    json_stream = PF.openJSONWriter(<my_dir_path>, <my_file_name>)
    LK.toFile(k, json_stream)
    PF.closeJSONStream(json_stream)
    
while to manage a `json_stream` for reading, consider the code fragment:

    json_stream = PF.openJSONReader(<my_dir_path>, <my_file_name>)
    k = LK.fromFile(LK.FiberKinematics, json_stream)
    PF.closeJSONStream(json_stream)

where `<my_dir_path>` is the path to your working directory wherein the file to be written, i.e., `<my_file_name>`, either exists or will be created. This file must have a .json extension.

```julia
advance!(k::FiberKinematics, dλ::PhysicalScalar)
```
Method `advance!` moves a solution from previous step `n-1` to current step `n` along a solution path with `N` solution nodes. It is assumed that stretch is controlled (and is therefore known as a function of time). Argument `dλ` denotes a differential, not a derivative. From these differentials, stretch rates are computed via third-order, finite-difference formulæ from which strain rates are then established.

This method updates counter `k.n`, plus those entries to its history arrays that are at the nᵗʰ array location in this `k` data structure; specifically: stretch `k.λ[n]` and its rate `k.λ′[n]`, plus strain `k.ε[n]` and its rate `k.ε′[n].` 

**Note**: Stretches are assigned at the end points of each solution interval; consequently, one is to send `dλ = λ(tₙ) - λ(tₙ₋₁)` for `n=1,2,…,N`.
    
```julia
update!(k::FiberKinematics, dλ::PhysicalScalar)
```
Method `update!` refines a solution at step `n` whenever an improvement can be made for the stretch differential `dλ` through an external algorithm, e.g., a finite element engine. There is no need to call `update!` unless `dλ` is being iteratively refined at a global step `n` via, say, some external optimization process.
"""
struct FiberKinematics
    # Properties of the arrays.
    dt::PF.PhysicalScalar          # time increment separating neighboring nodes
    N::Int                         # number of intervals along a solution path
    n::PF.MInteger                 # a counter that ratchets from 1 to N+1

    # Reference stretch (stretch at zero strain) of a 1D fiber element.
    λᵣ::PF.PhysicalScalar          # reference stretch, λᵣ = Lᵣ/L₀, L is length

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t::PF.ArrayOfPhysicalScalars   # time at the solution nodes

    # Arrays for the fiber stretch and its rate.
    λ::PF.ArrayOfPhysicalScalars   # stretches at the solution nodes
    λ′::PF.ArrayOfPhysicalScalars  # stretch rates at the solution nodes

    # Arrays for the thermodynamic (true) strain and its rate.
    ε::PF.ArrayOfPhysicalScalars   # strains at the solution nodes
    ε′::PF.ArrayOfPhysicalScalars  # strain rates at the solution nodes

    # Internal constructors.

    function FiberKinematics(dt::PF.PhysicalScalar, 
                             N::Int, 
                             λᵣ::PF.PhysicalScalar)

        # Convert all passed variables to CGS units.
        dt = PF.toCGS(dt)
        λᵣ = PF.toCGS(λᵣ)

        # Physical bounds:
        tₘᵢₙ = PF.PhysicalScalar(eps(Float64), TIME)
        λₘᵢₙ = PF.PhysicalScalar(eps(Float32), DIMENSIONLESS)

        # Verify inputs.
        if dt.units ≠ TIME
            error("The time increment dt does not have units of time.")
        end
        if dt < tₘᵢₙ
            error("The time increment dt must be positive valued.")
        end
        if N < 1
            error("Solution arrays must have a positive length.")
        end
        if λᵣ.units ≠ DIMENSIONLESS
            error("The reference stretch λᵣ must be dimensionless.")
        end
        if λᵣ < λₘᵢₙ
            error("The reference stretch λᵣ must be positive valued.")
        end

        # Initialize the counter.
        n = PF.MInteger(1)

        # Create and populate the array for nodal times.
        t = PF.ArrayOfPhysicalScalars(N+1, TIME)
        for n in 1:N
            t[n+1] = n * dt
        end

        # Create data arrays for the stretches and their rates.
        λ  = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        λ′ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)

        # Assign to these arrays their initial values.
        λ[1]  = PF.PhysicalScalar(1.0, DIMENSIONLESS)
        λ′[1] = PF.PhysicalScalar(0.0, TIME_RATE)

        # Create data arrays for thermodynamic strains and their rates: κᵣ ↦ κₙ.
        ε  = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ε′ = PF.ArrayOfPhysicalScalars(N+1, TIME_RATE)

        # Assign to these arrays their initial values.
        ε[1]  = PF.PhysicalScalar(log(λ[1]/λᵣ), DIMENSIONLESS)
        ε′[1] = λ′[1] / λ[1]

        # Return a new data structure for managing kinematics of a 1D fiber.
        new(dt, N, n, λᵣ, t, λ, λ′, ε, ε′)::FiberKinematics
    end

    # The internal constructor used by JSON3 and other external constructors.

    function FiberKinematics(dt::PF.PhysicalScalar,
                             N::Int, 
                             n::PF.MInteger, 
                             λᵣ::PF.PhysicalScalar, 
                             t::PF.ArrayOfPhysicalScalars,
                             λ::PF.ArrayOfPhysicalScalars,
                             λ′::PF.ArrayOfPhysicalScalars,
                             ε::PF.ArrayOfPhysicalScalars,
                             ε′::PF.ArrayOfPhysicalScalars)
        new(dt, N, n, λᵣ, t, λ, λ′, ε, ε′)::FiberKinematics
    end
end # FiberKinematics

# Methods

function Base.:(copy)(k::FiberKinematics)::FiberKinematics
    dt = copy(k.dt)
    N  = copy(k.N)
    n  = copy(k.n)
    λᵣ = copy(k.λᵣ)
    t  = copy(k.t)
    λ  = copy(k.λ)
    λ′ = copy(k.λ′)
    ε  = copy(k.ε)
    ε′ = copy(k.ε′)
    return FiberKinematics(dt, N, n, λᵣ, t, λ, λ′, ε, ε′)
end

# The histories of FiberKinematics are to be graphed, not printed, 
# so a toString method is not provided for objects of this type.

# Methods for storing and retrieving a FiberKinematics data structure 
# to and from a file.

StructTypes.StructType(::Type{FiberKinematics}) = StructTypes.Struct()

function toFile(k::FiberKinematics, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, k)
        write(json_stream, '\n')
    else
        error("The supplied JSON stream is not open.")
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{FiberKinematics}, json_stream::IOStream)::FiberKinematics
    if isopen(json_stream)
        k = JSON3.read(readline(json_stream), FiberKinematics)
    else
        error("The supplied JSON stream is not open.")
    end
    return k
end

# Methods that serve as a solver for objects of type FiberKinematics.

function advance!(k::FiberKinematics, dλ::PF.PhysicalScalar)
    # Advance the counter.
    if k.n < k.N+1
        PF.set!(k.n, PF.get(k.n)+1)
    else
        println("The data structure is full and cannot accept further data.")
        return nothing
    end
    n = PF.get(k.n)
    
    # Convert to CGS units.
    dλₙ = PF.toCGS(dλ)
    
    # Verify the input.
    if dλ.units ≠ DIMENSIONLESS
        error("An incremental change in stretch dλ must be dimensionless.")
    end
    
    # Update the stretch, noting that dλₙ = λₙ - λₙ₋₁,
    # with λ[1] = λ₀ = 1, and given that n = 2,3,…,N+1, then
    k.λ[n] = k.λ[n-1] + dλₙ
    
    # Approximate stretch rate using the following finite difference formulæ.
    if n == 2
        # Euler's first-order forward and backward difference formulæ.
        k.λ′[1] = (-k.λ[1] + k.λ[2]) / k.dt
        k.λ′[2] = ( k.λ[2] - k.λ[1]) / k.dt
    elseif n == 3
        # Second-order forward, central and backward difference formulæ.
        k.λ′[1] = (-3k.λ[1] + 4k.λ[2] - k.λ[3]) / (2k.dt)
        k.λ′[2] = (k.λ[3] - k.λ[1]) / (2k.dt)
        k.λ′[3] = ( 3k.λ[3] - 4k.λ[2] + k.λ[1]) / (2k.dt)
    elseif n == 4
        # Third-order forward and backward difference formulæ.
        k.λ′[1] = (-11k.λ[1] + 18k.λ[2] - 9k.λ[3] + 2k.λ[4]) / (6k.dt)
        k.λ′[4] = ( 11k.λ[4] - 18k.λ[3] + 9k.λ[2] - 2k.λ[1]) / (6k.dt)
    elseif n == 5
        # Third-order forward and backward difference formulæ.
        k.λ′[2] = (-11k.λ[2] + 18k.λ[3] - 9k.λ[4] + 2k.λ[5]) / (6k.dt)
        k.λ′[5] = ( 11k.λ[5] - 18k.λ[4] + 9k.λ[3] - 2k.λ[2]) / (6k.dt)
    elseif n == 6
        # Third-order forward and backward difference formulæ.
        k.λ′[3] = (-11k.λ[3] + 18k.λ[4] - 9k.λ[5] + 2k.λ[6]) / (6k.dt)
        k.λ′[6] = ( 11k.λ[6] - 18k.λ[5] + 9k.λ[4] - 2k.λ[3]) / (6k.dt)
    else
        # Third-order backward difference formula.
        k.λ′[n] = (11k.λ[n] - 18k.λ[n-1] + 9k.λ[n-2] - 2k.λ[n-3]) / (6k.dt)
    end

    # Compute the current strain and its rate (to third order).
    k.ε[n]  = PF.PhysicalScalar(log(k.λ[n]/k.λᵣ), DIMENSIONLESS)
    if n == 2
        k.ε′[1] = k.λ′[1] / k.λ[1]
    elseif n == 3
        k.ε′[1] = k.λ′[1] / k.λ[1]
        k.ε′[2] = k.λ′[2] / k.λ[2]
    elseif n == 4
        k.ε′[1] = k.λ′[1] / k.λ[1]
    elseif n == 5
        k.ε′[2] = k.λ′[2] / k.λ[2]
    elseif n == 6
        k.ε′[3] = k.λ′[3] / k.λ[3]
    else
        nothing
    end
    k.ε′[n] = k.λ′[n] / k.λ[n]
    
    return nothing
end # advance!

function update!(k::FiberKinematics, dλ::PF.PhysicalScalar)
    if k.n > 1
        PF.set!(k.n, PF.get(k.n)-1)
        advance!(k, dλ)
    end
    return nothing
end # update!

