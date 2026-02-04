"""
The fields for an object of type `FiberKinematics` are:

    # Properties of the arrays.
    dt  # time increment separating neighboring nodes
    N   # number of uniform intervals along a solution path
    n   # a nodal counter that ratchets from 1 to N+1

    # Reference (strain free) and initial (start of analysis) lengths.
    Lᵣ  # reference length of a fiber
    L₀  # initial length of a fiber, 0 < Lᵣ ≤ L₀

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array holding the independent variable, i.e., array of nodal times.
    t   # time at the solution nodes
    
    # Arrays for fiber stretch and its rate.
    λ   # stretch at solution nodes
    λ′  # stretch rate at solution nodes
    λ″  # rate of stretch rate at solution nodes
        
    # Arrays for thermodynamic (true) strain and its rate.
    ε   # strain at solution nodes
    ε′  # strain rate at solution nodes
    ε″  # rate of strain rate at solution nodes
    
*FiberKinematics* is a data structure that contains the physical fields needed to describe the kinematics of a 1D fiber. The arrays in this data structure allow for a history of these kinematic fields to be used, e.g., in a constitutive analysis, or for graphing, etc.

## Constructors

The constructor most likely to be used.
```julia
k = FiberKinematics(N::Int, 
                    dt::PF.PhysicalScalar,
                    Lᵣ::PF.PhysicalScalar,
                    L₀::PF.PhysicalScalar)
```
where

1. `N` is the total number of grid points or nodes where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays, e.g., t[1] = 0 s.
 
2. `dt` is a differential step in time that separates neighboring nodes, which themselves are taken to be uniformly spaced over time. *dt* has units of time.

3. `Lᵣ` is the length of a fiber in its reference, strain-free configuration κᵣ.  Lᵣ must have units of PF.CGS_LENGTH.

4. `L₀` is the initial length of a fiber in its initial configuration κ₀, i.e., a configuration at the start of an analysis. L₀ must have units of PF.CGS_LENGTH. Typically, Lᵣ ≤ L₀ implicating that the fiber is usually prestressed in its initial configuration.

The general constructor used by JSON3, etc.
```julia
k = FiberKinematics(dt::PF.PhysicalScalar, 
                    N::Int, 
                    n::PF.MInteger,
                    Lᵣ::PF.PhysicalScalar, 
                    L₀::PF.PhysicalScalar, 
                    t::PF.ArrayOfPhysicalScalars, 
                    λ::PF.ArrayOfPhysicalScalars,
                    λ′::PF.ArrayOfPhysicalScalars,
                    λ″::PF.ArrayOfPhysicalScalars,
                    ε::PF.ArrayOfPhysicalScalars,
                    ε′::PF.ArrayOfPhysicalScalars,
                    ε″::PF.ArrayOfPhysicalScalars)
```

## Methods

```julia
cc = copy(k::LK.FiberKinematics)
```
returns a copy `cc` for object `k` of type `FiberKinematics`.

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

where `<my_dir_path>` is the path to your working directory, wherein the file to be written, viz., `<my_file_name>`, either exists or will be created. This file must have a `.json` extension.

```julia
advance!(k::FiberKinematics, Lₙ::PhysicalScalar)
```
Method `advance!` moves a solution from previous step *n-1* to current step *n* along a solution path with *N* solution nodes. It is assumed that fiber length is controlled (and is therefore known as a function of time).  Argument `Lₙ` denotes fiber length in the current configuration κₙ, which must have physical units of *LENGTH*.  From this length all fields are updated, with rates and accelerations being approximated using third-order accurate, finite-difference formulæ.

This method updates counter *k.n*, plus it assigns entries to its history arrays at their nᵗʰ array location in this *k* data structure; specifically: position *k.x[n]*, velocity *k.v[n]*, acceleration *k.a[n]*, along with stretch *k.λ[n]* and its rate *k.λ′[n]*, plus strain *k.ε[n]* and its rate *k.ε′[n]*.


**Note**: It is not until step *n=6* that all tabulated derivatives in time become third-order accurate, or until step *n=8* that all tabulated second derivatives in time become third-order accurate. Derivative at the early steps are updated with more accurate approximations as more data becomes available, until they become third-order accurate.

```julia
update!(k::FiberKinematics, Lₙ::PhysicalScalar)
```
Method `update!` is to be called after `advance!`. It refines a solution at step *n* for an updated value of ` Lₙ` that is to be inserted into data structure `k`, e.g., updated by some external optimization process like a finite element engine.
"""

struct FiberKinematics
    # Properties of the arrays.
    dt::PF.PhysicalScalar          # time between neighboring solution nodes
    N::Int                         # number of intervals along a solution path
    n::PF.MInteger                 # a counter that ratchets from 1 to N+1

    # Reference (strain free) and initial (start of analysis) lengths.
    Lᵣ::PF.PhysicalScalar          # reference length of a fiber
    L₀::PF.PhysicalScalar          # initial length of a fiber, 0 < Lᵣ ≤ L₀

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array holding the independent variable, i.e., array of nodal times.
    t::PF.ArrayOfPhysicalScalars   # time at solution nodes
    
    # Arrays for the fiber stretch and its rates.
    λ::PF.ArrayOfPhysicalScalars   # stretch at solution nodes
    λ′::PF.ArrayOfPhysicalScalars  # stretch rate at solution nodes
    λ″::PF.ArrayOfPhysicalScalars  # rate of stretch rate at solution nodes

    # Arrays for the thermodynamic (true) strain and its rates.
    ε::PF.ArrayOfPhysicalScalars   # strain at solution nodes
    ε′::PF.ArrayOfPhysicalScalars  # strain rate at solution nodes
    ε″::PF.ArrayOfPhysicalScalars  # rate of strain rate at solution nodes

    # Internal constructors.

    function FiberKinematics(N::Int, 
                             dt::PF.PhysicalScalar,
                             Lᵣ::PF.PhysicalScalar,
                             L₀::PF.PhysicalScalar)

        # Convert all passed variables to CGS units.
        dt = PF.toCGS(dt)
        Lᵣ = PF.toCGS(Lᵣ)
        L₀ = PF.toCGS(L₀)

        # Physical bounds:
        tₘᵢₙ = PF.PhysicalScalar(eps(Float64), TIME)
        Lₘᵢₙ = PF.PhysicalScalar(eps(Float32), LENGTH)

        # Verify inputs.
        if N < 1
            error("Solution arrays must have a positive length.")
        end
        if dt.units ≠ TIME
            error("The time increment dt does not have units of time.")
        end
        if dt < tₘᵢₙ
            error("The time increment dt must be positive valued.")
        end
        if Lᵣ.units ≠ LENGTH
            error("The reference length Lᵣ must have units of length.")
        end
        if Lᵣ < Lₘᵢₙ
            error("The reference length Lᵣ must be positive valued.")
        end
        if L₀.units ≠ LENGTH
            error("The initial length L₀ must have units of length.")
        end
        if L₀ < Lᵣ
            error("Length L₀ must not be less than length Lᵣ.")
        end

        # Create and populate the array for nodal times.
        t = PF.ArrayOfPhysicalScalars(N+1, TIME)
        for n in 1:N
            t[n+1] = n * dt
        end
        
        # Create a mutable counter.
        n = PF.MInteger(1)

        # Create data arrays for the stretches and their rates.
        λ  = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        λ′ = PF.ArrayOfPhysicalScalars(N+1, RATE)
        λ″ = PF.ArrayOfPhysicalScalars(N+1, RATEofRATE)

        # Assign to these arrays their initial values.
        λ[1]  = PF.PhysicalScalar(1.0, DIMENSIONLESS)
        λ′[1] = PF.PhysicalScalar(0.0, RATE)
        λ″[1] = PF.PhysicalScalar(0.0, RATEofRATE)

        # Create data arrays for thermodynamic strain and its rate: κᵣ ↦ κₙ.
        ε  = PF.ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ε′ = PF.ArrayOfPhysicalScalars(N+1, RATE)
        ε″ = PF.ArrayOfPhysicalScalars(N+1, RATEofRATE)

        # Assign to these arrays their initial values.
        ε[1]  = PF.PhysicalScalar(log(L₀/Lᵣ), DIMENSIONLESS)
        ε′[1] = PF.PhysicalScalar(0.0, RATE)
        ε″[1] = PF.PhysicalScalar(0.0, RATEofRATE)

        # Return a new data structure for managing kinematics of a 1D fiber.
        new(dt, N, n, Lᵣ, L₀, t, λ, λ′, λ″, ε, ε′, ε″)::FiberKinematics
    end

    # The internal constructor used by JSON3 and other external constructors.

    function FiberKinematics(dt::PF.PhysicalScalar,
                             N::Int, 
                             n::PF.MInteger, 
                             Lᵣ::PF.PhysicalScalar, 
                             L₀::PF.PhysicalScalar,
                             t::PF.ArrayOfPhysicalScalars,
                             λ::PF.ArrayOfPhysicalScalars,
                             λ′::PF.ArrayOfPhysicalScalars,
                             λ″::PF.ArrayOfPhysicalScalars,
                             ε::PF.ArrayOfPhysicalScalars,
                             ε′::PF.ArrayOfPhysicalScalars,
                             ε″::PF.ArrayOfPhysicalScalars)
        new(dt, N, n, Lᵣ, L₀, t, λ, λ′, λ″, ε, ε′, ε″)::FiberKinematics
    end
end # FiberKinematics

# Methods

function Base.:(copy)(k::FiberKinematics)::FiberKinematics
    dt = copy(k.dt)
    N  = copy(k.N)
    n  = copy(k.n)
    Lᵣ = copy(k.Lᵣ)
    L₀ = copy(k.L₀)
    t  = copy(k.t)
    λ  = copy(k.λ)
    λ′ = copy(k.λ′)
    λ″ = copy(k.λ″)
    ε  = copy(k.ε)
    ε′ = copy(k.ε′)
    ε″ = copy(k.ε″)
    return FiberKinematics(dt, N, n, Lᵣ, L₀, t, λ, λ′, λ″, ε, ε′, ε″)
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
# Here Lₙ = L(k.t[n]), which is the length of a fiber evaluated at time tₙ.

function advance!(k::FiberKinematics, Lₙ::PF.PhysicalScalar)
    # Advance the counter.
    if k.n < k.N+1
        PF.set!(k.n, PF.get(k.n)+1)
    else
        println("The data structure is full and cannot accept further data.")
        return nothing
    end
    n = PF.get(k.n)
    
    # Convert to CGS units.
    Lₙ = PF.toCGS(Lₙ)
    
    # Verify the input.
    Lₘᵢₙ = PF.PhysicalScalar(eps(Float32), LENGTH)
    if Lₙ.units ≠ LENGTH
        error("Current length Lₙ must have units of length.")
    end
    if Lₙ < Lₘᵢₙ
        error("Current length Lₙ must be positive valued.")
    end
    
    # Finite differences are:
    #    Forward  difference:  ∆fₙ = fₙ₊₁ - fₙ
    #    Backward difference:  ∇fₙ = fₙ - fₙ₋₁
    
    # Advance the stretch and its rates.
    
    k.λ[n] = Lₙ / k.L₀
    
    # Stretch rates are accurate to a third-order error, viz., ∇⁴/4h.
    if n == 2
        # first-order estimates
        k.λ′[1] = (k.λ[2] - k.λ[1]) / k.dt   # forward  difference
        k.λ′[2] = (k.λ[2] - k.λ[1]) / k.dt   # backward difference
    elseif n == 3
        # second-order estimates
        k.λ′[1] = (-k.λ[3] + 4k.λ[2] - 3k.λ[1]) / (2k.dt)   # forward
        k.λ′[2] = ( k.λ[3]      -       k.λ[1]) / (2k.dt)   # central
        k.λ′[3] = (3k.λ[3] - 4k.λ[2] +  k.λ[1]) / (2k.dt)   # backward
    elseif n == 4
        # third-order estimates
        k.λ′[1] = ( 2k.λ[4] -  9k.λ[3] + 18k.λ[2] - 11k.λ[1]) / (6k.dt)
        k.λ′[4] = (11k.λ[4] - 18k.λ[3] +  9k.λ[2] -  2k.λ[1]) / (6k.dt)
    elseif n == 5
        k.λ′[2] = ( 2k.λ[5] -  9k.λ[4] + 18k.λ[3] - 11k.λ[2]) / (6k.dt)
        k.λ′[5] = (11k.λ[5] - 18k.λ[4] +  9k.λ[3] -  2k.λ[2]) / (6k.dt)
    elseif n == 6
        k.λ′[3] = ( 2k.λ[6] -  9k.λ[5] + 18k.λ[4] - 11k.λ[3]) / (6k.dt)
        k.λ′[6] = (11k.λ[6] - 18k.λ[5] +  9k.λ[4] -  2k.λ[3]) / (6k.dt)
    else
        k.λ′[n] = (11k.λ[n] - 18k.λ[n-1] + 9k.λ[n-2] - 2k.λ[n-3]) / (6k.dt)
    end
    
    # Rates of stretch rates are accurate to a third-order error, viz., 5∇⁵/6h².
    if n == 2
        k.λ″[2] = PF.PhysicalScalar(0.0, RATEofRATE)
    elseif n == 3
        # first-order estimate
        k.λ″[2] = (k.λ[3] - 2k.λ[2] + k.λ[1]) / (k.dt*k.dt)
        k.λ″[3] = PF.PhysicalScalar(0.0, RATEofRATE)
    elseif n == 4
        # first- and second-order estimates
        k.λ″[1] = (-k.λ[4] + 4k.λ[3] - 5k.λ[2] + 2k.λ[1]) / (k.dt*k.dt)
        k.λ″[3] = ( k.λ[4] - 2k.λ[3] +  k.λ[2])           / (k.dt*k.dt)
        k.λ″[4] = (2k.λ[4] - 5k.λ[3] + 4k.λ[2] -  k.λ[1]) / (k.dt*k.dt)
    elseif n == 5
        # second- and third-order estimates
        k.λ″[1] = (11k.λ[5] -  56k.λ[4] + 114k.λ[3] - 104k.λ[2] + 35k.λ[1]) / (k.dt*k.dt)
        k.λ″[2] = ( -k.λ[5] +   4k.λ[4] -   5k.λ[3] +   2k.λ[2])            / (k.dt*k.dt)
        k.λ″[5] = (35k.λ[5] - 104k.λ[4] + 114k.λ[3] -  56k.λ[2] + 11k.λ[1]) / (k.dt*k.dt)
    elseif n == 6
        k.λ″[2] = (11k.λ[6] -  56k.λ[5] + 114k.λ[4] - 104k.λ[3] + 35k.λ[2]) / (k.dt*k.dt)
        k.λ″[3] = ( -k.λ[6] +   4k.λ[5] -   5k.λ[4] +   2k.λ[3])            / (k.dt*k.dt)
        k.λ″[6] = (35k.λ[6] - 104k.λ[5] + 114k.λ[4] -  56k.λ[3] + 11k.λ[2]) / (k.dt*k.dt)
    elseif n == 7
        # third-order estimates
        k.λ″[3] = (11k.λ[7] -  56k.λ[6] + 114k.λ[5] - 104k.λ[4] + 35k.λ[3]) / (k.dt*k.dt)
        k.λ″[7] = (35k.λ[7] - 104k.λ[6] + 114k.λ[5] -  56k.λ[4] + 11k.λ[3]) / (k.dt*k.dt)
    elseif n == 8
        k.λ″[4] = (11k.λ[8] -  56k.λ[7] + 114k.λ[6] - 104k.λ[5] + 35k.λ[4]) / (k.dt*k.dt)
        k.λ″[8] = (35k.λ[8] - 104k.λ[7] + 114k.λ[6] -  56k.λ[5] + 11k.λ[4]) / (k.dt*k.dt)
    else   
        k.λ″[n] = (35k.λ[n] - 104k.λ[n-1] + 114k.λ[n-2] - 56k.λ[n-3] + 11k.λ[n-4]) / (k.dt*k.dt)
    end

    # Advance the strain and its rates.
    
    k.ε[n] = PF.PhysicalScalar(log(Lₙ/k.Lᵣ), DIMENSIONLESS)
    
    if n == 2
        k.ε′[1] = k.λ′[1] / k.λ[1]
        k.ε′[2] = k.λ′[2] / k.λ[2]
    elseif n == 3
        k.ε′[1] = k.λ′[1] / k.λ[1]
        k.ε′[2] = k.λ′[2] / k.λ[2]
        k.ε′[3] = k.λ′[3] / k.λ[3]
    elseif n == 4
        k.ε′[1] = k.λ′[1] / k.λ[1]
        k.ε′[4] = k.λ′[4] / k.λ[4]
    elseif n == 5
        k.ε′[2] = k.λ′[2] / k.λ[2]
        k.ε′[5] = k.λ′[5] / k.λ[5]
    elseif n == 6
        k.ε′[3] = k.λ′[3] / k.λ[3]
        k.ε′[6] = k.λ′[6] / k.λ[6]
    else
        k.ε′[n] = k.λ′[n] / k.λ[n]
    end
    
    if n == 2
        k.ε″[2] = (k.λ[2]*k.λ″[2] - k.λ′[2]*k.λ′[2]) / (k.λ[2]*k.λ[2])
    elseif n == 3
        k.ε″[2] = (k.λ[2]*k.λ″[2] - k.λ′[2]*k.λ′[2]) / (k.λ[2]*k.λ[2])
        k.ε″[3] = (k.λ[3]*k.λ″[3] - k.λ′[3]*k.λ′[3]) / (k.λ[3]*k.λ[3])
    elseif n == 4
        k.ε″[1] = (k.λ[1]*k.λ″[1] - k.λ′[1]*k.λ′[1]) / (k.λ[1]*k.λ[1])
        k.ε″[3] = (k.λ[3]*k.λ″[3] - k.λ′[3]*k.λ′[3]) / (k.λ[3]*k.λ[3])
        k.ε″[4] = (k.λ[4]*k.λ″[4] - k.λ′[4]*k.λ′[4]) / (k.λ[4]*k.λ[4])
    elseif n == 5
        k.ε″[1] = (k.λ[1]*k.λ″[1] - k.λ′[1]*k.λ′[1]) / (k.λ[1]*k.λ[1])
        k.ε″[2] = (k.λ[2]*k.λ″[2] - k.λ′[2]*k.λ′[2]) / (k.λ[2]*k.λ[2])
        k.ε″[5] = (k.λ[5]*k.λ″[5] - k.λ′[5]*k.λ′[5]) / (k.λ[5]*k.λ[5])
    elseif n == 6
        k.ε″[2] = (k.λ[2]*k.λ″[2] - k.λ′[2]*k.λ′[2]) / (k.λ[2]*k.λ[2])
        k.ε″[3] = (k.λ[3]*k.λ″[3] - k.λ′[3]*k.λ′[3]) / (k.λ[3]*k.λ[3])
        k.ε″[6] = (k.λ[6]*k.λ″[6] - k.λ′[6]*k.λ′[6]) / (k.λ[6]*k.λ[6])
    elseif n == 7
        k.ε″[3] = (k.λ[3]*k.λ″[3] - k.λ′[3]*k.λ′[3]) / (k.λ[3]*k.λ[3])
        k.ε″[7] = (k.λ[7]*k.λ″[7] - k.λ′[7]*k.λ′[7]) / (k.λ[7]*k.λ[7])
    elseif n == 8
        k.ε″[4] = (k.λ[4]*k.λ″[4] - k.λ′[4]*k.λ′[4]) / (k.λ[4]*k.λ[4])
        k.ε″[8] = (k.λ[8]*k.λ″[8] - k.λ′[8]*k.λ′[8]) / (k.λ[8]*k.λ[8])
    else
        k.ε″[n] = (k.λ[n]*k.λ″[n] - k.λ′[n]*k.λ′[n]) / (k.λ[n]*k.λ[n])
    end
    
    return nothing
end # advance!

function update!(k::FiberKinematics, Lₙ::PF.PhysicalScalar)
    if k.n > 1
        PF.set!(k.n, PF.get(k.n)-1)
        advance!(k, Lₙ)
    end
    return nothing
end # update!
