#=
Created on Mon 22 Nov 2021
Updated on Mon 27 Nov 2023
-------------------------------------------------------------------------------
This software, like the language it is written in, is published under the MIT
License, https://opensource.org/licenses/MIT.

Copyright (c) 2021-2023:
Alan Freed and John Clayton

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
-------------------------------------------------------------------------------
=#

"""
Type:\n
    FiberKinematics\n
        # Properties of the arrays.
        dt      # time increment separating neighboring nodes\n
        N       # total node count for traversing a solution path\n
        n       # a counter that ratchets from 1 to N+1\n

        # Array of nodal times.
        t       # array of current times of length N+1\n

        # Reference (strain free) values describing a 1D fiber element.\n
        Lᵣ      # reference length\n
        Aᵣ      # reference cross-sectional area\n

        # History arrays of length N+1 for holding the kinematic fields.\n
        # Initial values/conditions are stored in array location [1].\n
        L       # array of current lengths of length N+1\n
        L′      # array of current length rates of length N+1\n
        A       # array of current cross-sectional areas of length N+1\n
        A′      # array of current cross-sectional area rates of length N+1\n

        # Thermodynamic (true) strains and their rates.\n
        ϵ       # array of current strains of length N+1\n
        ϵ′      # array of current strain rates of length N+1\n
FiberKinematics is a data structure that contains the kinematic fields necessary to describe an isochoric 1D fiber. The arrays in this data structure allow for a history of these kinematic fields to be used, e.g., in a constitutive analysis, for graphing, etc.
"""
struct FiberKinematics
    # Properties of the arrays.
    dt::PhysicalScalar          # time increment separating neighboring nodes
    N::Integer                  # number of nodes to traverse a solution path
    n::MInteger                 # a counter that ratchets from 1 to N

    # Array of nodal times.
    t::ArrayOfPhysicalScalars   # times at the solution nodes

    # Reference (strain free) values describing an isochoric 1D fiber element.
    Lᵣ::PhysicalScalar          # reference length
    Aᵣ::PhysicalScalar          # reference cross-sectional area

    # History arrays of length N+1 for holding the kinematic fields.\n
    # Initial values/conditions are stored in array location [1].\n
    L::ArrayOfPhysicalScalars   # lengths at the solution nodes
    L′::ArrayOfPhysicalScalars  # length rates at the solution nodes
    A::ArrayOfPhysicalScalars   # areas at the solution nodes
    A′::ArrayOfPhysicalScalars  # area rates at the solution nodes

    # Thermodynamic (true) strains and their rates.
    ϵ::ArrayOfPhysicalScalars   # strains at the solution nodes
    ϵ′::ArrayOfPhysicalScalars  # strain rates at the solution nodes

    # Internal constructors.

"""
    Constructor:\n
        k = FiberKinematics(dt, N, Lᵣ, Aᵣ, L₀)\n
    Returns a new data structure `k` of type `FiberKinematics` that holds kinematic fields pertinent to the modeling of an isochoric 1D fiber. Arguments are: (i) A differential step in time `dt` that separates neighboring nodes, which themselves are taken to be uniformly spaced in time. (ii) The total number of grid points or nodes `N` where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays. (iii) The reference (or strain free) length `Lᵣ` and cross-sectional area `Aᵣ` of a fiber against which strains are to be measured. And (iv) a fiber's initial length `L₀` in some initial configuration selected for analysis κ₀, where typically L₀ ≥ Lᵣ. An isochoric (or constant volume) motion is assumed.
"""
    function FiberKinematics(dt::PhysicalScalar, N::Integer, Lᵣ::PhysicalScalar, Aᵣ::PhysicalScalar, L₀::PhysicalScalar)

        # Convert all passed variables to CGS units.
        d𝑡 = toCGS(dt)
        𝐿ᵣ = toCGS(Lᵣ)
        𝐴ᵣ = toCGS(Aᵣ)
        𝐿₀ = toCGS(L₀)

        # Physical bounds:
        tₘᵢₙ = PhysicalScalar(eps(Float64), TIME)
        Lₘᵢₙ = PhysicalScalar(eps(Float32), LENGTH)

        # Verify inputs.
        if d𝑡.units ≠ TIME
            msg = "The time increment dt does not have units of time."
            throw(ErrorException(msg))
        end
        if d𝑡 < tₘᵢₙ
            msg = "The time increment dt must be positive valued."
            throw(ErrorException(msg))
        end
        if N < 1
            msg = "Solution arrays must have a positive length."
            throw(ErrorException(msg))
        end
        if 𝐿ᵣ.units ≠ LENGTH
            msg = "The reference length Lᵣ does not have units of length."
            throw(ErrorException(msg))
        end
        if 𝐿ᵣ < Lₘᵢₙ
            msg = "The reference length Lᵣ must be positive valued."
            throw(ErrorException(msg))
        end
        if 𝐴ᵣ.units ≠ AREA
            msg = "The reference area Aᵣ does not have units of area."
            throw(ErrorException(msg))
        end
        if 𝐴ᵣ < Lₘᵢₙ*Lₘᵢₙ
            msg = "The reference area Aᵣ must be positive valued."
            throw(ErrorException(msg))
        end
        if 𝐿₀.units ≠ LENGTH
            msg = "The initial length L₀ does not have units of length."
            throw(ErrorException(msg))
        end
        if 𝐿₀ < Lₘᵢₙ
            msg = "The initial length L₀ must be positive valued."
            throw(ErrorException(msg))
        end

        # Initialize the counter.
        n = MInteger(1)

        # Create and populate an array for nodal times.
        t  = ArrayOfPhysicalScalars(N+1, TIME)
        for n in 1:N
            t[n+1] = n * d𝑡
        end

        # Create data arrays for the physical dimensions and their rates.
        L  = ArrayOfPhysicalScalars(N+1, LENGTH)
        L′ = ArrayOfPhysicalScalars(N+1, LENGTH_RATE)
        A  = ArrayOfPhysicalScalars(N+1, AREA)
        A′ = ArrayOfPhysicalScalars(N+1, AREA_RATE)

        # Assign to these arrays their initial values.
        L′₀ = PhysicalScalar(LENGTH_RATE)
        A₀  = 𝐴ᵣ * 𝐿ᵣ / 𝐿₀
        A′₀ = PhysicalScalar(AREA_RATE)
        L[1]  = 𝐿₀
        L′[1] = L′₀
        A[1]  = A₀
        A′[1] = A′₀

        # Create data arrays for thermodynamic strains and their rates: κᵣ ↦ κₙ.
        ϵ  = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ϵ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)

        # Assign to these arrays their initial values.
        ϵ[1]  = PhysicalScalar(log(𝐿₀/𝐿ᵣ), DIMENSIONLESS)
        ϵ′[1] = L′₀ / 𝐿₀

        # Return a new data structure for managing kinematics of a 1D fiber.
        new(d𝑡, N, n, t, 𝐿ᵣ, 𝐴ᵣ, L, L′, A, A′, ϵ, ϵ′)
    end

    # The internal constructor used by JSON3 and other external constructors.

    function FiberKinematics(dt::PhysicalScalar, N::Integer, n::MInteger, t::ArrayOfPhysicalScalars, Lᵣ::PhysicalScalar, Aᵣ::PhysicalScalar, L::ArrayOfPhysicalScalars, L′::ArrayOfPhysicalScalars, A::ArrayOfPhysicalScalars, A′::ArrayOfPhysicalScalars, ϵ::ArrayOfPhysicalScalars, ϵ′::ArrayOfPhysicalScalars)
        new(dt, N, n, t, Lᵣ, Aᵣ, L, L′, A, A′, ϵ, ϵ′)
    end
end # FiberKinematics

# Methods

function Base.:(copy)(k::FiberKinematics)::FiberKinematics
    dt = copy(k.dt)
    N  = copy(k.N)
    n  = copy(k.n)
    t  = copy(k.t)
    Lᵣ = copy(k.Lᵣ)
    Aᵣ = copy(k.Aᵣ)
    L  = copy(k.L)
    L′ = copy(k.L′)
    A  = copy(k.A)
    A′ = copy(k.A′)
    ϵ  = copy(k.ϵ)
    ϵ′ = copy(k.ϵ′)
    return FiberKinematics(dt, N, n, t, Lᵣ, Aᵣ, L, L′, A, A′, ϵ, ϵ′)
end

function Base.:(deepcopy)(k::FiberKinematics)::FiberKinematics
    dt = deepcopy(k.dt)
    N  = deepcopy(k.N)
    n  = deepcopy(k.n)
    t  = deepcopy(k.t)
    Lᵣ = deepcopy(k.Lᵣ)
    Aᵣ = deepcopy(k.Aᵣ)
    L  = deepcopy(k.L)
    L′ = deepcopy(k.L′)
    A  = deepcopy(k.A)
    A′ = deepcopy(k.A′)
    ϵ  = deepcopy(k.ϵ)
    ϵ′ = deepcopy(k.ϵ′)
    return FiberKinematics(dt, N, n, t, Lᵣ, Aᵣ, L, L′, A, A′, ϵ, ϵ′)
end

# The histories of FiberKinematics are to be graphed, not printed, so a toString method is not provided for objects of this type.

# Methods for storing and retrieving a FiberKinematics data structure to and from a file.

StructTypes.StructType(::Type{FiberKinematics}) = StructTypes.Struct()

"""
Method:\n
    toFile(k::LaplaceKinematics.FiberKinematics, json_stream::IOStream)\n
Writes data structure `k` to the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONWriter(<my_dir_path>, <my_file_name>)\n
    ...\n
    LaplaceKinematics.toFile(k, json_stream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream)\n
where <my_dir_path> is the path to your working directory wherein the file <my_file_name> to be written to either exists or will be created, and which must have a .json extension.
"""
function toFile(k::FiberKinematics, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, k)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    flush(json_stream)
    return nothing
end

"""
Method:\n
    fromFile(k::LaplaceKinematics.FiberKinematics, json_stream::IOStream)\n
Reads a FiberKinematics data structure from the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONReader(<my_dir_path>, <my_file_name>)\n
    ...\n
    k = LaplaceKinematics.fromFile(LaplaceKinematics.FiberKinematics, json_stream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream)\n
that returns `k,` which is an object of type FiberKinematics. Here <my_dir_path> is the path to your working directory wherein the file <my_file_name> to be read from must exist, and which is to have a .json extension.
"""
function fromFile(::Type{FiberKinematics}, json_stream::IOStream)::FiberKinematics
    if isopen(json_stream)
        k = JSON3.read(readline(json_stream), FiberKinematics)
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    return k
end

# Methods that serve as a solver for objects of type FiberKinematics.

"""
Method:\n
    advance!(k::FiberKinematics, L′::PhysicalScalar)\n
Method `advance!` moves a solution from previous step `n-1` to current step `n` along a solution path of N solution nodes by integrating its governing differential equation for length using a backward difference formula (BDF) when given the fiber's current time rate-of-change in length `L′`.

This method updates counter `k.n` and entries to its history arrays at the nᵗʰ array location in the `k` data structure; specifically: length `k.L[n]` and its rate `k.L′[n]`, area `k.A[n]` and its rate `k.A′[n]`, plus strain `k.ϵ[n]` and its rate `k.ϵ′[n]`.
"""
function advance!(k::FiberKinematics, L′::PhysicalScalar)
    # Advance the counter.
    if k.n < k.N+1
        set!(k.n, get(k.n)+1)
    else
        msg = "The data structure is full and cannot accept further data."
        throw(ErrorException(msg))
    end
    n = get(k.n)

    # Convert to CGS units.
    𝐿′ = toCGS(L′)

    # Verify the input.
    if 𝐿′.units ≠ LENGTH_RATE
        msg = "Fiber length rate L′ must have units of velocity."
        throw(ErrorException(msg))
    end
    k.L′[n] = 𝐿′

    # Integrate fiber length using a backward difference formula (BDF).
    if n == 2
        k.L[2] = k.L[1] + k.L′[2]*k.dt
    elseif n == 3
        k.L[3] = (4/3)*k.L[2] - (1/3)*k.L[1] + (2/3)*k.L′[3]*k.dt
    else
        k.L[n] = ((18/11)*k.L[n-1] - (9/11)*k.L[n-2] + (2/11)*k.L[n-3]
               + (6/11)*k.L′[n]*k.dt)
    end

    # Compute cross-sectional area and its rate assuming an isochoric response.
    k.A[n]  = k.Aᵣ * k.Lᵣ / k.L[n]
    k.A′[n] = -k.A[n] * k.L′[n] / k.L[n]

    # Compute the current strain and its rate.
    k.ϵ[n]  = PhysicalScalar(log(k.L[n]/k.Lᵣ), DIMENSIONLESS)
    k.ϵ′[n] = k.L′[n] / k.L[n]

    return nothing
end # advance!

"""
Method:\n
    update!(k::FiberKinematics, L′::PhysicalScalar)\n
Method `update!` refines a solution at step `n` by re-integrating its governing differential equation for fiber length, thereby allowing for iterative improvements to be made on length rate `L′` from an external algorithm, e.g., a finite element engine. There is no need to call `update!` unless `L′` is being iteratively refined at step `n`, e.g., by some external optimization process.
"""
function update!(k::FiberKinematics, L′::PhysicalScalar)
    if k.n > 1
        set!(k.n, get(k.n)-1)
        advance!(k, F′ₙ)
    end
    return nothing
end # update!
