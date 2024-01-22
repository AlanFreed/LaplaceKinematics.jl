"""
Type:\n
    FiberKinematics\n
        # Properties of the arrays.
        dt      # time increment separating neighboring nodes\n
        N       # total node count for traversing a solution path\n
        n       # a counter that ratchets from 1 to N+1\n

        # Reference (strain free) length of a 1D fiber element.\n
        Lᵣ      # reference length\n

        # History arrays of length N+1 for holding the kinematic fields.\n
        # Initial values/conditions are stored in array location [1].\n

        # Array of the independent variable, viz., array of nodal times.
        t       # time at the solution nodes\n

        # Arrays for the fiber length and its rate.\n
        L       # length at the solution nodes\n
        L′      # length rate at the solution nodes\n

        # Arrays for the thermodynamic (true) strains and their rates.\n
        ϵ       # strain at the solution nodes\n
        ϵ′      # strain rate at the solution nodes\n
FiberKinematics is a data structure that contains the physical fields necessary to describe kinematics of a 1D fiber. The arrays in this data structure allow for a history of these kinematic fields to be used, e.g., in a constitutive analysis, for graphing, etc.
"""
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

    # Internal constructors.

"""
    Constructor:\n
        k = FiberKinematics(dTime, N, midPtQuad, Lᵣ, L₀)\n
    Returns a new data structure `k` of type `FiberKinematics` that holds kinematic fields pertinent to the modeling of a 1D fiber. Arguments are: (i) A differential step in time `dTime` that separates neighboring nodes, which themselves are taken to be uniformly spaced in time. (ii) The total number of grid points or nodes `N` where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays. (iii) A boolean flag `midPtQuad` that, if true, implies the nodal spacing is for a mid-point quadrature rule and, if false, implies the nodal spacing is for an end-point quadrature rule. This determines how the array of independent times is populated. (iv) The reference (or strain free) length `Lᵣ` of a fiber against which strains are to be measured. And (v) a fiber's initial length `L₀` in some initial configuration selected for analysis κ₀ where, typically, L₀ ≥ Lᵣ.
"""
    function FiberKinematics(dTime::PhysicalScalar, N::Integer, midPtQuad::Bool, Lᵣ::PhysicalScalar, L₀::PhysicalScalar)

        # Convert all passed variables to CGS units.
        dt = toCGS(dTime)
        𝐿ᵣ = toCGS(Lᵣ)
        𝐿₀ = toCGS(L₀)

        # Physical bounds:
        tₘᵢₙ = PhysicalScalar(eps(Float64), TIME)
        Lₘᵢₙ = PhysicalScalar(eps(Float32), LENGTH)

        # Verify inputs.
        if dt.units ≠ TIME
            msg = "The time increment dTime does not have units of time."
            throw(ErrorException(msg))
        end
        if dt < tₘᵢₙ
            msg = "The time increment dTime must be positive valued."
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
        if midPtQuad
            # Assign times for a mid-point quadrature rule.
            for n in 1:N
                t[n+1] = 0.5dt + (n-1) * dt
            end
        else
            # Assign times for an end-point quadrature rule.
            for n in 1:N
                t[n+1] = n * dt
            end
        end

        # Create data arrays for the physical dimensions and their rates.
        L  = ArrayOfPhysicalScalars(N+1, LENGTH)
        L′ = ArrayOfPhysicalScalars(N+1, LENGTH_RATE)

        # Assign to these arrays their initial values.
        L′₀ = PhysicalScalar(LENGTH_RATE)
        L[1]  = 𝐿₀
        L′[1] = L′₀

        # Create data arrays for thermodynamic strains and their rates: κᵣ ↦ κₙ.
        ϵ  = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ϵ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)

        # Assign to these arrays their initial values.
        ϵ[1]  = PhysicalScalar(log(𝐿₀/𝐿ᵣ), DIMENSIONLESS)
        ϵ′[1] = L′₀ / 𝐿₀

        # Return a new data structure for managing kinematics of a 1D fiber.
        new(dt, N, n, 𝐿ᵣ, t, L, L′, ϵ, ϵ′)
    end

    # The internal constructor used by JSON3 and other external constructors.

    function FiberKinematics(dt::PhysicalScalar, N::Integer, n::MInteger, Lᵣ::PhysicalScalar, t::ArrayOfPhysicalScalars, L::ArrayOfPhysicalScalars, L′::ArrayOfPhysicalScalars,ϵ::ArrayOfPhysicalScalars, ϵ′::ArrayOfPhysicalScalars)
        new(dt, N, n, Lᵣ, t, L, L′, ϵ, ϵ′)
    end
end # FiberKinematics

# Methods

function Base.:(copy)(k::FiberKinematics)::FiberKinematics
    dt = copy(k.dt)
    N  = copy(k.N)
    n  = copy(k.n)
    Lᵣ = copy(k.Lᵣ)
    t  = copy(k.t)
    L  = copy(k.L)
    L′ = copy(k.L′)
    ϵ  = copy(k.ϵ)
    ϵ′ = copy(k.ϵ′)
    return FiberKinematics(dt, N, n, Lᵣ, t, L, L′, ϵ, ϵ′)
end

function Base.:(deepcopy)(k::FiberKinematics)::FiberKinematics
    dt = deepcopy(k.dt)
    N  = deepcopy(k.N)
    n  = deepcopy(k.n)
    Lᵣ = deepcopy(k.Lᵣ)
    t  = deepcopy(k.t)
    L  = deepcopy(k.L)
    L′ = deepcopy(k.L′)
    ϵ  = deepcopy(k.ϵ)
    ϵ′ = deepcopy(k.ϵ′)
    return FiberKinematics(dt, N, n, Lᵣ, t, L, L′, ϵ, ϵ′)
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
Method `advance!` moves a solution from previous step `n-1` to current step `n` along a solution path of N solution nodes by integrating its governing differential equation for length using a backward difference formula (BDF) when given the fiber's time rate-of-change in length `L′`. For a time-step interval of [tₙ₋₁, tₙ], L′ = dL/dt associates with either time tₙ when using end-point quadrature, or with time (tₙ₋₁ + tₙ)/2 when using mid-point quadrature.\n

This method updates counter `k.n` and entries to its history arrays at the nᵗʰ array location in the `k` data structure; specifically: length `k.L[n]` and its rate `k.L′[n]`, plus strain `k.ϵ[n]` and its rate `k.ϵ′[n].` These fields are evaluated at either the end-point, i.e. at time tₙ, or at the mid-point, i.e. at time (tₙ₋₁ + tₙ)/2, according to the argument `midPtQuad` supplied to its constructor.
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

    # Integrate fiber length rate using a backward difference formula (BDF).
    if k.t[2] ≈ 0.5k.dt
        # Integrated for nodes located at the mid-point of each time step.
        if n == 2
            k.L[2] = k.L[1] + 0.5k.L′[2]*k.dt
        elseif n == 3
            L₁ = k.L[1] - 0.5k.L′[2]*k.dt
            k.L[3] = (4/3)*k.L[2] - (1/3)*k.L[1] + (2/3)*k.L′[3]*k.dt
        elseif n == 4
            L₁ = k.L[1] - 0.5k.L′[2]*k.dt
            k.L[4] = ((18/11)*k.L[3] - (9/11)*k.L[2] + (2/11)*k.L[1]
                   + (6/11)*k.L′[4]*k.dt)
        else
            k.L[n] = ((18/11)*k.L[n-1] - (9/11)*k.L[n-2] + (2/11)*k.L[n-3]
                   + (6/11)*k.L′[n]*k.dt)
        end
    else
        # Integrated for nodes located at the end-point of each time step.
        if n == 2
            k.L[2] = k.L[1] + k.L′[2]*k.dt
        elseif n == 3
            k.L[3] = (4/3)*k.L[2] - (1/3)*k.L[1] + (2/3)*k.L′[3]*k.dt
        else
            k.L[n] = ((18/11)*k.L[n-1] - (9/11)*k.L[n-2] + (2/11)*k.L[n-3]
                   + (6/11)*k.L′[n]*k.dt)
        end
    end

    # Compute the current strain and its rate.
    k.ϵ[n]  = PhysicalScalar(log(k.L[n]/k.Lᵣ), DIMENSIONLESS)
    k.ϵ′[n] = k.L′[n] / k.L[n]

    return nothing
end # advance!

"""
Method:\n
    update!(k::FiberKinematics, L′::PhysicalScalar)\n
Method `update!` refines a solution at step `n` by re-integrating its governing differential equation for fiber length, thereby allowing for iterative improvements to be made on length rate `L′` from an external algorithm, e.g., a finite element engine. There is no need to call `update!` unless `L′` is being iteratively refined at step `n`, e.g., by some external optimization process. Here L′ = dL(k.t[n])/dt.
"""
function update!(k::FiberKinematics, L′::PhysicalScalar)
    if k.n > 1
        set!(k.n, get(k.n)-1)
        advance!(k, F′ₙ)
    end
    return nothing
end # update!
