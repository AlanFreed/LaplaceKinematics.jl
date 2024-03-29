"""
Type:\n
    FiberKinematics\n
        # Properties of the arrays.
        dt      # time increment separating neighboring nodes\n
        N       # number of intervals along a solution path\n
        n       # a counter that ratchets from 1 to N+1\n

        # Reference (strain free) stretch of a 1D fiber element.\n
        λᵣ      # reference stretch\n

        # History arrays are of length N+1 for holding the kinematic fields.\n
        # Initial values/conditions are stored in array location [1].\n

        # Array of the independent variable, viz., array of nodal times.\n
        t       # time at the solution nodes\n

        # Arrays for the fiber stretch and its rate.\n
        λ       # stretches at the solution nodes\n
        λ′      # stretch rates at the solution nodes\n

        # Arrays for the thermodynamic (true) strains and their rates.\n
        ϵ       # strains at the solution nodes\n
        ϵ′      # strain rates at the solution nodes\n
FiberKinematics is a data structure that contains the physical fields necessary to describe kinematics of a 1D fiber. The arrays in this data structure allow for a history of these kinematic fields to be used, e.g., in a constitutive analysis, for graphing, etc.
"""
struct FiberKinematics
    # Properties of the arrays.
    dt::PhysicalScalar          # time increment separating neighboring nodes
    N::Integer                  # number of intervals along a solution path
    n::MInteger                 # a counter that ratchets from 1 to N+1

    # Reference (strain free) length of a 1D fiber element.
    λᵣ::PhysicalScalar          # reference stretch, λᵣ = Lᵣ / L₀, L is length

    # History arrays are of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t::ArrayOfPhysicalScalars   # time at the solution nodes

    # Arrays for the fiber stretch and its rate.
    λ::ArrayOfPhysicalScalars   # stretches at the solution nodes
    λ′::ArrayOfPhysicalScalars  # stretch rates at the solution nodes

    # Arrays for the thermodynamic (true) strain and its rate.
    ϵ::ArrayOfPhysicalScalars   # strains at the solution nodes
    ϵ′::ArrayOfPhysicalScalars  # strain rates at the solution nodes

    # Internal constructors.

"""
    Constructor:\n
        k = FiberKinematics(dTime, N, midPtQuad, lambdaᵣ)\n
    Returns a new data structure `k` of type `FiberKinematics` that holds kinematic fields pertinent for the modeling of a 1D fiber. Arguments are: (i) A differential step in time `dTime` that separates neighboring nodes, which themselves are taken to be uniformly spaced over time. (ii) The total number of grid points or nodes `N` where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays, e.g., t[1] = 0. (iii) A boolean flag `midPtQuad` that, if true, implies the nodal spacing is for a mid-point quadrature rule, and if false, implies the nodal spacing is for an end-point quadrature rule. This determines how the array of independent times is to be populated. (iv) The reference (or strain free) stretch `lambdaᵣ` of a fiber against which strains are to be measured, viz., ϵ(λᵣ) = 0, with the fiber's initial stretch λ₀, associated with some initial configuration κ₀, being assigned a value of 1 with an outcome being that ϵ(λ₀) need not equal 0.
"""
    function FiberKinematics(dTime::PhysicalScalar, N::Integer, midPtQuad::Bool, lambdaᵣ::PhysicalScalar)

        # Convert all passed variables to CGS units.
        dt = toCGS(dTime)
        λᵣ = toCGS(lambdaᵣ)

        # Physical bounds:
        tₘᵢₙ = PhysicalScalar(eps(Float64), TIME)
        λₘᵢₙ = PhysicalScalar(eps(Float32), DIMENSIONLESS)

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
        if λᵣ.units ≠ DIMENSIONLESS
            msg = "The reference stretch lambdaᵣ must be dimensionless."
            throw(ErrorException(msg))
        end
        if λᵣ < λₘᵢₙ
            msg = "The reference stretch lambdaᵣ must be positive valued."
            throw(ErrorException(msg))
        end

        # Initialize the counter.
        n = MInteger(1)

        # Create and populate the array for nodal times.
        t  = ArrayOfPhysicalScalars(N+1, TIME)
        if midPtQuad
            # Assign times for a mid-point quadrature rule.
            for n in 1:N
                t[n+1] = (n-1)*dt + 0.5dt
            end
        else
            # Assign times for an end-point quadrature rule.
            for n in 1:N
                t[n+1] = n*dt
            end
        end

        # Create data arrays for the stretches and their rates.
        λ  = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        λ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)

        # Assign to these arrays their initial values.
        λ[1]  = PhysicalScalar(1.0, DIMENSIONLESS)
        λ′[1] = PhysicalScalar(0.0, TIME_RATE)

        # Create data arrays for thermodynamic strains and their rates: κᵣ ↦ κₙ.
        ϵ  = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ϵ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)

        # Assign to these arrays their initial values.
        ϵ[1]  = PhysicalScalar(log(λ[1]/λᵣ), DIMENSIONLESS)
        ϵ′[1] = λ′[1] / λ[1]

        # Return a new data structure for managing kinematics of a 1D fiber.
        new(dt, N, n, λᵣ, t, λ, λ′, ϵ, ϵ′)
    end

    # The internal constructor used by JSON3 and other external constructors.

    function FiberKinematics(dt::PhysicalScalar, N::Integer, n::MInteger, λᵣ::PhysicalScalar, t::ArrayOfPhysicalScalars, λ::ArrayOfPhysicalScalars, λ′::ArrayOfPhysicalScalars,ϵ::ArrayOfPhysicalScalars, ϵ′::ArrayOfPhysicalScalars)
        new(dt, N, n, λᵣ, t, λ, λ′, ϵ, ϵ′)
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
    ϵ  = copy(k.ϵ)
    ϵ′ = copy(k.ϵ′)
    return FiberKinematics(dt, N, n, λᵣ, t, λ, λ′, ϵ, ϵ′)
end

function Base.:(deepcopy)(k::FiberKinematics)::FiberKinematics
    dt = deepcopy(k.dt)
    N  = deepcopy(k.N)
    n  = deepcopy(k.n)
    λᵣ = deepcopy(k.λᵣ)
    t  = deepcopy(k.t)
    λ  = deepcopy(k.λ)
    λ′ = deepcopy(k.λ′)
    ϵ  = deepcopy(k.ϵ)
    ϵ′ = deepcopy(k.ϵ′)
    return FiberKinematics(dt, N, n, λᵣ, t, λ, λ′, ϵ, ϵ′)
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
    advance!(k::FiberKinematics, lambda′::PhysicalScalar)\n
Method `advance!` moves a solution from previous step `n-1` to current step `n` along a solution path of N solution nodes by integrating its governing differential equation for stretch using a backward difference formula (BDF) when given the fiber's time rate-of-change in stretch `lambda′`. For a time-step interval of [tₙ₋₁, tₙ], λ′ = dλ/dt associates with either time tₙ when using an end-point quadrature rule, or with time (tₙ₋₁ + tₙ)/2 when using a mid-point quadrature rule.\n

This method updates counter `k.n` and entries to its history arrays at the nᵗʰ array location in the `k` data structure; specifically: stretch `k.λ[n]` and its rate `k.λ′[n]`, plus strain `k.ϵ[n]` and its rate `k.ϵ′[n].` These fields are evaluated at either the end-point, i.e. at time tₙ, or at the mid-point, i.e. at time (tₙ₋₁ + tₙ)/2, according to the argument `midPtQuad` supplied to its constructor.
"""
function advance!(k::FiberKinematics, lambda′::PhysicalScalar)
    # Advance the counter.
    if k.n < k.N+1
        set!(k.n, get(k.n)+1)
    else
        msg = "The data structure is full and cannot accept further data."
        throw(ErrorException(msg))
    end
    n = get(k.n)

    # Convert to CGS units.
    λ′ = toCGS(lambda′)

    # Verify the input.
    if λ′.units ≠ TIME_RATE
        msg = "Fiber stretch rate lambda′ must have units of reciprocal time."
        throw(ErrorException(msg))
    end
    k.λ′[n] = λ′

    # Integrate fiber stretch rate using a backward difference formula (BDF).
    if k.t[2] ≈ 0.5k.dt
        # Integrated for nodes located at the mid-point of each time step.
        if n == 2
            k.λ[2] = k.λ[1] + 0.5k.λ′[2]*k.dt
        elseif n == 3
            λ₁ = k.λ[1] - 0.5k.λ′[2]*k.dt
            k.λ[3] = (4/3)*k.λ[2] - (1/3)*λ₁ + (2/3)*k.λ′[3]*k.dt
        elseif n == 4
            λ₁ = k.λ[1] - 0.5k.λ′[2]*k.dt
            k.λ[4] = ((18/11)*k.λ[3] - (9/11)*k.λ[2] + (2/11)*λ₁
                   + (6/11)*k.λ′[4]*k.dt)
        else
            k.λ[n] = ((18/11)*k.λ[n-1] - (9/11)*k.λ[n-2] + (2/11)*k.λ[n-3]
                   + (6/11)*k.λ′[n]*k.dt)
        end
    else
        # Integrated for nodes located at the end-point of each time step.
        if n == 2
            k.λ[2] = k.λ[1] + k.λ′[2]*k.dt
        elseif n == 3
            k.λ[3] = (4/3)*k.λ[2] - (1/3)*k.λ[1] + (2/3)*k.λ′[3]*k.dt
        else
            k.λ[n] = ((18/11)*k.λ[n-1] - (9/11)*k.λ[n-2] + (2/11)*k.λ[n-3]
                   + (6/11)*k.λ′[n]*k.dt)
        end
    end

    # Compute the current strain and its rate.
    k.ϵ[n]  = PhysicalScalar(log(k.λ[n]/k.λᵣ), DIMENSIONLESS)
    k.ϵ′[n] = k.λ′[n] / k.λ[n]

    return nothing
end # advance!

"""
Method:\n
    update!(k::FiberKinematics, lambda′::PhysicalScalar)\n
Method `update!` refines a solution at step `n` by re-integrating its governing differential equation for fiber stretch, thereby allowing for iterative improvements to be made on stretch rate `lambda′` from an external algorithm, e.g., a finite element engine. There is no need to call `update!` unless `λ′` is being iteratively refined at step `n`, e.g., by some external optimization process. Here λ′ = dλ(k.t[n])/dt.
"""
function update!(k::FiberKinematics, lambda′::PhysicalScalar)
    if k.n > 1
        set!(k.n, get(k.n)-1)
        advance!(k, lambda′)
    end
    return nothing
end # update!
