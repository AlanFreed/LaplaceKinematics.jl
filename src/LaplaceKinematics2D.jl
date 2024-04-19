# Matrices P2D₁ and P2D₂ are the two possible permutation matrices in 2-space.

one  = PhysicalScalar(1.0, DIMENSIONLESS)
P2D₁ = PhysicalTensor(2, 2, DIMENSIONLESS)  # (𝕚, 𝕛) ↦ (𝕖₁, 𝕖₂)
P2D₁[1,1] = one
P2D₁[2,2] = one
P2D₂ = PhysicalTensor(2, 2, DIMENSIONLESS)  # (𝕛, 𝕚) ↦ (𝕖₁, 𝕖₂)
P2D₂[1,2] = one
P2D₂[2,1] = one

# -----------------------------------------------------------------------------

"""
Type:\n
    MembraneKinematics\n
        # Properties of the arrays.
        dt      time increment separating neighboring nodes\n
        N       total node count for traversing a solution path\n
        n       a counter that ratchets from 1 to N+1\n

        # 2D Laplace stretch attributes for a reference deformation of κ₀ ↦ κᵣ.\n
        aᵣ      reference elongation (stretch) in 𝕚 direction\n
        bᵣ      reference elongation (stretch) in 𝕛 direction\n
        γᵣ      reference in-plane shear in (𝕚, 𝕛) plane in 𝕚 direction\n

        # History arrays of length N+1 for holding the kinematic fields.
        # Initial values/conditions are stored in array location [1].

        # Array of nodal times.
        t       times at the solution nodes, i.e., the tₙ\n

        # Unpivoted 2D deformation gradients for deformation κ₀ ↦ κₙ in (𝕚, 𝕛).
        F       deformation gradients at tₙ: Fₙ, κ₀ ↦ κₙ in (𝕚, 𝕛)\n
        F′      deformation gradient rates at tₙ: dFₙ/dtₙ, κₙ in (𝕚, 𝕛)\n
        motion  the motion case that applies at time tₙ:\n
                    1) with pure shear, no co-ordinate pivoting\n
                    2) with pure shear and co-ordinate pivoting\n
                    3) with rigid-body rotation, no pivoting\n
                    4) with rigid-body rotation and pivoting\n

        # Gram angles of rotation and their rates at tₙ, mapped to (𝕚, 𝕛)\n
        ωₙ      angular rotations ωₙ at tₙ:\n
                (𝕖₁, 𝕖₂) out of (𝕚, 𝕛) whenever P = P₁\n
                (𝕖₂, 𝕖₁) out of (𝕚, 𝕛) whenever P = P₂\n
        ω′ₙ     angular rates of rotation at tₙ, i.e., dωₙ/dtₙ\n

        # 2D Laplace stretch attributes for κ₀ ↦ κₙ, mapped to (𝕚, 𝕛)\n
        aₙ      elongations in 𝕚 direction at tₙ\n
        bₙ      elongations in 𝕛 direction at tₙ\n
        γₙ      in-plane shears in (𝕚, 𝕛) plane in 𝕚 direction at tₙ\n

        # 2D Laplace stretch-rate attributes at κₙ, mapped to (𝕚, 𝕛)\n
        a′ₙ     elongation rates in 𝕚 direction at tₙ: daₙ/dt\n
        b′ₙ     elongation rates in 𝕛 direction at tₙ: dbₙ/dt\n
        γ′ₙ     in-plane shear rates at tₙ in (𝕚, 𝕛) plane in 𝕚 direction: dγₙ/dt\n

        # 2D Laplace strain attributes for deformation κᵣ ↦ κₙ\n
        δ       strains of dilation at tₙ: δ\n
        ϵ       strains of squeeze at tₙ: ϵ\n
        γ       strains of shear at tₙ: γ\n

        # 2D Laplace strain-rate attributes at configuration κₙ\n
        δ′      strain rates of dilation at tₙ: dδ/dt\n
        ϵ′      strain rates of squeeze at tₙ: dϵ/dt\n
        γ′      strain rates of shear at tₙ: dγ/dt\n
`MembraneKinematics` is a data structure that contains the various Lagrangian fields associated with a Laplace (upper triangular) measure for stretch in an isochoric two-space. The arrays that comprise this data structure allow for a history of these kinematic variables to be compiled for later retrieval and use, e.g., for constitutive analysis, for graphing, etc. Fields of this type are evaluated in the user's co-ordinate frame (𝕚, 𝕛).
"""
struct MembraneKinematics
    # Properties of the arrays.
    dt::PhysicalScalar           # time step separating neighboring entries
    N::Int64                     # total number of steps or grid points
    n::MInteger                  # a counter that ratchets from 1 to N+1

    # 2D Laplace stretch attributes for a reference deformation of κ₀ ↦ κᵣ.
    aᵣ::PhysicalScalar           # reference elongation (stretch) in 𝕚 direction
    bᵣ::PhysicalScalar           # reference elongation (stretch) in 𝕛 direction
    γᵣ::PhysicalScalar           # reference in-plane shear in (𝕚,𝕛) plane

    # History arrays of length N+1 for holding the kinematic fields.
    # Initial values/conditions are stored in array location [1].

    # Array of the independent variable, viz., array of nodal times.
    t::ArrayOfPhysicalScalars    # time at the solution nodes, i.e., at tₙ

    # Unpivoted 2D deformation gradients for a deformation of κ₀ ↦ κₙ in (𝕚, 𝕛).
    F::ArrayOfPhysicalTensors    # deformation gradients at tₙ: Fₙ κ₀ ↦ κₙ
    F′::ArrayOfPhysicalTensors   # deformation gradient rates at tₙ: dFₙ/dtₙ
    motion::Vector{Int64}        # the motion case that applies at time tₙ:
                                 # 1) with pure shear, no co-ordinate pivoting
                                 # 2) with pure shear and co-ordinate pivoting
                                 # 3) with rigid-body rotation, no pivoting
                                 # 4) with rigid-body rotation and pivoting

    # Gram angles of rotation and their rates at tₙ, mapped to (𝕚, 𝕛)\n
    ωₙ::ArrayOfPhysicalScalars   # angular rotations at tₙ: ωₙ
                                 # (𝕖₁, 𝕖₂) out of (𝕚, 𝕛) whenever P = P₁
                                 # (𝕖₂, 𝕖₁) out of (𝕚, 𝕛) whenever P = P₂
    ω′ₙ::ArrayOfPhysicalScalars  # angular rates of rotation at tₙ: dωₙ/dtₙ

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
    ϵ::ArrayOfPhysicalScalars    # strains of squeeze at tₙ: ϵ
    γ::ArrayOfPhysicalScalars    # strains of shear at tₙ: γ

    # 2D Laplace strain-rate attributes at configuration κₙ, mapped to (𝕚, 𝕛)
    δ′::ArrayOfPhysicalScalars   # strain rates of dilation at tₙ: dδ/dt
    ϵ′::ArrayOfPhysicalScalars   # strain rates of squeeze at tₙ: dϵ/dt
    γ′::ArrayOfPhysicalScalars   # strain rates of shear at tₙ: dγ/dt

    # Internal constructors.

"""
    Constructor:\n
        k = MembraneKinematics(dTime, N, midPtQuad, aRef, bRef, γRef, PRef)\n
    Returns a new data structure `k` of type `MembraneKinematics` that holds a variety of kinematic fields. Arguments include: (i) A differential step in time `dTime` that separates neighboring nodes, which themselves are taken to be uniformly spaced over time. (ii) The number of grid points or nodes `N` where solutions are to be computed. The data arrays are of length N+1 with initial values/conditions being stored at location [1] in these arrays, e.g., t[1] = 0. (iii) A boolean flag `midPtQuad` that, if true, implies the nodal spacing is for a mid-point quadrature rule, and if false, implies the nodal spacing is for an end-point quadrature rule. This determines how the array of independent times is to be populated. (iv) The reference Laplace stretch attributes, viz., `aRef`, `bRef` and `γRef`, against which isochoric strains are to be established so that ϵ(aRef, bRef, γRef) = 0, with the membrane's initial deformation gradient F₀, associated with some initial configuration κ₀, being assigned the identity matrix I with an outcome being that ϵ(a₀, b₀, γ₀) need not equal 0. And (v) if γRef is to be a shearing in the 𝕚 direction then `PRef` is to equal 1, else if γᵣ is to be a shearing in the 𝕛 direction then `PRef` is to equal 2, where PRef denotes which permutation matrix it to be applied in the reference configuration.
"""
    function MembraneKinematics(dTime::PhysicalScalar, N::Int64, midPtQuad::Bool, aRef::PhysicalScalar, bRef::PhysicalScalar, γRef::PhysicalScalar, PRef::Int64)

        # Verify inputs.
        if (PRef < 1) || (PRef > 2)
            msg = "The permutation case for the reference configuration can be either 1 or 2."
            throw(ErrorException(msg))
        end

        # Convert all passed variables to CGS units.
        dt = toCGS(dTime)
        if PRef == 1
            aᵣ = toCGS(aRef)
            bᵣ = toCGS(bRef)
        else
            aᵣ = toCGS(bRef)
            bᵣ = toCGS(aRef)
        end
        gᵣ = toCGS(γRef)

        # Continue verification.
        if dt.units ≠ TIME
            msg = "The supplied time increment dt does not have units of time."
            throw(ErrorException(msg))
        end
        dtₘᵢₙ = PhysicalScalar(Float64(eps(Float64)), TIME)
        if dt < dtₘᵢₙ
            msg = "The supplied time increment dt must be positive valued."
            throw(ErrorException(msg))
        end
        if N < 1
            msg = string("Solution arrays must have a positive length.")
            throw(ErrorException(msg))
        end
        if !isDimensionless(aᵣ)
            msg = "The supplied reference stretch aRef is not dimensionless."
            throw(ErrorException(msg))
        end
        λₘᵢₙ = PhysicalScalar(Float64(eps(Float16)), DIMENSIONLESS)
        if aᵣ < λₘᵢₙ
            msg = "The supplied reference stretch aRef must be positive valued."
            throw(ErrorException(msg))
        end
        if !isDimensionless(bᵣ)
            msg = "The supplied reference stretch bRef is not dimensionless."
            throw(ErrorException(msg))
        end
        if bᵣ < λₘᵢₙ
            msg = "The supplied reference stretch bRef must be positive valued."
            throw(ErrorException(msg))
        end
        if !isDimensionless(gᵣ)
            msg = "The supplied reference in-plane shear γRef is not dimensionless."
            throw(ErrorException(msg))
        end

        # Establish the counter.
        n  = MInteger(1)

        # Create and populate an array for nodal times.
        t  = ArrayOfPhysicalScalars(N+1, TIME)
        if midPtQuad
            # Assign times for a mid-point quadrature rule.
            for n in 1:N
                t[n+1] = (n-1)*dt + 0.5*dt
            end
        else
            # Assign times for an end-point quadrature rule.
            for n in 1:N
                t[n+1] = n*dt
            end
        end

        # Create data arrays for the independent kinematic fields.
        F  = ArrayOfPhysicalTensors(N+1, 2, 2, DIMENSIONLESS)
        F′ = ArrayOfPhysicalTensors(N+1, 2, 2, TIME_RATE)

        # Assign values to deformation gradient in its initial configuration κ₀.
        F₀ = PhysicalTensor(2, 2, DIMENSIONLESS)
        F₀[1,1] = PhysicalScalar(1.0, DIMENSIONLESS)
        F₀[2,2] = PhysicalScalar(1.0, DIMENSIONLESS)
        F[1]  = F₀
        F′[1] = PhysicalTensor(2, 2, TIME_RATE)

        # Assign initial conditions to the Laplace stretch attributes.
        a₀ = PhysicalScalar(1.0, DIMENSIONLESS)
        b₀ = PhysicalScalar(1.0, DIMENSIONLESS)
        γ₀ = PhysicalScalar(DIMENSIONLESS)
        ω₀ = PhysicalScalar(DIMENSIONLESS)

        # Assign initial conditions to the Laplace stretch rate attributes.
        a′₀ = PhysicalScalar(TIME_RATE)
        b′₀ = PhysicalScalar(TIME_RATE)
        γ′₀ = PhysicalScalar(TIME_RATE)
        ω′₀ = PhysicalScalar(TIME_RATE)

        # Data array that holds the various cases of motion.
        motion = zeros(Int64, N+1)
        motion[1] = PRef

        # Data arrays that hold the Gram rotations and their rates: κ₀ ↦ κₙ.
        ωₙ  = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ω′ₙ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        ωₙ[1]  = ω₀
        ω′ₙ[1] = ω′₀

        # Data arrays for Laplace stretch attributes and their rates: κᵣ ↦ κₙ.
        aₙ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        bₙ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        γₙ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        aₙ[1] = a₀/aᵣ
        bₙ[1] = b₀/bᵣ
        γₙ[1] = (aᵣ/bᵣ)*(γ₀ - gᵣ)
        a′ₙ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        b′ₙ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        γ′ₙ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        a′ₙ[1] = a′₀/aᵣ
        b′ₙ[1] = b′₀/bᵣ
        γ′ₙ[1] = (aᵣ/bᵣ)*γ′₀

        # Data arrays for the thermodynamic strains and their rates: κᵣ ↦ κₙ.
        δ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ϵ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        γ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        δ[1] = PhysicalScalar(0.5*log((a₀/aᵣ)*(b₀/bᵣ)), DIMENSIONLESS)
        ϵ[1] = PhysicalScalar(0.5*log((a₀/aᵣ)*(bᵣ/b₀)), DIMENSIONLESS)
        γ[1] = γ₀ - gᵣ
        δ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        ϵ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        γ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        δ′[1] = 0.5*(a′₀/a₀ + b′₀/b₀)
        ϵ′[1] = 0.5*(a′₀/a₀ - b′₀/b₀)
        γ′[1] = γ′₀

        # Create and return a new data structure for Laplace kinematics in 2D.
        new(dt, N, n, aᵣ, bᵣ, gᵣ, t, F, F′, motion, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ϵ, γ, δ′, ϵ′, γ′)
    end

    # Internal constructor used by JSON3.

    function MembraneKinematics(dt::PhysicalScalar, N::Int64, n::MInteger, aᵣ::PhysicalScalar, bᵣ::PhysicalScalar, γᵣ::PhysicalScalar, t::ArrayOfPhysicalScalars, F::ArrayOfPhysicalTensors, F′::ArrayOfPhysicalTensors, motion::Vector{Int64}, ωₙ::ArrayOfPhysicalScalars, ω′ₙ::ArrayOfPhysicalScalars, aₙ::ArrayOfPhysicalScalars, bₙ::ArrayOfPhysicalScalars, γₙ::ArrayOfPhysicalScalars, a′ₙ::ArrayOfPhysicalScalars, b′ₙ::ArrayOfPhysicalScalars, γ′ₙ::ArrayOfPhysicalScalars, δ::ArrayOfPhysicalScalars, ϵ::ArrayOfPhysicalScalars, γ::ArrayOfPhysicalScalars, δ′::ArrayOfPhysicalScalars, ϵ′::ArrayOfPhysicalScalars, γ′::ArrayOfPhysicalScalars)

        new(dt, N, n, aᵣ, bᵣ, γᵣ, t, F, F′, motion, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ϵ, γ, δ′, ϵ′, γ′)
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
    ϵ   = copy(k.ϵ)
    ϵ′  = copy(k.ϵ′)
    γ   = copy(k.γ)
    γ′  = copy(k.γ′)
    return MembraneKinematics(dt, N, n, aᵣ, bᵣ, γᵣ, t, F, F′, motion, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ϵ, γ, δ′, ϵ′, γ′)
end

function Base.:(deepcopy)(k::MembraneKinematics)::MembraneKinematics
    dt  = deepcopy(k.dt)
    N   = deepcopy(k.N)
    n   = deepcopy(k.n)
    aᵣ  = deepcopy(k.aᵣ)
    bᵣ  = deepcopy(k.bᵣ)
    γᵣ  = deepcopy(k.γᵣ)
    t   = deepcopy(k.t)
    F   = deepcopy(k.F)
    F′  = deepcopy(k.F′)
    motion = deepcopy(k.motion)
    ωₙ  = deepcopy(k.ωₙ)
    ω′ₙ = deepcopy(k.ω′ₙ)
    aₙ  = deepcopy(k.aₙ)
    a′ₙ = deepcopy(k.a′ₙ)
    bₙ  = deepcopy(k.bₙ)
    b′ₙ = deepcopy(k.b′ₙ)
    γₙ  = deepcopy(k.γₙ)
    γ′ₙ = deepcopy(k.γ′ₙ)
    δ   = deepcopy(k.δ)
    δ′  = deepcopy(k.δ′)
    ϵ   = deepcopy(k.ϵ)
    ϵ′  = deepcopy(k.ϵ′)
    γ   = deepcopy(k.γ)
    γ′  = deepcopy(k.γ′)
    return MembraneKinematics(dt, N, n, aᵣ, bᵣ, γᵣ, t, F, F′, motion, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ϵ, γ, δ′, ϵ′, γ′)
end

# The histories of MembraneKinematics are to be graphed, not printed, so a toString method is not provided for objects of this type.

# Methods for storing and retrieving a MembraneKinematics data structure to and from a file.

StructTypes.StructType(::Type{MembraneKinematics}) = StructTypes.Struct()

"""
Method:\n
    toFile(k::LaplaceKinematics.MembraneKinematics, json_stream::IOStream)\n
Writes data structure `k` to the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONWriter(<my_dir_path>, <my_file_name>)\n
    ...\n
    toFile(k, json_stream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream)\n
where <my_dir_path> is the path to your working directory wherein the file <my_file_name> to be written to either exists or will be created, which must have a .json extension.
"""
function toFile(k::MembraneKinematics, json_stream::IOStream)
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
    fromFile(k::LaplaceKinematics.MembraneKinematics, json_stream::IOStream)\n
Reads a MembraneKinematics data structure from the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONReader(<my_dir_path>, <my_file_name>)\n
    ...\n
    k = fromFile(LaplaceKinematics.MembraneKinematics, json_stream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream)\n
which returns `k,` which is an object of type MembraneKinematics. Here <my_dir_path> is the path to your working directory wherein the file <my_file_name> to be read from must exist, which is to have a .json extension.
"""
function fromFile(::Type{MembraneKinematics}, json_stream::IOStream)::MembraneKinematics
    if isopen(json_stream)
        k = JSON3.read(readline(json_stream), MembraneKinematics)
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    return k
end

# Methods that serve as a solver for objects of type MembraneKinematics.

"""
Method:\n
    advance!(k, dF)\n
Method `advance!` moves a solution from previous step `n-1` to current step `n` along a solution path of N solution nodes by integrating its governing differential equation for the 2D deformation gradient using a backward difference formula (BDF). Supplied is the membrane's time rate-of-change in its deformation gradient `dF` evaluated in the user's co-ordinate system (𝕚, 𝕛). The initial condition of this differential equation is taken to be F₀ = I, i.e., k.F[1] = F₀ = I where I is the identity matrix. For a time-step interval of [tₙ₋₁, tₙ], dF = F′ = dF/dt associates with either time tₙ when using end-point quadrature, or with time (tₙ₋₁ + tₙ)/2 when using mid-point quadrature.\n

This method updates counter `k.n` and entries to its history arrays at the nᵗʰ array location in the `k` data structure; specifically: deformation gradient `k.F[n]` and its rate `k.F′[n]`, motion case `k.motion[n]`, Laplace stretch attributes `k.aₙ[n]`, `k.bₙ[n]`, `k.γₙ[n]` and `k.ωₙ[n]` and their rates `k.a′ₙ[n]`, `k.b′ₙ[n]`, `k.γ′ₙ[n]` and `k.ω′ₙ[n]`, plus Laplace strain attributes `k.δ[n]`, `k.ϵ[n]` and `k.γ[n]` and their rates `k.δ′[n]`, `k.ϵ′[n]` and `k.γ′[n]`, all mapped to the user's co-ordinate system whose base vectors are denoted as (𝕚, 𝕛). These fields are evaluated at either the end-point, i.e. at time tₙ, or at the mid-point, i.e. at time (tₙ₋₁ + tₙ)/2, according to the argument `midPtQuad` supplied to its constructor.
"""
function advance!(k::MembraneKinematics, dF::PhysicalTensor)
    # Advance the counter.
    if k.n < k.N+1
        set!(k.n, get(k.n)+1)
    else
        msg = "The data structure is full and cannot accept further data."
        throw(ErrorException(msg))
    end
    n = get(k.n)

    # Convert the passed variable to CGS units.
    F′ = toCGS(dF)

    # Verify inputs.
    if (F′.matrix.rows ≠ 2) || (F′.matrix.cols ≠ 2)
        msg = "Supplied deformation gradient rate dF must be a 2x2 matrix."
        throw(DimensionMismatch(msg))
    end
    if F′.units ≠ TIME_RATE
        msg = "Supplied deformation gradient rate dF must have units of reciprocal time."
        throw(ErrorException(msg))
    end
    k.F′[n] = F′

    # Advance the fields, i.e., insert new values into the data arrays.

    # Integrate F′ using a backward difference formula (BDF).
    Fₙ = PhysicalTensor(2, 2, DIMENSIONLESS)
    if k.t[2] ≈ 0.5k.dt # Use quadrature nodes that are located at mid span.
        if n == 2
            F₁ = k.F[1]
        elseif n == 3
            F₁ = PhysicalTensor(2, 2, DIMENSIONLESS)
            F₂ = k.F[2]
            F′₂ = k.F′[2]
            Fₙ₋₂ = k.F[1]
        elseif n == 4
            F₁ = PhysicalTensor(2, 2, DIMENSIONLESS)
            F₂ = k.F[2]
            F₃ = k.F[3]
            F′₂ = k.F′[2]
            Fₙ₋₃ = k.F[1]
        else
            Fₙ₋₁ = k.F[n-1]
            Fₙ₋₂ = k.F[n-2]
            Fₙ₋₃ = k.F[n-3]
        end
        for i in 1:2
            for j in 1:2
                if n == 2
                    Fₙ[i,j] = F₁[i,j] + 0.5*F′[i,j]*k.dt
                elseif n == 3
                    F₁[i,j] = Fₙ₋₂[i,j] - 0.5*F′₂[i,j]*k.dt
                    Fₙ[i,j] = (4/3)*F₂[i,j] - (1/3)*F₁[i,j] + (2/3)*F′[i,j]*k.dt
                elseif n == 4
                    F₁[i,j] = Fₙ₋₃[i,j] - 0.5*F′₂[i,j]*k.dt
                    Fₙ[i,j] = ((18/11)*F₃[i,j] - (9/11)*F₂[i,j]
                            + (2/11)*F₁[i,j] + (6/11)*F′[i,j]*k.dt)
                else
                    Fₙ[i,j] = ((18/11)*Fₙ₋₁[i,j] - (9/11)*Fₙ₋₂[i,j]
                            + (2/11)*Fₙ₋₃[i,j] + (6/11)*F′[i,j]*k.dt)
                end
            end
        end
    else # Use quadrature nodes that are located at the end of span.
        if n == 2
            F₁ = k.F[1]
        elseif n == 3
            F₂ = k.F[2]
            F₁ = k.F[1]
        else
            Fₙ₋₁ = k.F[n-1]
            Fₙ₋₂ = k.F[n-2]
            Fₙ₋₃ = k.F[n-3]
        end
        for i in 1:2
            for j in 1:2
                if n == 2
                    Fₙ[i,j] = F₁[i,j] + F′[i,j]*k.dt
                elseif n == 3
                    Fₙ[i,j] = (4/3)*F₂[i,j] - (1/3)*F₁[i,j] + (2/3)*F′[i,j]*k.dt
                else
                    Fₙ[i,j] = ((18/11)*Fₙ₋₁[i,j] - (9/11)*Fₙ₋₂[i,j]
                            + (2/11)*Fₙ₋₃[i,j] + (6/11)*F′[i,j]*k.dt)
                end
            end
        end
    end
    k.F[n]  = Fₙ
    k.F′[n] = F′

    # Get attributes for deformation F and rate of deformation F′ gradients.
    x  = Fₙ[1,1]
    y  = Fₙ[2,2]
    x′ = F′[1,1]
    y′ = F′[2,2]

    # Establish the Gram and Laplace attributes, and their rates.
    if (Fₙ[2,1] ≈ 0.0) || (sign(Fₙ[1,2]) == sign(Fₙ[2,1]))
        # The deformation has an attribute that is a pure shear.
        # g is this pure shear.
        if ((abs(Fₙ[1,2])/Fₙ[2,2] > abs(Fₙ[2,1])/Fₙ[1,1])
            || (abs(Fₙ[1,2])/Fₙ[2,2] ≈ abs(Fₙ[2,1])/Fₙ[1,1]))
            # Co-ordinates are indexed as supplied--the default condition.
            case = 1
            # Pure shear contributions.
            g  = Fₙ[2,1] / Fₙ[1,1]
            g′ = (Fₙ[1,1]*F′[2,1] - Fₙ[2,1]*F′[1,1]) / (Fₙ[1,1]*Fₙ[1,1])
            # Simple shear contributions.
            G  = (Fₙ[1,1]*Fₙ[1,2] - Fₙ[2,2]*Fₙ[2,1]) / (Fₙ[1,1]*Fₙ[2,2])
            G′ = (Fₙ[2,2]*F′[1,2] - Fₙ[1,2]*F′[2,2]) / (Fₙ[2,2]*Fₙ[2,2]) - g′
        else
            # The Gram co-ordinate system is left handed.
            case = 2
            # Pure shear contributions.
            g  = Fₙ[1,2] / Fₙ[2,2]
            g′ = (Fₙ[2,2]*F′[1,2] - Fₙ[1,2]*F′[2,2]) / (Fₙ[2,2]*Fₙ[2,2])
            # Simple shear contributions.
            G  = -(Fₙ[1,1]*Fₙ[1,2] - Fₙ[2,2]*Fₙ[2,1]) / (Fₙ[1,1]*Fₙ[2,2])
            G′ = (Fₙ[1,1]*F′[2,1] - Fₙ[2,1]*F′[1,1]) / (Fₙ[1,1]*Fₙ[1,1]) - g′
        end
        # Laplace stretch attributes in the (𝕚, 𝕛) co-ordinate frame.
        aₙ = x * sqrt(1+g*g)
        bₙ = y * (1 - g*(g+G)) / sqrt(1+g*g)
        γₙ = (y/x) * (2g+G) / (1+g*g)
        ωₙ = PhysicalScalar(atan(g), DIMENSIONLESS)

        # Rates of Laplace attributes in the (𝕚, 𝕛) co-ordinate frame.
        a′ₙ = aₙ*(x′/x + g*g′/(1+g*g))
        b′ₙ = bₙ*(y′/y - g*g′/(1+g*g)) - y*((2g+G)*g′ + g*G′)/sqrt(1+g*g)
        γ′ₙ = γₙ*(y′/y - x′/x - 2g*g′/(1+g*g)) + (y/x)*(2g′+G′)/(1+g*g)
        ω′ₙ = g′/(1+g*g)
    else
        # The deformation has an attribute that is a rigid-body rotation.
        # g is this rigid-body rotation.
        if ((abs(Fₙ[1,2])/Fₙ[2,2] > abs(Fₙ[2,1])/Fₙ[1,1])
            || (abs(Fₙ[1,2])/Fₙ[2,2] ≈ abs(Fₙ[2,1])/Fₙ[1,1]))
            # The Gram co-ordinate system is right handed.
            case = 3
            # Rigid-body rotation contributions.
            g  = -Fₙ[2,1] / Fₙ[1,1]
            g′ = -(Fₙ[1,1]*F′[2,1] - Fₙ[2,1]*F′[1,1]) / (Fₙ[1,1]*Fₙ[1,1])
            # Angle of rigid-body rotation in the (𝕚, 𝕛) co-ordinate frame.
            ωₙ  = PhysicalScalar(-atan(g), DIMENSIONLESS)
            ω′ₙ = -g′/(1+g*g)
        else
            # The Gram co-ordinate system is left handed.
            case = 4
            # Rigid-body rotation contributions.
            g  = -Fₙ[1,2] / Fₙ[2,2]
            g′ = -(Fₙ[2,2]*F′[1,2] - Fₙ[1,2]*F′[2,2]) / (Fₙ[2,2]*Fₙ[2,2])
            # Angle of rigid-body rotation in the (𝕚, 𝕛) co-ordinate frame.
            ωₙ  = PhysicalScalar(atan(g), DIMENSIONLESS)
            ω′ₙ = g′/(1+g*g)
        end
        # G is a simple shear.
        G  = (Fₙ[1,1]*Fₙ[1,2] + Fₙ[2,2]*Fₙ[2,1]) / (Fₙ[1,1]*Fₙ[2,2])
        G′ = ((Fₙ[2,2]*F′[1,2] - Fₙ[1,2]*F′[2,2]) / (Fₙ[2,2]*Fₙ[2,2]) +
            (Fₙ[1,1]*F′[2,1] - Fₙ[2,1]*F′[1,1]) / (Fₙ[1,1]*Fₙ[1,1]))

        # Laplace attributes in the (𝕚, 𝕛) co-ordinate frame.
        aₙ = x * sqrt(1+g*g)
        bₙ = y * (1 + g*(g+G)) / sqrt(1+g*g)
        γₙ = (y/x) * G / (1+g*g)

        # Rates of Laplace attributes in the (𝕚, 𝕛) co-ordinate frame.
        a′ₙ = aₙ*(x′/x + g*g′/(1+g*g))
        b′ₙ = bₙ*(y′/y - g*g′/(1+g*g)) + y*((2g+G)*g′ + g*G′)/sqrt(1+g*g)
        γ′ₙ = γₙ*(y′/y - x′/x - 2g*g′/(1+g*g)) + (y/x)*G′/(1+g*g)
    end

    # Advance the data array that holds pivoting cases.
    k.motion[n] = case

    # Advance the Laplace stretch attributes and their rates for κᵣ ↦ κₙ.
    k.aₙ[n]  = aₙ/k.aᵣ
    k.bₙ[n]  = bₙ/k.bᵣ
    k.γₙ[n]  = (k.aᵣ/k.bᵣ)*(γₙ - k.γᵣ)
    k.ωₙ[n]  = ωₙ
    k.a′ₙ[n] = a′ₙ/k.aᵣ
    k.b′ₙ[n] = b′ₙ/k.bᵣ
    k.γ′ₙ[n] = (k.aᵣ/k.bᵣ)*γ′ₙ
    k.ω′ₙ[n] = ω′ₙ

    # Advance the thermodynamic Laplace strains and their rates for κᵣ ↦ κₙ.
    k.δ[n]  = PhysicalScalar(0.5*log(k.aₙ[n]*k.bₙ[n]), DIMENSIONLESS)
    k.ϵ[n]  = PhysicalScalar(0.5*log(k.aₙ[n]/k.bₙ[n]), DIMENSIONLESS)
    k.γ[n]  = γₙ - k.γᵣ
    k.δ′[n] = 0.5*(a′ₙ/aₙ + b′ₙ/bₙ)
    k.ϵ′[n] = 0.5*(a′ₙ/aₙ - b′ₙ/bₙ)
    k.γ′[n] = γ′ₙ

    return nothing
end # advance!

"""
Method:\n
    update!(k, dF)\n
Method `update!` provides a capability to refine a solution at step `n` by re-integrating the deformation gradient rate `dF`, thereby allowing for iterative improvements to be made on the deformation rate `dF` from an external algorithm, e.g., a finite element engine. There is no need to call `update!` unless `dF` is being iteratively refined at step `n` by some external optimization process.
"""
function update!(k::MembraneKinematics, dF::PhysicalTensor)
    if k.n > 1
        set!(k.n, get(k.n)-1)
        advance!(k, dF)
    end
    return nothing
end # update!
