#=
Created on Mon 22 Nov 2021
updated on Mon 27 Nov 2023
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

        # Array of nodal times.
        t       times at the solution nodes, i.e., the tₙ\n

        # Unpivoted 2D deformation gradients for a deformation of κ₀ ↦ κₙ in (𝕚, 𝕛),
        # where F₃₃, the third (thickness) direction, makes deformation isochoric.
        F       deformation gradients at tₙ: Fₙ, κ₀ ↦ κₙ in (𝕚, 𝕛)\n
        F′      deformation gradient rates at tₙ: dFₙ/dtₙ, κₙ in (𝕚, 𝕛)\n
        P       permutation case at tₙ: i.e., i in Pᵢ, i ∈ {1, 2}\n
                where {𝕖₁ 𝕖₂} = {𝕚 𝕛}[Pᵢ], i ∈ {1, 2}\n

        # 2D Laplace stretch attributes for reference deformation κ₀ ↦ κᵣ\n
        aᵣ      reference elongation in 𝕚 direction\n
        bᵣ      reference elongation in 𝕛 direction\n
        γᵣ      reference in-plane shear in (𝕚, 𝕛) plane\n

        # Gram angles of rotation and their rates at tₙ mapped to (𝕚, 𝕛)\n
        ωₙ      angular rotations ωₙ at tₙ:\n
                (𝕖₁, 𝕖₂) out of (𝕚, 𝕛) whenever P = P₁\n
                (𝕖₂, 𝕖₁) out of (𝕚, 𝕛) whenever P = P₂\n
        ω′ₙ     angular rates of rotation at tₙ, i.e., dωₙ/dtₙ\n

        # 2D Laplace stretch attributes for κ₀ ↦ κₙ, mapped to (𝕚, 𝕛)\n
        aₙ      elongations in 𝕚 direction at tₙ\n
        bₙ      elongations in 𝕛 direction at tₙ\n
        γₙ      in-plane shears in (𝕚, 𝕛) plane at tₙ\n

        # 2D Laplace stretch-rate attributes at κₙ, mapped to (𝕚, 𝕛)\n
        a′ₙ     elongation rates in 𝕚 direction at tₙ: daₙ/dt\n
        b′ₙ     elongation rates in 𝕛 direction at tₙ: dbₙ/dt\n
        γ′ₙ     in-plane shear rates at tₙ in (𝕚, 𝕛) plane: dγₙ/dt\n

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
    N::Integer                   # total number of steps or grid points
    n::MInteger                  # a counter that ratchets from 1 to N

    # Array of nodal times.
    t::ArrayOfPhysicalScalars    # times at the solution nodes, i.e., the tₙ

    # Unpivoted 2D deformation gradients for a deformation of κ₀ ↦ κₙ in (𝕚, 𝕛),
    # where F₃₃, the third (thickness) direction, makes deformation isochoric.
    F::ArrayOfPhysicalTensors    # deformation gradients at tₙ: Fₙ κ₀ ↦ κₙ
    F′::ArrayOfPhysicalTensors   # deformation gradient rates at tₙ: dFₙ/dtₙ
    P::Vector                    # permutation case at tₙ: i.e., i in Pᵢ,
                                 # where {𝕖₁ 𝕖₂} = {𝕚 𝕛}[Pᵢ], i ∈ {1, 2}

    # 2D Laplace stretch attributes for reference deformation κ₀ ↦ κᵣ
    aᵣ::PhysicalScalar           # reference elongation in 𝕚 direction
    bᵣ::PhysicalScalar           # reference elongation in 𝕛 direction
    γᵣ::PhysicalScalar           # reference in-plane shear in (𝕚,𝕛) plane

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

    # Internal constructors.

"""
    Constructor:\n
        k = MembraneKinematics(dt, N, aᵣ, bᵣ, γᵣ, F₀)\n
    Returns a new data structure `k` of type `MembraneKinematics` that holds a variety of kinematic fields. Arguments include: (i) a differential step in time `dt` that separates neighboring nodes, uniformly spaced in time, that belong to the data arrays that are to be populated, (ii) the number of grid points or nodes `N` where solutions are to be computed, (iii) the reference Laplace stretch attributes, viz., `aᵣ`, `bᵣ` and `γᵣ`, against which isochoric strains are to be established, and (iv) a deformation gradient `F₀` belonging to some initial configuration κ₀ evaluated in an user specified co-ordinate system with base vectors (𝕚, 𝕛). Laplace tensors are evaluated in a frame-indifferent co-ordinate system (𝕖₁, 𝕖₂) that are then mapped to the user's co-ordinate system (𝕚, 𝕛). It is in the (𝕚, 𝕛) co-ordinate system that the kinematic fields of this data structure are quantified in.
"""
    function MembraneKinematics(dt::PhysicalScalar, N::Integer, aᵣ::PhysicalScalar, bᵣ::PhysicalScalar, γᵣ::PhysicalScalar, F₀::PhysicalTensor)

        # Convert all passed variables to CGS units.
        d𝑡 = toCGS(dt)
        𝑎ᵣ = toCGS(aᵣ)
        𝑏ᵣ = toCGS(bᵣ)
        𝑔ᵣ = toCGS(γᵣ)
        𝐹₀ = toCGS(F₀)

        # Verify inputs.
        if d𝑡.units ≠ TIME
            msg = "The supplied time increment dt does not have units of time."
            throw(ErrorException(msg))
        end
        dtₘᵢₙ = PhysicalScalar(Float64(eps(Float64)), TIME)
        if d𝑡 < dtₘᵢₙ
            msg = "The supplied time increment dt must be positive valued."
            throw(ErrorException(msg))
        end
        if N < 1
            msg = string("Solution arrays must have a positive length.")
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑎ᵣ)
            msg = "The supplied reference stretch aᵣ is not dimensionless."
            throw(ErrorException(msg))
        end
        λₘᵢₙ = PhysicalScalar(Float64(eps(Float16)), DIMENSIONLESS)
        if 𝑎ᵣ < λₘᵢₙ
            msg = "The supplied reference stretch aᵣ must be positive valued."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑏ᵣ)
            msg = "The supplied reference stretch bᵣ is not dimensionless."
            throw(ErrorException(msg))
        end
        if 𝑏ᵣ < λₘᵢₙ
            msg = "The supplied reference stretch bᵣ must be positive valued."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑔ᵣ)
            msg = "The supplied reference in-plane shear γᵣ is not dimensionless."
            throw(ErrorException(msg))
        end
        if (𝐹₀.matrix.rows ≠ 2) || (𝐹₀.matrix.cols ≠ 2)
            msg = "The supplied deformation gradient F₀ must be a 2x2 matrix."
            throw(DimensionMismatch(msg))
        end
        if !isDimensionless(𝐹₀)
            msg = "The supplied deformation gradient F₀ is not dimensionless."
            throw(ErrorException(msg))
        end

        # Establish the counter.
        n  = MInteger(1)

        # Create data arrays for the independent kinematic fields.
        t  = ArrayOfPhysicalScalars(N+1, TIME)
        F  = ArrayOfPhysicalTensors(N+1, 2, 2, DIMENSIONLESS)
        F′ = ArrayOfPhysicalTensors(N+1, 2, 2, TIME_RATE)
        t[1]  = PhysicalScalar(TIME)
        F[1]  = 𝐹₀
        F′[1] = PhysicalTensor(2, 2, TIME_RATE)

        # Assign attributes for the deformation 𝐹 and rate of deformation 𝐹′ gradients.
        x  = 𝐹₀[1,1]
        y  = 𝐹₀[2,2]
        x′ = PhysicalScalar(TIME_RATE)
        y′ = PhysicalScalar(TIME_RATE)

        # Establish the Gram and Laplace attributes, and their rates.
        if sign(𝐹₀[1,2]) == sign(𝐹₀[2,1])
            # The deformation has an attribute that is a pure shear.
            # g is this pure shear.
            if ((abs(𝐹₀[1,2])/𝐹₀[2,2] > abs(𝐹₀[2,1])/𝐹₀[1,1])
                || (abs(𝐹₀[1,2])/𝐹₀[2,2] ≈ abs(𝐹₀[2,1])/𝐹₀[1,1]))
                # Co-ordinates are indexed as supplied--the default condition.
                case = 1
                # Pure shear contributions.
                g  = 𝐹₀[2,1] / 𝐹₀[1,1]
                g′ = PhysicalScalar(TIME_RATE)
                # Simple shear contributions.
                G  = (𝐹₀[1,1]*𝐹₀[1,2] - 𝐹₀[2,2]*𝐹₀[2,1]) / (𝐹₀[1,1]*𝐹₀[2,2])
                G′ = PhysicalScalar(TIME_RATE)
            else
                # The Gram co-ordinate system is left handed.
                case = 2
                # Pure shear contributions.
                g  = 𝐹₀[1,2] / 𝐹₀[2,2]
                g′ = PhysicalScalar(TIME_RATE)
                # Simple shear contributions.
                G  = -(𝐹₀[1,1]*𝐹₀[1,2] - 𝐹₀[2,2]*𝐹₀[2,1]) / (𝐹₀[1,1]*𝐹₀[2,2])
                G′ = PhysicalScalar(TIME_RATE)
            end
            # Laplace stretch attributes in the co-ordinate frame 𝕚 × 𝕛.
            a₀ = x * sqrt(1+g*g)
            b₀ = y * (1 - g*(g+G)) / sqrt(1+g*g)
            γ₀ = (y/x) * (2g+G) / (1+g*g)
            ω₀ = PhysicalScalar(atan(g), DIMENSIONLESS)

            # Rates of Laplace attributes in the co-ordinate frame 𝕚 × 𝕛.
            a′₀ = a₀*(x′/x + g*g′/(1+g*g))
            b′₀ = b₀*(y′/y - g*g′/(1+g*g)) - y*((2g+G)*g′ + g*G′)/sqrt(1+g*g)
            γ′₀ = γ₀*(y′/y - x′/x - 2g*g′/(1+g*g)) + (y/x)*(2g′+G′)/(1+g*g)
            ω′₀ = g′/(1+g*g)
        else
            # The deformation has an attribute that is a rigid-body rotation.
            # g is this rigid-body rotation.
            if ((abs(𝐹₀[1,2])/𝐹₀[2,2] > abs(𝐹₀[2,1])/𝐹₀[1,1])
                || (abs(𝐹₀[1,2])/𝐹₀[2,2] ≈ abs(𝐹₀[2,1])/𝐹₀[1,1]))
                # The Gram co-ordinate system is right handed.
                case = 3
                # Rigid-body rotation contributions.
                g  = -𝐹₀[2,1] / 𝐹₀[1,1]
                g′ = PhysicalScalar(TIME_RATE)
                # Angle of rigid-body rotation in the co-ordinate frame 𝕚 × 𝕛.
                ω₀  = PhysicalScalar(-atan(g), DIMENSIONLESS)
                ω′₀ = -g′/(1+g*g)
            else
                # The Gram co-ordinate system is left handed.
                case = 4
                # Rigid-body rotation contributions.
                g  = -𝐹₀[1,2] / 𝐹₀[2,2]
                g′ = PhysicalScalar(TIME_RATE)
                # Angle of rigid-body rotation in the co-ordinate frame 𝕚 × 𝕛.
                ω₀  = PhysicalScalar(atan(g), DIMENSIONLESS)
                ω′₀ = g′/(1+g*g)
            end
            # G is a simple shear.
            G  = (𝐹₀[1,1]*𝐹₀[1,2] + 𝐹₀[2,2]*𝐹₀[2,1]) / (𝐹₀[1,1]*𝐹₀[2,2])
            G′ = PhysicalScalar(TIME_RATE)

            # Laplace attributes in the co-ordinate frame 𝕖₁ × 𝕖₂.
            a₀ = x * sqrt(1+g*g)
            b₀ = y * (1 + g*(g+G)) / sqrt(1+g*g)
            γ₀ = (y/x) * G / (1+g*g)

            # Rates of Laplace attributes in the co-ordinate frame 𝕚 × 𝕛.
            a′₀ = a₀*(x′/x + g*g′/(1+g*g))
            b′₀ = b₀*(y′/y - g*g′/(1+g*g)) + y*((2g+G)*g′ + g*G′)/sqrt(1+g*g)
            γ′₀ = γ₀*(y′/y - x′/x - 2g*g′/(1+g*g)) + (y/x)*G′/(1+g*g)
        end

        # Data array that holds the pivoting cases.
        P = zeros(Int64, N+1)
        P[1] = case

        # Data arrays that hold the Gram rotations and their rates: κ₀ ↦ κₙ.
        ωₙ  = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ω′ₙ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        ωₙ[1]  = ω₀
        ω′ₙ[1] = ω′₀

        # Data arrays for Laplace stretch attributes and their rates: κᵣ ↦ κₙ.
        aₙ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        bₙ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        γₙ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        aₙ[1] = a₀/𝑎ᵣ
        bₙ[1] = b₀/𝑏ᵣ
        γₙ[1] = (𝑎ᵣ/𝑏ᵣ)*(γ₀ - 𝑔ᵣ)
        a′ₙ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        b′ₙ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        γ′ₙ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        a′ₙ[1] = a′₀/𝑎ᵣ
        b′ₙ[1] = b′₀/𝑏ᵣ
        γ′ₙ[1] = (𝑎ᵣ/𝑏ᵣ)*γ′₀

        # Data arrays for the thermodynamic strains and their rates: κᵣ ↦ κₙ.
        δ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        ϵ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        γ = ArrayOfPhysicalScalars(N+1, DIMENSIONLESS)
        δ[1] = PhysicalScalar(0.5log((a₀/𝑎ᵣ)*(b₀/𝑏ᵣ)), DIMENSIONLESS)
        ϵ[1] = PhysicalScalar(0.5log((a₀/𝑎ᵣ)*(𝑏ᵣ/b₀)), DIMENSIONLESS)
        γ[1] = γ₀ - 𝑔ᵣ
        δ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        ϵ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        γ′ = ArrayOfPhysicalScalars(N+1, TIME_RATE)
        δ′[1] = 0.5(a′₀/a₀ + b′₀/b₀)
        ϵ′[1] = 0.5(a′₀/a₀ - b′₀/b₀)
        γ′[1] = γ′₀

        # Create and return a new data structure for Laplace kinematics in 2D.
        new(d𝑡, N, n, t, F, F′, P, 𝑎ᵣ, 𝑏ᵣ, 𝑔ᵣ, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ϵ, γ, δ′, ϵ′, γ′)
    end

    # Internal constructor used by JSON3.

    function MembraneKinematics(dt::PhysicalScalar, N::Integer, n::MInteger, t::ArrayOfPhysicalScalars, F::ArrayOfPhysicalTensors, F′::ArrayOfPhysicalTensors, P::Vector, aᵣ::PhysicalScalar, bᵣ::PhysicalScalar, γᵣ::PhysicalScalar, ωₙ::ArrayOfPhysicalScalars, ω′ₙ::ArrayOfPhysicalScalars, aₙ::ArrayOfPhysicalScalars, bₙ::ArrayOfPhysicalScalars, γₙ::ArrayOfPhysicalScalars, a′ₙ::ArrayOfPhysicalScalars, b′ₙ::ArrayOfPhysicalScalars, γ′ₙ::ArrayOfPhysicalScalars, δ::ArrayOfPhysicalScalars, ϵ::ArrayOfPhysicalScalars, γ::ArrayOfPhysicalScalars, δ′::ArrayOfPhysicalScalars, ϵ′::ArrayOfPhysicalScalars, γ′::ArrayOfPhysicalScalars)

        new(dt, N, n, t, F, F′, P, aᵣ, bᵣ, γᵣ, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ϵ, γ, δ′, ϵ′, γ′)
    end
end # MembraneKinematics

# Methods

function Base.:(copy)(k::MembraneKinematics)::MembraneKinematics
    dt  = copy(k.dt)
    N   = copy(k.N)
    n   = copy(k.n)
    t   = copy(k.t)
    F   = copy(k.F)
    F′  = copy(k.F′)
    P   = copy(k.P)
    aᵣ  = copy(k.aᵣ)
    bᵣ  = copy(k.bᵣ)
    γᵣ  = copy(k.γᵣ)
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
    return MembraneKinematics(dt, N, n, t, F, F′, P, aᵣ, bᵣ, γᵣ, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ϵ, γ, δ′, ϵ′, γ′)
end

function Base.:(deepcopy)(k::MembraneKinematics)::MembraneKinematics
    dt  = deepcopy(k.dt)
    N   = deepcopy(k.N)
    n   = deepcopy(k.n)
    t   = deepcopy(k.t)
    F   = deepcopy(k.F)
    F′  = deepcopy(k.F′)
    P   = deepcopy(k.P)
    aᵣ  = deepcopy(k.aᵣ)
    bᵣ  = deepcopy(k.bᵣ)
    γᵣ  = deepcopy(k.γᵣ)
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
    return MembraneKinematics(dt, N, n, t, F, F′, P, aᵣ, bᵣ, γᵣ, ωₙ, ω′ₙ, aₙ, bₙ, γₙ, a′ₙ, b′ₙ, γ′ₙ, δ, ϵ, γ, δ′, ϵ′, γ′)
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
    advance!(k, F′ₙ)\n

Method `advance!` moves a solution from previous step `n-1` to current step `n` along a solution path comprised of N solution nodes overall. This is accomplished by integrating the supplied rate of deformation gradient `F′ₙ` using a backward difference formula (BDF) from which all other terms are then derived. This method updates counter `k.n` plus those entries to its history arrays located at the nᵗʰ array entries in the `k` data structure. Specifically: deformation gradient `k.F[n]` and its rate `k.F′[n]`, pivot case `k.P[n]`, Laplace attributes `k.aₙ[n]`, `k.bₙ[n]`, `k.γₙ[n]` and `k.ωₙ[n]` and their rates `k.a′ₙ[n]`, `k.b′ₙ[n]`, `k.γ′ₙ[n]` and `k.ω′ₙ[n]`, plus Laplace strains `k.δ[n]`, `k.ϵ[n]` and `k.γ[n]` and their rates `k.δ′[n]`, `k.ϵ′[n]` and `k.γ′[n]`, all evaluated in the user's co-ordinate system whose base vectors are denoted as (𝕚, 𝕛).
"""
function advance!(k::MembraneKinematics, F′ₙ::PhysicalTensor)
    # Advance the counter.
    if k.n < k.N+1
        set!(k.n, get(k.n)+1)
    else
        msg = "The data structure is full and cannot accept further data."
        throw(ErrorException(msg))
    end
    n = get(k.n)

    # Convert the passed variable to CGS units.
    𝐹′ₙ = toCGS(F′ₙ)

    # Verify inputs.
    if (𝐹′ₙ.matrix.rows ≠ 2) || (𝐹′ₙ.matrix.cols ≠ 2)
        msg = "Supplied deformation gradient rate F′ₙ must be a 2x2 matrix."
        throw(DimensionMismatch(msg))
    end
    if 𝐹′ₙ.units ≠ TIME_RATE
        msg = "Supplied deformation gradient rate F′ₙ must have units of reciprocal time."
        throw(ErrorException(msg))
    end

    # Advance the fields, i.e., insert new values into the data arrays.

    # Integrate F′ₙ using a backward difference formula (BDF).
    if n == 2
        F₁ = k.F[1]
    elseif n == 3
        F₁ = k.F[1]
        F₂ = k.F[2]
    else
        Fₙ₋₁ = k.F[n-1]
        Fₙ₋₂ = k.F[n-2]
        Fₙ₋₃ = k.F[n-3]
    end
    Fₙ = PhysicalTensor(2, 2, DIMENSIONLESS)
    for i in 1:2
        for j in 1:2
            if n == 2
                Fₙ[i,j] = F₁[i,j] + 𝐹′ₙ[i,j]*k.dt
            elseif n == 3
                Fₙ[i,j] = (4/3)*F₂[i,j] - (1/3)*F₁[i,j] + (2/3)*𝐹′ₙ[i,j]*k.dt
            else
                Fₙ[i,j] = ((18/11)*Fₙ₋₁[i,j] - (9/11)*Fₙ₋₂[i,j]
                        + (2/11)*Fₙ₋₃[i,j] + (6/11)*𝐹′ₙ[i,j]*k.dt)
            end
        end
    end
    k.F[n]  = Fₙ
    k.F′[n] = 𝐹′ₙ

    # Get attributes for deformation F and rate of deformation 𝐹′ gradients.
    x  = Fₙ[1,1]
    y  = Fₙ[2,2]
    x′ = 𝐹′ₙ[1,1]
    y′ = 𝐹′ₙ[2,2]

    # Establish the Gram and Laplace attributes, and their rates.
    if sign(Fₙ[1,2]) == sign(Fₙ[2,1])
        # The deformation has an attribute that is a pure shear.
        # g is this pure shear.
        if ((abs(Fₙ[1,2])/Fₙ[2,2] > abs(Fₙ[2,1])/Fₙ[1,1])
            || (abs(Fₙ[1,2])/Fₙ[2,2] ≈ abs(Fₙ[2,1])/Fₙ[1,1]))
            # Co-ordinates are indexed as supplied--the default condition.
            case = 1
            # Pure shear contributions.
            g  = Fₙ[2,1] / Fₙ[1,1]
            g′ = (Fₙ[1,1]*𝐹′ₙ[2,1] - Fₙ[2,1]*𝐹′ₙ[1,1]) / (Fₙ[1,1]*Fₙ[1,1])
            # Simple shear contributions.
            G  = (Fₙ[1,1]*Fₙ[1,2] - Fₙ[2,2]*Fₙ[2,1]) / (Fₙ[1,1]*Fₙ[2,2])
            G′ = (Fₙ[2,2]*𝐹′ₙ[1,2] - Fₙ[1,2]*𝐹′ₙ[2,2]) / (Fₙ[2,2]*Fₙ[2,2]) - g′
        else
            # The Gram co-ordinate system is left handed.
            case = 2
            # Pure shear contributions.
            g  = Fₙ[1,2] / Fₙ[2,2]
            g′ = (Fₙ[2,2]*𝐹′ₙ[1,2] - Fₙ[1,2]*𝐹′ₙ[2,2]) / (Fₙ[2,2]*Fₙ[2,2])
            # Simple shear contributions.
            G  = -(Fₙ[1,1]*Fₙ[1,2] - Fₙ[2,2]*Fₙ[2,1]) / (Fₙ[1,1]*Fₙ[2,2])
            G′ = (Fₙ[1,1]*𝐹′ₙ[2,1] - Fₙ[2,1]*𝐹′ₙ[1,1]) / (Fₙ[1,1]*Fₙ[1,1]) - g′
        end
        # Laplace stretch attributes in the co-ordinate frame 𝕚 × 𝕛.
        aₙ = x * sqrt(1+g*g)
        bₙ = y * (1 - g*(g+G)) / sqrt(1+g*g)
        γₙ = (y/x) * (2g+G) / (1+g*g)
        ωₙ = PhysicalScalar(atan(g), DIMENSIONLESS)

        # Rates of Laplace attributes in the co-ordinate frame 𝕚 × 𝕛.
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
            g′ = -(Fₙ[1,1]*𝐹′ₙ[2,1] - Fₙ[2,1]*𝐹′ₙ[1,1]) / (Fₙ[1,1]*Fₙ[1,1])
            # Angle of rigid-body rotation in the co-ordinate frame 𝕚 × 𝕛.
            ωₙ  = PhysicalScalar(-atan(g), DIMENSIONLESS)
            ω′ₙ = -g′/(1+g*g)
        else
            # The Gram co-ordinate system is left handed.
            case = 4
            # Rigid-body rotation contributions.
            g  = -Fₙ[1,2] / Fₙ[2,2]
            g′ = -(Fₙ[2,2]*𝐹′ₙ[1,2] - Fₙ[1,2]*𝐹′ₙ[2,2]) / (Fₙ[2,2]*Fₙ[2,2])
            # Angle of rigid-body rotation in the co-ordinate frame 𝕚 × 𝕛.
            ωₙ  = PhysicalScalar(atan(g), DIMENSIONLESS)
            ω′ₙ = g′/(1+g*g)
        end
        # G is a simple shear.
        G  = (Fₙ[1,1]*Fₙ[1,2] + Fₙ[2,2]*Fₙ[2,1]) / (Fₙ[1,1]*Fₙ[2,2])
        G′ = ((Fₙ[2,2]*𝐹′ₙ[1,2] - Fₙ[1,2]*𝐹′ₙ[2,2]) / (Fₙ[2,2]*Fₙ[2,2]) +
            (Fₙ[1,1]*𝐹′ₙ[2,1] - Fₙ[2,1]*𝐹′ₙ[1,1]) / (Fₙ[1,1]*Fₙ[1,1]))

        # Laplace attributes in the co-ordinate frame 𝕚 × 𝕛.
        aₙ = x * sqrt(1+g*g)
        bₙ = y * (1 + g*(g+G)) / sqrt(1+g*g)
        γₙ = (y/x) * G / (1+g*g)

        # Rates of Laplace attributes in the co-ordinate frame 𝕚 × 𝕛.
        a′ₙ = aₙ*(x′/x + g*g′/(1+g*g))
        b′ₙ = bₙ*(y′/y - g*g′/(1+g*g)) + y*((2g+G)*g′ + g*G′)/sqrt(1+g*g)
        γ′ₙ = γₙ*(y′/y - x′/x - 2g*g′/(1+g*g)) + (y/x)*G′/(1+g*g)
    end

    # Advance the data array that holds pivoting cases.
    k.P[n] = case

    # Advance the Laplace stretch attributes and their rates for κᵣ ↦ κₙ.
    k.aₙ[n]  = aₙ/k.aᵣ
    k.bₙ[n]  = bₙ/k.bᵣ
    k.γₙ[n]  = (k.aᵣ/k.bᵣ)*(γₙ - k.γᵣ)
    k.ωₙ[n]  = ωₙ
    k.a′ₙ[n] = a′ₙ/k.aᵣ
    k.b′ₙ[n] = b′ₙ/k.bᵣ
    k.γ′ₙ[n] = (k.aᵣ/k.bᵣ)*γ′ₙ
    k.ω′ₙ[n] = ω′ₙ

    # Advance the thermodynamic strains and their rates for κᵣ ↦ κₙ.
    k.δ[n]  = PhysicalScalar(0.5log(k.aₙ[n]*k.bₙ[n]), DIMENSIONLESS)
    k.ϵ[n]  = PhysicalScalar(0.5log(k.aₙ[n]/k.bₙ[n]), DIMENSIONLESS)
    k.γ[n]  = γₙ - k.γᵣ
    k.δ′[n] = 0.5(a′ₙ/aₙ + b′ₙ/bₙ)
    k.ϵ′[n] = 0.5(a′ₙ/aₙ - b′ₙ/bₙ)
    k.γ′[n] = γ′ₙ

    return nothing
end # advance!

"""
Method:\n
    update!(k, F′ₙ)\n
Method `update!` provides a capability to refine a solution at step `n` by re-integrating the deformation gradient rate `F′ₙ`, thereby allowing for iterative improvements to be made on the deformation rate `F′ₙ` from an external algorithm, e.g., a finite element engine. There is no need to call `update!` unless `F′ₙ` is being iteratively refined at step `n` by some external optimization process.
"""
function update!(k::MembraneKinematics, F′ₙ::PhysicalTensor)
    if k.n > 1
        set!(k.n, get(k.n)-1)
        advance!(k, F′ₙ)
    end
    return nothing
end # update!
