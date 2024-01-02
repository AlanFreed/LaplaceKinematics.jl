"""
Matrices P3D₁, …, P3D₆ are the six possible permutation matrices in 3 space.
"""
one = newPhysicalScalar(1.0, DIMENSIONLESS)
P3D₁ = newPhysicalTensor(3, 3, DIMENSIONLESS)  # (𝕚, 𝕛, 𝕜) ↦ (𝔼₁, 𝔼₂, 𝔼₃)
P3D₁[1,1] = one
P3D₁[2,2] = one
P3D₁[3,3] = one
P3D₂ = newPhysicalTensor(3, 3, DIMENSIONLESS)  # (𝕚, 𝕜, 𝕛) ↦ (𝔼₁, 𝔼₂, 𝔼₃)
P3D₂[1,1] = one
P3D₂[2,3] = one
P3D₂[3,2] = one
P3D₃ = newPhysicalTensor(3, 3, DIMENSIONLESS)  # (𝕛, 𝕜, 𝕚) ↦ (𝔼₁, 𝔼₂, 𝔼₃)
P3D₃[1,3] = one
P3D₃[2,1] = one
P3D₃[3,2] = one
P3D₄ = newPhysicalTensor(3, 3, DIMENSIONLESS)  # (𝕛, 𝕚, 𝕜) ↦ (𝔼₁, 𝔼₂, 𝔼₃)
P3D₄[1,2] = one
P3D₄[2,1] = one
P3D₄[3,3] = one
P3D₅ = newPhysicalTensor(3, 3, DIMENSIONLESS)  # (𝕜, 𝕚, 𝕛) ↦ (𝔼₁, 𝔼₂, 𝔼₃)
P3D₅[1,2] = one
P3D₅[2,3] = one
P3D₅[3,1] = one
P3D₆ = newPhysicalTensor(3, 3, DIMENSIONLESS)  # (𝕜, 𝕛, 𝕚) ↦ (𝔼₁, 𝔼₂, 𝔼₃)
P3D₆[1,3] = one
P3D₆[2,2] = one
P3D₆[3,1] = one
# Cases 1, 3, 5 have a determinant of +1.
# Cases 2, 4, 6 have a determinant of -1.
#=
--------------------------------------------------------------------------------
=#
"""
Function:\n
    (case, 𝐹) = pivotF(F)\n
Given a deformation gradient 'F' evaluated in some user selected co-ordinate system, either (𝕚, 𝕛) for 2D analyses or (𝕚, 𝕛, 𝕜) for 3D analyses, this function returns a tuple comprised of: (i) the 'case' of permutation, i.e., either 1 or 2 for 2D analyses or 1, 2, … or 6 for 3D analyses, and (ii) a deformation gradient '𝐹' of appropriate dimension that is now evaluated in a permuted co-ordinate system of either (𝔼₁, 𝔼₂) or (𝔼₁, 𝔼₂, 𝔼₃), respectively, to ensure frame indifference. This co-ordinate frame is suitable for constitutive construction.
"""
function pivotF(F::PhysicalTensor)::Tuple
    # Verify input.
    if ((isCGS(F) && (F.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(F) && (F.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied deformation gradient F is not dimensionless."
        throw(ErrorException(msg))
    end
    # Establish the appropriate co-ordinate permutation for frame indifference.
    if (F.r == 2) && (F.c == 2)
        # Determine the appropriate permutation to apply to (𝕚, 𝕛).
        Σ₁ = abs(F[1,2]) / F[2,2]
        Σ₂ = abs(F[2,1]) / F[1,1]
        # Establish the permutation case.
        if Σ₁ > Σ₂ || Σ₁ ≈ Σ₂
            case = 1
            P = P2D₁
        else
            case = 2
            P = P2D₂
        end
    elseif (F.r == 3) && (F.c == 3)
        # Determine the appropriate permutation to apply to (𝕚, 𝕛, 𝕜).
        f₁ = newPhysicalVector(3, DIMENSIONLESS)
        f₂ = newPhysicalVector(3, DIMENSIONLESS)
        f₃ = newPhysicalVector(3, DIMENSIONLESS)
        if isEulerian
            # Compute the three row vectors needed to establish 𝔼₃.
            for i in 1:3
                f₁[i] = F[1,i]
                f₂[i] = F[2,i]
                f₃[i] = F[3,i]
            end
        else # Lagrangian
            # Compute the three column vectors needed to establish 𝔼₃.
            for i in 1:3
                f₁[i] = F[i,1]
                f₂[i] = F[i,2]
                f₃[i] = F[i,3]
            end
        end
        # Establish the permutation case.
        Σ₁ = norm(f₁)
        Σ₂ = norm(f₂)
        Σ₃ = norm(f₃)
        if ((Σ₁ < Σ₂) || (Σ₁ ≈ Σ₂)) && ((Σ₁ < Σ₃) || (Σ₁ ≈ Σ₃))
            if (Σ₃ > Σ₂) || (Σ₃ ≈ Σ₂)
                case = 1
                P = P3D₁
            else
                case = 2
                P = P3D₂
            end
        elseif ((Σ₂ < Σ₁) || (Σ₂ ≈ Σ₁)) && ((Σ₂ < Σ₃) || (Σ₂ ≈ Σ₃))
            if (Σ₃ > Σ₁) || (Σ₃ ≈ Σ₁)
                case = 4
                P = P3D₄
            else
                case = 3
                P = P3D₃
            end
        else
            if (Σ₂ > Σ₁) || (Σ₂ ≈ Σ₁)
                case = 5
                P = P3D₅
            else
                case = 6
                P = P3D₆
            end
        end
    else
        # Verify input.
        msg = "The deformation gradient F must be either a 2x2 or 3x3 matrix."
        throw(ErrorException(msg))
    end
    # Permute the deformation gradient.
    𝐹 = P * F * transpose(P)
    return (case, 𝐹)
end

"""
Function:\n
    𝑣 = permuteIn(v, case)\n
This function maps the components of a vector 'v' evaluated in the user's co-ordinate system of either (𝕚, 𝕛) or (𝕚, 𝕛, 𝕜) into its associated components '𝑣' evaluated in a permuted co-ordinate system of either (𝔼₁, 𝔼₂) or (𝔼₁, 𝔼₂, 𝔼₃) for 2D or 3D vectors, respectively. There are two possible permutations that can arise in 2-space, and six possible permutations that can arise in 3-space. Which permutation applies is specified by argument 'case', viz., 'case' ∈ {1, 2} for 2D vectors, or 'case' ∈ {1, 2, …, 6} for 3D vectors.
"""
function permuteIn(v::PhysicalVector, case::Integer)::PhysicalVector
    # Map the vector.
    if v.l == 2
        if case == 1
            P = P2D₁
        elseif case == 2
            P = P2D₂
        else
            msg = "2D permutation cases must come from the set {1, 2}."
            throw(ErrorException(msg))
        end
    elseif v.l == 3
        if case == 1
            P = P3D₁
        elseif case == 2
            P = P3D₂
        elseif case == 3
            P = P3D₃
        elseif case == 4
            P = P3D₄
        elseif case == 5
            P = P3D₅
        elseif case == 6
            P = P3D₆
        else
            msg = "3D permutation cases must come from the set {1, 2, …, 6}."
            throw(ErrorException(msg))
        end
    else
        msg = "The supplied vector has an inadmissible length."
        throw(ErrorException(msg))
    end
    # Permute the vector.
    𝑣 = P * v
    return 𝑣
end

"""
Function:\n
    v = permuteOut(𝑣, case)\n
This function maps the components of a vector '𝑣' evaluated in a permuted co-ordinate system of either (𝔼₁, 𝔼₂) or (𝔼₁, 𝔼₂, 𝔼₃) out-to its components 'v' evaluated in the user's co-ordinate system of either (𝕚, 𝕛) or (𝕚, 𝕛, 𝕜) for 2D or 3D vectors, respectively. There are two possible permutations that can arise in 2-space, and six possible permutations that can arise in 3-space. Which permutation applies is specified by argument 'case', viz., 'case' ∈ {1, 2} for 2D vectors, or 'case' ∈ {1, 2, …, 6} for 3D vectors.
"""
function permuteOut(𝑣::PhysicalVector, case::Integer)::PhysicalVector
    # Map the vector.
    if 𝑣.l == 2
        if case == 1
            P = P2D₁
        elseif case == 2
            P = P2D₂
        else
            msg = "2D permutation cases must come from the set {1, 2}."
            throw(ErrorException(msg))
        end
    elseif 𝑣.l == 3
        if case == 1
            P = P3D₁
        elseif case == 2
            P = P3D₂
        elseif case == 3
            P = P3D₃
        elseif case == 4
            P = P3D₄
        elseif case == 5
            P = P3D₅
        elseif case == 6
            P = P3D₆
        else
            msg = "3D permutation cases must come from the set {1, 2, …, 6}."
            throw(ErrorException(msg))
        end
    else
        msg = "The supplied vector has an inadmissible length."
        throw(ErrorException(msg))
    end
    # Un-permute the vector.
    v = transpose(P) * 𝑣
    return v
end

"""
Function:\n
    𝑚 = permuteIn(m, case)\n
This function maps the components of a matrix 'm' evaluated in the user's co-ordinate system of either (𝕚, 𝕛) or (𝕚, 𝕛, 𝕜) into its associated components '𝑚' evaluated in a permuted co-ordinate system of either (𝔼₁, 𝔼₂) or (𝔼₁, 𝔼₂, 𝔼₃) for 2D or 3D matrices, respectively. There are two possible permutations that can arise in 2-space, and six possible permutations that can arise in 3-space. Which permutation applies is specified by argument 'case', viz., 'case' ∈ {1, 2} for 2D matrices, or 'case' ∈ {1, 2, …, 6} for 3D matrices.
"""
function permuteIn(m::PhysicalTensor, case::Integer)::PhysicalTensor
    # Map the matrix.
    if (m.r == 2) && (m.c == 2)
        if case == 1
            P = P2D₁
        elseif case == 2
            P = P2D₂
        else
            msg = "2D permutation cases must come from the set {1, 2}."
            throw(ErrorException(msg))
        end
    elseif (m.r == 3) && (m.c == 3)
        if case == 1
            P = P3D₁
        elseif case == 2
            P = P3D₂
        elseif case == 3
            P = P3D₃
        elseif case == 4
            P = P3D₄
        elseif case == 5
            P = P3D₅
        elseif case == 6
            P = P3D₆
        else
            msg = "3D permutation cases must come from the set {1, 2, …, 6}."
            throw(ErrorException(msg))
        end
    else
        msg = "The supplied matrix has an inadmissible dimension."
        throw(ErrorException(msg))
    end
    # Permute the matrix.
    𝑚 = P * m * transpose(P)
    return 𝑚
end

"""
Function:\n
    m = permuteOut(𝑚, case)\n
This function maps the components of a matrix '𝑚' evaluated in a permuted co-ordinate system of either (𝔼₁, 𝔼₂) or (𝔼₁, 𝔼₂, 𝔼₃) out-to its components 'm' evaluated in the user's co-ordinate system of either (𝕚, 𝕛) or (𝕚, 𝕛, 𝕜) for 2D or 3D matrices, respectively. There are two possible permutations that can arise in 2-space, and six possible permutations that can arise in 3-space. Which permutation applies is specified by argument 'case', viz., 'case' ∈ {1, 2} for 2D matrices, or 'case' ∈ {1, 2, …, 6} for 3D matrices.
"""
function permuteOut(𝑚::PhysicalTensor, case::Integer)::PhysicalTensor
    # Map the matrix.
    if (𝑚.r == 2) && (𝑚.c == 2)
        if case == 1
            P = P2D₁
        elseif case == 2
            P = P2D₂
        else
            msg = "2D permutation cases must come from the set {1, 2}."
            throw(ErrorException(msg))
        end
    elseif (𝑚.r == 3) && (𝑚.c == 3)
        if case == 1
            P = P3D₁
        elseif case == 2
            P = P3D₂
        elseif case == 3
            P = P3D₃
        elseif case == 4
            P = P3D₄
        elseif case == 5
            P = P3D₅
        elseif case == 6
            P = P3D₆
        else
            msg = "3D permutation cases must come from the set {1, 2, …, 6}."
            throw(ErrorException(msg))
        end
    else
        msg = "The supplied matrix has an inadmissible dimension."
        throw(ErrorException(msg))
    end
    # Un-permute the matrix.
    m = transpose(P) * 𝑚 * P
    return m
end

"""
Function:\n
    (q₀, q₁, q₂, q₃) = quaternion(G)\n
Uses Spurrier's algorithm to compute the Euler-Rodrigues quaternion, i.e., the set {q₀, q₁, q₂, q₃}, for computing the axis and angle of rotation from a given orthogonal, 3D, Gram, rotation matrix 'G'. (See function spin below.)
"""
function quaternion(G::PhysicalTensor)::Tuple
    # Verify input.
    if (G.r ≠ 3) || (G.c ≠ 3)
        msg = "The orthogonal matrix G must be a 3x3 matrix."
        throw(ErrorException(msg))
    end
    if ((isCGS(G) && (G.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(G) && (G.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied orthogonal matrix G is not dimensionless."
        throw(ErrorException(msg))
    end
    # Determine the Euler-Rodrigues quaternion.
    G₀ = tr(G)
    G₁ = G[1,1]
    G₂ = G[2,2]
    G₃ = G[3,3]
    if (G₀ ≥ G₁) && (G₀ ≥ G₂) && (G₀ ≥ G₃)
        q₀ = √(1 + G₀) / 2
        q₁ = (G[3,2] - G[2,3]) / (4q₀)
        q₂ = (G[1,3] - G[3,1]) / (4q₀)
        q₃ = (G[2,1] - G[1,2]) / (4q₀)
    elseif (G₁ ≥ G₀) && (G₁ ≥ G₂) && (G₁ ≥ G₃)
        q₁ = √(2G₁ + (1 - G₀)) / 2
        q₀ = (G[3,2] - G[2,3]) / (4q₁)
        q₂ = (G[2,1] + G[1,2]) / (4q₁)
        q₃ = (G[3,1] + G[1,3]) / (4q₁)
    elseif (G₂ ≥ G₀) && (G₂ ≥ G₁) && (G₂ ≥ G₃)
        q₂ = √(2G₂ + (1 - G₀)) / 2
        q₀ = (G[1,3] - G[3,1]) / (4q₂)
        q₁ = (G[1,2] + G[2,1]) / (4q₂)
        q₃ = (G[3,2] + G[2,3]) / (4q₂)
    else
        q₃ = √(2G₃ + (1 - G₀)) / 2
        q₀ = (G[2,1] - G[1,2]) / (4q₃)
        q₁ = (G[1,3] + G[3,1]) / (4q₃)
        q₂ = (G[2,3] + G[3,2]) / (4q₃)
    end
    return (q₀, q₁, q₂, q₃)
end

"""
function:\n
    Ωₙ = spin(Gₙ, Gₙ₋₁, Gₙ₋₂, dt)\n
Given the current 'Gₙ', previous 'Gₙ₋₁', and previous-previous 'Gₙ₋₂' 3D Gram rotation matrices, which are orthogonal, and are obtained from a deformation gradient history whose entries are separated in time by an interval 'dt', this function returns the spin of this Gram rotation, viz., Ωₙ = (dGₙ/dt)⋅Gₙᵀ.
"""
function spin(Gₙ::PhysicalTensor, Gₙ₋₁::PhysicalTensor, Gₙ₋₂::PhysicalTensor, dt::PhysicalScalar)::PhysicalTensor
    # Determine angles θᵢ and axes 𝕘ᵢ of Gram rotations Gᵢ.
    # At current step n:
    (q₀ₙ, q₁ₙ, q₂ₙ, q₃ₙ) = quaternion(Gₙ)
    θₙ = newPhysicalScalar(2acos(q₀ₙ), DIMENSIONLESS)
    if θₙ ≈ 1
        axisExistsₙ = false
    else
        axisExistsₙ = true
        𝕘ₙ = newPhysicalVector(3, DIMENSIONLESS)
        𝕘ₙ[1] = q₁ₙ / sin(θₙ/2)
        𝕘ₙ[2] = q₂ₙ / sin(θₙ/2)
        𝕘ₙ[3] = q₃ₙ / sin(θₙ/2)
    end
    # At previous step n-1:
    (q₀ₙ₋₁, q₁ₙ₋₁, q₂ₙ₋₁, q₃ₙ₋₁) = quaternion(Gₙ₋₁)
    θₙ₋₁ = newPhysicalScalar(2acos(q₀ₙ₋₁), DIMENSIONLESS)
    if θₙ₋₁ ≈ 1
        axisExistsₙ₋₁ = false
    else
        axisExistsₙ₋₁ = true
        𝕘ₙ₋₁ = newPhysicalVector(3, DIMENSIONLESS)
        𝕘ₙ₋₁[1] = q₁ₙ₋₁ / sin(θₙ₋₁/2)
        𝕘ₙ₋₁[2] = q₂ₙ₋₁ / sin(θₙ₋₁/2)
        𝕘ₙ₋₁[3] = q₃ₙ₋₁ / sin(θₙ₋₁/2)
    end
    # At previous-previous step n-2:
    (q₀ₙ₋₂, q₁ₙ₋₂, q₂ₙ₋₂, q₃ₙ₋₂) = quaternion(Gₙ₋₂)
    θₙ₋₂ = newPhysicalScalar(2acos(q₀ₙ₋₂), DIMENSIONLESS)
    if θₙ₋₂ ≈ 1
        axisExistsₙ₋₂ = false
    else
        axisExistsₙ₋₂ = true
        𝕘ₙ₋₂ = newPhysicalVector(3, DIMENSIONLESS)
        𝕘ₙ₋₂[1] = q₁ₙ₋₂ / sin(θₙ₋₂/2)
        𝕘ₙ₋₂[2] = q₂ₙ₋₂ / sin(θₙ₋₂/2)
        𝕘ₙ₋₂[3] = q₃ₙ₋₂ / sin(θₙ₋₂/2)
    end
    # Determine the axis of spin ω.
    if axisExistsₙ
        if axisExistsₙ₋₁
            if axisExistsₙ₋₂
                # Use the second-order backward difference formula.
                dθₙ = (3θₙ - 4θₙ₋₁ + θₙ₋₂) / (2dt)
                d𝕘ₙ = (3𝕘ₙ - 4𝕘ₙ₋₁ + 𝕘ₙ₋₂) / (2dt)
            else
                # Use the first-order backward difference formula.
                dθₙ = (θₙ - θₙ₋₁) / dt
                d𝕘ₙ = (𝕘ₙ - 𝕘ₙ₋₁) / dt
            end
            ω = dθₙ * 𝕘ₙ + sin(θₙ) * d𝕘ₙ + (1 - cos(θₙ)) * cross(𝕘ₙ, d𝕘ₙ)
        else
            # There is an axis of rotation, but no axis of spin.
            ω = newPhysicalVector(3, RATE)
        end
    else
        # There are no axes of rotation or spin.
        ω = newPhysicalVector(3, RATE)
    end
    # Create the 3D skew-symmetric spin tensor.
    Ω = newPhysicalTensor(3, 3, RATE)
    Ω[1,2] = -ω[3]
    Ω[1,3] =  ω[2]
    Ω[2,1] =  ω[3]
    Ω[2,3] = -ω[1]
    Ω[3,1] = -ω[2]
    Ω[3,2] =  ω[1]
    return Ω
end
#=
-------------------------------------------------------------------------------
=#
"""
Type:\n
    Kinematics3D\n
        isEulerian::    Eulerian → true; Lagrangian → false\n
        dt::    size of the time step separating neighboring entries\n
        N::     total number of steps or grid points to be considered\n
        n::     steps since the initial state/condition whereat n = 1\n
        m::     steps since last loss in differentiability in deformation\n
        F::     array of deformation gradients evaluated in (𝔼₁, 𝔼₂, 𝔼₃)\n
        P::     array of permutation cases, viz., i in Pᵢ, i ∈ {1, 2, ..., 6}\n
        L::     array of Laplace stretches evaluated in (𝔼₁, 𝔼₂, 𝔼₃)\n
        L⁻¹::   array of inverse Laplace stretches in (𝔼₁, 𝔼₂, 𝔼₃)\n
        G::     array of Gram rotations out of co-ordinates (𝔼₁, 𝔼₂, 𝔼₃)\n
        Ω::     array of Gram spins (dG/dt)Gᵀ evaluated in (𝔼₁, 𝔼₂, 𝔼₃)\n
        Lᵣ::    reference Laplace stretch evaluated in (𝔼₁, 𝔼₂, 𝔼₃)\n
        Lᵣ⁻¹::  inverse of reference Laplace stretch in (𝔼₁, 𝔼₂, 𝔼₃)\n
        aᵣ::    reference in-plane elongation for a\n
        bᵣ::    reference in-plane elongation for b\n
        cᵣ::    reference out-of-plane elongation for c\n
        αᵣ::    reference out-of-plane shear α\n
        βᵣ::    reference out-of-plane shear β\n
        γᵣ::    reference in-plane shear γ\n
        a::     array of in-plane elongations for a\n
        b::     array of in-plane elongations for b\n
        c::     array of out-of-plane elongations for c\n
        α::     array of out-of-plane shears α\n
        β::     array of out-of-plane shears β\n
        γ::     array of in-plane shears γ\n
        da::    array of in-plane elongation rates da/dt\n
        db::    array of in-plane elongation rates db/dt\n
        dc::    array of out-of-plane elongation rates dc/dt\n
        dα::    array of out-of-plane shear rates dα/dt\n
        dβ::    array of out-of-plane shear rates dβ/dt\n
        dγ::    array of in-plane shear rates dγ/dt\n
        δ::     array of dilatation strains δ\n
        ϵ₁::    array of squeeze strains ϵ₁\n
        ϵ₂::    array of squeeze strains ϵ₂\n
        ϵ₃::    array of squeeze strains ϵ₃\n
        γ₁::    array of shear strains γ₁\n
        γ₂::    array of shear strains γ₂\n
        γ₃::    array of shear strains γ₃\n
        dδ::    array of dilatation strain rates dδ/dt\n
        dϵ₁::   array of squeeze strain rates dϵ₁/dt\n
        dϵ₂::   array of squeeze strain rates dϵ₂/dt\n
        dϵ₃::   array of squeeze strain rates dϵ₃/dt\n
        dγ₁::   array of shear strain rates dγ₁/dt\n
        dγ₂::   array of shear strain rates dγ₂/dt\n
        dγ₃::   array of shear strain rates dγ₃/dt\n
Kinematics3D is a data structure that contains the various fields associated with a Laplace (triangular) measure for stretch in three space. The arrays that comprise this data structure allow for a history of kinematic variables to be created for later retrieval and use, e.g., in a post analysis or for graphing. All fields are evaluated in a frame-indifferent co-ordinate system (𝔼₁, 𝔼₂, 𝔼₃).
"""
struct Kinematics3D
    isEulerian::Bool    # true for Eulerian; false for Lagrangian
    dt::    PhysicalScalar   # differential in time between grid points
    N::     UInt32      # total number of steps or grid points considered
    n::     MInteger    # steps taken since the initial state whereat n = 1
    m::     MInteger    # steps taken since the last loss in differentiability
    # kinematic tensor fields located at the nodal points of the grid
    F::     ArrMatrix   # pivoted or frame-indifferent deformation gradients
    P::     Array       # permutation case used for co-ordinate pivoting
    L::     ArrMatrix   # lower or upper-triangular Laplace stretch
    L⁻¹::   ArrMatrix   # inverse of triangular Laplace stretch
    G::     ArrMatrix   # orthogonal rotations of Gram co-ordinates
    Ω::     ArrMatrix   # skew-symmetric spins of Gram co-ordinates
    # Laplace stretch attributes for deformation κ₀ ↦ κᵣ
    Lᵣ::    PhysicalTensor   # reference triangular Laplace stretch
    Lᵣ⁻¹::  PhysicalTensor   # inverse of reference triangular Laplace stretch
    aᵣ::    PhysicalScalar   # reference in-plane elongation for a
    bᵣ::    PhysicalScalar   # reference in-plane elongation for b
    cᵣ::    PhysicalScalar   # reference out-of-plane elongation for c
    αᵣ::    PhysicalScalar   # reference out-of-plane shear for α
    βᵣ::    PhysicalScalar   # reference out-of-plane shear for β
    γᵣ::    PhysicalScalar   # reference in-plane shear for γ
    # Laplace stretch attributes for deformation κ₀ ↦ κₙ
    a::     ArrScalar   # current in-plane elongation for a
    b::     ArrScalar   # current in-plane elongation for b
    c::     ArrScalar   # current out-of-plane elongation for c
    α::     ArrScalar   # current out-of-plane shear for α
    β::     ArrScalar   # current out-of-plane shear for β
    γ::     ArrScalar   # current in-plane shear for γ
    # Rates of Laplace stretch attributes at configuration κₙ
    da::    ArrScalar   # current in-plane elongation rate for da/dt
    db::    ArrScalar   # current in-plane elongation rate for db/dt
    dc::    ArrScalar   # current out-of-plane elongation rate for dc/dt
    dα::    ArrScalar   # current out-of-plane shear rate for dα/dt
    dβ::    ArrScalar   # current out-of-plane shear rate for dβ/dt
    dγ::    ArrScalar   # current in-plane shear rate for dγ/dt
    # Laplace strain attributes for deformation κᵣ ↦ κₙ
    δ::     ArrScalar   # dilatation strain δ
    ϵ₁::    ArrScalar   # first squeeze strain ϵ₁
    ϵ₂::    ArrScalar   # second squeeze strain ϵ₂
    ϵ₃::    ArrScalar   # third squeeze strain ϵ₃
    γ₁::    ArrScalar   # first out-of-plane shear strain γ₁
    γ₂::    ArrScalar   # second out-of-plane shear strain γ₂
    γ₃::    ArrScalar   # third in-plane shear strain γ₃
    # Laplace strain-rate attributes
    dδ::    ArrScalar   # dilatation strain rate dδ/dt
    dϵ₁::   ArrScalar   # first squeeze strain rate dϵ₁/dt
    dϵ₂::   ArrScalar   # second squeeze strain rate dϵ₂/dt
    dϵ₃::   ArrScalar   # third squeeze strain rate dϵ₃/dt
    dγ₁::   ArrScalar   # first out-of-plane shear strain rate dγ₁/dt
    dγ₂::   ArrScalar   # second out-of-plane shear strain rate dγ₂/dt
    dγ₃::   ArrScalar   # third in-plane shear strain rate dγ₃/dt
end

"""
Constructor:\n 
    k = newKinematics3D(N, dt, aᵣ, bᵣ, cᵣ, αᵣ, βᵣ, γᵣ, F₀, isEulerian)\n
Returns a new data structure 'k' of type 'Kinematics3D' that holds a variety of kinematic fields. An Eulerian Laplace stretch is determined whenever 'isEulerian' is true; otherwise, a Lagrangian Laplace stretch is determined. Arguments include: (i) the number of grid points or nodes 'N' where solutions are to be computed, (ii) an uniform differential step in time 'dt' separating neighboring nodes in the data arrays to be populated, (iii) the reference Laplace stretch attributes, viz., 'aᵣ', 'bᵣ', 'cᵣ', 'αᵣ', 'βᵣ' and 'γᵣ', against which strains are to be defined, (iv) an initial deformation gradient 'F₀' evaluated in the user's co-ordinate system with base vectors (𝕚, 𝕛, 𝕜), and (v) a flag 'isEulerian' that is used to select an Eulerian vs. Lagrangian construction for the supplied fields. The fields in this data structure are evaluated in a frame-indifferent co-ordinate system (𝔼₁, 𝔼₂, 𝔼₃), which is a permutation of (𝕚, 𝕛, 𝕜).
"""
function newKinematics3D(N::Integer, dt::PhysicalScalar, aᵣ::PhysicalScalar, bᵣ::PhysicalScalar, cᵣ::PhysicalScalar, αᵣ::PhysicalScalar, βᵣ::PhysicalScalar, γᵣ::PhysicalScalar, F₀::PhysicalTensor, isEulerian::Bool=true)::Kinematics3D
    # Verify inputs.
    if (N < 1) || (N > typemax(UInt32))
        msg = string("Arrays can only have lengths ∈ [1…4,294,967,295].")
        throw(ErrorException(msg))
    end
    if ((isCGS(dt) && (dt.u ≠ PhysicalFields.CGS_SECOND)) ||
        (isSI(dt) && (dt.u ≠ PhysicalFields.SI_SECOND)))
        msg = "The supplied time increment dt does not have units of time."
        throw(ErrorException(msg))
    end
    if dt < newPhysicalScalar(Float64(eps(Float32)), TIME)
        msg = "The supplied time increment dt must be positive valued."
        throw(ErrorException(msg))
    end
    if ((isCGS(aᵣ) && (aᵣ.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(aᵣ) && (aᵣ.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied reference stretch aᵣ is not dimensionless."
        throw(ErrorException(msg))
    end
    if aᵣ < Float64(eps(Float16))
        msg = "The supplied reference strecth aᵣ must be positive valued."
        throw(ErrorException(msg))
    end
    if ((isCGS(bᵣ) && (bᵣ.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(bᵣ) && (bᵣ.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied reference stretch bᵣ is not dimensionless."
        throw(ErrorException(msg))
    end
    if bᵣ < Float64(eps(Float16))
        msg = "The supplied reference strecth bᵣ must be positive valued."
        throw(ErrorException(msg))
    end
    if ((isCGS(cᵣ) && (cᵣ.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(cᵣ) && (cᵣ.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied reference stretch cᵣ is not dimensionless."
        throw(ErrorException(msg))
    end
    if cᵣ < Float64(eps(Float16))
        msg = "The supplied reference strecth cᵣ must be positive valued."
        throw(ErrorException(msg))
    end
    if ((isCGS(αᵣ) && (αᵣ.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(αᵣ) && (αᵣ.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied reference shear αᵣ is not dimensionless."
        throw(ErrorException(msg))
    end
    if ((isCGS(βᵣ) && (βᵣ.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(βᵣ) && (βᵣ.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied reference shear βᵣ is not dimensionless."
        throw(ErrorException(msg))
    end
    if ((isCGS(γᵣ) && (γᵣ.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(γᵣ) && (γᵣ.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied reference shear γᵣ is not dimensionless."
        throw(ErrorException(msg))
    end
    if (F₀.r ≠ 3) || (F₀.c ≠ 3)
        msg = "The supplied deformation gradient F₀ must be a 3x3 matrix."
        throw(DimensionMismatch(msg))
    end
    if ((isCGS(F₀) && (F₀.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(F₀) && (F₀.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied deformation gradient F₀ is not dimensionless."
        throw(ErrorException(msg))
    end
    m = 1
    n = 1
    # Reference Laplace stretch and its inverse.
    Lᵣ = newPhysicalTensor(3, 3, STRETCH)
    if isEulerian
        Lᵣ[1,1] = aᵣ
        Lᵣ[2,1] = aᵣ * γᵣ
        Lᵣ[2,2] = bᵣ
        Lᵣ[3,1] = aᵣ * βᵣ
        Lᵣ[3,2] = bᵣ * αᵣ
        Lᵣ[3,3] = cᵣ
    else # Lagrangian
        Lᵣ[1,1] = aᵣ
        Lᵣ[1,2] = aᵣ * γᵣ
        Lᵣ[1,3] = aᵣ * βᵣ
        Lᵣ[2,2] = bᵣ
        Lᵣ[2,3] = bᵣ * αᵣ
        Lᵣ[3,3] = cᵣ
    end
    Lᵣ⁻¹ = newPhysicalTensor(3, 3, STRETCH)
    if isEulerian
        Lᵣ⁻¹[1,1] = 1 / aᵣ
        Lᵣ⁻¹[2,1] = -γᵣ / bᵣ
        Lᵣ⁻¹[2,2] = 1 / bᵣ
        Lᵣ⁻¹[3,1] = -(βᵣ - αᵣ * γᵣ) / cᵣ
        Lᵣ⁻¹[3,2] = -αᵣ / cᵣ
        Lᵣ⁻¹[3,3] = 1 / cᵣ
    else # Lagrangian
        Lᵣ⁻¹[1,1] = 1 / aᵣ
        Lᵣ⁻¹[1,2] = -γᵣ / bᵣ
        Lᵣ⁻¹[1,3] = -(βᵣ - αᵣ * γᵣ) / cᵣ
        Lᵣ⁻¹[2,2] = 1 / bᵣ
        Lᵣ⁻¹[2,3] = -αᵣ / cᵣ
        Lᵣ⁻¹[3,3] = 1 / cᵣ
    end
    # Create and insert the reference values in their associated arrays.
    (case, 𝐹₀) = pivotF(F₀, isEulerian)
    F = newArrMtx(N, 𝐹₀)
    P = zeros(UInt8, N)
    P[1] = UInt8(case)
    # 𝐹ₙ = LₙGₙ; Lₙ is a triangular Laplace stretch
    # 𝐹ₙ = GₙLₙ; Gₙ is a Gram co-ordinate rotation
    if isEulerian
        (L₀, G₀) = lq(𝐹₀)    # L₀ is lower (left) triangular
    else  # Lagrangian
        (G₀, L₀) = qr(𝐹₀)    # L₀ is upper (right) triangular
    end
    Ω₀ = newPhysicalTensor(3, 3, RATE)
    G = newArrMtx(N, G₀)
    L = newArrMtx(N, L₀)
    Ω = newArrMtx(N, Ω₀)
    if isEulerian
        a₀ = L₀[1,1]
        γ₀ = L₀[2,1] / a₀
        b₀ = L₀[2,2]
        β₀ = L₀[3,1] / a₀
        α₀ = L₀[3,2] / b₀
        c₀ = L₀[3,3]
    else # Lagrangian
        a₀ = L₀[1,1]
        γ₀ = L₀[1,2] / a₀
        β₀ = L₀[1,3] / a₀
        b₀ = L₀[2,2]
        α₀ = L₀[2,3] / b₀
        c₀ = L₀[3,3]
    end
    a = newArrSca(N, a₀)
    b = newArrSca(N, b₀)
    c = newArrSca(N, c₀)
    α = newArrSca(N, α₀)
    β = newArrSca(N, β₀)
    γ = newArrSca(N, γ₀)
    L₀⁻¹ = newPhysicalTensor(3, 3, STRETCH)
    if isEulerian
        L₀⁻¹[1,1] = 1 / a₀
        L₀⁻¹[2,1] = -γ₀ / b₀
        L₀⁻¹[2,2] = 1 / b₀
        L₀⁻¹[3,1] = -(β₀ - α₀ * γ₀) / c₀
        L₀⁻¹[3,2] = -α₀ / c₀
        L₀⁻¹[3,3] = 1 / c₀
    else # Lagrangian
        L₀⁻¹[1,1] = 1 / a₀
        L₀⁻¹[1,2] = -γ₀ / b₀
        L₀⁻¹[1,3] = -(β₀ - α₀ * γ₀) / c₀
        L₀⁻¹[2,2] = 1 / b₀
        L₀⁻¹[2,3] = -α₀ / c₀
        L₀⁻¹[3,3] = 1 / c₀
    end
    L⁻¹ = newArrMtx(N, L₀⁻¹)
    da₀ = newPhysicalScalar(0.0, RATE)
    da  = newArrSca(N, da₀)
    db₀ = newPhysicalScalar(0.0, RATE)
    db  = newArrSca(N, db₀)
    dc₀ = newPhysicalScalar(0.0, RATE)
    dc  = newArrSca(N, dc₀)
    dα₀ = newPhysicalScalar(0.0, RATE)
    dα  = newArrSca(N, dα₀)
    dβ₀ = newPhysicalScalar(0.0, RATE)
    dβ  = newArrSca(N, dβ₀)
    dγ₀ = newPhysicalScalar(0.0, RATE)
    dγ  = newArrSca(N, dγ₀)
    δ₀  = newPhysicalScalar(log(∛((a₀ * b₀ * c₀) / (aᵣ * bᵣ * cᵣ))), STRETCH)
    δ   = newArrSca(N, δ₀)
    ϵ₁₀ = newPhysicalScalar(log(∛((a₀ / aᵣ) * (bᵣ / b₀))), STRETCH)
    ϵ₁  = newArrSca(N, ϵ₁₀)
    ϵ₂₀ = newPhysicalScalar(log(∛((b₀ / bᵣ) * (cᵣ / c₀))), STRETCH)
    ϵ₂  = newArrSca(N, ϵ₂₀)
    ϵ₃₀ = newPhysicalScalar(log(∛((c₀ / cᵣ) * (aᵣ / a₀))), STRETCH)
    ϵ₃  = newArrSca(N, ϵ₃₀)
    γ₁₀ = α₀ - αᵣ
    γ₁  = newArrSca(N, γ₁₀)
    γ₂₀ = β₀ - βᵣ
    γ₂  = newArrSca(N, γ₂₀)
    γ₃₀ = γ₀ - γᵣ
    γ₃  = newArrSca(N, γ₃₀)
    dδ₀  = newPhysicalScalar(0.0, RATE)
    dδ   = newArrSca(N, dδ₀)
    dϵ₁₀ = newPhysicalScalar(0.0, RATE)
    dϵ₁  = newArrSca(N, dϵ₁₀)
    dϵ₂₀ = newPhysicalScalar(0.0, RATE)
    dϵ₂  = newArrSca(N, dϵ₂₀)
    dϵ₃₀ = newPhysicalScalar(0.0, RATE)
    dϵ₃  = newArrSca(N, dϵ₃₀)
    dγ₁₀ = newPhysicalScalar(0.0, RATE)
    dγ₁  = newArrSca(N, dγ₁₀)
    dγ₂₀ = newPhysicalScalar(0.0, RATE)
    dγ₂  = newArrSca(N, dγ₂₀)
    dγ₃₀ = newPhysicalScalar(0.0, RATE)
    dγ₃  = newArrSca(N, dγ₃₀)
    # Create and return a data structure for Laplace kinematics.
    return Kinematics3D(isEulerian, dt, UInt32(N), MInteger(n), MInteger(m), F, P, L, L⁻¹, G, Ω, Lᵣ, Lᵣ⁻¹, aᵣ, bᵣ, cᵣ, αᵣ, βᵣ, γᵣ, a, b, c, α, β, γ, da, db, dc, dα, dβ, dγ, δ, ϵ₁, ϵ₂, ϵ₃, γ₁, γ₂, γ₃, dδ, dϵ₁, dϵ₂, dϵ₃, dγ₁, dγ₂, dγ₃)
end

"""
Function:\n
    advance!(k, Fₙ, restart)\n
This function inserts new entries into the various data arrays supplied by a kinematic data structure 'k' of type 'Kinematics3D'. These new entries associate with a configuration κₙ, wherein the most recent data held by these arrays associate with a configuration κₙ₋₁. These new data are computed from an user supplied deformation gradient 'Fₙ' evaluated in the user's co-ordinate system with base vectors (𝕚, 𝕛, 𝕜). The data arrays have components evaluated in a permuted, frame-indifferent, co-ordinate system with base vectors (𝔼₁, 𝔼₂, 𝔼₃). This permuted co-ordinate frame is suitable for constitutive construction. Whenever a lose in differentiability exists in a deformation history, e.g., a wavefront passes through a mass point, the 'restart' flag is to be toggled to true. By default it is set to false. This flag affects how derivatives are evaluated.
"""
function advance!(k::Kinematics3D, Fₙ::PhysicalTensor, restart::Bool=false)
    # Verify inputs.
    if (Fₙ.r ≠ 3) || (Fₙ.c ≠ 3)
        msg = "The supplied deformation gradient Fₙ must be a 3x3 matrix."
        throw(DimensionMismatch(msg))
    end
    if ((isCGS(Fₙ) && (Fₙ.u ≠ PhysicalFields.CGS_DIMENSIONLESS)) ||
        (isSI(Fₙ) && (Fₙ.u ≠ PhysicalFields.SI_DIMENSIONLESS)))
        msg = "The supplied deformation gradient Fₙ is not dimensionless."
        throw(ErrorException(msg))
    end
    # Advance the fields, i.e., insert new values onto the arrays.
    k.n.n = k.n.n + 1
    n = toInteger(k.n)
    (Pₙ, 𝐹ₙ) = pivotF(Fₙ, k.isEulerian)
    if restart || ((n > 1) && (Pₙ ≠ k.P[n-1]))
        k.m.n = 1
    else
        k.m.n = k.m.n + 1
    end
    m = toInteger(k.m)
    # Insert new values for the nᵗʰ entry of the kinematic arrays.
    k.F[n] = 𝐹ₙ
    k.P[n] = UInt8(Pₙ)
    # 𝐹ₙ = LₙGₙ; Lₙ is a triangular Laplace stretch
    # 𝐹ₙ = GₙLₙ; Gₙ is a Gram co-ordinate rotation
    if k.isEulerian
        (Lₙ, Gₙ) = lq(𝐹ₙ)    # Lₙ is lower (left) triangular
    else  # Lagrangian
        (Gₙ, Lₙ) = qr(𝐹ₙ)    # Lₙ is upper (right) triangular
    end
    k.G[n] = Gₙ
    k.L[n] = Lₙ
    if k.isEulerian
        aₙ = Lₙ[1,1]
        γₙ = Lₙ[2,1] / aₙ
        bₙ = Lₙ[2,2]
        βₙ = Lₙ[3,1] / aₙ
        αₙ = Lₙ[3,2] / bₙ
        cₙ = Lₙ[3,3]
    else # Lagrangian
        aₙ = Lₙ[1,1]
        γₙ = Lₙ[1,2] / aₙ
        βₙ = Lₙ[1,3] / aₙ
        bₙ = Lₙ[2,2]
        αₙ = Lₙ[2,3] / bₙ
        cₙ = Lₙ[3,3]
    end
    k.a[n] = aₙ
    k.b[n] = bₙ
    k.c[n] = cₙ
    k.α[n] = αₙ
    k.β[n] = βₙ
    k.γ[n] = γₙ
    Lₙ⁻¹ = newPhysicalTensor(3, 3, STRETCH)
    if k.isEulerian
        Lₙ⁻¹[1,1] = 1 / aₙ
        Lₙ⁻¹[2,1] = -γₙ / bₙ
        Lₙ⁻¹[2,2] = 1 / bₙ
        Lₙ⁻¹[3,1] = -(βₙ - αₙ * γₙ) / cₙ
        Lₙ⁻¹[3,2] = -αₙ / cₙ
        Lₙ⁻¹[3,3] = 1 / cₙ
    else # Lagrangian
        Lₙ⁻¹[1,1] = 1 / aₙ
        Lₙ⁻¹[1,2] = -γₙ / bₙ
        Lₙ⁻¹[1,3] = -(βₙ - αₙ * γₙ) / cₙ
        Lₙ⁻¹[2,2] = 1 / bₙ
        Lₙ⁻¹[2,3] = -αₙ / cₙ
        Lₙ⁻¹[3,3] = 1 / cₙ
    end
    k.L⁻¹[n] = Lₙ⁻¹
    # Spin is computed very differently in 3-space than in 2-space.
    if m == 1
        # There is no spin across a non-differentiability in deformation.
        Ωₙ = newPhysicalTensor(3, 3, RATE)
    elseif m == 2
        # Use first-order finite difference formula to estimate spin Ω.
        I = newPhysicalTensor(3, 3, DIMENSIONLESS)
        for i in 1:3
            I[i,i] = newPhysicalScalar(1, DIMENSIONLESS)
        end
        Ωₙ₋₁ = spin(k.G[n], k.G[n-1], I, k.dt)
        k.Ω[n-1] = Ωₙ₋₁
        Ωₙ = deepcopy(Ωₙ₋₁)
    else
        # Use second-order finite difference formula to estimate spin Ω.
        Ωₙ = spin(k.G[n], k.G[n-1], k.G[n-2], k.dt)
    end
    k.Ω[n] = Ωₙ
    # Rate for elongation a.
    if m == 1
        # There is no velocity da at the IC or at a wavefront.
        daₙ = newPhysicalScalar(0.0, RATE)
        k.da[n] = daₙ
    elseif m == 2
        # Use first-order forward difference to update previous velocity da.
        daₙ₋₁ = (k.a[n] - k.a[n-1]) / k.dt
        k.da[n-1] = daₙ₋₁
        # Use first-order backward difference to compute current velocity da.
        daₙ = (k.a[n] - k.a[n-1]) / k.dt
        k.da[n] = daₙ
    elseif m == 3
        # Use second-order forward difference to update prevprev velocity da.
        daₙ₋₂ = (-k.a[n] + 4k.a[n-1] - 3k.a[n-2]) / (2k.dt)
        k.da[n-2] = daₙ₋₂
        # Use second-order central difference to update previous velocity da.
        daₙ₋₁ = (k.a[n] - k.a[n-2]) / (2k.dt)
        k.da[n-1] = daₙ₋₁
        # Use second-order backward difference to compute current velocity da.
        daₙ = (3k.a[n] - 4k.a[n-1] + k.a[n-2]) / (2k.dt)
        k.da[n] = daₙ
    else
        # Use second-order backward difference to compute current velocity da.
        daₙ = (3k.a[n] - 4k.a[n-1] + k.a[n-2]) / (2k.dt)
        k.da[n] = daₙ
    end
    # Rate for elongation b.
    if m == 1
        # There is no velocity db at the IC or at a wavefront.
        dbₙ = newPhysicalScalar(0.0, RATE)
        k.db[n] = dbₙ
    elseif m == 2
        # Use first-order forward difference to update previous velocity db.
        dbₙ₋₁ = (k.b[n] - k.b[n-1]) / k.dt
        k.db[n-1] = dbₙ₋₁
        # Use first-order backward difference to compute current velocity db.
        dbₙ = (k.b[n] - k.b[n-1]) / k.dt
        k.db[n] = dbₙ
    elseif m == 3
        # Use second-order forward difference to update prevprev velocity db.
        dbₙ₋₂ = (-k.b[n] + 4k.b[n-1] - 3k.b[n-2]) / (2k.dt)
        k.db[n-2] = dbₙ₋₂
        # Use second-order central difference to update previous velocity db.
        dbₙ₋₁ = (k.b[n] - k.b[n-2]) / (2k.dt)
        k.db[n-1] = dbₙ₋₁
        # Use second-order backward difference to compute current velocity db.
        dbₙ = (3k.b[n] - 4k.b[n-1] + k.b[n-2]) / (2k.dt)
        k.db[n] = dbₙ
    else
        # Use second-order backward difference to compute current velocity db.
        dbₙ = (3k.b[n] - 4k.b[n-1] + k.b[n-2]) / (2k.dt)
        k.db[n] = dbₙ
    end
    # Rate for elongation c.
    if m == 1
        # There is no velocity dc at the IC or at a wavefront.
        dcₙ = newPhysicalScalar(0.0, RATE)
        k.dc[n] = dcₙ
    elseif m == 2
        # Use first-order forward difference to update previous velocity dc.
        dcₙ₋₁ = (k.c[n] - k.c[n-1]) / k.dt
        k.dc[n-1] = dcₙ₋₁
        # Use first-order backward difference to compute current velocity dc.
        dcₙ = (k.c[n] - k.c[n-1]) / k.dt
        k.dc[n] = dcₙ
    elseif m == 3
        # Use second-order forward difference to update prevprev velocity dc.
        dcₙ₋₂ = (-k.c[n] + 4k.c[n-1] - 3k.c[n-2]) / (2k.dt)
        k.dc[n-2] = dcₙ₋₂
        # Use second-order central difference to update previous velocity dc.
        dcₙ₋₁ = (k.c[n] - k.c[n-2]) / (2k.dt)
        k.dc[n-1] = dcₙ₋₁
        # Use second-order backward difference to compute current velocity dc.
        dcₙ = (3k.c[n] - 4k.c[n-1] + k.c[n-2]) / (2k.dt)
        k.dc[n] = dcₙ
    else
        # Use second-order backward difference to compute current velocity dc.
        dcₙ = (3k.c[n] - 4k.c[n-1] + k.c[n-2]) / (2k.dt)
        k.dc[n] = dcₙ
    end
    # Rate for shear α.
    if m == 1
        # There is no velocity dα at the IC or at a wavefront.
        dαₙ = newPhysicalScalar(0.0, RATE)
        k.dα[n] = dαₙ
    elseif m == 2
        # Use first-order forward difference to update previous velocity dα.
        dαₙ₋₁ = (k.b[n] / k.b[n-1]) * (k.α[n] - k.α[n-1]) / k.dt
        k.dα[n-1] = dαₙ₋₁
        # Use first-order backward difference to compute current velocity dα.
        dαₙ = (k.b[n-1] / k.b[n]) * (k.α[n] - k.α[n-1]) / k.dt
        k.dα[n] = dαₙ
    elseif m == 3
        # Use second-order forward difference to update prevprev velocity dα.
        dαₙ₋₂ = (2(k.b[n-1] / k.b[n-2]) * (k.α[n-1] - k.α[n-2]) / k.dt
            - (k.b[n] / k.b[n-2]) * (k.α[n] - k.α[n-2]) / (2k.dt))
        k.dα[n-2] = dαₙ₋₂
        # Use second-order central difference to update previous velocity dα.
        dαₙ₋₁ = ((k.b[n] / k.b[n-1]) * (k.α[n] - k.α[n-1]) / (2k.dt)
            + (k.b[n-2] / k.b[n-1]) * (k.α[n-1] - k.α[n-2]) / (2k.dt))
        k.dα[n-1] = dαₙ₋₁
        # Use second-order backward difference to compute current velocity dα.
        dαₙ = (2(k.b[n-1] / k.b[n]) * (k.α[n] - k.α[n-1]) / k.dt
            - (k.b[n-2] / k.b[n]) * (k.α[n] - k.α[n-2]) / (2k.dt))
        k.dα[n] = dαₙ
    else
        # Use second-order backward difference to compute current velocity dα.
        dαₙ = (2(k.b[n-1] / k.b[n]) * (k.α[n] - k.α[n-1]) / k.dt
            - (k.b[n-2] / k.b[n]) * (k.α[n] - k.α[n-2]) / (2k.dt))
        k.dα[n] = dαₙ
    end
    # Rate for shear β.
    if m == 1
        # There is no velocity dβ at the IC or at a wavefront.
        dβₙ = newPhysicalScalar(0.0, RATE)
        k.dβ[n] = dβₙ
    elseif m == 2
        # Use first-order forward difference to update previous velocity dβ.
        dβₙ₋₁ = (k.a[n] / k.a[n-1]) * (k.β[n] - k.β[n-1]) / k.dt
        k.dβ[n-1] = dβₙ₋₁
        # Use first-order backward difference to compute current velocity dβ.
        dβₙ = (k.a[n-1] / k.a[n]) * (k.β[n] - k.β[n-1]) / k.dt
        k.dβ[n] = dβₙ
    elseif m == 3
        # Use second-order forward difference to update prevprev velocity dβ.
        dβₙ₋₂ = (2(k.a[n-1] / k.a[n-2]) * (k.β[n-1] - k.β[n-2]) / k.dt
            - (k.a[n] / k.a[n-2]) * (k.β[n] - k.β[n-2]) / (2k.dt))
        k.dβ[n-2] = dβₙ₋₂
        # Use second-order central difference to update previous velocity dβ.
        dβₙ₋₁ = ((k.a[n] / k.a[n-1]) * (k.β[n] - k.β[n-1]) / (2k.dt)
            + (k.a[n-2] / k.a[n-1]) * (k.β[n-1] - k.β[n-2]) / (2k.dt))
        k.dβ[n-1] = dβₙ₋₁
        # Use second-order backward difference to compute current velocity dβ.
        dβₙ = (2(k.a[n-1] / k.a[n]) * (k.β[n] - k.β[n-1]) / k.dt
            - (k.a[n-2] / k.a[n]) * (k.β[n] - k.β[n-2]) / (2k.dt))
        k.dβ[n] = dβₙ
    else
        # Use second-order backward difference to compute current velocity dβ.
        dβₙ = (2(k.a[n-1] / k.a[n]) * (k.β[n] - k.β[n-1]) / k.dt
            - (k.a[n-2] / k.a[n]) * (k.β[n] - k.β[n-2]) / (2k.dt))
        k.dβ[n] = dβₙ
    end
    # Rate for shear γ.
    if m == 1
        # There is no velocity dγ at the IC or at a wavefront.
        dγₙ = newPhysicalScalar(0.0, RATE)
        k.dγ[n] = dγₙ
    elseif m == 2
        # Use first-order forward difference to update previous velocity dγ.
        dγₙ₋₁ = (k.a[n] / k.a[n-1]) * (k.γ[n] - k.γ[n-1]) / k.dt
        k.dγ[n-1] = dγₙ₋₁
        # Use first-order backward difference to compute current velocity dγ.
        dγₙ = (k.a[n-1] / k.a[n]) * (k.γ[n] - k.γ[n-1]) / k.dt
        k.dγ[n] = dγₙ
    elseif m == 3
        # Use second-order forward difference to update prevprev velocity dγ.
        dγₙ₋₂ = (2(k.a[n-1] / k.a[n-2]) * (k.γ[n-1] - k.γ[n-2]) / k.dt
            - (k.a[n] / k.a[n-2]) * (k.γ[n] - k.γ[n-2]) / (2k.dt))
        k.dγ[n-2] = dγₙ₋₂
        # Use second-order central difference to update previous velocity dγ.
        dγₙ₋₁ = ((k.a[n] / k.a[n-1]) * (k.γ[n] - k.γ[n-1]) / (2k.dt)
            + (k.a[n-2] / k.a[n-1]) * (k.γ[n-1] - k.γ[n-2]) / (2k.dt))
        k.dγ[n-1] = dγₙ₋₁
        # Use second-order backward difference to compute current velocity dγ.
        dγₙ = (2(k.a[n-1] / k.a[n]) * (k.γ[n] - k.γ[n-1]) / k.dt
            - (k.a[n-2] / k.a[n]) * (k.γ[n] - k.γ[n-2]) / (2k.dt))
        k.dγ[n] = dγₙ
    else
        # Use second-order backward difference to compute current velocity dγ.
        dγₙ = (2(k.a[n-1] / k.a[n]) * (k.γ[n] - k.γ[n-1]) / k.dt
            - (k.a[n-2] / k.a[n]) * (k.γ[n] - k.γ[n-2]) / (2k.dt))
        k.dγ[n] = dγₙ
    end
    δₙ = newPhysicalScalar(log(∛((aₙ * bₙ * cₙ) / (k.aᵣ * k.bᵣ * k.cᵣ))), STRETCH)
    k.δ[n] = δₙ
    ϵ₁ₙ = newPhysicalScalar(log(∛((aₙ / k.aᵣ) * (k.bᵣ / bₙ))), STRETCH)
    k.ϵ₁[n] = ϵ₁ₙ
    ϵ₂ₙ = newPhysicalScalar(log(∛((bₙ / k.bᵣ) * (k.cᵣ / cₙ))), STRETCH)
    k.ϵ₂[n] = ϵ₂ₙ
    ϵ₃ₙ = newPhysicalScalar(log(∛((cₙ / k.cᵣ) * (k.aᵣ / aₙ))), STRETCH)
    k.ϵ₃[n] = ϵ₃ₙ
    γ₁ₙ = αₙ - k.αᵣ
    k.γ₁[n] = γ₁ₙ
    γ₂ₙ = βₙ - k.βᵣ
    k.γ₂[n] = γ₂ₙ
    γ₃ₙ = γₙ - k.γᵣ
    k.γ₃[n] = γ₃ₙ
    dδₙ = (daₙ / aₙ + dbₙ / bₙ + dcₙ / cₙ) / 3
    k.dδ[n] = dδₙ
    dϵ₁ₙ = (daₙ / aₙ - dbₙ / bₙ) / 3
    k.dϵ₁[n] = dϵ₁ₙ
    dϵ₂ₙ = (dbₙ / bₙ - dcₙ / cₙ) / 3
    k.dϵ₂[n] = dϵ₂ₙ
    dϵ₃ₙ = (dcₙ / cₙ - daₙ / aₙ) / 3
    k.dϵ₃[n] = dϵ₃ₙ
    dγ₁ₙ = dαₙ
    k.dγ₁[n] = dγ₁ₙ
    dγ₂ₙ = dβₙ
    k.dγ₂[n] = dγ₂ₙ
    dγ₃ₙ = dγₙ
    k.dγ₃[n] = dγ₃ₙ
    return nothing
end # advance!

"""
Function:\n
    update!(k, Fₙ)\n
This function updates (overwrites) the most recent entries found in the various arrays supplied by a kinematic data structure 'k' of type 'Kinematics3D'. These entries associate with a configuration κₙ. The updated data are computed from an user supplied, updated, deformation gradient 'Fₙ' evaluated in the user's co-ordinate system with base vectors (𝕚, 𝕛, 𝕜). All data arrays in 'k' have components evaluated in a frame-indifferent co-ordinate system with base vectors (𝔼₁, 𝔼₂, 𝔼₃). This co-ordinate frame is suitable for constitutive construction. Function update! can be called multiple times between calls to advance!, the latter of which propagates a solution along its path. Whenever a lose in differentiability occurs in a deformation history, e.g., a wavefront passes through a mass point, the flag 'restart' is to be toggled to true. By default it is set to false. This flag affects how derivatives are evaluated.
"""
function update!(k::Kinematics3D, Fₙ::PhysicalTensor, restart::Bool=false)
    if k.m > 1
        k.m.n = k.m.n - 1
    end
    if k.n > 1
        k.n.n = k.n.n - 1
    end
    advance!(k, Fₙ, restart)
    return nothing
end # update!
