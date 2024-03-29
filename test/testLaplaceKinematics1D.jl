module testLaplaceKinematics1D

using
    CairoMakie,       # Pixel based figure construction.
    FijLung,
    LaplaceKinematics,
    PhysicalFields

export
    persistence,
    figures1D
#=
--------------------------------------------------------------------------------
=#

"""
Function:\n
    persistence(midPtQuad, myDirPath)\n
where\n
    midPtQuad  is true if a mid-point quadrature is to be used.\n
    myDirPath  is the path where the user wants his file to be written.\n
This function tests writing and reading a `FiberKinematics` object to and from a file for its ability to recreate the object from file.
"""
function persistence(midPtQuad::Bool, myDirPath::String)
    N = 3 # Considered so the JSON file would not be too big.

    arrayOfTimes = FijLung.t_loc1()
    Nₛ = arrayOfTimes.array.len
    if midPtQuad
        splineF = FijLung.splineAtMidPoints(1, Nₛ)
    else
        splineF = FijLung.splineAtEndPoints(1, Nₛ)
    end

    # Consider the reference and initial states to be the same, i.e., κᵣ = κ₀.
    λᵣ = PhysicalScalar(1.0, CGS_DIMENSIONLESS)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = splineF.t[N] - splineF.t[N-1]
    k = LaplaceKinematics.FiberKinematics(dt, N, midPtQuad, λᵣ)

    # Populate this data structure.
    for n in 1:N
        F′ₙ = splineF.F′[n+1]
        λ′ₙ = F′ₙ[2,2]
        LaplaceKinematics.advance!(k, λ′ₙ)
    end

    # Convert this data structure to a JSON stream.
    json_stream = PhysicalFields.openJSONWriter(myDirPath, "test1D.json")
    LaplaceKinematics.toFile(k, json_stream)
    PhysicalFields.closeJSONStream(json_stream)

    # Retrieve this data structure from a JSON stream.
    json_stream = PhysicalFields.openJSONReader(myDirPath, "test1D.json")
    k1 = LaplaceKinematics.fromFile(LaplaceKinematics.FiberKinematics, json_stream)
    PhysicalFields.closeJSONStream(json_stream)

    # Verify what was read in is equivalent to what was written to.
    equal = true
    if k.dt ≠ k1.dt
        equal = false
    end
    if k.N ≠ k1.N
        equal = false
    end
    if k.λᵣ ≠ k1.λᵣ
        equal = false
    end
    for i in 1:N+1
        if k.t[i] ≠ k1.t[i]
            equal = false
        end
        if k.λ[i] ≠ k1.λ[i]
            equal = false
        end
        if k.λ′[i] ≠ k1.λ′[i]
            equal = false
        end
        if k.ϵ[i] ≠ k1.ϵ[i]
            equal = false
        end
        if k.ϵ′[i] ≠ k1.ϵ′[i]
            equal = false
        end
    end
    if equal
        println("PASSED: The retrieved object equals the saved object.")
    else
        println("FAILED: The retrieved object did not equal the saved object.")
    end
end # persistence

"""
Function:\n
    figures1D(N, midPtQuad, myDirPath)\n
where\n
    N          is the number of nodes or knots in the B-spline.\n
    midPtQuad  is true if a mid-point quadrature is to be used.\n
    myDirPath  is the path where the user wants his file to be written.\n
This function tests the exported functions of `LaplaceKinematics` for the 1D case; specifically, in the 2 direction, viz., in the direction of the spine.
"""
function figures1D(N::Integer, midPtQuad::Bool, myDirPath::String)

    CairoMakie.activate!(type = "png")
    println("For these figures, the number of intervals is N = ", string(N), ".")

    println("Creating B-splines for the deformation gradient data.")
    if midPtQuad
        splineF1 = FijLung.splineAtMidPoints(1, N) # near visceral pleura
        splineF2 = FijLung.splineAtMidPoints(2, N) # deep in the lung
        splineF3 = FijLung.splineAtMidPoints(3, N) # near bronchiole tube
    else
        splineF1 = FijLung.splineAtEndPoints(1, N) # near visceral pleura
        splineF2 = FijLung.splineAtEndPoints(2, N) # deep in the lung
        splineF3 = FijLung.splineAtEndPoints(3, N) # near bronchiole tube
    end

    # The last splined node is not included in the plots because it may not be
    # 'natural', i.e., actual derivatives at the last node may not be accurate.

    # Times at the solution nodal points. Times t[1] hold the initial time.
    t1 = zeros(Float64, N)
    t2 = zeros(Float64, N)
    t3 = zeros(Float64, N)
    for n in 1:N
        t1[n] = get(splineF1.t[n])
        t2[n] = get(splineF2.t[n])
        t3[n] = get(splineF3.t[n])
    end

    # These deformation gradient components associate with the 2 direction.
    println("Deformations are in the direction of the spine.")

    # Create data arrays for the fiber stretches λ at the solution nodes.
    # Stretches λ[1] hold their initial values.
    λ1 = zeros(Float64, N)
    λ2 = zeros(Float64, N)
    λ3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        λ1[n] = get(Fᵢⱼ1[2,2])
        Fᵢⱼ2 = splineF2.F[n]
        λ2[n] = get(Fᵢⱼ2[2,2])
        Fᵢⱼ3 = splineF3.F[n]
        λ3[n] = get(Fᵢⱼ3[2,2])
    end
    # Create data arrays for the fiber stretch rates dλ/dt at the nodes.
    λ′1 = zeros(Float64, N)
    λ′2 = zeros(Float64, N)
    λ′3 = zeros(Float64, N)
    for n in 1:N
        F′ᵢⱼ1 = splineF1.F′[n]
        λ′1[n] = get(F′ᵢⱼ1[2,2])
        F′ᵢⱼ2 = splineF2.F′[n]
        λ′2[n] = get(F′ᵢⱼ2[2,2])
        F′ᵢⱼ3 = splineF3.F′[n]
        λ′3[n] = get(F′ᵢⱼ3[2,2])
    end

    println("Building the Laplace kinematics data structures.")
    # Consider the reference and initial states to be the same, i.e., κᵣ = κ₀.
    λᵣ = PhysicalScalar(1.0, CGS_LENGTH)
    # Create the variable to hold stretch rates.
    TIME_RATE = PhysicalFields.PhysicalUnits("CGS", 0, 0, 0, -1, 0, 0, 0)
    λ′ₙ = PhysicalScalar(TIME_RATE)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = PhysicalScalar(t1[N]-t1[N-1], CGS_SECOND)
    k1 = LaplaceKinematics.FiberKinematics(dt, N, midPtQuad, λᵣ)
    # Populate this data structure.
    for n in 2:N
        set!(λ′ₙ, λ′1[n])
        LaplaceKinematics.advance!(k1, λ′ₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 2.
    set!(dt, t2[N]-t2[N-1])
    k2 = LaplaceKinematics.FiberKinematics(dt, N, midPtQuad, λᵣ)
    # Populate this data structure.
    for n in 2:N
        set!(λ′ₙ, λ′2[n])
        LaplaceKinematics.advance!(k2, λ′ₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 3.
    set!(dt, t3[N]-t3[N-1])
    k3 = LaplaceKinematics.FiberKinematics(dt, N, midPtQuad, λᵣ)
    # Populate this data structure.
    for n in 2:N
        set!(λ′ₙ, λ′3[n])
        LaplaceKinematics.advance!(k3, λ′ₙ)
    end

    # Create the figures.
    println("Working on figures for ϵ and dϵ/dt.")
    ϵ1 = zeros(Float64, N)
    ϵ′1 = zeros(Float64, N)
    ϵ2 = zeros(Float64, N)
    ϵ′2 = zeros(Float64, N)
    ϵ3 = zeros(Float64, N)
    ϵ′3 = zeros(Float64, N)
    for n in 1:N
        ϵ1[n] = get(k1.ϵ[n])
        ϵ′1[n] = get(k1.ϵ′[n])
        ϵ2[n] = get(k2.ϵ[n])
        ϵ′2[n] = get(k2.ϵ′[n])
        ϵ3[n] = get(k3.ϵ[n])
        ϵ′3[n] = get(k3.ϵ′[n])
    end
    fig = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "strain ϵ",
        title = "Strain ϵ and its rate dϵ/dt.",
        titlesize = 24,
        ylabelsize = 20)
    lines!(ax, t1, ϵ1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, ϵ2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, ϵ3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "1DStrainAtMidPoints.png"), fig)
    else
        save(string(myDirPath, "1DStrainAtEndPoints.png"), fig)
    end

    fig = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "strain rate dϵ/dt (s⁻¹)",
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, ϵ′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, ϵ′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, ϵ′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "1DStrainRateAtMidPoints.png"), fig)
    else
        save(string(myDirPath, "1DStrainRateAtEndPoints.png"), fig)
    end
end # figures1D

end # testLaplaceKinematics1D.jl