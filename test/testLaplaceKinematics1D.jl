module testLaplaceKinematics1D

using
    CairoMakie,       # Pixel based figure construction.
    FijLung,
    PhysicalFields,
    ..LaplaceKinematics

export
    persistence,
    figures1D
#=
--------------------------------------------------------------------------------
=#

function persistence(myDirPath::String)
    arrayOfTimes = FijLung.t_loc1()
    N = arrayOfTimes.array.len
    splineF = FijLung.SplineF(1, N)

    # Consider the reference and initial states to be the same, i.e., κᵣ = κ₁.
    Lᵣ = PhysicalScalar(1.0, CGS_LENGTH)

    # Build a data structure for Laplace kinematics at lung location 1.
    N = 3 # Considered so the JSON file would not be to big.
    dt = splineF.t[2] - splineF.t[1]
    L₀ = PhysicalScalar(1.0, CGS_LENGTH)
    F′₀ = splineF.F′[1]
    midPtQuad = false
    k = LaplaceKinematics.FiberKinematics(dt, N, midPtQuad, Lᵣ, L₀)

    # Populate this data structure.
    for n in 2:N
        F′ₙ = splineF.F′[n]
        L′ₙ = L₀ * F′ₙ[2,2]
        LaplaceKinematics.advance!(k, L′ₙ)
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
    if k.Lᵣ ≠ k1.Lᵣ
        equal = false
    end
    for i in 1:3
        if k.t[i] ≠ k1.t[i]
            equal = false
        end
        if k.L[i] ≠ k1.L[i]
            equal = false
        end
        if k.L′[i] ≠ k1.L′[i]
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
    figures1D(N)\n
where\n
    N is the number of nodes in the B-spline.\n
This function tests the exported functions of `LaplaceKinematics` for the 1D case; specifically, in the 1 direction, viz., from shoulder to shoulder.
"""
function figures1D(N::Integer, myDirPath::String)
    CairoMakie.activate!(type = "png")
    println("For these figures, N = ", string(N), ".")

    println("Creating B-splines of the deformation gradient data.")
    splineF1 = FijLung.SplineF(1, N) # lung location 1: near visceral pleura
    splineF2 = FijLung.SplineF(2, N) # lung location 2: deep lung
    splineF3 = FijLung.SplineF(3, N) # lung location 3: near bronchiole tube

    t1 = zeros(Float64, N)
    t2 = zeros(Float64, N)
    t3 = zeros(Float64, N)
    for n in 1:N
        t1[n] = get(splineF1.t[n])
        t2[n] = get(splineF2.t[n])
        t3[n] = get(splineF3.t[n])
    end
    midPtQuad = false

    # These deformation gradient components associate with the 1 direction.
    println("Deformations are in the direction of the spine.")
    # Create data arrays for a fiber of unit length.
    L1 = zeros(Float64, N)
    L2 = zeros(Float64, N)
    L3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        L1[n] = get(Fᵢⱼ1[2,2])
        Fᵢⱼ2 = splineF2.F[n]
        L2[n] = get(Fᵢⱼ2[2,2])
        Fᵢⱼ3 = splineF3.F[n]
        L3[n] = get(Fᵢⱼ3[2,2])
    end
    # Create the data arrays for dL/dt.
    L′1 = zeros(Float64, N)
    L′2 = zeros(Float64, N)
    L′3 = zeros(Float64, N)
    for n in 1:N
        F′ᵢⱼ1 = splineF1.F′[n]
        L′1[n] = get(F′ᵢⱼ1[2,2])
        F′ᵢⱼ2 = splineF2.F′[n]
        L′2[n] = get(F′ᵢⱼ2[2,2])
        F′ᵢⱼ3 = splineF3.F′[n]
        L′3[n] = get(F′ᵢⱼ3[2,2])
    end

    println("Building the Laplace kinematics data structures.")
    # Consider the reference and initial states to be the same, i.e., κᵣ = κ₁.
    Lᵣ = PhysicalScalar(1.0, CGS_LENGTH)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = PhysicalScalar(t1[2]-t1[1], CGS_SECOND)
    L₀ = PhysicalScalar(L1[1], CGS_LENGTH)
    k1 = LaplaceKinematics.FiberKinematics(dt, N, midPtQuad, Lᵣ, L₀)
    # Populate this data structure.
    L′ₙ = PhysicalScalar(CGS_VELOCITY)
    for n in 2:N
        set!(L′ₙ, L′1[n])
        LaplaceKinematics.advance!(k1, L′ₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 2.
    set!(dt, t2[2]-t2[1])
    set!(L₀, L2[1])
    k2 = LaplaceKinematics.FiberKinematics(dt, N, midPtQuad, Lᵣ, L₀)
    # Populate this data structure.
    for n in 2:N
        set!(L′ₙ, L′2[n])
        LaplaceKinematics.advance!(k2, L′ₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 3.
    set!(dt, t3[2]-t3[1])
    set!(L₀, L3[1])
    k3 = LaplaceKinematics.FiberKinematics(dt, N, midPtQuad, Lᵣ, L₀)
    # Populate this data structure.
    for n in 2:N
        set!(L′ₙ, L′3[n])
        LaplaceKinematics.advance!(k3, L′ₙ)
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
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
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
    save(string(myDirPath, "1Dstrain.png"), fig)
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
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
    save(string(myDirPath, "1Dstrainrate.png"), fig)
end # figures1D

end # testLaplaceKinematics1D.jl