module testLaplaceKinematics2D

using
    CairoMakie,      # Pixel based figure construction.
    FijLung,
    PhysicalFields,
    ..LaplaceKinematics

export
    figures2D,
    persistence
#=
-------------------------------------------------------------------------------
=#

function persistence(myDirPath::String)
    arrayOfTimes = FijLung.t_loc1()
    N = arrayOfTimes.array.len
    splineF = FijLung.SplineF(1, N)

    # Consider the reference and initial states to be the same, i.e., κᵣ = κ₁.
    aᵣ = PhysicalScalar(1.0, CGS_STRETCH)
    bᵣ = PhysicalScalar(1.0, CGS_STRETCH)
    γᵣ = PhysicalScalar(CGS_DIMENSIONLESS)

    # Build a data structure for Laplace kinematics at lung location 1.
    N = 3 # Considered so the JSON file would not be to big.
    dt = splineF.t[2] - splineF.t[1]
    F = splineF.F[1]
    F₀ = PhysicalTensor(2, 2, CGS_DIMENSIONLESS)
    F₀[1,1] = F[1,1]
    F₀[1,2] = F[1,2]
    F₀[2,1] = F[2,1]
    F₀[2,2] = F[2,2]
    F′ = splineF.F′[1]
    F′₀ = PhysicalTensor(2, 2, CGS_STRAIN_RATE)
    F′₀[1,1] = F′[1,1]
    F′₀[1,2] = F′[1,2]
    F′₀[2,1] = F′[2,1]
    F′₀[2,2] = F′[2,2]
    k = LaplaceKinematics.MembraneKinematics(dt, N, aᵣ, bᵣ, γᵣ, F₀)

    # Populate this data structure.
    for n in 2:N
        F′ = splineF.F′[n]
        F′ₙ = PhysicalTensor(2, 2, CGS_STRAIN_RATE)
        F′ₙ[1,1] = F′[1,1]
        F′ₙ[1,2] = F′[1,2]
        F′ₙ[2,1] = F′[2,1]
        F′ₙ[2,2] = F′[2,2]
        LaplaceKinematics.advance!(k, F′ₙ)
    end

    # Convert this data structure to a JSON stream.
    json_stream = PhysicalFields.openJSONWriter(myDirPath, "test2D.json")
    LaplaceKinematics.toFile(k, json_stream)
    PhysicalFields.closeJSONStream(json_stream)

    # Retrieve this data structure from a JSON stream.
    json_stream = PhysicalFields.openJSONReader(myDirPath, "test2D.json")
    k1 = LaplaceKinematics.fromFile(LaplaceKinematics.MembraneKinematics, json_stream)
    PhysicalFields.closeJSONStream(json_stream)

    # Verify what was read in is equivalent to what was written to.
    equal = true
    if k.dt ≠ k1.dt
        equal = false
    end
    if k.N ≠ k1.N
        equal = false
    end
    if k.aᵣ ≠ k1.aᵣ
        equal = false
    end
    if k.bᵣ ≠ k1.bᵣ
        equal = false
    end
    if k.γᵣ ≠ k1.γᵣ
        equal = false
    end
    for i in 1:3
        if k.t[i] ≠ k1.t[i]
            equal = false
        end
        if k.F[i] ≠ k1.F[i]
            equal = false
        end
        if k.F′[i] ≠ k1.F′[i]
            equal = false
        end
        if k.P[i] ≠ k1.P[i]
            equal = false
        end
        if k.ωₙ[i] ≠ k1.ωₙ[i]
            equal = false
        end
        if k.ω′ₙ[i] ≠ k1.ω′ₙ[i]
            equal = false
        end
        if k.aₙ[i] ≠ k1.aₙ[i]
            equal = false
        end
        if k.bₙ[i] ≠ k1.bₙ[i]
            equal = false
        end
        if k.γₙ[i] ≠ k1.γₙ[i]
            equal = false
        end
        if k.a′ₙ[i] ≠ k1.a′ₙ[i]
            equal = false
        end
        if k.b′ₙ[i] ≠ k1.b′ₙ[i]
            equal = false
        end
        if k.γ′ₙ[i] ≠ k1.γ′ₙ[i]
            equal = false
        end
        if k.δ[i] ≠ k1.δ[i]
            equal = false
        end
        if k.ϵ[i] ≠ k1.ϵ[i]
            equal = false
        end
        if k.γ[i] ≠ k1.γ[i]
            equal = false
        end
        if k.δ′[i] ≠ k1.δ′[i]
            equal = false
        end
        if k.ϵ′[i] ≠ k1.ϵ′[i]
            equal = false
        end
        if k.γ′[i] ≠ k1.γ′[i]
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
    figures2D(N, myDirPath)\n
where\n
    N is the number of nodes in the B-spline, and\n
    myDirPath specifies the directory into which the figures are saved.\n 
This function tests the exported functions of `LaplaceKinematics` for the 2D case; specifically, for the 23 plane for Fᵢⱼ exported by FijLung.
"""
function figures2D(N::Integer, myDirPath::String)
    CairoMakie.activate!(type = "png")
    println("For these figures, N = ", string(N), ".")

    println("Creating B Splines for the deformation gradient data.")
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

    # These deformation gradient components associate with the 23 plane.
    println("Deformations are for the spine/breast-bone plane.")
    # Create the data arrays for F₁₁.
    F₁₁1 = zeros(Float64, N)
    F₁₁2 = zeros(Float64, N)
    F₁₁3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₁₁1[n] = get(Fᵢⱼ1[2,2])
        Fᵢⱼ2 = splineF2.F[n]
        F₁₁2[n] = get(Fᵢⱼ2[2,2])
        Fᵢⱼ3 = splineF3.F[n]
        F₁₁3[n] = get(Fᵢⱼ3[2,2])
    end
    # Create the data arrays for dF₁₁/dt.
    F′₁₁1 = zeros(Float64, N)
    F′₁₁2 = zeros(Float64, N)
    F′₁₁3 = zeros(Float64, N)
    for n in 1:N
        F′ᵢⱼ1 = splineF1.F′[n]
        F′₁₁1[n] = get(F′ᵢⱼ1[2,2])
        F′ᵢⱼ2 = splineF2.F′[n]
        F′₁₁2[n] = get(F′ᵢⱼ2[2,2])
        F′ᵢⱼ3 = splineF3.F′[n]
        F′₁₁3[n] = get(F′ᵢⱼ3[2,2])
    end
    # Create the data arrays for F₁₂
    F₁₂1 = zeros(Float64, N)
    F₁₂2 = zeros(Float64, N)
    F₁₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₁₂1[n] = get(Fᵢⱼ1[2,3])
        Fᵢⱼ2 = splineF2.F[n]
        F₁₂2[n] = get(Fᵢⱼ2[2,3])
        Fᵢⱼ3 = splineF3.F[n]
        F₁₂3[n] = get(Fᵢⱼ3[2,3])
    end
    # Create the data arrays dF₁₂/dt.
    F′₁₂1 = zeros(Float64, N)
    F′₁₂2 = zeros(Float64, N)
    F′₁₂3 = zeros(Float64, N)
    for n in 1:N
        F′ᵢⱼ1 = splineF1.F′[n]
        F′₁₂1[n] = get(F′ᵢⱼ1[2,3])
        F′ᵢⱼ2 = splineF2.F′[n]
        F′₁₂2[n] = get(F′ᵢⱼ2[2,3])
        F′ᵢⱼ3 = splineF3.F′[n]
        F′₁₂3[n] = get(F′ᵢⱼ3[2,3])
    end
    # Create the data arrays for F₂₁.
    F₂₁1 = zeros(Float64, N)
    F₂₁2 = zeros(Float64, N)
    F₂₁3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₂₁1[n] = get(Fᵢⱼ1[3,2])
        Fᵢⱼ2 = splineF2.F[n]
        F₂₁2[n] = get(Fᵢⱼ2[3,2])
        Fᵢⱼ3 = splineF3.F[n]
        F₂₁3[n] = get(Fᵢⱼ3[3,2])
    end
    # Create the data arrays for dF₂₁/dt.
    F′₂₁1 = zeros(Float64, N)
    F′₂₁2 = zeros(Float64, N)
    F′₂₁3 = zeros(Float64, N)
    for n in 1:N
        F′ᵢⱼ1 = splineF1.F′[n]
        F′₂₁1[n] = get(F′ᵢⱼ1[3,2])
        F′ᵢⱼ2 = splineF2.F′[n]
        F′₂₁2[n] = get(F′ᵢⱼ2[3,2])
        F′ᵢⱼ3 = splineF3.F′[n]
        F′₂₁3[n] = get(F′ᵢⱼ3[3,2])
    end
    # Create the data arrays for F₂₂.
    F₂₂1 = zeros(Float64, N)
    F₂₂2 = zeros(Float64, N)
    F₂₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₂₂1[n] = get(Fᵢⱼ1[3,3])
        Fᵢⱼ2 = splineF2.F[n]
        F₂₂2[n] = get(Fᵢⱼ2[3,3])
        Fᵢⱼ3 = splineF3.F[n]
        F₂₂3[n] = get(Fᵢⱼ3[3,3])
    end
    # Create the data arrays for dF₂₂/dt.
    F′₂₂1 = zeros(Float64, N)
    F′₂₂2 = zeros(Float64, N)
    F′₂₂3 = zeros(Float64, N)
    for n in 1:N
        F′ᵢⱼ1 = splineF1.F′[n]
        F′₂₂1[n] = get(F′ᵢⱼ1[3,3])
        F′ᵢⱼ2 = splineF2.F′[n]
        F′₂₂2[n] = get(F′ᵢⱼ2[3,3])
        F′ᵢⱼ3 = splineF3.F′[n]
        F′₂₂3[n] = get(F′ᵢⱼ3[3,3])
    end

    println("Building the Laplace kinematics data structures.")
    # Consider the reference and initial states to be the same, i.e., κᵣ = κ₁.
    aᵣ = PhysicalScalar(1.0, CGS_STRETCH)
    bᵣ = PhysicalScalar(1.0, CGS_STRETCH)
    γᵣ = PhysicalScalar(CGS_DIMENSIONLESS)

    Fᵢⱼ  = PhysicalScalar(CGS_STRETCH)
    F′ᵢⱼ = PhysicalScalar(CGS_STRETCH_RATE)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt1 = PhysicalScalar(t1[2]-t1[1], CGS_SECOND)
    F₁1 = PhysicalTensor(2, 2, CGS_STRETCH)
    set!(Fᵢⱼ, F₁₁1[1])
    F₁1[1,1] = Fᵢⱼ
    set!(Fᵢⱼ, F₁₂1[1])
    F₁1[1,2] = Fᵢⱼ
    set!(Fᵢⱼ, F₂₁1[1])
    F₁1[2,1] = Fᵢⱼ
    set!(Fᵢⱼ, F₂₂1[1])
    F₁1[2,2] = Fᵢⱼ
    k1 = LaplaceKinematics.MembraneKinematics(dt1, N, midPtQuad, aᵣ, bᵣ, γᵣ, F₁1)
    # Populate this data structure.
    F′ₙ1 = PhysicalTensor(2, 2, CGS_STRETCH_RATE)
    for n in 2:N
        set!(F′ᵢⱼ, F′₁₁1[n])
        F′ₙ1[1,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₁₂1[n])
        F′ₙ1[1,2] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₁1[n])
        F′ₙ1[2,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₂1[n])
        F′ₙ1[2,2] = F′ᵢⱼ
        LaplaceKinematics.advance!(k1, F′ₙ1)
    end

    # Build a data structure for Laplace kinematics at lung location 2.
    dt2 = PhysicalScalar(t2[2]-t2[1], CGS_SECOND)
    F₁2 = PhysicalTensor(2, 2, CGS_STRETCH)
    set!(Fᵢⱼ, F₁₁2[1])
    F₁2[1,1] = Fᵢⱼ
    set!(Fᵢⱼ, F₁₂2[1])
    F₁2[1,2] = Fᵢⱼ
    set!(Fᵢⱼ, F₂₁2[1])
    F₁2[2,1] = Fᵢⱼ
    set!(Fᵢⱼ, F₂₂2[1])
    F₁2[2,2] = Fᵢⱼ
    k2 = LaplaceKinematics.MembraneKinematics(dt2, N, midPtQuad, aᵣ, bᵣ, γᵣ, F₁2)
    F′ₙ2 = PhysicalTensor(2, 2, CGS_STRETCH_RATE)
    for n in 2:N
        set!(F′ᵢⱼ, F′₁₁2[n])
        F′ₙ2[1,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₁₂2[n])
        F′ₙ2[1,2] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₁2[n])
        F′ₙ2[2,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₂2[n])
        F′ₙ2[2,2] = F′ᵢⱼ
        LaplaceKinematics.advance!(k2, F′ₙ2)
    end

    # Build a data structure for Laplace kinematics at lung location 3.
    dt3 = PhysicalScalar(t3[2]-t3[1], CGS_SECOND)
    F₁3 = PhysicalTensor(2, 2, CGS_STRETCH)
    set!(Fᵢⱼ, F₁₁3[1])
    F₁3[1,1] = Fᵢⱼ
    set!(Fᵢⱼ, F₁₂3[1])
    F₁3[1,2] = Fᵢⱼ
    set!(Fᵢⱼ, F₂₁3[1])
    F₁3[2,1] = Fᵢⱼ
    set!(Fᵢⱼ, F₂₂3[1])
    F₁3[2,2] = Fᵢⱼ
    k3 = LaplaceKinematics.MembraneKinematics(dt3, N, midPtQuad, aᵣ, bᵣ, γᵣ, F₁3)
    F′ₙ3 = PhysicalTensor(2, 2, CGS_STRETCH_RATE)
    for n in 2:N
        set!(F′ᵢⱼ, F′₁₁3[n])
        F′ₙ3[1,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₁₂3[n])
        F′ₙ3[1,2] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₁3[n])
        F′ₙ3[2,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₂3[n])
        F′ₙ3[2,2] = F′ᵢⱼ
        LaplaceKinematics.advance!(k3, F′ₙ3)
    end

    println("Working on figures for a and a′ = da/dt.")
    a1  = zeros(Float64, N)
    a′1 = zeros(Float64, N)
    a2  = zeros(Float64, N)
    a′2 = zeros(Float64, N)
    a3  = zeros(Float64, N)
    a′3 = zeros(Float64, N)
    for n in 1:N
        a1[n]  = get(k1.aₙ[n])
        a′1[n] = get(k1.a′ₙ[n])
        a2[n]  = get(k2.aₙ[n])
        a′2[n] = get(k2.a′ₙ[n])
        a3[n]  = get(k3.aₙ[n])
        a′3[n] = get(k3.a′ₙ[n])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "elongation a",
        title = "Elongation a and its rate da/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, a1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, a2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, a3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2Da.png"), fig)
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "elongation rate da/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, a′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, a′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, a′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2Dda.png"), fig)

    println("Working on figures for b and db/dt.")
    b1  = zeros(Float64, N)
    b′1 = zeros(Float64, N)
    b2  = zeros(Float64, N)
    b′2 = zeros(Float64, N)
    b3  = zeros(Float64, N)
    b′3 = zeros(Float64, N)
    for n in 1:N
        b1[n]  = get(k1.bₙ[n])
        b′1[n] = get(k1.b′ₙ[n])
        b2[n]  = get(k2.bₙ[n])
        b′2[n] = get(k2.b′ₙ[n])
        b3[n]  = get(k3.bₙ[n])
        b′3[n] = get(k3.b′ₙ[n])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "elongation b",
        title = "Elongation b and its rate db/dt",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, b1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, b2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, b3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2Db.png"), fig)
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "elongation rate db/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, b′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, b′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, b′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2Ddb.png"), fig)

    println("Working on figures for γ and dγ/dt.")
    γ1  = zeros(Float64, N)
    γ′1 = zeros(Float64, N)
    γ2  = zeros(Float64, N)
    γ′2 = zeros(Float64, N)
    γ3  = zeros(Float64, N)
    γ′3 = zeros(Float64, N)
    for n in 1:N
        γ1[n]  = get(k1.γₙ[n])
        γ′1[n] = get(k1.γ′ₙ[n])
        γ2[n]  = get(k2.γₙ[n])
        γ′2[n] = get(k2.γ′ₙ[n])
        γ3[n]  = get(k3.γₙ[n])
        γ′3[n] = get(k3.γ′ₙ[n])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "shear γ",
        title = "In-plane shear γ and its rate dγ/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, γ1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, γ2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, γ3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2Dgamma.png"), fig)
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "shear rate dγ/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, γ′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, γ′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, γ′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rc)
    save(string(myDirPath, "2Ddgamma.png"), fig)

    println("Working on figures for δ and dδ/dt.")
    δ1  = zeros(Float64, N)
    δ′1 = zeros(Float64, N)
    δ2  = zeros(Float64, N)
    δ′2 = zeros(Float64, N)
    δ3  = zeros(Float64, N)
    δ′3 = zeros(Float64, N)
    for n in 1:N
        δ1[n]  = get(k1.δ[n])
        δ′1[n] = get(k1.δ′[n])
        δ2[n]  = get(k2.δ[n])
        δ′2[n] = get(k2.δ′[n])
        δ3[n]  = get(k3.δ[n])
        δ′3[n] = get(k3.δ′[n])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "dilation δ",
        title = "Dilation δ and its rate dδ/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, δ1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, δ2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, δ3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2Ddelta.png"), fig)
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "dilation rate dδ/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, δ′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, δ′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, δ′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2DdDelta.png"), fig)

    println("Working on figures for ϵ and dϵ/dt.")
    ϵ1  = zeros(Float64, N)
    ϵ′1 = zeros(Float64, N)
    ϵ2  = zeros(Float64, N)
    ϵ′2 = zeros(Float64, N)
    ϵ3  = zeros(Float64, N)
    ϵ′3 = zeros(Float64, N)
    for n in 1:N
        ϵ1[n]  = get(k1.ϵ[n])
        ϵ′1[n] = get(k1.ϵ′[n])
        ϵ2[n]  = get(k2.ϵ[n])
        ϵ′2[n] = get(k2.ϵ′[n])
        ϵ3[n]  = get(k3.ϵ[n])
        ϵ′3[n] = get(k3.ϵ′[n])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "squeeze ϵ",
        title = "Squeeze ϵ and its rate dϵ/dt.",
        titlesize = 24,
        xlabelsize = 20,
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
    save(string(myDirPath, "2Depsilon.png"), fig)
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "squeeze rate dϵ/dt (s⁻¹)",
        titlesize = 24,
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
    save(string(myDirPath, "2DdEpsilon.png"), fig)

    println("Working on figures for Laplace γ and dγ/dt.")
    γ1  = zeros(Float64, N)
    γ′1 = zeros(Float64, N)
    γ2  = zeros(Float64, N)
    γ′2 = zeros(Float64, N)
    γ3  = zeros(Float64, N)
    γ′3 = zeros(Float64, N)
    for n in 1:N
        γ1[n]  = get(k1.γ[n])
        γ′1[n] = get(k1.γ′[n])
        γ2[n]  = get(k2.γ[n])
        γ′2[n] = get(k2.γ′[n])
        γ3[n]  = get(k3.γ[n])
        γ′3[n] = get(k3.γ′[n])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "shear strain γ",
        title = "Shear strain γ and its rate dγ/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, γ1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, γ2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, γ3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2DLaplaceGamma.png"), fig)
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "shear strain rate dγ/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, γ′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, γ′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, γ′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rc)
    save(string(myDirPath, "2DdLaplaceGamma.png"), fig)

    println("Working on figures for Gram rotation ω and spin dω/dt.")
    ω1  = zeros(Float64, N)
    ω′1 = zeros(Float64, N)
    ω2  = zeros(Float64, N)
    ω′2 = zeros(Float64, N)
    ω3  = zeros(Float64, N)
    ω′3 = zeros(Float64, N)
    for n in 1:N
        ω1[n]  = (180 / π) * get(k1.ωₙ[n])
        ω′1[n] = (180 / π) * get(k1.ω′ₙ[n])
        ω2[n]  = (180 / π) * get(k2.ωₙ[n])
        ω′2[n] = (180 / π) * get(k2.ω′ₙ[n])
        ω3[n]  = (180 / π) * get(k3.ωₙ[n])
        ω′3[n] = (180 / π) * get(k3.ω′ₙ[n])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title = "Gram Rotation: ω.",
        xlabel = "time (s)",
        ylabel = "Gram rotation ω (°)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, ω1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, ω2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, ω3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(myDirPath, "2Domega.png"), fig)
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title = "Gram Spin: dω/dt.",
        xlabel = "time (s)",
        ylabel = "Gram spin dω / dt (°⋅s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, ω′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, ω′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, ω′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rb)
    save(string(myDirPath, "2DdOmega.png"), fig)

    println("Working on figure of pivoting for frame indifference.")
    P1 = zeros(UInt8, N)
    P2 = zeros(UInt8, N)
    P3 = zeros(UInt8, N)
    for n in 1:N
        P1[n] = k1.P[n]
        P2[n] = k2.P[n]
        P3[n] = k3.P[n]
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "pivot case",
        yticks = [1, 2, 3, 4],
        title = "Pivot Case for Frame Indifference",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, P1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, P2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, P3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rc)
    save(string(myDirPath, "2Dpivot.png"), fig)
    countP1_1 = 0
    countP1_2 = 0
    countP1_3 = 0
    countP1_4 = 0
    countP2_1 = 0
    countP2_2 = 0
    countP2_3 = 0
    countP2_4 = 0
    countP3_1 = 0
    countP3_2 = 0
    countP3_3 = 0
    countP3_4 = 0
    for n in 1:N
        if P1[n] == 1
            countP1_1 += 1
        elseif P1[n] == 2
            countP1_2 += 1
        elseif P1[n] == 3
            countP1_3 += 1
        else
            countP1_4 += 1
        end
        if P2[n] == 1
            countP2_1 += 1
        elseif P2[n] == 2
            countP2_2 += 1
        elseif P2[n] == 3
            countP2_3 += 1
        else
            countP2_4 += 1
        end
        if P3[n] == 1
            countP3_1 += 1
        elseif P3[n] == 2
            countP3_2 += 1
        elseif P3[n] == 3
            countP3_3 += 1
        else
            countP3_4 += 1
        end
    end
    println("At location 1:")
    println("   case 1 occurred ", string(countP1_1), ", 2 occurred ", string(countP1_2),", 3 occurred ", string(countP1_3),", and 4 occurred ", string(countP1_4), " times.")
    println("At location 2:")
    println("   case 1 occurred ", string(countP2_1), ", 2 occurred ", string(countP2_2),", 3 occurred ", string(countP2_3),", and 4 occurred ", string(countP2_4), " times.")
    println("At location 3:")
    println("   case 1 occurred ", string(countP3_1), ", 2 occurred ", string(countP3_2),", 3 occurred ", string(countP3_3),", and 4 occurred ", string(countP3_4), " times.")
end # figures2D

end # testLaplaceKinematics2D.jl