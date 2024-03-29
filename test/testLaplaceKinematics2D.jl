module testLaplaceKinematics2D

using
    CairoMakie,      # Pixel based figure construction.
    FijLung,
    LaplaceKinematics,
    PhysicalFields

export
    figures2D,
    persistence
#=
-------------------------------------------------------------------------------
=#

"""
Function:\n
    persistence(midPtQuad, myDirPath)\n
where\n
    midPtQuad  is true if a mid-point quadrature is to be used.\n
    myDirPath  is the path where the user wants his file to be written.\n
This function tests writing and reading a `MembraneKinematics` object to and from a file for its ability to recreate the object from file.
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
    aᵣ = PhysicalScalar(1.0, CGS_STRETCH)
    bᵣ = PhysicalScalar(1.0, CGS_STRETCH)
    γᵣ = PhysicalScalar(CGS_STRETCH)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = splineF.t[N] - splineF.t[N-1]
    Pᵣ = 1
    k = LaplaceKinematics.MembraneKinematics(dt, N, midPtQuad, aᵣ, bᵣ, γᵣ, Pᵣ)

    # Populate this data structure.
    for n in 1:N
        F′ = splineF.F′[n+1]
        F′ₙ = PhysicalTensor(2, 2, CGS_STRETCH_RATE)
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
    for i in 1:N+1
        if k.t[i] ≠ k1.t[i]
            equal = false
        end
        if k.F[i] ≠ k1.F[i]
            equal = false
        end
        if k.F′[i] ≠ k1.F′[i]
            equal = false
        end
        if k.motion[i] ≠ k1.motion[i]
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
    figures2D(N, midPtQuad, myDirPath)\n
where\n
    N          is the number of nodes or knots in the B-spline.\n
    midPtQuad  is true if a mid-point quadrature is to be used.\n
    myDirPath  is the path where the user wants his file to be written.\n
This function tests the exported functions of `LaplaceKinematics` for the 2D case; specifically, for the 23 plane for Fᵢⱼ exported by FijLung.
"""
function figures2D(N::Integer, midPtQuad::Bool, myDirPath::String)

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
    # Consider reference and initial states to be the same, i.e., κᵣ = κ₁.
    aᵣ = PhysicalScalar(1.0, CGS_STRETCH)
    bᵣ = PhysicalScalar(1.0, CGS_STRETCH)
    γᵣ = PhysicalScalar(CGS_STRETCH)
    # This shearing is in the 1 direction.
    Pᵣ = 1
    # Create the variables to hold the rate of deformation gradient.
    F′ᵢⱼ = PhysicalScalar(CGS_STRETCH_RATE)
    F′ₙ  = PhysicalTensor(2, 2, CGS_STRETCH_RATE)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = PhysicalScalar(t1[N]-t1[N-1], CGS_SECOND)
    k1 = LaplaceKinematics.MembraneKinematics(dt, N, midPtQuad, aᵣ, bᵣ, γᵣ, Pᵣ)
    # Populate this data structure.
    for n in 2:N
        set!(F′ᵢⱼ, F′₁₁1[n])
        F′ₙ[1,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₁₂1[n])
        F′ₙ[1,2] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₁1[n])
        F′ₙ[2,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₂1[n])
        F′ₙ[2,2] = F′ᵢⱼ
        LaplaceKinematics.advance!(k1, F′ₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 2.
    set!(dt, t2[N]-t2[N-1])
    k2 = LaplaceKinematics.MembraneKinematics(dt, N, midPtQuad, aᵣ, bᵣ, γᵣ, Pᵣ)
    # Populate this data structure.
    for n in 2:N
        set!(F′ᵢⱼ, F′₁₁2[n])
        F′ₙ[1,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₁₂2[n])
        F′ₙ[1,2] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₁2[n])
        F′ₙ[2,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₂2[n])
        F′ₙ[2,2] = F′ᵢⱼ
        LaplaceKinematics.advance!(k2, F′ₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 3.
    set!(dt, t3[N]-t3[N-1])
    k3 = LaplaceKinematics.MembraneKinematics(dt, N, midPtQuad, aᵣ, bᵣ, γᵣ, Pᵣ)
    # Populate this data structure.
    for n in 2:N
        set!(F′ᵢⱼ, F′₁₁3[n])
        F′ₙ[1,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₁₂3[n])
        F′ₙ[1,2] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₁3[n])
        F′ₙ[2,1] = F′ᵢⱼ
        set!(F′ᵢⱼ, F′₂₂3[n])
        F′ₙ[2,2] = F′ᵢⱼ
        LaplaceKinematics.advance!(k3, F′ₙ)
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
    fig1 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax1 = Axis(fig1[1, 1];
        xlabel = "time (s)",
        ylabel = "elongation a",
        title = "Elongation a and its rate da/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax1, t1, a1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax1, t2, a2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax1, t3, a3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DaAtMidPoints.png"), fig1)
    else
        save(string(myDirPath, "2DaAtEndPoints.png"), fig1)
    end

    fig2 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax2 = Axis(fig2[1, 1];
        xlabel = "time (s)",
        ylabel = "elongation rate da/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax2, t1, a′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax2, t2, a′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax2, t3, a′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DdaAtMidPoints.png"), fig2)
    else
        save(string(myDirPath, "2DdaAtEndPoints.png"), fig2)
    end

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
    fig3 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax3 = Axis(fig3[1, 1];
        xlabel = "time (s)",
        ylabel = "elongation b",
        title = "Elongation b and its rate db/dt",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax3, t1, b1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax3, t2, b2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax3, t3, b3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DbAtMidPoints.png"), fig3)
    else
        save(string(myDirPath, "2DbAtEndPoints.png"), fig3)
    end

    fig4 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax4 = Axis(fig4[1, 1];
        xlabel = "time (s)",
        ylabel = "elongation rate db/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax4, t1, b′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax4, t2, b′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax4, t3, b′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DdbAtMidPoints.png"), fig4)
    else
        save(string(myDirPath, "2DdbAtEndPoints.png"), fig4)
    end

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
    fig5 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax5 = Axis(fig5[1, 1];
        xlabel = "time (s)",
        ylabel = "shear γ",
        title = "In-plane shear γ and its rate dγ/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax5, t1, γ1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax5, t2, γ2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax5, t3, γ3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DgammaAtMidPoints.png"), fig5)
    else
        save(string(myDirPath, "2DgammaAtEndPoints.png"), fig5)
    end

    fig6 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax6 = Axis(fig6[1, 1];
        xlabel = "time (s)",
        ylabel = "shear rate dγ/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax6, t1, γ′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax6, t2, γ′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax6, t3, γ′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rc)
    if midPtQuad
        save(string(myDirPath, "2DdgammaAtMidPoints.png"), fig6)
    else
        save(string(myDirPath, "2DdgammaAtEndPoints.png"), fig6)
    end

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
    fig7 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax7 = Axis(fig7[1, 1];
        xlabel = "time (s)",
        ylabel = "dilation δ",
        title = "Dilation δ and its rate dδ/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax7, t1, δ1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax7, t2, δ2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax7, t3, δ3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DdeltaAtMidPoints.png"), fig7)
    else
        save(string(myDirPath, "2DdeltaAtEndPoints.png"), fig7)
    end

    fig8 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax8 = Axis(fig8[1, 1];
        xlabel = "time (s)",
        ylabel = "dilation rate dδ/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax8, t1, δ′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax8, t2, δ′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax8, t3, δ′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DdDeltaAtMidPoints.png"), fig8)
    else
        save(string(myDirPath, "2DdDeltaAtEndPoints.png"), fig8)
    end

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

    fig9 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax9 = Axis(fig9[1, 1];
        xlabel = "time (s)",
        ylabel = "squeeze ϵ",
        title = "Squeeze ϵ and its rate dϵ/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax9, t1, ϵ1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax9, t2, ϵ2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax9, t3, ϵ3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DepsilonAtMidPoints.png"), fig9)
    else
        save(string(myDirPath, "2DepsilonAtEndPoints.png"), fig9)
    end

    fig10 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax10 = Axis(fig10[1, 1];
        xlabel = "time (s)",
        ylabel = "squeeze rate dϵ/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax10, t1, ϵ′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax10, t2, ϵ′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax10, t3, ϵ′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DdEpsilonAtMidPoints.png"), fig10)
    else
        save(string(myDirPath, "2DdEpsilonAtEndPoints.png"), fig10)
    end

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

    fig11 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax11 = Axis(fig11[1, 1];
        xlabel = "time (s)",
        ylabel = "shear strain γ",
        title = "Shear strain γ and its rate dγ/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax11, t1, γ1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax11, t2, γ2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax11, t3, γ3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DLaplaceGammaAtMidPoints.png"), fig11)
    else
        save(string(myDirPath, "2DLaplaceGammaAtEndPoints.png"), fig11)
    end

    fig12 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax12 = Axis(fig12[1, 1];
        xlabel = "time (s)",
        ylabel = "shear strain rate dγ/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax12, t1, γ′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax12, t2, γ′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax12, t3, γ′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rc)
    if midPtQuad
        save(string(myDirPath, "2DdLaplaceGammaAtMidPoints.png"), fig12)
    else
        save(string(myDirPath, "2DdLaplaceGammaAtEndPoints.png"), fig12)
    end

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

    fig13 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax13 = Axis(fig13[1, 1];
        title = "Gram Rotation: ω.",
        xlabel = "time (s)",
        ylabel = "Gram rotation ω (°)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax13, t1, ω1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax13, t2, ω2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax13, t3, ω3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    if midPtQuad
        save(string(myDirPath, "2DomegaAtMidPoints.png"), fig13)
    else
        save(string(myDirPath, "2DomegaAtEndPoints.png"), fig13)
    end

    fig14 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax14 = Axis(fig14[1, 1];
        title = "Gram Spin: dω/dt.",
        xlabel = "time (s)",
        ylabel = "Gram spin dω / dt (°⋅s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax14, t1, ω′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax14, t2, ω′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax14, t3, ω′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rb)
    if midPtQuad
        save(string(myDirPath, "2DdOmegaAtMidPoints.png"), fig14)
    else
        save(string(myDirPath, "2DdOmegaAtEndPoints.png"), fig14)
    end

    println("Working on figure of pivoted motions for frame indifference.")
    motion1 = zeros(UInt8, N)
    motion2 = zeros(UInt8, N)
    motion3 = zeros(UInt8, N)
    for n in 1:N
        motion1[n] = k1.motion[n]
        motion2[n] = k2.motion[n]
        motion3[n] = k3.motion[n]
    end

    fig15 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax15 = Axis(fig15[1, 1];
        xlabel = "time (s)",
        ylabel = "motion case",
        yticks = [1, 2, 3, 4],
        title = "Motion Case for Frame Indifference",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax15, t1, motion1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax15, t2, motion2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax15, t3, motion3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rc)
    if midPtQuad
        save(string(myDirPath, "2DmotionAtMidPoints.png"), fig15)
    else
        save(string(myDirPath, "2DmotionAtEndPoints.png"), fig15)
    end

    countMotion1_1 = 0
    countMotion1_2 = 0
    countMotion1_3 = 0
    countMotion1_4 = 0
    countMotion2_1 = 0
    countMotion2_2 = 0
    countMotion2_3 = 0
    countMotion2_4 = 0
    countMotion3_1 = 0
    countMotion3_2 = 0
    countMotion3_3 = 0
    countMotion3_4 = 0
    for n in 1:N
        if motion1[n] == 1
            countMotion1_1 += 1
        elseif motion1[n] == 2
            countMotion1_2 += 1
        elseif motion1[n] == 3
            countMotion1_3 += 1
        else
            countMotion1_4 += 1
        end
        if motion2[n] == 1
            countMotion2_1 += 1
        elseif motion2[n] == 2
            countMotion2_2 += 1
        elseif motion2[n] == 3
            countMotion2_3 += 1
        else
            countMotion2_4 += 1
        end
        if motion3[n] == 1
            countMotion3_1 += 1
        elseif motion3[n] == 2
            countMotion3_2 += 1
        elseif motion3[n] == 3
            countMotion3_3 += 1
        else
            countMotion3_4 += 1
        end
    end
    if midPtQuad
        println("Assigning nodes to the mid-points of each solution interval:")
    else
        println("Assigning nodes to the end-points of each solution interval:")
    end
    println("At location 1:")
    println("   motion 1 occurred ", string(countMotion1_1), ", 2 occurred ", string(countMotion1_2),", 3 occurred ", string(countMotion1_3),", and 4 occurred ", string(countMotion1_4), " times.")
    println("At location 2:")
    println("   motion 1 occurred ", string(countMotion2_1), ", 2 occurred ", string(countMotion2_2),", 3 occurred ", string(countMotion2_3),", and 4 occurred ", string(countMotion2_4), " times.")
    println("At location 3:")
    println("   motion 1 occurred ", string(countMotion3_1), ", 2 occurred ", string(countMotion3_2),", 3 occurred ", string(countMotion3_3),", and 4 occurred ", string(countMotion3_4), " times.")
end # figures2D

end # testLaplaceKinematics2D.jl