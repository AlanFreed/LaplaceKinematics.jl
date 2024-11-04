module testLaplaceKinematics2D

using
    CairoMakie,      # Pixel based figure construction.
    FijLung,
    LaplaceKinematics,
    PhysicalFields

import
    LaplaceKinematics as LK,
    PhysicalFields    as PF

export
    figures2D,
    persistence
#=
-------------------------------------------------------------------------------
=#

"""
```julia
persistence()
```
This function tests writing and reading a *MembraneKinematics* object to and from a file for its ability to recreate the object from file.
"""
function persistence()
    my_dir_path = string(pwd(), "/test/files/")
    if !isdir(my_dir_path)
        mkdir(my_dir_path)
    end
    
    N = 3 # Considered so the JSON file would not be too big.

    array_of_times = FijLung.t_loc1()
    Nₛ = Int(array_of_times.array.len)
    location = 1  # Next to the visceral pleura.
    splineF = FijLung.splineAtEndPoints(location, Nₛ)

    # The motion from the initial to reference state, i.e., κ₀ ↦ κᵣ.
    aᵣ = PF.PhysicalScalar(0.95,  PF.CGS_DIMENSIONLESS)
    bᵣ = PF.PhysicalScalar(0.95,  PF.CGS_DIMENSIONLESS)
    γᵣ = PF.PhysicalScalar(-0.05, PF.CGS_DIMENSIONLESS)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = splineF.t[N] - splineF.t[N-1]
    Pᵣ = 1
    k = LK.MembraneKinematics(dt, N, aᵣ, bᵣ, γᵣ, Pᵣ)

    # Populate this data structure.
    for n in 1:N
        dF  = splineF.F[n+1] - splineF.F[n]
        dFₙ = PF.PhysicalTensor(2, 2, PF.CGS_DIMENSIONLESS)
        dFₙ[1,1] = dF[1,1]
        dFₙ[1,2] = dF[1,2]
        dFₙ[2,1] = dF[2,1]
        dFₙ[2,2] = dF[2,2]
        LK.advance!(k, dFₙ)
    end

    # Convert this data structure to a JSON stream.
    json_stream = PF.openJSONWriter(my_dir_path, "test2D.json")
    LK.toFile(k, json_stream)
    PF.closeJSONStream(json_stream)

    # Retrieve this data structure from a JSON stream.
    json_stream = PF.openJSONReader(my_dir_path, "test2D.json")
    k1 = LK.fromFile(LK.MembraneKinematics, json_stream)
    PF.closeJSONStream(json_stream)

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
        if k.ε[i] ≠ k1.ε[i]
            equal = false
        end
        if k.γ[i] ≠ k1.γ[i]
            equal = false
        end
        if k.δ′[i] ≠ k1.δ′[i]
            equal = false
        end
        if k.ε′[i] ≠ k1.ε′[i]
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
```julia
figures2D(N)
```
where

    N  is the number of nodes or knots in the B-spline.

This function tests the exported functions of *LaplaceKinematics* for the 2D case; specifically, for the 23 plane for Fᵢⱼ exported by FijLung.
"""
function figures2D(N::Int)
    my_dir_path = string(pwd(), "/test/figures/")
    if !isdir(my_dir_path)
        mkdir(my_dir_path)
    end

    CairoMakie.activate!(type = "png")
    println("For these figures, the number of intervals is N = ", string(N), ".")

    println("Creating B-splines for the deformation gradient data.")
    splineF1 = FijLung.splineAtEndPoints(1, N) # near visceral pleura
    splineF2 = FijLung.splineAtEndPoints(2, N) # deep in the lung
    splineF3 = FijLung.splineAtEndPoints(3, N) # near bronchiole tube

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
        F₁₁1[n] = PF.get(Fᵢⱼ1[2,2])
        Fᵢⱼ2 = splineF2.F[n]
        F₁₁2[n] = PF.get(Fᵢⱼ2[2,2])
        Fᵢⱼ3 = splineF3.F[n]
        F₁₁3[n] = PF.get(Fᵢⱼ3[2,2])
    end
    # Create the data arrays for dF₁₁.
    dF₁₁1 = zeros(Float64, N)
    dF₁₁2 = zeros(Float64, N)
    dF₁₁3 = zeros(Float64, N)
    for n in 2:N
        dFᵢⱼ1 = splineF1.F[n] - splineF1.F[n-1]
        dF₁₁1[n] = PF.get(dFᵢⱼ1[2,2])
        dFᵢⱼ2 = splineF2.F[n] - splineF2.F[n-1]
        dF₁₁2[n] = PF.get(dFᵢⱼ2[2,2])
        dFᵢⱼ3 = splineF3.F[n] - splineF3.F[n-1]
        dF₁₁3[n] = PF.get(dFᵢⱼ3[2,2])
    end
    # Create the data arrays for F₁₂
    F₁₂1 = zeros(Float64, N)
    F₁₂2 = zeros(Float64, N)
    F₁₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₁₂1[n] = PF.get(Fᵢⱼ1[2,3])
        Fᵢⱼ2 = splineF2.F[n]
        F₁₂2[n] = PF.get(Fᵢⱼ2[2,3])
        Fᵢⱼ3 = splineF3.F[n]
        F₁₂3[n] = PF.get(Fᵢⱼ3[2,3])
    end
    # Create the data arrays dF₁₂.
    dF₁₂1 = zeros(Float64, N)
    dF₁₂2 = zeros(Float64, N)
    dF₁₂3 = zeros(Float64, N)
    for n in 2:N
        dFᵢⱼ1 = splineF1.F[n] - splineF1.F[n-1]
        dF₁₂1[n] = PF.get(dFᵢⱼ1[2,3])
        dFᵢⱼ2 = splineF2.F[n] - splineF2.F[n-1]
        dF₁₂2[n] = PF.get(dFᵢⱼ2[2,3])
        dFᵢⱼ3 = splineF3.F[n] - splineF3.F[n-1]
        dF₁₂3[n] = PF.get(dFᵢⱼ3[2,3])
    end
    # Create the data arrays for F₂₁.
    F₂₁1 = zeros(Float64, N)
    F₂₁2 = zeros(Float64, N)
    F₂₁3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₂₁1[n] = PF.get(Fᵢⱼ1[3,2])
        Fᵢⱼ2 = splineF2.F[n]
        F₂₁2[n] = PF.get(Fᵢⱼ2[3,2])
        Fᵢⱼ3 = splineF3.F[n]
        F₂₁3[n] = PF.get(Fᵢⱼ3[3,2])
    end
    # Create the data arrays for dF₂₁.
    dF₂₁1 = zeros(Float64, N)
    dF₂₁2 = zeros(Float64, N)
    dF₂₁3 = zeros(Float64, N)
    for n in 2:N
        dFᵢⱼ1 = splineF1.F[n] - splineF1.F[n-1]
        dF₂₁1[n] = PF.get(dFᵢⱼ1[3,2])
        dFᵢⱼ2 = splineF2.F[n] - splineF2.F[n-1]
        dF₂₁2[n] = PF.get(dFᵢⱼ2[3,2])
        dFᵢⱼ3 = splineF3.F[n] - splineF3.F[n-1]
        dF₂₁3[n] = PF.get(dFᵢⱼ3[3,2])
    end
    # Create the data arrays for F₂₂.
    F₂₂1 = zeros(Float64, N)
    F₂₂2 = zeros(Float64, N)
    F₂₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₂₂1[n] = PF.get(Fᵢⱼ1[3,3])
        Fᵢⱼ2 = splineF2.F[n]
        F₂₂2[n] = PF.get(Fᵢⱼ2[3,3])
        Fᵢⱼ3 = splineF3.F[n]
        F₂₂3[n] = PF.get(Fᵢⱼ3[3,3])
    end
    # Create the data arrays for dF₂₂.
    dF₂₂1 = zeros(Float64, N)
    dF₂₂2 = zeros(Float64, N)
    dF₂₂3 = zeros(Float64, N)
    for n in 2:N
        dFᵢⱼ1 = splineF1.F[n] - splineF1.F[n-1]
        dF₂₂1[n] = PF.get(dFᵢⱼ1[3,3])
        dFᵢⱼ2 = splineF2.F[n] - splineF2.F[n-1]
        dF₂₂2[n] = PF.get(dFᵢⱼ2[3,3])
        dFᵢⱼ3 = splineF3.F[n] - splineF3.F[n-1]
        dF₂₂3[n] = PF.get(dFᵢⱼ3[3,3])
    end

    println("Building the Laplace kinematics data structures.")
    # The motion from the initial to reference state, i.e., κ₀ ↦ κᵣ.
    aᵣ = PF.PhysicalScalar(0.95,  PF.CGS_DIMENSIONLESS)
    bᵣ = PF.PhysicalScalar(0.95,  PF.CGS_DIMENSIONLESS)
    γᵣ = PF.PhysicalScalar(-0.05, PF.CGS_DIMENSIONLESS)
    # This shearing is in the 1 direction.
    Pᵣ = 1
    
    s = "The Laplace stretch attributes for motion κ₀ ↦ κᵣ are:\n"
    s = string(s, "   aᵣ = ", PF.toString(aᵣ), ",\n")
    s = string(s, "   bᵣ = ", PF.toString(bᵣ), ",\n")
    s = string(s, "   γᵣ = ", PF.toString(γᵣ), " in the ", Pᵣ, "-direction.")
    println(s)
    
    # Create the variable to hold increments of the deformation gradient.
    dFᵢⱼ = PF.PhysicalScalar(PF.CGS_DIMENSIONLESS)
    dFₙ  = PF.PhysicalTensor(2, 2, PF.CGS_DIMENSIONLESS)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = PF.PhysicalScalar(t1[N]-t1[N-1], PF.CGS_SECOND)
    k1 = LK.MembraneKinematics(dt, N, aᵣ, bᵣ, γᵣ, Pᵣ)
    # Populate this data structure.
    for n in 2:N
        PF.set!(dFᵢⱼ, dF₁₁1[n])
        dFₙ[1,1] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₁₂1[n])
        dFₙ[1,2] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₂₁1[n])
        dFₙ[2,1] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₂₂1[n])
        dFₙ[2,2] = dFᵢⱼ
        LK.advance!(k1, dFₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 2.
    PF.set!(dt, t2[N]-t2[N-1])
    k2 = LK.MembraneKinematics(dt, N, aᵣ, bᵣ, γᵣ, Pᵣ)
    # Populate this data structure.
    for n in 2:N
        PF.set!(dFᵢⱼ, dF₁₁2[n])
        dFₙ[1,1] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₁₂2[n])
        dFₙ[1,2] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₂₁2[n])
        dFₙ[2,1] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₂₂2[n])
        dFₙ[2,2] = dFᵢⱼ
        LK.advance!(k2, dFₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 3.
    PF.set!(dt, t3[N]-t3[N-1])
    k3 = LK.MembraneKinematics(dt, N, aᵣ, bᵣ, γᵣ, Pᵣ)
    # Populate this data structure.
    for n in 2:N
        PF.set!(dFᵢⱼ, dF₁₁3[n])
        dFₙ[1,1] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₁₂3[n])
        dFₙ[1,2] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₂₁3[n])
        dFₙ[2,1] = dFᵢⱼ
        PF.set!(dFᵢⱼ, dF₂₂3[n])
        dFₙ[2,2] = dFᵢⱼ
        LK.advance!(k3, dFₙ)
    end

    println("Working on figures for a and da/dt.")
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
    save(string(my_dir_path, "LaplaceAttribute_a_2D.png"), fig1)

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
    save(string(my_dir_path, "LaplaceAttribute_da_2D.png"), fig2)

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
    save(string(my_dir_path, "LaplaceAttribute_b_2D.png"), fig3)
    
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
    save(string(my_dir_path, "LaplaceAttribute_db_2D.png"), fig4)

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
    save(string(my_dir_path, "LaplaceAttribute_γ_2D.png"), fig5)

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
    save(string(my_dir_path, "LaplaceAttribute_dγ_2D.png"), fig6)

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
    save(string(my_dir_path, "LaplaceStrain_δ_2D.png"), fig7)

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
    save(string(my_dir_path, "LaplaceStrain_dδ_2D.png"), fig8)

    println("Working on figures for ε and dε/dt.")
    ε1  = zeros(Float64, N)
    ε′1 = zeros(Float64, N)
    ε2  = zeros(Float64, N)
    ε′2 = zeros(Float64, N)
    ε3  = zeros(Float64, N)
    ε′3 = zeros(Float64, N)
    for n in 1:N
        ε1[n]  = get(k1.ε[n])
        ε′1[n] = get(k1.ε′[n])
        ε2[n]  = get(k2.ε[n])
        ε′2[n] = get(k2.ε′[n])
        ε3[n]  = get(k3.ε[n])
        ε′3[n] = get(k3.ε′[n])
    end

    fig9 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax9 = Axis(fig9[1, 1];
        xlabel = "time (s)",
        ylabel = "squeeze ε",
        title = "Squeeze ε and its rate dε/dt.",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax9, t1, ε1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax9, t2, ε2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax9, t3, ε3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(my_dir_path, "LaplaceStrain_ε_2D.png"), fig9)
    
    fig10 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax10 = Axis(fig10[1, 1];
        xlabel = "time (s)",
        ylabel = "squeeze rate dε/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax10, t1, ε′1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax10, t2, ε′2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax10, t3, ε′3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save(string(my_dir_path, "LaplaceStrain_dε_2D.png"), fig10)
    
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
    save(string(my_dir_path, "LaplaceStrain_γ_2D.png"), fig11)

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
    save(string(my_dir_path, "LaplaceStrain_dγ_2D.png"), fig12)

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
    save(string(my_dir_path, "GramRotation_ω_2D.png"), fig13)

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
    save(string(my_dir_path, "GramRotation_dω_2D.png"), fig14)

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
        title = "Motion Case for Frame Indifference.",
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
    save(string(my_dir_path, "motion_2D.png"), fig15)

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

    println()
    println("At location 1:")
    println("   motion 1 occurred ", string(countMotion1_1), ", 2 occurred ", string(countMotion1_2),", 3 occurred ", string(countMotion1_3),", and 4 occurred ", string(countMotion1_4), " times.")
    println("At location 2:")
    println("   motion 1 occurred ", string(countMotion2_1), ", 2 occurred ", string(countMotion2_2),", 3 occurred ", string(countMotion2_3),", and 4 occurred ", string(countMotion2_4), " times.")
    println("At location 3:")
    println("   motion 1 occurred ", string(countMotion3_1), ", 2 occurred ", string(countMotion3_2),", 3 occurred ", string(countMotion3_3),", and 4 occurred ", string(countMotion3_4), " times.")
end # figures2D

end # testLaplaceKinematics2D.jl
