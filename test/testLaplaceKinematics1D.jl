#=
Created on Sat 19 Feb 2022
Updated on Fri 23 Jan 2026
=#

module testLaplaceKinematics1D

using
    FijLung,
    LaplaceKinematics,
    Measures,        # needed to pad white space (margins) around a plot
    PhysicalFields,
    Plots
    
import
    LaplaceKinematics as LK,
    PhysicalFields    as PF

export
    figures1D,
    persistence
#=
--------------------------------------------------------------------------------
=#

"""
```julia
persistence()
```
This function tests writing and reading a *FiberKinematics* object to and from a file for its ability to recreate the object from file.
"""
function persistence()
    my_dir_path = string(pwd(), "/files/")
    if !isdir(my_dir_path)
        mkdir(my_dir_path)
    end
    
    N = 3 # Considered so the JSON file would not be too big.

    array_of_times = FijLung.t_loc1()
    Nₛ = Int(array_of_times.array.len)
    location = 1  # Next to the visceral pleura.
    splineF = FijLung.splineAtEndPoints(location, Nₛ)

    # The stretch from initial to reference state, i.e., κ₀ ↦ κᵣ  .
    Lᵣ = PF.PhysicalScalar(0.98, CGS_LENGTH)
    L₀ = PF.PhysicalScalar(1.00, CGS_LENGTH)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = splineF.t[N+1] - splineF.t[N]
    k = LK.FiberKinematics(N, dt, Lᵣ, L₀)

    # Populate this data structure.
    for n in 2:N+1
        Fₙ = splineF.F[n]
        Lₙ = Fₙ[2,2] * L₀
        LK.advance!(k, Lₙ)
    end

    # Convert this data structure to a JSON stream.
    json_stream = PF.openJSONWriter(my_dir_path, "test1D.json")
    LK.toFile(k, json_stream)
    PF.closeJSONStream(json_stream)

    # Retrieve this data structure from a JSON stream.
    json_stream = PF.openJSONReader(my_dir_path, "test1D.json")
    k1 = LK.fromFile(LK.FiberKinematics, json_stream)
    PF.closeJSONStream(json_stream)

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
    if k.L₀ ≠ k1.L₀
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
        if k.λ″[i] ≠ k1.λ″[i]
            equal = false
        end
        if k.ε[i] ≠ k1.ε[i]
            equal = false
        end
        if k.ε′[i] ≠ k1.ε′[i]
            equal = false
        end
        if k.ε″[i] ≠ k1.ε″[i]
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
figures1D(N)
```
where

    N  is the number of nodes or knots in the B-spline.

This function tests the exported functions of *LaplaceKinematics* for the 1D case; specifically, in the 2 direction, viz., in the direction of the spine.
"""
function figures1D(N::Int)
    my_dir_path = string(pwd(), "/test/figures/")
    if !isdir(my_dir_path)
        mkdir(my_dir_path)
    end

    # set the graphics backend to GR
    ENV["QT_QPA_PLATFORM"] = "wayland"
    gr()

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
        t1[n] = PF.get(splineF1.t[n])
        t2[n] = PF.get(splineF2.t[n])
        t3[n] = PF.get(splineF3.t[n])
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
        λ1[n] = PF.get(Fᵢⱼ1[2,2])
        Fᵢⱼ2 = splineF2.F[n]
        λ2[n] = PF.get(Fᵢⱼ2[2,2])
        Fᵢⱼ3 = splineF3.F[n]
        λ3[n] = PF.get(Fᵢⱼ3[2,2])
    end
    
    println("Building the Laplace kinematics data structures.")
    # The stretch from initial to reference state, i.e., κ₀ ↦ κᵣ.
    Lᵣ = PF.PhysicalScalar(0.98, PF.CGS_LENGTH)
    L₀ = PF.PhysicalScalar(1.00, PF.CGS_LENGTH)
    Lₙ = PF.PhysicalScalar(PF.CGS_LENGTH)

    # Build a data structure for Laplace kinematics at lung location 1.
    dt = PF.PhysicalScalar(t1[N]-t1[N-1], PF.CGS_SECOND)
    k1 = LK.FiberKinematics(N, dt, Lᵣ, L₀)
    # Populate this data structure.
    for n in 2:N
        Lₙ = L₀ * λ1[n]
        LK.advance!(k1, Lₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 2.
    PF.set!(dt, t2[N]-t2[N-1])
    k2 = LK.FiberKinematics(N, dt, Lᵣ, L₀)
    # Populate this data structure.
    for n in 2:N
        Lₙ = L₀ * λ2[n]
        LK.advance!(k2, Lₙ)
    end

    # Build a data structure for Laplace kinematics at lung location 3.
    PF.set!(dt, t3[N]-t3[N-1])
    k3 = LK.FiberKinematics(N, dt, Lᵣ, L₀)
    # Populate this data structure.
    for n in 2:N
        Lₙ = L₀ * λ3[n]
        LK.advance!(k3, Lₙ)
    end

    # Create the figures.
    
    println("Working on figures for ε and dε/dt.")
    ε1 = zeros(Float64, N)
    ε′1 = zeros(Float64, N)
    ε″1 = zeros(Float64, N)
    ε2 = zeros(Float64, N)
    ε′2 = zeros(Float64, N)
    ε″2 = zeros(Float64, N)
    ε3 = zeros(Float64, N)
    ε′3 = zeros(Float64, N)
    ε″3 = zeros(Float64, N)
    for n in 1:N
        ε1[n] = PF.get(k1.ε[n])
        ε′1[n] = PF.get(k1.ε′[n])
        ε″1[n] = PF.get(k1.ε″[n])
        ε2[n] = PF.get(k2.ε[n])
        ε′2[n] = PF.get(k2.ε′[n])
        ε″2[n] = PF.get(k2.ε″[n])
        ε3[n] = PF.get(k3.ε[n])
        ε′3[n] = PF.get(k3.ε′[n])
        ε″3[n] = PF.get(k3.ε″[n])
    end
    
    plot(t1, ε1, label="1: pleural", linecolor=:black, linewidth=3)
    plot!(t2, ε2, label="2: interior", linecolor=:blue, linewidth=3)
    plot!(t3, ε3, label="3: bronchiole", linecolor=:red, linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!("Strain ε and its Rate dε/dt")
    xlabel!("time (s)")
    ylabel!("strain ε")
    
    figName = string("1DStrain.png")
    figPath = string(my_dir_path, figName)
    savefig(figPath)
    
    plot(t1, ε′1, label="1: pleural", linecolor=:black, linewidth=3)
    plot!(t2, ε′2, label="2: interior", linecolor=:blue, linewidth=3)
    plot!(t3, ε′3, label="3: bronchiole", linecolor=:red, linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    xlabel!("time (s)")
    ylabel!("strain rate dε/dt (s⁻¹)")
    
    figName = string("1DStrainRate.png")
    figPath = string(my_dir_path, figName)
    savefig(figPath)
    
    plot(t1, ε″1, label="1: pleural", linecolor=:black, linewidth=3)
    plot!(t2, ε″2, label="2: interior", linecolor=:blue, linewidth=3)
    plot!(t3, ε″3, label="3: bronchiole", linecolor=:red, linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    xlabel!("time (s)")
    ylabel!("rate of strain rate d²ε/dt² (s⁻²)")
    
    figName = string("1DRateOfStrainRate.png")
    figPath = string(my_dir_path, figName)
    savefig(figPath)
    
end # figures1D

end # testLaplaceKinematics1D.jl

