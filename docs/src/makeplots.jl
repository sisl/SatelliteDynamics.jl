using Plots, LaTeXStrings

Plots.gr()

function example_orbit_propagation_kepler(plot_dir)
    # Simulate Keplerian Orbit
    plot_name = "keplerian_orbit"
    println("Generating: $plot_name.svg")

    epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0) 
    oe0  = [R_EARTH + 500e3, 0.01, 75.0, 45.0, 30.0, 0.0]
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # Simulate orbit for one orbit
    T    = orbit_period(oe0[1])
    epcf = epc0 + T
    orb  = EarthInertialState(epc0, eci0, dt=1.0,
            mass=1.0, n_grav=0, m_grav=0,
            drag=false, srp=false,
            moon=false, sun=false,
            relativity=false
    )

    t, epc, eci = sim!(orb, epcf)

    # Create plot with one series
    AX_LIM = R_EARTH + 1000e3
    plt = plot3d(1, xlim=(-AX_LIM,AX_LIM), 
                    ylim=(-AX_LIM,AX_LIM), 
                    zlim=(-AX_LIM, AX_LIM), 
                    title = "Keplerian Orbit",
                    xformatter=:scientific,
                    yformatter=:scientific,
                    zformatter=:scientific,
                    linewidth=2,
                    legend=false)

    for i in 1:length(t)
        push!(plt, eci[1, i], eci[2, i], eci[3, i])
    end
    Plots.savefig(plot_dir * "/$plot_name.svg")

    # Plot Keplerian Elements
    plot_name = "keplerian_elements"
    println("Generating: $plot_name.svg")

    # Compute orbital Elements
    elements = zeros(Float64, 6, length(t))
    plots = []
    for (i, e) in enumerate(epc)
        elements[:, i] = sCARTtoOSC(eci[:, i], use_degrees=true)
    end

    ylabels = [L"a \; [km]", L"e", L"i \; [deg]", L"\Omega \; [deg]", L"\omega \; [deg]", L"M \; [deg]"]
    for i in 1:6
        y = elements[i, 1]
        if i == 6
            push!(plots, plot(t/T, elements[i, :],
                        xlim=(0.0, 1.0),
                        ylim=(0.0, 360.0), 
                        ylabel=ylabels[i],
                        legend=false, 
                        linewidth=2))
        else
            push!(plots, plot(t/T, elements[i, :],
                        xlim=(0.0, 1.0),
                        ylim=(y-0.1*y, y+0.1*y), 
                        ylabel=ylabels[i], 
                        legend=false,
                        linewidth=2))
        end
    end

    plot(plots..., layout=(3,2), xlabel=L"T \; [Orbits]")

    Plots.savefig(plot_dir * "/$plot_name.svg")
end

function example_orbit_propagation_fullforce(plot_dir)
    # Simulate Keplerian Orbit
    plot_name = "fullforce_orbit"
    println("Generating: $plot_name.svg")

    epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0) 
    oe0  = [R_EARTH + 500e3, 0.01, 75.0, 45.0, 30.0, 0.0]
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # Simulate orbit for one orbit
    T    = orbit_period(oe0[1])
    epcf = epc0 + T
    orb  = EarthInertialState(epc0, eci0, dt=1.0,
            mass=100.0, n_grav=20, m_grav=20,
            drag=true, srp=true,
            moon=true, sun=true,
            relativity=true
    )

    t, epc, eci = sim!(orb, epcf)

    # Create plot with one series
    AX_LIM = R_EARTH + 1000e3
    plt = plot3d(1, xlim=(-AX_LIM,AX_LIM), 
                    ylim=(-AX_LIM,AX_LIM), 
                    zlim=(-AX_LIM, AX_LIM), 
                    title = "Keplerian Orbit",
                    xformatter=:scientific,
                    yformatter=:scientific,
                    zformatter=:scientific,
                    linewidth=2,
                    legend=false)

    for i in 1:length(t)
        push!(plt, eci[1, i], eci[2, i], eci[3, i])
    end
    Plots.savefig(plot_dir * "/$plot_name.svg")

    # Plot Keplerian Elements
    plot_name = "fullforce_elements"
    println("Generating: $plot_name.svg")

    # Compute orbital Elements
    elements = zeros(Float64, 6, length(t))
    plots = []
    for (i, e) in enumerate(epc)
        elements[:, i] = sCARTtoOSC(eci[:, i], use_degrees=true)
    end

    # ylabels = [L"a \; [km]", L"e", L"i \; [deg]", L"\Omega \; [deg]", L"\omega \; [deg]", L"M \; [deg]"]
    ylabels = ["a [km]", "e", "i [deg]", "\\Omega [deg]", "\\omega [deg]", "M [deg]"]
    for i in 1:6
        y_min = minimum(elements[i, :])
        y_max = maximum(elements[i, :])
        if i == 6
            push!(plots, plot(t/T, elements[i, :],
                        xlim=(0.0, 1.0),
                        ylim=(0.0, 360.0), 
                        ylabel=ylabels[i],
                        legend=false, 
                        linewidth=2))
        else
            push!(plots, plot(t/T, elements[i, :],
                        xlim=(0.0, 1.0),
                        ylim=(y_min-0.01*y_min, y_max+0.01*y_max), 
                        ylabel=ylabels[i], 
                        legend=false,
                        linewidth=2))
        end
    end

    # plot(plots..., layout=(3,2), xlabel=L"T \; [Orbits]")
    plot(plots..., layout=(3,2), xlabel="T [Orbits]")

    Plots.savefig(plot_dir * "/$plot_name.svg")
end

# Generate all plots
function makeplots()
    # Start building plots
    println("Generating plots")

    # Make output directory
    plot_dir = (pwd()[end-3:end] == "docs") ? "build/plots" : "docs/build/plots"
    if isdir(plot_dir)
        rm(plot_dir, recursive=true)
    end
    mkdir(plot_dir)

    example_orbit_propagation_kepler(plot_dir)
    example_orbit_propagation_fullforce(plot_dir)
end