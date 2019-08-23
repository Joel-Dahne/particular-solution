"""
Creates two plots, one showing the estimated convergence of the
eigenvalues and one showing the convergence of the radius of the
computed enclosures.
"""
function plotconvergence(enclosure::ArbReal, g::Geometry, opt::Options)
    enclosures = getindex.(particularsolution(enclosure, g, opt), 2)

    Ns = opt.N_beg:opt.N_step:opt.N_end
    νs = midpoint.(enclosures)

    # Plot the estimated convergence
    p1 = plot(Ns[1:end-1], abs.(νs[1:end-1] .- νs[end]),
              xaxis = ("N", Ns),
              yaxis = ("Error", :log10),
              marker = :circle,
              legend = :none)

    # Plot the convergence of the radius
    p2 = plot(Ns, radius.(enclosures),
              xaxis = ("N", Ns),
              yaxis = ("Width", :log10),
              marker = :circle,
              legend = :none)

    return (p1, p2)
end

function plotconvergence(::Type{T}, triangle::Int) where T <: Geometry
    (g, enclosure, opt) = getdomain(T, triangle)

    return plotconvergence(enclosure, g, opt)
end

function plotconvergence(::Type{T}, triangles, Ns) where T <: Geometry
    p1 = plot()
    p2 = plot()

    for triangle in triangles
        (g, enclosure, opt) = getdomain(T, triangle)
        opt = Options(opt, Ns)

        enclosures = getindex.(particularsolution(enclosure, g, opt), 2)

        νs = midpoint.(enclosures)

        plot!(p1, Ns[1:end-1], abs.(νs[1:end-1] .- νs[end]),
              xaxis = ("N", Ns),
              yaxis = ("Error", :log10),
              marker = :circle,
              label = "Triangle $triangle")

        plot!(p2, Ns, radius.(enclosures),
              xaxis = ("N", Ns),
              yaxis = ("Width", :log10),
              marker = :circle,
              label = "Triangle $triangle")
    end

    return (p1, p2)
end

"""
Plot the σ-value for different values of ν.
"""
function plotsigma(interval::Tuple{ArbReal, ArbReal}, g::Geometry, p::Points,
                   opt::Options, N::Int)
    (start, stop) = interval
    νs = Array{ArbReal}(undef, opt.plot_n)
    σs = Array{ArbReal}(undef, opt.plot_n)

    @showprogress for i in 1:opt.plot_n
        νs[i] = start + (i - 1)*(stop - start)/(opt.plot_n - 1)
        σs[i] = sigma(νs[i], g, p, N)
    end

    p1 = plot(νs, σs,
              xaxis = L"\nu",
              yaxis = L"\sigma(\nu)",
              legend = :none)

    return p1
end

function plotsigma(::Type{T}, triangle::Int, N::Int, numpoints::Int = 0) where T <: Geometry
    (g, enclosure, opt) = getdomain(T, triangle)
    p = Points(g, 2N, 2N)

    if numpoints > 0
        opt.plot_n = numpoints
    end

    plotsigma(interval(enclosure), g, p, opt, N)
end

"""
Plot the eigenfunction on the boundary.
"""
function ploteigenfunctionboundary(ν::ArbReal, coeffs::Array{ArbReal},
                                   g::Geometry, numpoints::Int, vertex::Int = 1)
    points = boundary(Points(g, numpoints, 0))
    values = map(p -> eigenfunction(p, ν, coeffs, g, vertex), points)
    ϕs = map(p -> p.ϕ, points)

    p1 = plot(ϕs, values,
              xaxis = (L"\phi"),
              yaxis = ("Eigenfunction"),
              legend = :none)

    return p1
end

"""
Plot the eigenfunction in the interior of the spherical triangle. This
plot is a bit hard to get good. We give two different versions, the
first is 3d scatter plot in Cartesian coordinates whereas the second
is a 2d scatter plot in spherical coordinates. The precise options
depends a lot on which plotting backend is used, currently it's
optimized for pyplot which seems to handle the 3d color plots the
best.
"""
function ploteigenfunctioninterior(ν::ArbReal, coeffs::Array{ArbReal},
                                   g::Geometry, numpoints::Int = 1000,
                                   vertex::Int = 1)
    p = Points(g, 0, numpoints)
    setinterior2!(p, g)
    ps = interior(p)
    points = Array{SVector{3}{ArbReal}}(undef, numpoints)
    values = zeros(ArbReal, numpoints)
    @showprogress for i in 1:numpoints
        points[i] = CartesianFromSpherical()(ps[i])
        values[i] = eigenfunction(ps[i], ν, coeffs, g, vertex)
    end

    p1 = scatter3d(getindex.(points, 1), getindex.(points, 2),
                   getindex.(points, 3), zcolor = Float64.(values),
                   ms = 10)

    p2 = scatter(map(x -> x.θ, ps), map(x -> x.ϕ, ps), zcolor = Float64.(values),
                 ms = 10)

    return (p1, p2)
end
