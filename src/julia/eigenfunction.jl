"""
Compute the basis eigenfunctions, the approximate eigenfunction is
given by a linear combination of these. It takes as input a point and
two parameters.

For the spherical case the point is given by (θ, ϕ) and the parameters
by ν and μ, where ν corresponds to the eigenvalue and μ is the index
of the function.

For the Cartesian case the point is given by (r, θ) and the parameters
by λ and ν, where λ is the eigenvalue and ν the index of the function.
"""
function eigenfunction_basis(θ::ArbReal, ϕ::ArbReal, ν::ArbReal, μ::ArbReal)
    res = ArbReal(0)

    ccall((:eigenfunction_basis, "build/particular_solution"), Cvoid,
          (Ref{ArbReal}, Ref{ArbReal}, Ref{ArbReal}, Ref{ArbReal}, Ref{ArbReal}, Int),
          res, θ, ϕ, ν, μ, workingprecision(ArbReal))

    return res
end


function eigenfunction_basis(point::ISOSpherical, ν, μ)
    eigenfunction_basis(point.θ, point.ϕ, ν, μ)
end

function eigenfunction_basis(point::Polar, λ, ν)
    eigenfunction_basis(point.r, point.θ, λ, ν)
end

"""
Compute the Taylor expansion of the basis eigenfunctions. It takes as
input the Taylor expansion of a point and two parameters. The
parameters are the same as for eigenfunction_basis.

The point is either (θ, ϕ) for the spherical case or (r, θ) for the
Cartesian case. In both cases they should be the n-th order Taylor
expansions of the respective values. The output is the computed to the
same order.
"""
function eigenfunction_basis_series(θ::Vector{T}, ϕ::Vector{T},
                                    ν::ArbReal, μ::ArbReal) where T <: ArbReal
    @assert length(θ) == length(ϕ)
    n = length(θ)

    θPtr = _arb_vec_init(n)
    ϕPtr = _arb_vec_init(n)
    resPtr = _arb_vec_init(n)
    unsafe_store_ArbRealPtr!(θPtr, θ)
    unsafe_store_ArbRealPtr!(ϕPtr, ϕ)

    ccall((:eigenfunction_basis_series, "build/particular_solution"), Cvoid,
          (Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}, Ref{ArbReal}, Ref{ArbReal}, Int, Int),
          resPtr, θPtr, ϕPtr, ν, μ, n, workingprecision(ν))

    res = unsafe_load_ArbRealPtr(resPtr, n)

    _arb_vec_clear(θPtr, n)
    _arb_vec_clear(ϕPtr, n)
    _arb_vec_clear(resPtr, n)

    return res
end

function eigenfunction_basis_series(point::Vector{T}, ν, μ) where T <: ISOSpherical
    eigenfunction_basis_series(map(p -> p.θ, point), map(p -> p.ϕ, point), ν, μ)
end

function eigenfunction_basis_series(point::Vector{T}, ν, μ) where T <: Polar
    eigenfunction_basis_series(map(p -> p.r, point), map(p -> p.θ, point), ν, μ)
end

"""
Compute the approximate eigenfunction. The eigenfunction is given by a
linear combination of the basis eigenfunction with the given
coefficients. It takes as input a point, a parameter, the coefficients
and the geometrical data needed to determine the function.

For the spherical case the point is given by (θ, ϕ) and the parameter
by ν which corresponds to the eigenvalue. The geometrical information
is a GeometrySpherical as well as which vertex it should be evaluated
from.

For the Cartesian case the point is given by (r, θ) and the parameter
by λ which is the eigenvalue. The geometrical information is a
GeometryCartesian, the vertex should always be 1 since support for
other vertices is not implemented.
"""
function eigenfunction(θ::ArbReal, ϕ::ArbReal, ν::ArbReal, coeffs::Vector{T},
                       g::Geometry, vertex::Int = 1) where T <: ArbReal
    res = ArbReal(0)

    N = length(coeffs)
    coeffsPtr = _arb_vec_init(N)
    unsafe_store_ArbRealPtr!(coeffsPtr, coeffs)

    ccall((:eigenfunction, "build/particular_solution"), Cvoid,
          (Ref{ArbReal}, Ref{Geometry}, Ptr{ArbReal}, Int, Ref{ArbReal},
           Ref{ArbReal}, Ref{ArbReal}, Int, Int),
          res, g, coeffsPtr, N, θ, ϕ, ν, vertex - 1, workingprecision(ν))

    _arb_vec_clear(coeffsPtr, N)

    return res
end

function eigenfunction(point::ISOSpherical, ν, coeffs, g::GeometrySpherical,
                       vertex = 1)
    eigenfunction(point.θ, point.ϕ, ν, coeffs, g, vertex)
end

function eigenfunction(point::Polar, λ, coeffs, g::GeometryCartesian, vertex = 1)
    eigenfunction(point.r, point.θ, λ, coeffs, g, 1)
end

"""
Compute the Taylor expansion of the approximate eigenfunction. It takes as
input the Taylor expansion of a point, a parameter, the coefficients
and the geometrical data needed to determine the function. The
parameter, the coefficients and geometrical data is the same as for
eigenfunction.

The point is either (θ, ϕ) for the spherical case or (r, θ) for the
Cartesian case. In both cases they should be the n-th order Taylor
expansions of the respective values. The output is the computed to the
same order.
"""
function eigenfunction_series(θ::Vector{T}, ϕ::Vector{T}, ν::ArbReal,
                              coeffs::Vector{S} where S <: ArbReal, g::Geometry,
                              vertex::Int = 1) where T <: ArbReal
    @assert length(θ) == length(ϕ)
    n = length(θ)

    θPtr = _arb_vec_init(n)
    ϕPtr = _arb_vec_init(n)
    resPtr = _arb_vec_init(n)
    unsafe_store_ArbRealPtr!(θPtr, θ)
    unsafe_store_ArbRealPtr!(ϕPtr, ϕ)

    N = length(coeffs)
    coeffsPtr = _arb_vec_init(N)
    unsafe_store_ArbRealPtr!(coeffsPtr, coeffs)

    ccall((:eigenfunction_series, "build/particular_solution"), Cvoid,
          (Ptr{ArbReal}, Ref{Geometry}, Ptr{ArbReal}, Int, Ptr{ArbReal},
           Ptr{ArbReal}, Ref{ArbReal}, Int, Int, Int),
          resPtr, g, coeffsPtr, N, θPtr, ϕPtr, ν, vertex - 1, n, workingprecision(ν))

    res = unsafe_load_ArbRealPtr(resPtr, n)

    _arb_vec_clear(θPtr, n)
    _arb_vec_clear(ϕPtr, n)
    _arb_vec_clear(resPtr, n)

    _arb_vec_clear(coeffsPtr, N)

    return res
end

function eigenfunction_series(point::Vector{T}, ν, coeffs, g::GeometrySpherical,
                       vertex = 1) where T <: ISOSpherical
    eigenfunction_series(map(p -> p.θ, point), map(p -> p.ϕ, point), ν, coeffs, g, vertex)
end

function eigenfunction_series(point::Vector{T}, λ, coeffs,
                              g::GeometryCartesian) where T <: Polar
    eigenfunction_series(map(p -> p.r, point), map(p -> p.θ, point), λ, coeffs, g, 1)
end

"""
Enclose the eigenfunction using Taylor expansions.
"""
function eigenfunction_series_enclosure(t::ArbReal, ν::ArbReal,
                                        coeffs::Vector{T}, g::Geometry,
                                        order::Int, vertex::Int = 1) where T <: ArbReal
    # Compute the Taylor polynomial of the eigenfunction at the
    # midpoint of t to the given order
    (θ, ϕ) = parameterization(midpoint(t), g, order, vertex)
    midpointseries = eigenfunction_series(θ, ϕ, ν, coeffs, g, vertex)

    # Enclose the Taylor polynomial evaluated at t - midpoint(t)
    enclosure = ArbReal(0)
    midpointseriesPtr = _arb_vec_init(order + 1)
    unsafe_store_ArbRealPtr!(midpointseriesPtr, midpointseries)

    ccall((:_arb_poly_evaluate, :libarb), Cvoid,
          (Ref{ArbReal}, Ptr{ArbReal}, Int, Ref{ArbReal}, Int),
          enclosure, midpointseriesPtr, order + 1, t - midpoint(t), workingprecision(ν))

    _arb_vec_clear(midpointseriesPtr, order + 1)
    #@show enclosure
    #@show radius(enclosure)

    # Compute the Taylor polynomial of the eigenfunction on the
    # interval t of order one higher than the given one
    (θ, ϕ) = parameterization(t, g, order + 1, vertex)
    series = eigenfunction_series(θ, ϕ, ν, coeffs, g, vertex)
    #@show series[end]
    #@show radius(series[end])


    # Compute the rest term
    restterm = series[end]*pow(t - midpoint(t), order)
    #@show restterm
    #@show radius(restterm)

    return enclosure + restterm
end
