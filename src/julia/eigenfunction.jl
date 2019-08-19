function eigenfunction_basis(θ::ArbReal, ϕ::ArbReal, ν::ArbReal, μ::ArbReal)
    res = ArbReal(0)

    ccall((:eigenfunction_basis, "build/particular_solution"), Cvoid,
          (Ref{ArbReal}, Ref{ArbReal}, Ref{ArbReal}, Ref{ArbReal}, Ref{ArbReal}, Int),
          res, θ, ϕ, ν, μ, workingprecision(ArbReal))

    return res
end

function eigenfunction_basis(point::ISOSpherical{ArbReal{T}}, ν::ArbReal, μ::ArbReal) where T
    eigenfunction_basis(point.θ, point.ϕ, ν, μ)
end

"""
Compute the basis eigenfunctions for the plane. The input is the same
as for the spherical case but with different naming conventions. The
argument names here correspond to the naming in the plane, but it then
calls a method where the naming corresponds to the sphere.
"""
function eigenfunction_basis(point::Polar{ArbReal{T}}, λ::ArbReal, ν::ArbReal) where T
    eigenfunction_basis(point.r, point.θ, λ, ν)
end

function eigenfunction(θ::ArbReal, ϕ::ArbReal, ν::ArbReal,
                       coeffs::Array{ArbReal, 1}, g::GeometrySpherical, vertex::Int = 1)
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

function eigenfunction(point::ISOSpherical{ArbReal{T}}, ν::ArbReal,
                       coeffs::Array{ArbReal, 1}, g::GeometrySpherical, vertex::Int = 1) where T
    eigenfunction(point.θ, point.ϕ, ν, coeffs, g, vertex)
end

function eigenfunction(r::ArbReal, θ::ArbReal, λ::ArbReal,
                       coeffs::Array{ArbReal, 1}, g::GeometryCartesian)
    res = ArbReal(0)
    N = length(coeffs)
    coeffsPtr = _arb_vec_init(N)
    unsafe_store_ArbRealPtr!(coeffsPtr, coeffs)

    ccall((:eigenfunction, "build/particular_solution"), Cvoid,
          (Ref{ArbReal}, Ref{Geometry}, Ptr{ArbReal}, Int, Ref{ArbReal},
           Ref{ArbReal}, Ref{ArbReal}, Int, Int),
          res, g, coeffsPtr, N, r, θ, λ, 0, workingprecision(λ))

    _arb_vec_clear(coeffsPtr, N)

    return res
end

function eigenfunction(point::Polar{ArbReal{T}}, λ::ArbReal,
                       coeffs::Array{ArbReal, 1}, g::GeometryCartesian) where T
    eigenfunction(point.r, point.θ, λ, coeffs, g)
end
