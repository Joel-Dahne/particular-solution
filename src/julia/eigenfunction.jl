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

function eigenfunction(θ::ArbReal, ϕ::ArbReal, ν::ArbReal,
                       coeffs::Array{ArbReal, 1}, g::Geometry, vertex::Int = 1)
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
                       coeffs::Array{ArbReal, 1}, g::Geometry, vertex::Int = 1) where T
    eigenfunction(point.θ, point.ϕ, ν, coeffs, g, vertex)
end
