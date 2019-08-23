function norm(g::Geometry, coeffs::Array{ArbReal{P}}, ν::ArbReal) where P
    res = ArbReal(0)
    N = length(coeffs)
    coeffsPtr = _arb_vec_init(N)
    unsafe_store_ArbRealPtr!(coeffsPtr, coeffs)

    ccall((:integral_norm, "build/particular_solution"), Cvoid,
          (Ref{ArbReal}, Ref{Geometry}, Ptr{ArbReal}, Int, Ref{ArbReal},
           Int, Int), res, g, coeffsPtr, N, ν, 0, workingprecision(ν))

    _arb_vec_clear(coeffsPtr, N)

    return res
end
