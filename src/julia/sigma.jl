function sigma(ν::ArbReal, g::Geometry, p::Points, N::Int)
    res = typeof(ν)(0)
    ccall((:sigma, "build/particular_solution"), Cvoid,
          (Ref{ArbReal}, Ref{Geometry}, Ref{Points}, Int, Ref{ArbReal}, Int),
          res, g, p, N, ν, workingprecision(ν))

    return res
end

function minimizesigma(enclosure::ArbReal, g::Geometry, p::Points, N::Int,
                       tol::ArbReal)
    res = typeof(enclosure)(0)

    ccall((:minimize_sigma, "build/particular_solution"), Cvoid,
          (Ref{ArbReal}, Ref{Geometry}, Ref{Points}, Int, Ref{ArbReal},
           Ref{ArbReal}, Int),
          res, g, p, N, enclosure, tol, workingprecision(enclosure))

    return res
end

function coefficientssigma(ν::ArbReal, g::Geometry, p::Points, N::Int)
    v = [_arb_vec_init(N), _arb_vec_init(N), _arb_vec_init(N)]
    ccall((:coefs_sigma, "build/particular_solution"), Cvoid,
          (Ref{Ptr{ArbReal}}, Ref{Geometry}, Ref{Points}, Int, Ref{ArbReal}, Int),
          v, g, p, N, ν, workingprecision(ν))

    res = (unsafe_load_ArbRealPtr(v[1], N),
           unsafe_load_ArbRealPtr(v[2], N),
           unsafe_load_ArbRealPtr(v[3], N))

    _arb_vec_clear(v[1], N)
    _arb_vec_clear(v[2], N)
    _arb_vec_clear(v[3], N)

    return res
end
