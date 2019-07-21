function enclose(enclosure::ArbReal, ν::ArbReal, coeffs::Array{ArbReal, 1},
                 g::Geometry, vertex::Int)
    N = length(coeffs)
    coeffsPtr = [_arb_vec_init(N), _arb_vec_init(N), _arb_vec_init(N)]
    unsafe_store_ArbRealPtr!(coeffsPtr[vertex], coeffs)

    ccall((:enclose, "build/particular_solution"), Cvoid,
          (Ref{ArbReal}, Ref{Geometry}, Ref{Ptr{ArbReal}}, Int, Ref{ArbReal}, Int, Int),
          enclosure, g, coeffsPtr, N, ν, vertex - 1, workingprecision(ν))

    for i in 1:3
        _arb_vec_clear(coeffsPtr[i], N)
    end

    return enclosure
end
