abstract type Geometry end

function geometry_clear(g::Geometry)
    ccall((:geom_clear, "build/particular_solution"), Nothing, (Ref{Geometry},), g)
end

function recompute(g::Geometry, prec::Int = workingprecision(ArbReal))
    ccall((:geom_compute, "build/particular_solution"), Nothing,
          (Ref{Geometry}, Int), g, prec)
end

"""
Returns the Taylor polynomial of the parameterization of the boundary.
"""
function parameterization(t::ArbReal, g::Geometry, order::Int, vertex::Int = 1)
    n = order + 1
    r = _arb_vec_init(n)
    θ = _arb_vec_init(n)
    ccall((:parametrization, "build/particular_solution"), Cvoid,
          (Ptr{ArbReal}, Ptr{ArbReal}, Ref{Geometry}, Ref{ArbReal},
           Int, Int, Int), r, θ, g, t, n, vertex - 1, workingprecision(t))

    res = (unsafe_load_ArbRealPtr(r, n), unsafe_load_ArbRealPtr(θ, n))

    _arb_vec_clear(r, n)
    _arb_vec_clear(θ, n)

    return res
end
