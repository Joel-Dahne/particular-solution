abstract type Geometry end

function geometry_clear(g::Geometry)
    ccall((:geom_clear, "build/particular_solution"), Nothing, (Ref{Geometry},), g)
end

function recompute(g::Geometry, prec::Int = workingprecision(ArbReal))
    ccall((:geom_compute, "build/particular_solution"), Nothing,
          (Ref{Geometry}, Int), g, prec)
end
