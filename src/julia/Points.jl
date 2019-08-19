abstract type Points end

function points_clear(p::Points)
    ccall((:points_clear, "build/particular_solution"), Cvoid, (Ref{Points},), p)
end

Points(g::GeometrySpherical, args...) = PointsSpherical(g, args...)

Points(g::GeometryCartesian, args...) = PointsCartesian(g, args...)
