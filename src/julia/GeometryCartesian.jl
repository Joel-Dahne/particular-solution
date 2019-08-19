mutable struct GeometryCartesian <: Geometry
    angle::Ptr{fmpq}
    side_length::Ptr{ArbReal}
    x::Ptr{ArbReal}
    y::Ptr{ArbReal}
    vertices::Tuple{Int32, Int32, Int32}
    half_edge::Tuple{Int32, Int32, Int32}

    function GeometryCartesian()
        g = new()
        ccall((:geom_init, "build/particular_solution"), Cvoid,
              (Ref{GeometryCartesian},), g)
        finalizer(geometry_clear, g)
        return g
    end
end

function GeometryCartesian(angle, side_length)
    g = GeometryCartesian()
    g.vertices = (1, 0, 0)
    g.half_edge = (0, 0, 0)
    ccall((:fmpq_set, :libarb), Cvoid, (Ref{fmpq}, Ref{fmpq}),
          g.angle, fmpq(angle))
    ccall((:arb_set, :libarb), Cvoid, (Ref{ArbReal}, Ref{ArbReal}),
          g.side_length, side_length)
    recompute(g)
    g
end

function angle(g::GeometryCartesian)
    unsafe_load(g.angle, 1)
end

function sidelength(g::GeometryCartesian)
    unsafe_load_ArbRealPtr(g.side_length, 1)[1]
end

function xy(g::GeometryCartesian)
    x = unsafe_load_ArbRealPtr(g.x, 1)[1]
    y = unsafe_load_ArbRealPtr(g.y, 1)[1]
    SVector{2}([x, y])
end
