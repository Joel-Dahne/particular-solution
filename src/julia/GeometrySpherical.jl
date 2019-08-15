mutable struct GeometrySpherical
    angles::Ptr{fmpq}
    v1::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}
    v2::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}
    v3::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}
    theta_lower::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}
    theta_upper::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}
    vertices::Tuple{Int32, Int32, Int32}
    half_edge::Tuple{Int32, Int32, Int32}

    function GeometrySpherical()
        g = new()
        ccall((:geom_init, "build/particular_solution"), Cvoid,
              (Ref{GeometrySpherical},), g)
        finalizer(geometry_clear, g)
        return g
    end
end

function geometry_clear(g::GeometrySpherical)
    ccall((:geom_clear, "build/particular_solution"), Nothing, (Ref{GeometrySpherical},), g)
end

function recompute(g::GeometrySpherical, prec::Int = workingprecision(ArbReal))
    ccall((:geom_compute, "build/particular_solution"), Nothing,
          (Ref{GeometrySpherical}, Int), g, prec)
end

function activevertices(g::GeometrySpherical)
    findall(x -> x == 1, g.vertices)
end

function numactivevertices(g::GeometrySpherical, before::Int = 4)
    res = Int(0)
    for i in 1:min(before-1, length(g.vertices))
        res += g.vertices[i]
    end
    res
end

function angles(g::GeometrySpherical)
    (unsafe_load(g.angles, 1),
     unsafe_load(g.angles, 2),
     unsafe_load(g.angles, 3))
end

function v1(g::GeometrySpherical, vertex::Int = 1)
    v = unsafe_load_ArbRealPtr(g.v1[vertex], 3)

    SVector{3}(v)
end

function v2(g::GeometrySpherical, vertex::Int = 1)
    v = unsafe_load_ArbRealPtr(g.v2[vertex], 3)

    SVector{3}(v)
end

function v3(g::GeometrySpherical, vertex::Int = 1)
    v = unsafe_load_ArbRealPtr(g.v3[vertex], 3)

    SVector{3}(v)
end

function vectors(g::GeometrySpherical, vertex::Int = 1)
    return [v1(g, vertex), v2(g, vertex), v3(g, vertex)]
end

function θ(g::GeometrySpherical, vertex::Int = 1)
    (unsafe_load_ArbRealPtr(g.theta_lower[vertex], 1)[1],
     unsafe_load_ArbRealPtr(g.theta_upper[vertex], 1)[1])
end
