mutable struct GeometrySpherical <: Geometry
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

"""
Return coefficients a, b, c giving the plane ax + by + cz = 0 which is
parallel to the great circle that defines the bottom boundary of the
spherical triangle.
"""
function greatcircleplane(g::GeometrySpherical, vertex::Int = 1)
    # Compute the normal vector of the plane
    w = cross(v2(g, vertex), v3(g, vertex))
    # w = [a, b, c]
    return (w[1], w[2], w[3])
end

"""
Compute the θ value corresponding to the give ϕ value on the great
circle defined by the plane ax + by + cz = 0.
"""
function greatcircle(ϕ, a, b, c)
    θ = atan(-c, a*cos(ϕ) + b*sin(ϕ))
    if θ < 0
        θ += ArbReal(π)
    end
    return θ
end
