mutable struct PointsSpherical <: Points
    boundary::Int
    interior::Int
    total::Int
    thetas::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}
    phis::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}

    function PointsSpherical(g::GeometrySpherical, perboundary::Int, interior::Int)
        p = new()
        ccall((:points_init, "build/particular_solution"), Cvoid,
              (Ref{Points}, Ref{GeometrySpherical}, Int, Int), p, g, perboundary, interior)
        finalizer(points_clear, p)
        ccall((:boundary, "build/particular_solution"), Cvoid,
              (Ref{Points}, Ref{GeometrySpherical}, Int), p, g, workingprecision(ArbReal))
        ccall((:interior, "build/particular_solution"), Cvoid,
              (Ref{Points}, Ref{GeometrySpherical}, Int), p, g, workingprecision(ArbReal))

        return p
    end
end

function PointsSpherical(g::GeometrySpherical,
                boundary::Tuple{Array{ISOSpherical{ArbReal}, 1},
                                Array{ISOSpherical{ArbReal}, 1},
                                Array{ISOSpherical{ArbReal}, 1}},
                interior::Tuple{Array{ISOSpherical{ArbReal}, 1},
                                Array{ISOSpherical{ArbReal}, 1},
                                Array{ISOSpherical{ArbReal}, 1}})
    p = Points(g, length(boundary[1])/sum(g.vertices), length(interior[1]))

    setboundary!(p, boundary)
    setinterior!(p, interior)

    return p
end

function boundary(p::PointsSpherical, vertex::Int = 1)
    #start = Int(sum(collect(g.vertices)[1:vertex-1]))*p.boundary + 1
    θ = unsafe_load_ArbRealPtr(p.thetas[vertex], p.boundary)
    ϕ = unsafe_load_ArbRealPtr(p.phis[vertex], p.boundary)

    res = Array{ISOSpherical{ArbReal{workingprecision(ArbReal)}}}(undef, p.boundary)
    for i in 1:p.boundary
        res[i] = ISOSpherical(one(ArbReal), θ[i], ϕ[i])
    end

    return res
end

function interior(p::PointsSpherical, vertex::Int = 1)
    θ = unsafe_load_ArbRealPtr(p.thetas[vertex], p.total)[p.boundary + 1:end]
    ϕ = unsafe_load_ArbRealPtr(p.phis[vertex], p.total)[p.boundary + 1:end]

    res = Array{ISOSpherical{ArbReal{workingprecision(ArbReal)}}}(undef, p.interior)
    for i in 1:p.interior
        res[i] = ISOSpherical(one(ArbReal), θ[i], ϕ[i])
    end

    return res
end

function setboundary!(p::PointsSpherical, values::Array{ISOSpherical{T}, 1},
                      vertex::Int = 1) where T<:ArbReal
    @assert length(values) == p.boundary
    θ = map(x -> x.θ, values)
    ϕ = map(x -> x.ϕ, values)

    unsafe_store_ArbRealPtr!(p.thetas[vertex], θ)
    unsafe_store_ArbRealPtr!(p.phis[vertex], ϕ)
end

function setboundary!(p::PointsSpherical, values::Tuple{Array{ISOSpherical{T}, 1},
                                               Array{ISOSpherical{T}, 1},
                                               Array{ISOSpherical{T}, 1}}) where T <: ArbReal
    for vertex in 1:3
        setboundary!(p, values[vertex], vertex)
    end
end

"""
Compute points on the boundary.
"""
function setboundary!(p::PointsSpherical, g::GeometrySpherical)
    n = Int(p.boundary/numactivevertices(g))
    P = workingprecision(ArbReal)
    values = (Array{ISOSpherical{ArbReal{P}}}(undef, p.boundary),
              Array{ISOSpherical{ArbReal{P}}}(undef, p.boundary),
              Array{ISOSpherical{ArbReal{P}}}(undef, p.boundary))

    for vertex = 1:3
        for i in 1:p.boundary
            values[vertex][i] = ISOSpherical(indeterminate(), indeterminate(),
                                             indeterminate())
        end
    end

    for vertex in activevertices(g)
        v = v2(g, vertex)
        w = v3(g, vertex)

        start = numactivevertices(g, vertex)*n

        for i in 1:n
            # This points does represent one on the current
            # boundary and should thus be computed.
            if g.half_edge[vertex] != 0
                t = ArbReal(i)/(2n)
            else
                t = ArbReal(i)/(n + 1)
            end

            xyz = v + t.*(w - v)
            values[vertex][start + i] = ISOSphericalFromCartesian()(xyz)
        end
    end

    setboundary!(p, values)
end

function setinterior!(p::PointsSpherical, values::Array{ISOSpherical{T}, 1},
                      vertex::Int = 1) where T <: ArbReal
    @assert length(values) == p.interior
    θ = map(x -> x.θ, values)
    ϕ = map(x -> x.ϕ, values)

    unsafe_store_ArbRealPtr!(p.thetas[vertex], θ, p.boundary + 1)
    unsafe_store_ArbRealPtr!(p.phis[vertex], ϕ, p.boundary + 1)
end

function setinterior!(p::PointsSpherical, values::Tuple{Array{ISOSpherical{T}, 1},
                                               Array{ISOSpherical{T}, 1},
                                               Array{ISOSpherical{T}, 1}}) where T <: ArbReal
    for vertex in 1:3
        setinterior!(p, values[vertex], vertex)
    end
end

function setinterior!(p::PointsSpherical, g::GeometrySpherical, seed::Int = 42)
    P = workingprecision(ArbReal)
    Random.seed!(seed)
    n = p.interior
    values = (Array{ISOSpherical{ArbReal{P}}}(undef, n),
              Array{ISOSpherical{ArbReal{P}}}(undef, n),
              Array{ISOSpherical{ArbReal{P}}}(undef, n))

    # We first generate random points with the parameterization
    # originating from the first vertex. Then we compute coordinates
    # for these points in the other parameterizations as well.
    v2mv1 = v2(g, 1) - v1(g, 1)
    v3mv1 = v3(g, 1) - v1(g, 1)

    # Transformation from the first coordinate system to the second.
    # Given by rotation by c around the y-axis followed by rotation by
    # B around the z-axis. Here c is the angle between v_1 and v_2 and
    # B is the angle of the triangle at the vertex at the north pole.
    # The names c and B corresponds with the naming used on Wikipedia
    # if we place v_1 on the top
    # https://en.wikipedia.org/wiki/File:Spherical_trigonometry_basic_triangle.svg
    c = -acos(v1(g)'v2(g))
    B = (ArbReal(Rational(angles(g)[2])) - 1)*ArbReal(π)
    T1 = LinearMap(RotZ(B)) ∘ LinearMap(RotY(c))

    # Transformation from the second coordinate system to the third.
    # Given by rotation by a around the y-axis followed by rotation by
    # C around the z-axis. Here a is the angle between v_2 and v_3 and
    # C is the angle of the triangle at the vertex in the bottom left.
    # The names a and C corresponds with the naming used on Wikipedia
    # if we place v_1 on the top
    # https://en.wikipedia.org/wiki/File:Spherical_trigonometry_basic_triangle.svg
    c = -acos(v2(g)'v3(g))
    B = (ArbReal(Rational(angles(g)[3])) - 1)*ArbReal(π)
    T2 = LinearMap(RotZ(B)) ∘ LinearMap(RotY(c))

    for i in 1:n
        s = rand()
        t = rand()

        # If s >= 1 - t set s = 1 - t, t = 1 - s to make them satisfy
        # s < 1 - t.
        if s >= 1 - t
            s, t = 1 - t, 1 - s
        end

        xyz = v1(g, 1) + v2mv1*s + v3mv1*t

        values[1][i] = ISOSphericalFromCartesian()(xyz)

        xyz = T1(xyz)
        values[2][i] = ISOSphericalFromCartesian()(xyz)

        xyz = T2(xyz)
        values[3][i] = ISOSphericalFromCartesian()(xyz)
    end

    setinterior!(p, values)
end

"""
Similar to setinterior! but randomly chooses the points in a way that
gives a uniform distribution on the spherical triangle.
"""
function setinterior2!(p::PointsSpherical, g::GeometrySpherical, seed::Int = 42)
    P = workingprecision(ArbReal)
    Random.seed!(seed)
    n = p.interior
    values = (Array{ISOSpherical{ArbReal{P}}}(undef, n),
              Array{ISOSpherical{ArbReal{P}}}(undef, n),
              Array{ISOSpherical{ArbReal{P}}}(undef, n))

    # The points are generated by uniformly choosing a value for ϕ and
    # for θ we uniformly choose an u ∈ [-1, 1] and let θ = acos(1 -
    # u), this gives a uniform distribution on the sphere. To ensure
    # that the generated point is indeed inside the triangle we
    # restrict the possible values for ϕ, given by the angle at the
    # north pole. For θ we generate it and then check that it's inside
    # the triangle by checking that it's above the great circle
    # determining the lower boundary of the spherical triangle, if
    # it's outside we generate ϕ and θ again until we get something
    # which is inside.

    # Transformation from the first coordinate system to the second.
    # Given by rotation by c around the y-axis followed by rotation by
    # B-π around the z-axis. Here c is the angle between v_1 and v_2 and
    # B is the angle of the triangle at the vertex at the north pole.
    # The names c and B corresponds with the naming used on Wikipedia
    # if we place the first angle on the top
    # https://en.wikipedia.org/wiki/File:Spherical_trigonometry_basic_triangle.svg
    c = -acos(v1(g)'v2(g))
    B = (ArbReal(Rational(angles(g)[2])) - 1)*ArbReal(π)
    T1 = LinearMap(RotZ(B)) ∘ LinearMap(RotY(c))

    # Transformation from the second coordinate system to the third.
    # Given by rotation by a around the y-axis followed by rotation by
    # C-π around the z-axis. Here a is the angle between v_2 and v_3 and
    # C is the angle of the triangle at the vertex in the bottom left.
    # The names a and C corresponds with the naming used on Wikipedia
    # if we place the first angle on the top
    # https://en.wikipedia.org/wiki/File:Spherical_trigonometry_basic_triangle.svg
    a = -acos(v2(g)'v3(g))
    C = (ArbReal(Rational(angles(g)[3])) - 1)*ArbReal(π)
    T2 = LinearMap(RotZ(C)) ∘ LinearMap(RotY(a))

    # Maximum value of ϕ
    A = Float64(Rational(angles(g)[1]))*π

    # Compute coefficients for plane determining the great circle
    # giving the lower boundary of the spherical triangle.
    (d, e, f) = greatcircleplane(g)

    for i in 1:n
        # Randomly pick θ and ϕ
        θ = acos(1 - 2rand())
        ϕ = rand()*A

        # To determine if the point (θ, ϕ) is inside the triangle we
        # check if θ is above the θ value of great circle going
        # between the two vertices of the triangle not on the north
        # pole at the value ϕ.
        while θ > greatcircle(ϕ, d, e, f)
            θ = acos(1 - 2rand())
            ϕ = rand()*A
        end

        xyz = CartesianFromSpherical()(ISOSpherical(ArbReal(1), ArbReal(θ), ArbReal(ϕ)))

        values[1][i] = ISOSpherical(ArbReal(1), ArbReal(θ), ArbReal(ϕ))

        xyz = T1(xyz)
        values[2][i] = ISOSphericalFromCartesian()(xyz)

        xyz = T2(xyz)
        values[3][i] = ISOSphericalFromCartesian()(xyz)
    end

    setinterior!(p, values)
end
