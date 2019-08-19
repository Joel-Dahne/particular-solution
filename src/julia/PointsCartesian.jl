mutable struct PointsCartesian <: Points
    boundary::Int
    interior::Int
    total::Int
    r::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}
    θ::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}

    function PointsCartesian(g::GeometryCartesian, boundary::Int, interior::Int)
        p = new()
        ccall((:points_init, "build/particular_solution"), Cvoid,
              (Ref{Points}, Ref{GeometryCartesian}, Int, Int), p, g, boundary, interior)
        finalizer(points_clear, p)
        setboundary!(p, g)
        setinterior!(p, g)

        return p
    end
end

function Points(g::GeometryCartesian,
                boundary::Array{Polar{ArbReal}, 1},
                interior::Array{Polar{ArbReal}, 1})
    p = Points(g, length(boundary), length(interior))

    setboundary!(p, boundary)
    setinterior!(p, interior)

    return p
end

function boundary(p::PointsCartesian)
    #start = Int(sum(collect(g.vertices)[1:vertex-1]))*p.boundary + 1
    r = unsafe_load_ArbRealPtr(p.r[1], p.boundary)
    θ = unsafe_load_ArbRealPtr(p.θ[1], p.boundary)

    res = Array{Polar{ArbReal{workingprecision(ArbReal)}}}(undef, p.boundary)
    for i in 1:p.boundary
        res[i] = Polar(r[i], θ[i])
    end

    return res
end

function interior(p::PointsCartesian)
    r = unsafe_load_ArbRealPtr(p.r[1], p.total)[p.boundary + 1:end]
    θ = unsafe_load_ArbRealPtr(p.θ[1], p.total)[p.boundary + 1:end]

    res = Array{Polar{ArbReal{workingprecision(ArbReal)}}}(undef, p.interior)
    for i in 1:p.interior
        res[i] = Polar(r[i], θ[i])
    end

    return res
end

function setboundary!(p::PointsCartesian, values::Array{Polar{T}, 1}) where T<:ArbReal
    @assert length(values) == p.boundary
    r = map(x -> x.r, values)
    θ = map(x -> x.θ, values)

    unsafe_store_ArbRealPtr!(p.r[1], r)
    unsafe_store_ArbRealPtr!(p.θ[1], θ)
end

function setboundary!(p::PointsCartesian, g::GeometryCartesian) where T<:ArbReal
    n = p.boundary
    P = workingprecision(ArbReal)
    values = Array{Polar{ArbReal{P}}}(undef, n)

    v = xy(g)
    w = SVector{2}(ArbReal.([1, 0]))
    for i in 1:n
        t = ArbReal(i)/(n + 1)

        values[i] = PolarFromCartesian()(v + t.*(w - v))
    end

    setboundary!(p, values)
end

function setinterior!(p::PointsCartesian, values::Array{Polar{T}, 1}) where T <: ArbReal
    @assert length(values) == p.interior
    r = map(x -> x.r, values)
    θ = map(x -> x.θ, values)

    unsafe_store_ArbRealPtr!(p.r[1], r, p.boundary + 1)
    unsafe_store_ArbRealPtr!(p.θ[1], θ, p.boundary + 1)
end

function setinterior!(p::PointsCartesian, g::GeometryCartesian) where T <: ArbReal
    n = p.interior
    P = workingprecision(ArbReal)
    values = Array{Polar{ArbReal{P}}}(undef, n)

    v = xy(g)
    w = SVector{2}(ArbReal.([1, 0]))
    for i in 1:n
        (c1, c2, c3) = ArbReal.((rand(), rand(), rand()))
        s = c1 + c2 + c3
        c1 /= s
        c2 /= s
        c3 /= s

        values[i] = PolarFromCartesian()(c2.*w .+ c3.*v)
    end

    setinterior!(p, values)
end
