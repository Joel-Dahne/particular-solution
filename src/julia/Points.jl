mutable struct Points
    boundary::Int
    interior::Int
    total::Int
    thetas::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}
    phis::Tuple{Ptr{ArbReal}, Ptr{ArbReal}, Ptr{ArbReal}}

    function Points(g::Geometry, perboundary::Int, interior::Int)
        p = new()
        ccall((:points_init, "build/particular_solution"), Cvoid,
              (Ref{Points}, Ref{Geometry}, Int, Int), p, g, perboundary, interior)
        finalizer(points_clear, p)
        ccall((:boundary, "build/particular_solution"), Cvoid,
              (Ref{Points}, Ref{Geometry}, Int), p, g, workingprecision(ArbReal))
        ccall((:interior, "build/particular_solution"), Cvoid,
              (Ref{Points}, Ref{Geometry}, Int), p, g, workingprecision(ArbReal))

        return p
    end
end

function Points(g::GeometrySpherical,
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

function points_clear(p::Points)
    ccall((:points_clear, "build/particular_solution"), Cvoid, (Ref{Points},), p)
end

function setboundary!(p::Points, values::Array{ISOSpherical{ArbReal}, 1},
                      vertex::Int = 1)
    @assert length(values) == p.boundary
    θ = map(x -> x.θ, values)
    ϕ = map(x -> x.ϕ, values)

    unsafe_store_ArbRealPtr!(p.thetas[vertex], θ)
    unsafe_store_ArbRealPtr!(p.phis[vertex], ϕ)
end

function setboundary!(p::Points, values::Tuple{Array{ISOSpherical{ArbReal}, 1},
                                               Array{ISOSpherical{ArbReal}, 1},
                                               Array{ISOSpherical{ArbReal}, 1}})
    for vertex in 1:3
        setboundary!(p, values[vertex], vertex)
    end
end

function setinterior!(p::Points, values::Array{ISOSpherical{ArbReal}, 1},
                      vertex::Int = 1)
    @assert length(values) == p.interior
    θ = map(x -> x.θ, values)
    ϕ = map(x -> x.ϕ, values)

    unsafe_store_ArbRealPtr!(p.thetas[vertex], θ, p.boundary + 1)
    unsafe_store_ArbRealPtr!(p.phis[vertex], ϕ, p.boundary + 1)
end

function setinterior!(p::Points, values::Tuple{Array{ISOSpherical{ArbReal}, 1},
                                               Array{ISOSpherical{ArbReal}, 1},
                                               Array{ISOSpherical{ArbReal}, 1}})
    for vertex in 1:3
        setinterior!(p, values[vertex], vertex)
    end
end

function boundary(p::Points, vertex::Int = 1)
    θ = unsafe_load_ArbRealPtr(p.thetas[vertex], p.boundary)
    ϕ = unsafe_load_ArbRealPtr(p.phis[vertex], p.boundary)

    res = Array{ISOSpherical{ArbReal{workingprecision(ArbReal)}}}(undef, p.boundary)
    for i in 1:p.boundary
        res[i] = ISOSpherical(one(ArbReal), θ[i], ϕ[i])
    end

    return res
end

function interior(p::Points, vertex::Int = 1)
    θ = unsafe_load_ArbRealPtr(p.thetas[vertex], p.total)[p.boundary + 1:end]
    ϕ = unsafe_load_ArbRealPtr(p.phis[vertex], p.total)[p.boundary + 1:end]

    res = Array{ISOSpherical{ArbReal{workingprecision(ArbReal)}}}(undef, p.interior)
    for i in 1:p.interior
        res[i] = ISOSpherical(one(ArbReal), θ[i], ϕ[i])
    end

    return res
end
