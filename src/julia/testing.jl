using ArbNumerics
using ArbNumerics: radius, midpoint, ball
using Nemo
using StaticArrays
using CoordinateTransformations
using ProgressMeter

include("arb.jl")
include("GeometrySpherical.jl")
Geometry = GeometrySpherical
include("Options.jl")
include("Points.jl")

include("sigma.jl")
include("enclose.jl")
include("particularsolution.jl")

function getdomain(i::Int)
    g = Geometry()
    enclosure = ArbReal(0)
    opt = Options()

    ccall((:get_domain, "build/particular_solution"), Cvoid,
          (Ref{Geometry}, Ref{ArbReal}, Ref{Options}, Int32),
          g, enclosure, opt, i)
    ccall((:geom_compute, "build/particular_solution"), Cvoid,
          (Ref{Geometry}, Int), g, workingprecision(ArbReal))

    return (g, enclosure, opt)
end

function runall()
    numtriangles = 10
    res = Array{ArbReal}(undef, numtriangles)
    @showprogress for i in 1:numtriangles
        (g, enclosure, opt) = getdomain(i - 1)
        res[i] = particularsolution_c(enclosure, g, opt)
    end

    return res
end

setworkingprecision(ArbReal, 64)

(g, enclosure, opt) = getdomain(0)
N = 4
#p = Points(g, 2N, 2N)
#coeffs = coefficientssigma(midpoint(enclosure), g, p, N)[1]
#res = enclose(enclosure, midpoint(enclosure), coeffs, g, 1)
#res = particularsolution(enclosure, g, opt)
