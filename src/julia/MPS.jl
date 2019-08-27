using ArbNumerics
import ArbNumerics
using ArbNumerics: radius, midpoint, ball
using Nemo
using StaticArrays
using CoordinateTransformations
using ProgressMeter
using LaTeXStrings
using Random
using LinearAlgebra
using TimerOutputs

include("arb.jl")
include("Geometry.jl")
include("GeometrySpherical.jl")
include("GeometryCartesian.jl")
include("Points.jl")
include("PointsSpherical.jl")
include("PointsCartesian.jl")
include("Options.jl")

include("sigma.jl")
include("enclose.jl")
include("particularsolution.jl")
include("eigenfunction.jl")
include("norm.jl")

include("examples.jl")
include("timings.jl")

function getdomain(::Type{T}, i::Int) where T <: Geometry
    g = T()
    enclosure = ArbReal(0)
    opt = Options()

    ccall((:get_domain, "build/particular_solution"), Cvoid,
          (Ref{T}, Ref{ArbReal}, Ref{Options}, Int32),
          g, enclosure, opt, i)
    ccall((:geom_compute, "build/particular_solution"), Cvoid,
          (Ref{T}, Int), g, workingprecision(ArbReal))

    return (g, enclosure, opt)
end

setworkingprecision(ArbReal, 64)

(g, enclosure, opt) = getdomain(GeometryCartesian, 0)

N = 8
p = Points(g, 2N, 2N)

#coefs = coefficientssigma(midpoint(enclosure), g, p, N)[1]
#newenclosure = enclose(enclosure, midpoint(enclosure), coefs, g, 1)
#res = particularsolution(enclosure, g, opt)
