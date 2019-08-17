function particularsolution_c(enclosure::ArbReal, g::Geometry, opt::Options)
    res = copy(enclosure)

    ccall((:particular_solution_enclosure, "build/particular_solution"), Nothing,
          (Ref{ArbReal}, Ref{Geometry}, Ref{Options}, Int), res, g, opt,
          workingprecision(enclosure))

    return res
end

function particularsolution(enclosure::ArbReal, g::Geometry, opt::Options)
    enclosurelocal = copy(enclosure)
    recompute(g, 2*workingprecision(ArbReal))

    Ns = opt.N_beg:opt.N_step:opt.N_end
    res = Array{Tuple{Array{ArbReal, 1}, ArbReal}}(undef, length(Ns))
    @showprogress for i in 1:length(Ns)
        N = Ns[i]
        # Compute tolerance and precision to use
        tol = radius(enclosurelocal)*opt.tol_relative

        prec = -Int(ceil(log(2, tol)*opt.prec_factor))
        if prec > workingprecision(enclosurelocal)
            setworkingprecision(ArbReal, prec)
            enclosurelocal = ArbReal{prec}(enclosurelocal)
            recompute(g, 2*prec)
        end

        #@show tol
        #@show prec

        p = Points(g, 2N, 2N)
        #setboundary!(p, g)
        #setinterior!(p, g)

        ν = minimizesigma(enclosurelocal, g, p, N, tol)

        #println("Minimum: $ν")

        coeffs = coefficientssigma(ν, g, p, N)[1]

        @show coeffs

        enclosurelocal = enclose(enclosurelocal, ν, coeffs, g, 1)

        @show enclosurelocal
        res[i] = (coeffs, copy(enclosurelocal))
    end

    return res
end
