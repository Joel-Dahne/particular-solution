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

    to = TimerOutput()

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

        @timeit to "$N" begin
            p = @timeit to "Points" Points(g, 2N, 2N)

            ν = @timeit to "Minimum" minimizesigma(enclosurelocal, g, p, N, tol)
            #@show ν

            coeffs = @timeit to "Coefficients" coefficientssigma(ν, g, p, N)[1]
            #@show coeffs

            enclosurelocal = @timeit to "Enclosure" enclose(enclosurelocal, ν, coeffs, g, 1)
            #@show enclosurelocal
        end

        res[i] = (coeffs, copy(enclosurelocal))
    end

    return res
    # To also return timing information
    #return (res, to)
end
