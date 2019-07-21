mutable struct Options
    tol_relative::Float64
    prec_factor::Float64
    N_beg::Int
    N_end::Int
    N_step::Int
    plot_n::Int
    output::Int32
    output_final::Int32
    output_time::Int32

    function Options()
        opt = new()
        ccall((:options_default, "build/particular_solution"),
              Nothing, (Ref{Options},), opt)
        return opt
    end
end

function Options(opt::Options, Ns::StepRange)
    opt.N_beg = Ns.start
    opt.N_end = Ns.stop
    opt.N_step = Ns.step
    return opt
end

function Options(Ns::StepRange)
    opt = Options()
    options(opt, Ns)
end
