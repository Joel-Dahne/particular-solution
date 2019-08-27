function plottimeoutput(to, opt::Options)
    Ns = opt.N_beg:opt.N_step:opt.N_end
    minimumtime = map(N -> TimerOutputs.time(to["$N"]["Minimum"])*1e-9, Ns)
    enclosuretime = map(N -> TimerOutputs.time(to["$N"]["Enclosure"])*1e-9, Ns)

    p = scatter([Ns Ns], [minimumtime enclosuretime],
                label = ["Minimum" "Enclosure"],
                xticks = Ns,
                xlabel = "N",
                ylabel = "seconds")
end
