function indeterminate!(x)
    ccall((:arb_indeterminate, :libarb), Cvoid, (ArbReal,), x)
end

function indeterminate()
    x = ArbReal(0)
    indeterminate!(x)
    return x
end

function unsafe_load_ArbRealPtr(ptr::Ptr{ArbReal}, size::Int)
    res = Array{ArbReal}(undef, size)

    for i in 1:size
        res[i] = ArbReal(0)
        ccall((:arb_set, :libarb), Cvoid, (Ref{ArbReal}, Ref{ArbReal}),
              res[i], ptr + (i - 1)*sizeof(ArbReal))
    end

    return res
end

function unsafe_store_ArbRealPtr!(ptr::Ptr{ArbReal}, values::Array{ArbReal, 1},
                                  start::Int = 1)
    for i in 1:length(values)
        ccall((:arb_set, :libarb), Cvoid, (Ref{ArbReal}, Ref{ArbReal}),
              ptr + ((start - 1) + (i - 1))*sizeof(ArbReal), values[i])
    end
end

function _arb_vec_init(n::Int)
    ccall((:_arb_vec_init, :libarb), Ptr{ArbReal}, (Int,), n)
end

function _arb_vec_clear(v::Ptr{ArbReal}, n::Int)
    ccall((:_arb_vec_clear, :libarb), Cvoid, (Ptr{ArbReal}, Int), v, n)
end
