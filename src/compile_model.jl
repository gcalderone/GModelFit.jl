function convert_patch_expression(fp::FitProblem, fd::FunctDesc, cur::String)
    @assert length(fd.args) in [1,2]
    out = deepcopy(fd.display)
    i = findfirst("->", out)
    @assert !isnothing(i)
    out = out[(maximum(i)+1):end]

    sym = string(fd.args[1])
    if length(fp.mevals) > 1
        for i in 1:length(fp.mevals)
            meval = fp.mevals[i]
            for (cname, comp) in meval.model.comps
                for (pname, par) in getparams(comp)
                    from = "(($(sym)[$i])[:$(cname)]).$(pname)"
                    to = "m$(i)_$(cname)_$(pname)"
                    out = replace(out, from => to)
                end
            end
        end
    else
        meval = fp.mevals[1]
        for (cname, comp) in meval.model.comps
            for (pname, par) in getparams(comp)
                from = "($(sym)[:$(cname)]).$(pname)"
                to = "m1_$(cname)_$(pname)"
                out = replace(out, from => to)
            end
        end
    end

    if length(fd.args) == 2
        sym = fd.args[2]
        out = "($(sym) -> ($(out)))($(cur))"
    end

    return out
end


function compile_model(fp::FitProblem{TFitStat}) where TFitStat
    tmp = Vector{String}()

    accum = Vector{Any}()
    push!(accum, :fp => fp)

    # Domains
    for i in 1:length(fp.mevals)
        push!(accum, Symbol("m$(i)domain") => fp.mevals[i].domain)
    end

    # Components
    empty!(tmp)
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        for (cname, comp) in meval.model.comps
            push!(accum, Symbol("m$(i)_$(cname)") => fp.mevals[i].model[cname])
        end
    end

    # Parameters
    fip = 0
    guess = Vector{Float64}()
    lowb  = Vector{Float64}()
    highb = Vector{Float64}()
    params = Vector{String}()
    actual = Vector{String}()
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        ip = 0
        for (cname, comp) in meval.model.comps
            for (pname, par) in getparams(comp)
                ip += 1
                spname = "m$(i)_$(cname)_$(pname)"
                if ip in meval.ifree
                    fip += 1
                    push!(guess, getproperty(fp.mevals[i].model[cname], pname).val)
                    push!(lowb , getproperty(fp.mevals[i].model[cname], pname).low)
                    push!(highb, getproperty(fp.mevals[i].model[cname], pname).high)
                    push!(params, "$(spname) = params[$(fip)]")

                    if !isnothing(par.patch)
                        @assert length(par.patch.args) == 2
                        push!(actual, "$(spname) = " * convert_patch_expression(fp, par.patch, spname))
                    elseif !isnothing(par.mpatch)
                        @assert length(par.mpatch.args) == 2
                        push!(actual, "$(spname) = " * convert_patch_expression(fp, par.mpatch, spname))
                    end
                else
                    if !isnothing(par.patch)
                        if isa(par.patch, Symbol)
                            push!(actual, "$(spname) = m$(i)_$(par.patch)_$(pname)")
                        else
                            @assert length(par.patch.args) == 1
                            push!(actual, "$(spname) = " * convert_patch_expression(fp, par.patch, spname))
                        end
                    elseif !isnothing(par.mpatch)
                        @assert length(par.mpatch.args) == 1
                        push!(actual, "$(spname) = " * convert_patch_expression(fp, par.mpatch, spname))
                    else
                        push!(accum, Symbol(spname) => getproperty(fp.mevals[i].model[cname], pname).val)  # fixed, non-patched value
                        push!(params, "$(spname) = shared.$(spname)")
                    end
                end
            end
        end
    end
    push!(accum, :guess => guess)
    push!(accum, :lowb  => lowb)
    push!(accum, :highb => highb)
    push!(accum, :prealloc => Dict{DataType, Vector{Vector}}())

    # evaluate! function
    io = IOBuffer()
    println(io, "function (output::AbstractArray{T, N}, params, shared) where {T, N}")

    println(io, "    # Buffers to store component evaluation")
    println(io, "    if !(T in keys(shared.prealloc))")
    println(io, "        shared.prealloc[T] = Vector{Vector{T}}()")
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        for (cname, ceval) in meval.cevals
            # Avoid using SizedVector here: an error may occur when ndata is large because of the attempt to unroll loops
            println(io, "        push!(shared.prealloc[T], Vector{T}(undef, " * string(length(meval.domain)) * "))")
        end
    end
    println(io, "    end")
    println(io, "    prealloc = shared.prealloc[T]")

    ip = 0
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        for (cname, ceval) in meval.cevals
            ip += 1
            println(io, "    m$(i)_$(cname)::Vector{T} = prealloc[$(ip)]")
        end
    end

    println(io, "\n    # Parameter values")
    for s in params
        println(io, "    ", s)
    end

    if length(actual) > 0
        println(io, "\n    # Patched parameters")
        for s in actual
            println(io, "    ", s)
        end
    end

    println(io, "\n    # Component evaluation")
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        for cname in meval.seq
            ceval = fp.mevals[i].cevals[cname]
            print(io, "    evaluate!(shared.m$(i)_$(cname), shared.m$(i)domain, ")
            print(io, "m$(i)_$(cname)")

            empty!(tmp)
            id = 1
            for d in dependencies(meval.model, cname, select_domain=true)
                push!(tmp, "coords(shared.m$(i)domain, $id)")
                id += 1
            end
            for d in dependencies(meval.model, cname, select_domain=false)
                push!(tmp, "m$(i)_$(d)")
            end
            if length(tmp) > 0
                print(io, ", [", join(tmp, ", "), "]")
            end

            empty!(tmp)
            for (pname, par) in getparams(fp.mevals[i].model[cname])
                push!(tmp, "m$(i)_$(cname)_$(pname)")
            end
            if length(tmp) > 0
                print(io, ", ", join(tmp, ", "))
            end
            println(io, ")")
        end
    end

    println(io, "\n    # Populate residuals")
    print(io, "    return populate_residuals!(shared.fp, ")
    empty!(tmp)
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        push!(tmp, "m$(i)_$(meval.maincomp)")
    end
    println(io, "[", join(tmp, ", "), "], output)")
    println(io, "end")

    funcdef = String(take!(io))
    # println(funcdef)
    return NamedTuple(accum), eval(Meta.parse(funcdef))
end
