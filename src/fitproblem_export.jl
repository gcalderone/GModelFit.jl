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

function export_model(name::String, fp::FitProblem{TFitStat}) where TFitStat
    @assert !isfile(name * ".jl")
    tmp = Vector{String}()
    
    f = stdout # open(name * ".jl", "w")
    println(f)
    println(f, "using GModelFit, StaticArrays")
    println(f)
    println(f, "function $(name)_init(fp::GModelFit.FitProblem)")

    println(f, "    return (")
    println(f, "        fp = fp,")

    # Domains
    println(f, join(["        m$(i)domain = fp.mevals[$i].domain" for i in 1:length(fp.mevals)], ",\n"), ",")

    # Components
    empty!(tmp)
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        for (cname, comp) in meval.model.comps
            push!(tmp, "m$(i)_$(cname) = fp.mevals[$i].model[:$(cname)]")
        end
    end
    println(f, join("        " .* tmp, ",\n"), ",")

    # Parameters
    fip = 0
    empty!(tmp)
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
                    push!(tmp, "fp.mevals[$i].model[:$(cname)].$(pname).val")  # this is the guess value
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
                        println(f, "        $(spname) = fp.mevals[$i].model[:$(cname)].$(pname).val,")  # fixed, non-patched value
                        push!(params, "$(spname) = fp.$(spname)")
                    end
                end
            end
        end
    end
    println(f, "        guess = [")
    for i in 1:length(tmp)
        println(f, "            ", tmp[i], (i < length(tmp))  ?  ","  :  "", "  #  $i")
    end
    println(f, "        ])")
    println(f, "end")

    # evaluate! function
    println(f)
    println(f)
    println(f, "function $(name)_evaluate!(output::AbstractArray{T, N}, params, fp) where {T, N}") 

    println(f, "    # Buffers to store component evaluation")
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        for (cname, ceval) in meval.cevals
            println(f, "    m$(i)_$(cname) = SizedVector{" * string(length(meval.domain)) * ", T}(undef)")
        end
    end

    println(f, "\n    # Parameter values")
    for s in params
        println(f, "    ", s)
    end

    if length(actual) > 0
        println(f, "\n    # Patched parameters")
        for s in actual
            println(f, "    ", s)
        end
    end

    println(f, "\n    # Component evaluation")
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        for cname in meval.seq
            ceval = fp.mevals[i].cevals[cname]
            print(f, "    GModelFit.evaluate!(fp.m$(i)_$(cname), fp.m$(i)domain, ")
            print(f, "m$(i)_$(cname)")
            
            empty!(tmp)
            id = 1
            for d in dependencies(meval.model, cname, select_domain=true)
                push!(tmp, "coords(fp.m$(i)domain, $id)")
                id += 1
            end
            for d in dependencies(meval.model, cname, select_domain=false)
                push!(tmp, "m$(i)_$(d)")
            end
            if length(tmp) > 0
                print(f, ", [", join(tmp, ", "), "]")
            end

            empty!(tmp)
            for (pname, par) in getparams(fp.mevals[i].model[cname])
                push!(tmp, "m$(i)_$(cname)_$(pname)")
            end
            if length(tmp) > 0
                print(f, ", ", join(tmp, ", "))
            end
            println(f, ")")
        end
    end

    println(f, "\n    # Populate residuals")
    print(f, "    return GModelFit.populate_residuals!(fp.fp, ")
    empty!(tmp)
    for i in 1:length(fp.mevals)
        meval = fp.mevals[i]
        push!(tmp, "m$(i)_$(meval.maincomp)")
    end
    println("[", join(tmp, ", "), "], output)")
    println(f, "end")


    # close(f)
end


#=

fp = GModelFit.FitProblem(model, data);

GModelFit.export_model("aa", fp)

aa = aa_init(fp)
aa_evaluate!(fp.buffer, aa.guess, aa)


result = NonlinearSolve.solve(NonlinearSolve.NonlinearLeastSquaresProblem(
        NonlinearSolve.NonlinearFunction(aa_evaluate!, resid_prototype = zeros(GModelFit.ndata(fp))),
        aa.guess, aa))
                                       wrap.solver)


=#
