function evaluate!(comp::CompSlotEval{T, TComp, TDomain}, pvalues::AbstractVector) where {TComp, TDomain, T}
    doeval = (comp.counter == 0)                                  ||
        any(comp.lastparvalues   .!= pvalues)                     ||
        length(comp.lastdepscounter) != length(comp.deps)         ||
        any(comp.lastdepscounter .!= getfield.(comp.deps, :counter))
    if doeval
        if length(comp.deps) > 0
            evaluate!(comp.comp, comp.domain, comp.evalbuf, getfield.(comp.deps, :evalbuf), pvalues...)
            if length(  comp.lastdepscounter) != length(comp.deps)
                empty!( comp.lastdepscounter)
                append!(comp.lastdepscounter, getfield.(comp.deps, :counter))
            else
                comp.lastdepscounter .= getfield.(comp.deps, :counter)
            end
        else
            evaluate!(comp.comp, comp.domain, comp.evalbuf, pvalues...)
        end
        comp.lastparvalues .= pvalues
        comp.counter += 1
    end
    return comp.evalbuf
end

function evaluate!(model::AbstractModelEval)
    isdirty(model)  ||  return nothing

    empty!(model.ifree)
    for (gkey, group) in getgroups(model)
        empty!(group.seq)
        empty!(group.patched)
        append!(group.seq, find_sequence(model, gkey))
    end
    for (ckey, comp) in getcomps(model)
        gkey = group_key(model, ckey)
        empty!(comp.deps)
        for d in dependencies(comp.comp)
            @assert haskey(getcomps(model), (gkey..., d)) "No component has name $((gkey..., d))"
            push!(comp.deps, getcomps(model)[(gkey..., d)])
        end
    end
    i = 0
    for (pkey, par) in getparams(model)
        gkey  = group_key(model, pkey)
        ckey  = comp_key( model, pkey)
        cpkey = copa_key( model, pkey)
        cname = comp_name(model, pkey)
        i += 1
        model.pvalues[pkey] = par.val
        model.pactual[pkey] = par.val
        if par.fixed  ||  !isnothing(par.patch)  ||  isfrozen(getcomps(model)[ckey])
            setfield!(par, :actually_fixed, true)
        else
            setfield!(par, :actually_fixed, false)
            push!(model.ifree, i)
        end
        @assert isnothing(par.patch)  ||  isnothing(par.reparam) "A parameter can either be patch-ed or reparametrized, not both!"
        if !isnothing(par.patch)  ||  !isnothing(par.reparam)
            getgroups(model)[gkey].patched[cpkey] = getparams(model)[pkey]
        end
    end
    for (ckey, comp) in getcomps(model)
        haspath(model.pactual, ckey)  &&  (comp.pactual = values_view(view(model.pactual, ckey)))
    end
    setdirty!(model, false)
    evaluate!(model, getproperty.(collect(values(getfreepars(model))), :val))
    return nothing
end

function evaluate!(model::AbstractModelEval, freepvals::AbstractVector)
    i = model.ifree
    model.pvalues.values[i] .= freepvals
    model.pactual.values[i] .= freepvals

    for (gkey, group) in getgroups(model)
        pvalues = view(model.pvalues   , gkey)
        pactual = view(model.pactual   , gkey)
        for (cpkey, par) in group.patched
            if !isnothing(par.patch)
                if isa(par.patch, Symbol)  # use same param. value from a different component
                    pactual[cpkey] = pvalues[(par.patch, cpkey[2])]
                else                       # invoke a patch function
                    pactual[cpkey] = par.patch(pvalues)
                end
            else
                pactual[cpkey] = par.reparam(pvalues, pvalues[cpkey])
            end
        end

        comps = getcomps(model[gkey])
        for cname in group.seq
            comp = comps[cname]
            evaluate!(comp, comp.pactual)
        end
        main_comp = getcomps(model)[(gkey..., group.seq[end])]
        fold_model!(group.resp, main_comp.evalbuf, group.evalbuf)
    end

    if isa(model, Likelihood)
        i1 = 1
        for (gkey, group) in getgroups(model)
            nn = length(group.evalbuf)
            if nn > 0
                evaluate_resid!(group)
                i2 = i1 + nn - 1
                model.evalbuf[i1:i2] .= group.residuals
                i1 += nn
            end
        end
    end
    return nothing
end

function evaluate_resid!(lh::Likelihood, freepvals::Vector{T}) where {T}
    evaluate!(lh, freepvals)
    return convert(Vector{T}, lh.evalbuf)
end

function _set_bestfit!(model::AbstractModelEval, freepvals::Vector{Float64}, freepuncerts::Vector{Float64})
    evaluate!(model, freepvals)
    i = 1
    for (pkey, par) in getparams(model)
        if !par.fixed  &&  !par.actually_fixed
            par.val = freepvals[i]
            setfield!(par, :unc, freepuncerts[i])
            iszero(par.unc)  &&  setfield!(par, :unc, NaN)
            i += 1
        else
            setfield!(par, :unc, NaN)
        end
        setfield!(par, :mval, model.pactual[pkey])
    end
    evaluate!(model)
end

evaluate!(::LikelihoodSnapshot) = nothing
