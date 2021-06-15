# ====================================================================
# MultiModel
#
mutable struct MultiModel
    models::Vector{Model}
    patchfuncts::Vector{ExprFunction}
    patchcomps::Vector{OrderedDict{Symbol, PatchComp}}

    function MultiModel(v::Vararg{Model})
        model = new([v...], Vector{ExprFunction}(),
                    Vector{OrderedDict{Symbol, PatchComp}}())
        evaluate!(model)
        return model
    end
end

Base.getindex(m::MultiModel, id::Int) = m.models[id]
Base.length(m::MultiModel) = length(m.models)


function evaluate!(multi::MultiModel)
    for i in 1:length(multi.models)
        evaluate!(multi.models[i])
    end
    empty!(multi.patchcomps)
    for i in 1:length(multi.models)
        push!(multi.patchcomps, multi.models[i].peval.patchcomps)
    end
    patch_params(multi)
    quick_evaluate(multi)
end

function patch_params(multi::MultiModel)
    for i in 1:length(multi.models)
        patch_params(multi.models[i])
    end
    for pf in multi.patchfuncts
        pf.funct(multi.patchcomps)
    end
    nothing
end

function quick_evaluate(multi::MultiModel)
    for i in 1:length(multi.models)
        quick_evaluate(multi.models[i])
    end
end

function push!(multi::MultiModel, model::Model)
    push!(multi.models, model)
    evaluate!(multi)
end

function patch!(multi::MultiModel, exfunc::ExprFunction)
    push!(multi.patchfuncts, exfunc)
    evaluate!(multi)
    return multi
end


function fit!(multi::MultiModel, data::Vector{Measures{N}};
              minimizer=lsqfit()) where N
    timestamp = now()
    evaluate!(multi)

    lparams = Vector{Parameter}()
    free = Vector{Bool}()
    source = Vector{NTuple{2,Int}}()
    for id in 1:length(multi)
        model = multi[id]
        i = 1
        for (cname, ceval) in model.cevals
            for (pid, par) in ceval.params
                if !(par.low <= par.val <= par.high)
                    s = "Value outside limits for param [$(cname)].$(pid.name):\n" * string(par)
                    error(s)
                end
                push!(lparams, par)
                push!(free, (!par.fixed)  &&  (model.cevals[cname].cfixed == 0))
                push!(source, (id, i))
                i += 1
            end
        end
    end
    ifree = findall(free)
    @assert length(ifree) > 0 "No free parameter in the model"

    # Flatten empirical data
    data1d = [flatten(data[id], multi[id].domain) for id in 1:length(multi)]

    prog = ProgressUnknown("Minimizer iteration:", dt=0.5, showspeed=true)
    ndata = [length(multi[id]()) for id in 1:length(multi)]
    dof = sum(ndata) - length(ifree)
    @assert dof >= 1

    # Evaluate normalized residuals starting from free parameter values
    residuals = fill(NaN, sum(ndata))
    function pval2resid(pvalues_free::Vector{Float64})
        for i in 1:length(ifree)
            j = ifree[i]
            multi[source[j][1]].peval.pvalues[source[j][2]] = pvalues_free[i]  # update parameter values
        end
        patch_params(multi)
        quick_evaluate(multi)
        i1 = 1
        for id in 1:length(multi)
            i2 = i1 + ndata[id] - 1
            residuals[i1:i2] .= (multi[id]() .- data1d[id].val) ./ data1d[id].unc
            i1 += ndata[id]
        end
        ret = residuals
        evaluate_showvalues(ret) = () -> [(:red_chisq, sum(abs2.(ret)) / dof)]
        ProgressMeter.next!(prog; showvalues = evaluate_showvalues(ret))
        return ret
    end

    (status, best_val, best_unc) = minimize(minimizer, pval2resid, lparams[ifree])
    pval2resid(best_val)
    ProgressMeter.finish!(prog)

    # Copy best fit values back into components.  This is needed since
    # the evaluated components are stored in the Model (rather than in
    # BestFitResult), hence I do this to maintain a coherent status.
    for id in 1:length(multi)
        model = multi[id]
        i = 1
        for (cname, ceval) in model.cevals
            for (pid, par) in ceval.params
                par.val = model.peval.pvalues[i]
                i += 1
            end
        end
    end

    # Prepare output
    bfmodels = Vector{OrderedDict{Symbol, BestFitComp}}()
    for id in 1:length(multi)
        model = multi[id]
        comps = OrderedDict{Symbol, BestFitComp}()
        i = 1
        for (cname, ceval) in model.cevals
            comps[cname] = BestFitComp()
            for (pid, par) in ceval.params
                if (!par.fixed)  &&  (model.cevals[cname].cfixed == 0)
                    bfpar = BestFitParam(model.peval.pvalues[i], popfirst!(best_unc), false, model.peval.patched[i])
                else
                    bfpar = BestFitParam(model.peval.pvalues[i], NaN                , true , model.peval.patched[i])
                end
                push!(comps[cname], pid.name, bfpar)
                i += 1
            end
        end
        push!(bfmodels, comps)
    end
    
    cost = sum(abs2, residuals)
    elapsed = now() - timestamp
    @assert isa(elapsed, Millisecond)
    result = BestFitMultiResult(bfmodels, length(residuals), dof, cost, status,
                                logccdf(Chisq(dof), cost) * log10(exp(1)),
                                timestamp, elapsed.value / 1.e3)
    return result
end
