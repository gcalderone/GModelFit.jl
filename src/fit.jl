function fit!(model::Model, data::Measures{N};
              minimizer=lsqfit()) where N
    timestamp = now()
    evaluate!(model)
    mdc = MDComparison(model, data)
    @assert length(mdc.ifree) > 0 "No free parameter in the model"
    @assert mdc.dof >= 1

    prog = ProgressUnknown("Model evaluations:", dt=0.5, showspeed=true)
    function pval2resid(pvalues::Vector{Float64})
        out = try_pvalues(mdc, model, pvalues)
        evaluate_showvalues(x) = () -> [(:red_chisq, sum(abs2.(x)) / mdc.dof)]
        ProgressMeter.next!(prog; showvalues = evaluate_showvalues(out))
        return out
    end
    result = minimize(minimizer, pval2resid, mdc.params)
    pval2resid(result.best)
    ProgressMeter.finish!(prog)
    save_bestfit!(mdc, model)

    # Prepare output
    comps = OrderedDict{Symbol, BestFitComp}()
    i = 1
    unc = deepcopy(result.unc)
    for (cname, ceval) in model.cevals
        comps[cname] = BestFitComp()
        for (pid, par) in ceval.params
            if (!par.fixed)  &&  (model.cevals[cname].cfixed == 0)
                bfpar = BestFitParam(model.peval.pvalues[i], popfirst!(unc), false, model.peval.patched[i])
            else
                bfpar = BestFitParam(model.peval.pvalues[i], NaN           , true , model.peval.patched[i])
            end
            push!(comps[cname], pid.name, bfpar)
            i += 1
        end
    end

    elapsed = now() - timestamp
    @assert isa(elapsed, Millisecond)
    return BestFitResult(timestamp, elapsed.value / 1.e3, result,
                         comps, mdc)
end


function fit!(multi::MultiModel, data::Vector{Measures{N}};
              minimizer=lsqfit()) where N
    timestamp = now()
    evaluate!(multi)
    mdc = MDMultiComparison(multi, data)
    @assert length(mdc.ifree) > 0 "No free parameter in the model"
    @assert mdc.dof >= 1

    prog = ProgressUnknown("Model evaluations:", dt=0.5, showspeed=true)
    function pval2resid(pvalues::Vector{Float64})
        out = try_pvalues(mdc, multi, pvalues)
        evaluate_showvalues(x) = () -> [(:red_chisq, sum(abs2.(x)) / mdc.dof)]
        ProgressMeter.next!(prog; showvalues = evaluate_showvalues(mdc.residuals))
        return mdc.residuals
    end
    result = minimize(minimizer, pval2resid, mdc.params)
    pval2resid(result.best)
    ProgressMeter.finish!(prog)
    save_bestfit!(mdc, multi)

    # Prepare output
    bfmodels = Vector{OrderedDict{Symbol, BestFitComp}}()
    unc = deepcopy(result.unc)
    for id in 1:length(multi)
        model = multi[id]
        comps = OrderedDict{Symbol, BestFitComp}()
        i = 1
        for (cname, ceval) in model.cevals
            comps[cname] = BestFitComp()
            for (pid, par) in ceval.params
                if (!par.fixed)  &&  (model.cevals[cname].cfixed == 0)
                    bfpar = BestFitParam(model.peval.pvalues[i], popfirst!(unc), false, model.peval.patched[i])
                else
                    bfpar = BestFitParam(model.peval.pvalues[i], NaN           , true , model.peval.patched[i])
                end
                push!(comps[cname], pid.name, bfpar)
                i += 1
            end
        end
        push!(bfmodels, comps)
    end

    elapsed = now() - timestamp
    @assert isa(elapsed, Millisecond)
    return BestFitMultiResult(timestamp, elapsed.value / 1.e3, result,
                              bfmodels, mdc)
end
