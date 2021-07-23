# --------------------------------------------------------------------
struct FitResult
    timestamp::DateTime
    elapsed::Float64
    ndata::Int
    nfree::Int
    dof::Int
    fitstat::Float64
    gofstat::Float64
    log10testprob::Float64
    resid::Vector{Float64}
    mzer::Union{Nothing, AbstractMinimizerStatus}
end


function fit!(model::Model, data::Measures{N};
              minimizer=lsqfit(), dry=false) where N
    timestamp = now()
    evaluate!(model)

    data1d = flatten(data)
    resid1d = fill(NaN, length(data1d))
    params = model.meval.params[model.meval.ifree]
    @assert length(params) > 0 "No free parameter in the model"

    if !dry
        dof = length(resid1d) - length(params)
        prog = ProgressUnknown("Model (dof=$dof) evaluations:", dt=0.5, showspeed=true)
        function private_func(pvalues::Vector{Float64})
            model.meval.pvalues[model.meval.ifree] .= pvalues
            patch_params(model)
            quick_evaluate(model)
            resid1d .= (model() .- data1d.val) ./ data1d.unc

            evaluate_showvalues(x) = () -> begin
                dof = (length(resid1d) - length(pvalues))
                [(:fit_stat, sum(abs2, x) / dof)]
            end
            ProgressMeter.next!(prog; showvalues=evaluate_showvalues(resid1d))
            return resid1d
        end
        result = minimize(minimizer, private_func, model.meval.params[model.meval.ifree])
        private_func(result.best)
        update_params!(model, result.unc)
        ProgressMeter.finish!(prog)
    else
        private_func(getfield.(params, :val))
    end

    # Prepare output
    ndata = length(resid1d)
    nfree = length(params)
    dof = ndata - nfree
    fitstat = sum(abs2, resid1d)
    gofstat = sum(abs2, resid1d)
    elapsed = now() - timestamp
    @assert isa(elapsed, Millisecond)
    return FitResult(timestamp, elapsed.value / 1.e3,
                     ndata, nfree, dof, fitstat, gofstat,
                     logccdf(Chisq(dof), gofstat) * log10(exp(1)),
                     resid1d, (dry  ?  nothing  :  result))
end


function fit!(multi::MultiModel, data::Vector{Measures{N}};
              minimizer=lsqfit(), dry=false) where N
    timestamp = now()
    evaluate!(multi)

    data1d = [flatten(data[id]) for id in 1:length(multi)]
    resid1d = fill(NaN, sum(length.(data1d)))
    params = Vector{Parameter}()
    for id in 1:length(multi)
        append!(params, multi[id].meval.params[multi[id].meval.ifree])
    end
    @assert length(params) > 0 "No free parameter in the model"

    if !dry
        dof = length(resid1d) - length(params)
        prog = ProgressUnknown("Model (dof=$dof) evaluations:", dt=0.5, showspeed=true)
        function private_func(pvalues::Vector{Float64})
            i1 = 1
            for id in 1:length(multi)
                nn = length(multi[id].meval.ifree)
                i2 = i1 + nn - 1
                multi[id].meval.pvalues[multi[id].meval.ifree] = pvalues[i1:i2]
                i1 += nn
            end
            patch_params(multi)
            quick_evaluate(multi)

            i1 = 1
            for id in 1:length(multi)
                nn = length(data1d[id])
                i2 = i1 + nn - 1
                resid1d[i1:i2] .= (multi[id]() .- data1d[id].val) ./ data1d[id].unc
                i1 += nn
            end
            evaluate_showvalues(x) = () -> begin
                dof = (length(resid1d) - length(pvalues))
                [(:fit_stat, sum(abs2, x) / dof)]
            end
            ProgressMeter.next!(prog; showvalues=evaluate_showvalues(resid1d))
            return resid1d
        end
        result = minimize(minimizer, private_func, params)
        private_func(result.best)
        i1 = 1
        for id in 1:length(multi)
            nn = length(multi[id].meval.ifree)
            i2 = i1 + nn - 1
            update_params!(multi[id], result.unc[i1:i2])
            i1 += nn
        end
        ProgressMeter.finish!(prog)        
    else
        private_func(getfield.(params, :val))
    end

    # Prepare output
    ndata = length(resid1d)
    nfree = length(params)
    dof = ndata - nfree
    fitstat = sum(abs2, resid1d)
    gofstat = sum(abs2, resid1d)
    elapsed = now() - timestamp
    @assert isa(elapsed, Millisecond)
    return FitResult(timestamp, elapsed.value / 1.e3,
                     ndata, nfree, dof, fitstat, gofstat,
                     logccdf(Chisq(dof), gofstat) * log10(exp(1)),
                     resid1d, (dry  ?  nothing  :  result))
end



function free_param_names(model::Model)
    out = Vector{String}()
    for (cname, ceval) in model.cevals
        for (pid, par) in ceval.params
            if (!par.fixed)  &&  (model.cevals[cname].cfixed == 0)
                if pid.index == 0
                    push!(out, "[:$(cname)].$(pid.name)")
                else
                    push!(out, "[:$(cname)].$(pid.name)[$(pid.index)]")
                end
            end
        end
    end
    return out
end

function free_param_names(multi::MultiModel)
    out = Vector{String}()
    for id in 1:length(multi)
        append!(out, "[$(string(id))]" .* free_param_names(multi[id]))
    end
    return out
end


function print_param_covariance(model::Union{Model, MultiModel}, fitres::FitResult;
                                param::Union{Nothing, String}=nothing, sort=false, thresh=0.)
    @assert isa(fitres.mzer.specific, CMPFit.Result)
    if isa(model, Model)
        @assert isnothing(model.parent)
    end
    names = free_param_names(model)
    @assert length(names)^2 == length(fitres.mzer.specific.covar)
    ii = Vector{Int}()
    jj = Vector{Int}()
    covar = Vector{Float64}()
    for i in 1:length(names)
        for j in i+1:length(names)
            push!(covar, fitres.mzer.specific.covar[i, j])
            push!(ii, i)
            push!(jj, j)
        end
    end
    if sort
        ii    = ii[   sortperm(abs.(covar))]
        jj    = jj[   sortperm(abs.(covar))]
        covar = covar[sortperm(abs.(covar))]
    end
    for i in 1:length(ii)
        if !isnothing(param)
            (names[ii[i]] != param)  &&  continue
        end
        (abs(covar[i]) < thresh)  &&  continue
        @printf "%-30s  %-30s  %10.4f\n" names[ii[i]] names[jj[i]] covar[i]
    end
end
