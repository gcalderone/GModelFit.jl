# --------------------------------------------------------------------
mutable struct MDComparison
    params::Vector{Parameter}
    ifree::Vector{Int}
    data::Measures{1}
    ndata::Int
    nfree::Int
    dof::Int
    residuals::Vector{Float64}
    fitstat::Float64
    log10testprob::Float64

    function MDComparison(model::Model, data::Measures{N}) where N
        evaluate!(model)
        params = Vector{Parameter}()
        ifree = Vector{Int}()
        i = 1
        for (cname, ceval) in model.cevals
            for (pid, par) in ceval.params
                if !(par.low <= par.val <= par.high)
                    s = "Value outside limits for param [$(cname)].$(pid.name):\n" * string(par)
                    error(s)
                end
                if (!par.fixed)  &&  (model.cevals[cname].cfixed == 0)
                    push!(params, par)
                    push!(ifree, i)
                end
                i += 1
            end
        end
        
        data1d = flatten(data)
        mdc = new(params, ifree, data1d, length(data1d), length(ifree),
                  length(data1d) - length(ifree),
                  fill(NaN, length(data1d)), NaN, NaN)
        update_residuals!(mdc, model)
        update_stats!(mdc)
        return mdc
    end
end

function update_residuals!(mdc::MDComparison, model::Model)
    mdc.residuals .= (model() .- mdc.data.val) ./ mdc.data.unc
end

function try_pvalues(mdc::MDComparison, model::Model, pvalues::Vector{Float64})
    model.peval.pvalues[mdc.ifree] .= pvalues
    patch_params(model)
    quick_evaluate(model)
    update_residuals!(mdc, model)
    return mdc.residuals
end

function save_bestfit!(mdc::MDComparison, model::Model)
    # Copy best fit values back into components
    i = 1
    for (cname, ceval) in model.cevals
        for (pid, par) in ceval.params
            par.val = model.peval.pvalues[i]
            i += 1
        end
    end
    update_stats!(mdc)
end


# --------------------------------------------------------------------
mutable struct MDMultiComparison
    params::Vector{Parameter}
    modelid::Vector{Int}
    ifree::Vector{Int}
    data::Vector{Measures{1}}
    ndata::Int
    nfree::Int
    dof::Int
    residuals::Vector{Float64}
    fitstat::Float64
    log10testprob::Float64

    function MDMultiComparison(multi::MultiModel, data::Vector{Measures{N}}) where N
        evaluate!(multi)
        params = Vector{Parameter}()
        modelid = Vector{Int}()
        ifree = Vector{Int}()
        for id in 1:length(multi)
            model = multi[id]
            i = 1
            for (cname, ceval) in model.cevals
                for (pid, par) in ceval.params
                    if !(par.low <= par.val <= par.high)
                        s = "Value outside limits for param [$(cname)].$(pid.name):\n" * string(par)
                        error(s)
                    end
                    if (!par.fixed)  &&  (model.cevals[cname].cfixed == 0)
                        push!(params, par)
                        push!(modelid, id)
                        push!(ifree, i)
                    end
                    i += 1
                end
            end
        end

        data1d = [flatten(data[id]) for id in 1:length(multi)]
        ndata = sum(length.(data1d))
        mdc = new(params, modelid, ifree, data1d, ndata, length(ifree),
                  ndata - length(ifree),
                  fill(NaN, ndata), NaN, NaN)
        update_residuals!(mdc, multi)
        update_stats!(mdc)
    end
end

function update_residuals!(mdc::MDMultiComparison, multi::MultiModel)
    i1 = 1
    for id in 1:length(multi)
        nn = length(mdc.data[id])
        i2 = i1 + nn - 1
        mdc.residuals[i1:i2] .= (multi[id]() .- mdc.data[id].val) ./ mdc.data[id].unc
        i1 += nn
    end
    return mdc
end

function try_pvalues(mdc::MDMultiComparison, multi::MultiModel, pvalues::Vector{Float64})
    for i in 1:length(mdc.ifree)
        multi[mdc.modelid[i]].peval.pvalues[mdc.ifree[i]] = pvalues[i]
    end
    patch_params(multi)
    quick_evaluate(multi)
    update_residuals!(mdc, multi)
    return mdc.residuals
end

function save_bestfit!(mdc::MDMultiComparison, multi::MultiModel)
    # Copy best fit values back into components
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
    update_stats!(mdc)
end


function update_stats!(mdc::Union{MDComparison, MDMultiComparison})
    mdc.fitstat = sum(abs2, mdc.residuals)
    mdc.log10testprob = logccdf(Chisq(mdc.dof), mdc.fitstat) * log10(exp(1))
    return mdc
end
