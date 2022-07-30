
# ====================================================================
struct FitProblem{T <: AbstractMeasures} <: AbstractFitProblem
    model::Model
    measures::AbstractMeasures
    resid::Vector{Float64}
    nfree::Int
    dof::Int

    function FitProblem(model::Model, data::T) where T <: AbstractMeasures
        evaluate(model)
        resid = fill(NaN, length(data))
        nfree = length(free_params(model))
        @assert nfree > 0 "No free parameter in the model"
        return new{T}(model, data, resid, nfree, length(resid) - nfree)
    end
end

free_params(fp::FitProblem) = free_params(fp.model)
residuals(fp::FitProblem) = fp.resid

function update!(fp::FitProblem, pvalues::Vector{Float64})
    eval_step1(fp.model, pvalues)
    eval_step2(fp.model)
    eval_step3(fp.model)
    update_residuals!(fp)
    return fp.resid
end

function update_residuals!(fp::FitProblem{Measures{N}}) where N
    fp.resid .= (fp.model() .- values(fp.measures)) ./ uncerts(fp.measures)
end
fit_stat(fp::FitProblem{Measures{N}}) where N =
    sum(abs2, fp.resid) / fp.dof

function finalize!(fp::FitProblem, best::Vector{Float64}, uncerts::Vector{Float64})
    @assert fp.nfree == length(best) == length(uncerts)
    update!(fp, best)
    eval_step4(fp.model, uncerts)
end

function error!(fp::FitProblem)
    eval_step4(fp.model, fill(NaN, fp.nfree))
end



# ====================================================================
struct MultiFitProblem <: AbstractFitProblem
    multi::MultiModel
    fp::Vector{FitProblem}
    resid::Vector{Float64}
    nfree::Int
    dof::Int

    function MultiFitProblem(multi::MultiModel, datasets::Vector{T}) where T <: AbstractMeasures
        @assert length(multi) == length(datasets)
        evaluate(multi)
        fp = [FitProblem(multi[id], datasets[id]) for id in 1:length(multi)]
        resid = fill(NaN, sum(length.(getfield.(fp, :resid))))
        nfree = sum(getfield.(fp, :nfree))
        @assert nfree > 0 "No free parameter in the model"
        return new(multi, fp, resid, nfree, length(resid) - nfree)
    end
end

free_params(fp::MultiFitProblem) = free_params(fp.multi)
residuals(fp::MultiFitProblem) = fp.resid

function update!(fp::MultiFitProblem, pvalues::Vector{Float64})
    # We need to copy all parameter values before evaluation to ensure
    # all patch functions use the current parameter values
    for (id, i1, i2) in free_params_indices(fp.multi)
        eval_step1(fp.multi[id], pvalues[i1:i2])
    end
    eval_step2(fp.multi)
    eval_step3(fp.multi)

    # Populate resid vector
    i1 = 1
    for id in 1:length(fp.multi)
        update_residuals!(fp.fp[id])
        nn = length(fp.fp[id].resid)
        if nn > 0
            i2 = i1 + nn - 1
            fp.resid[i1:i2] .= fp.fp[id].resid
            i1 += nn
        end
    end
    return fp.resid
end

# TODO: Handle the case where one (or more) dataset is not a `Measures`
fit_stat(fp::MultiFitProblem) =
    sum(abs2, fp.resid) / fp.dof

function finalize!(fp::MultiFitProblem, best::Vector{Float64}, unc::Vector{Float64})
    for (id, i1, i2) in free_params_indices(fp.multi)
        finalize!(fp.fp[id], best[i1:i2], unc[i1:i2])
    end
end

function error!(fp::MultiFitProblem)
    for id in 1:length(fp.multi)
        error!(fp.fp[id])
    end
end


# ====================================================================
struct FitResult
    timestamp::DateTime
    elapsed::Float64
    ndata::Int
    nfree::Int
    dof::Int
    fitstat::Float64
    gofstat::Float64
    log10testprob::Float64
    status::AbstractMinimizerStatus

    function FitResult(timestamp::DateTime, fp::AbstractFitProblem, status::AbstractMinimizerStatus)
        gof_stat = sum(abs2, residuals(fp))
        tp = NaN
        try
            tp = logccdf(Chisq(fp.dof), gof_stat) * log10(exp(1))
        catch; end
        new(timestamp, (now() - timestamp).value / 1e3,
            length(residuals(fp)), length(residuals(fp)) - fp.dof, fp.dof,
            gof_stat, gof_stat, tp, status)
    end
end



# ====================================================================
fit!(model::Model, data::Measures) =
    fit!(lsqfit(), model, data)

function fit!(mzer::AbstractMinimizer, model::Model, data::Measures)
    ts = now()
    fp = FitProblem(model, data)
    status = fit!(mzer, fp)
    return FitResult(ts, fp, status)
end


fit!(multi::MultiModel, data::Vector{Measures{N}}) where N =
    fit!(lsqfit(), multi, data)

function fit!(mzer::AbstractMinimizer, multi::MultiModel, data::Vector{Measures{N}}) where N
    ts = now()
    fp = MultiFitProblem(multi, data)
    status = fit!(mzer, fp)
    return FitResult(ts, fp, status)
end


# ====================================================================
function print_param_covariance(model::Union{Model, MultiModel}, fitres::FitResult;
                                param::Union{Nothing, String}=nothing, sort=false, thresh=0.)
    @assert isa(fitres.mzer.internal, CMPFit.Result)
    if isa(model, Model)
        @assert isnothing(model.parent)
    end
    names = free_param_names(model)
    @assert length(names)^2 == length(fitres.mzer.internal.covar)
    ii = Vector{Int}()
    jj = Vector{Int}()
    covar = Vector{Float64}()
    for i in 1:length(names)
        for j in i+1:length(names)
            push!(covar, fitres.mzer.internal.covar[i, j])
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
