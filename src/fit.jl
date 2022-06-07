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
    status::AbstractMinimizerStatus

    function FitResult(fd::AbstractFitData, status::AbstractMinimizerStatus)
        gof_stat = sum(abs2, residuals(fd))
        tp = NaN
        try
            tp = logccdf(Chisq(fd.dof), gof_stat) * log10(exp(1))
        catch; end
        new(fd.timestamp, (now() - fd.timestamp).value / 1e3,
            length(residuals(fd)), length(residuals(fd)) - fd.dof, fd.dof,
            gof_stat, gof_stat, tp, status)
    end
end



fit!(model::Model, data::Measures) =
    fit!(lsqfit(), model, data)

function fit!(mzer::AbstractMinimizer, model::Model, data::Measures)
    fd = FitData(model, data)
    status = fit!(mzer, fd)
    return FitResult(fd, status)
end


fit!(multi::MultiModel, data::Vector{Measures{N}}) where N =
    fit!(lsqfit(), multi, data)

function fit!(mzer::AbstractMinimizer, multi::MultiModel, data::Vector{Measures{N}}) where N
    fd = MultiFitData(multi, data)
    status = fit!(mzer, fd)
    return FitResult(fd, status)
end


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
