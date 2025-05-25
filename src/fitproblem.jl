# ====================================================================
struct MultiModelEval{N}
    v::Vector{ModelEval}

    function MultiModelEval(models::Vector{Model}, domains::Vector{<: AbstractDomain})
        @assert length(models) == length(domains)
        return MultiModelEval([ModelEval(models[i], domains[i]) for i in 1:length(models)])
    end

    function MultiModelEval(mevals::Vector{ModelEval})
        out = new{length(mevals)}(mevals)
        scan_model!(out)
        return out
    end
end
length(mevals::MultiModelEval) = length(mevals.v)
getindex(mevals::MultiModelEval, i) = mevals.v[i]


function free_params_indices(mevals::MultiModelEval)
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(mevals)
        nn = length(mevals[id].ifree)
        if nn > 0
            i2 = i1 + nn - 1
            push!(out, (id, i1, i2))
            i1 += nn
        end
    end
    return out
end


scan_model!(mevals::MultiModelEval{1}) = scan_model!(mevals[1])
function scan_model!(mevals::MultiModelEval)
    for i in 1:length(mevals)
        scan_model!(mevals[i])
    end
    pv1 = [mevals[i].tpar.pvalues for i in 1:length(mevals)]
    for i in 1:length(mevals)
        empty!( mevals[i].tpar.pvmulti)
        append!(mevals[i].tpar.pvmulti, pv1)
    end
    pv2 = [mevals[i].tparad.pvalues for i in 1:length(mevals)]
    for i in 1:length(mevals)
        empty!( mevals[i].tparad.pvmulti)
        append!(mevals[i].tparad.pvmulti, pv2)
    end
end


function free_params(mevals::MultiModelEval)
    out = Vector{Parameter}()
    for id in 1:length(mevals)
        append!(out, free_params(mevals[id]))
    end
    return out
end
free_params_val(mevals::MultiModelEval) = getfield.(free_params(mevals), :val)
nfree(mevals::MultiModelEval) = sum(nfree.(mevals.v))


function set_pvalues!(mevals::MultiModelEval{N}, pvalues::AbstractVector) where N
    if N == 1
        set_pvalues!(mevals.v[1], pvalues)
    else
        for (id, i1, i2) in free_params_indices(mevals)
            set_pvalues!(mevals[id], pvalues[i1:i2])
        end
    end
end

function update_eval!(mevals::MultiModelEval, pvalues::AbstractVector{Float64})
    set_pvalues!(mevals, pvalues)
    return update_eval!.(mevals.v)
end

function update_eval!(mevals::MultiModelEval{N}, pvalues::AbstractVector) where N
    set_pvalues!(mevals, pvalues)
    return update_eval_ad!.(mevals.v)
end

update_eval!(mevals::MultiModelEval) = update_eval!.(mevals.v)


# ====================================================================
# Fit statistics
abstract type AbstractFitStat end
abstract type ChiSquared <: AbstractFitStat end

# ====================================================================
struct FitProblem{T <: AbstractFitStat}
    mevals::MultiModelEval
    data::Vector{<: AbstractMeasures}
    bestfit::Vector{PVModel{Parameter}}
    buffer::Vector{Float64}

    FitProblem(models::Vector{Model}, datasets::Vector{Measures{N}}) where N = FitProblem(MultiModelEval(models, getfield.(datasets, :domain)), datasets)

    function FitProblem(mevals::MultiModelEval, datasets::Vector{Measures{N}}) where N
        @assert length(mevals) == length(datasets)
        fp = new{ChiSquared}(mevals, datasets, Vector{PVModel{Parameter}}(), Vector{Float64}(undef, sum(length.(datasets))))
        update_eval!(fp, fp.buffer, free_params_val(fp.mevals))
        return fp
    end
end

free_params(fitprob::FitProblem) = free_params(fitprob.mevals)
nfree(fitprob::FitProblem) = nfree(fitprob.mevals)
ndata(fitprob::FitProblem) = sum(length.(fitprob.data))

function set_bestfit!(fitprob::FitProblem, pvalues::Vector{Float64}, puncerts::Vector{Float64})
    update_eval!(fitprob.mevals, pvalues)

    empty!(fitprob.bestfit)
    for id in 1:length(fitprob.mevals)
        push!(fitprob.bestfit, PVModel{Parameter}())
    end

    for (id, i1, i2) in free_params_indices(fitprob.mevals)
        meval = fitprob.mevals[id]
        i = 1
        for (cname, comp) in meval.model.comps
            for (pname, _par) in getparams(comp)
                par = deepcopy(_par)
                par.val    = meval.tpar.pvalues[cname][pname]
                par.actual = meval.tpar.pactual[cname][pname]
                push!(fitprob.bestfit[id], cname, pname, par)
                if length(fitprob.bestfit[id].data) in meval.ifree
                    par.unc = puncerts[i1:i2][i]
                    par.fixed = false
                    i += 1
                else
                    par.fixed = true
                end
            end
        end
    end
    nothing
end


# FitProblem{ChiSquared} specific methods
dof(fitprob::FitProblem{ChiSquared}) = ndata(fitprob) - nfree(fitprob)
fitstat(fitprob::FitProblem{ChiSquared}) = sum(abs2, fitprob.buffer) / dof(fitprob)

function update_eval!(fitprob::FitProblem{ChiSquared}, output::Vector{T}, pvalues::Vector{T}) where T
    evals = update_eval!(fitprob.mevals, pvalues)
    i1 = 1
    for i in 1:length(fitprob.mevals)
        nn = length(evals[i])
        if nn > 0
            i2 = i1 + nn - 1
            output[i1:i2] .= reshape((reshape(fitprob.mevals[i].domain, evals[i]) .- values(fitprob.data[i])) ./ uncerts(fitprob.data[i]), :)
            i1 += nn
        end
    end

    (T == Float64)  &&  (fitprob.buffer .= output)
    return output
end
